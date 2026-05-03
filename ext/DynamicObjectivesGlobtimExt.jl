"""
DynamicObjectivesGlobtimExt — TOML-driven experiment pipeline

This package extension is loaded when both `Globtim` and `GlobtimPostProcessing`
are available. It provides the `run_experiment_from_config()` orchestrator that
ties together TOML config parsing, model resolution (ODE catalogue or analytical
benchmark functions), polynomial approximation, HC critical point solving,
optional NelderMead/BFGS refinement, and optional Newton-based CP validation.

Triggered by: using DynamicObjectives; using Globtim; using GlobtimPostProcessing
"""
module DynamicObjectivesGlobtimExt

using DynamicObjectives
using Globtim
using GlobtimPostProcessing
using Optim: NelderMead, BFGS, LBFGS
using Printf

import DynamicObjectives: run_experiment_from_config, _resolve_solver
import DynamicObjectives:
    Tsit5, Vern7, Vern9, Rodas5, Rosenbrock23, TRBDF2, KenCarp4, AutoTsit5

# ═══════════════════════════════════════════════════════════════════════════════
# run_experiment_from_config — TOML-driven experiment orchestrator
# ═══════════════════════════════════════════════════════════════════════════════

"""
    run_experiment_from_config(path::String; io::IO=stdout) -> Dict{Symbol, Any}

Run a complete experiment pipeline from a TOML configuration file.

This is the single entry point for all TOML-driven experiments. It:
1. Parses and validates the TOML config via `Globtim.load_experiment_config`
2. Resolves the model (catalogue entry OR analytical function from FUNCTION_REGISTRY)
3. Builds the objective function and domain bounds
4. Resolves output directory
5. Copies the TOML config into output dir for reproducibility
6. Runs `Globtim.run_standard_experiment` for polynomial approximation + HC solving
7. Optionally runs NelderMead/BFGS refinement via GlobtimPostProcessing
8. Optionally validates raw CPs via Newton on ∇f=0, deduplicates, classifies via Hessian

The output directory contains all artifacts needed for standalone post-processing:
`experiment_config.toml`, `results_summary.jld2`, `results_summary.json`, and per-degree CSVs.

Requires `using Globtim` and `using GlobtimPostProcessing` to be loaded.

# Arguments
- `path::String`: Path to the TOML config file
- `io::IO`: Output stream for progress messages (default: stdout)

# Returns
A Dict with keys:
- `:config` — the parsed `ExperimentPipelineConfig`
- `:experiment_result` — return value from `run_standard_experiment`
- `:degree_results` — vector of `DegreeResult` from the experiment
- `:output_dir` — path to the output directory
- `:entry` — the `CatalogueEntry` (if catalogue mode, else `nothing`)
- `:refinement_results` — refinement results (if enabled, else `nothing`)
- `:known_cps` — `KnownCriticalPoints` from Newton analysis (if enabled, else `nothing`)
- `:objective` — the objective function (for downstream analysis)
- `:bounds` — the domain bounds used

# Example
```julia
using DynamicObjectives, Globtim, GlobtimPostProcessing

result = run_experiment_from_config("examples/configs/levy_3d.toml")
result[:degree_results]  # Vector of DegreeResult
result[:known_cps]       # KnownCriticalPoints from Newton analysis
```
"""
function DynamicObjectives.run_experiment_from_config(path::String; io::IO = stdout)
    # ── 1. Load and validate config ──
    config = Globtim.load_experiment_config(path)
    println(io, "Loaded config: $(config.name)")
    config.description != "" && println(io, "  $(config.description)")

    # ── 2. Resolve model ──
    entry = nothing
    objective = nothing
    bounds = nothing
    p_true = nothing
    obj_name = config.name

    if config.catalogue_path !== nothing
        # ── 2a. Catalogue mode: load entry, build TolerantObjective ──
        entries = load_catalogue(config.catalogue_path; model_name = config.entry_name)
        matching = filter(e -> e.name == config.entry_name, entries)
        isempty(matching) && error(
            "Entry '$(config.entry_name)' not found in catalogue '$(config.catalogue_path)'. " *
            "Available: $(join([e.name for e in entries], ", "))",
        )
        entry = first(matching)

        # p_true: config overrides catalogue
        p_true = config.p_true !== nothing ? config.p_true : entry.p_true

        println(io, "  Model: $(config.entry_name) from $(config.catalogue_path)")
        if config.p_true !== nothing
            println(
                io,
                "  p_true: $p_true (from config, overriding catalogue $(entry.p_true))",
            )
        else
            println(io, "  p_true: $p_true")
        end
        if config.time_interval !== nothing
            println(
                io,
                "  time_interval: $(config.time_interval) (overriding catalogue default $(entry.time_interval))",
            )
        end
        if config.sample_times !== nothing
            println(
                io,
                "  sample_times: $(length(config.sample_times)) explicit time points [$(config.sample_times[1]), ..., $(config.sample_times[end])]",
            )
        end

        # Build TolerantObjective from catalogue entry
        model, _, _, outputs = entry.model_fn()

        # Resolve solver settings (config overrides, then defaults)
        solver =
            config.solver_method !== nothing ? _resolve_solver(config.solver_method) :
            Tsit5()
        abstol = config.solver_abstol !== nothing ? config.solver_abstol : 1e-4
        reltol = config.solver_reltol !== nothing ? config.solver_reltol : 1e-4

        # sample_times overrides both time_interval and numpoints
        if config.sample_times !== nothing
            time_interval = [config.sample_times[1], config.sample_times[end]]
            numpoints = length(config.sample_times)
            uneven_sampling_times = config.sample_times
        else
            # numpoints: TOML config overrides catalogue entry default
            numpoints =
                config.solver_numpoints !== nothing ? config.solver_numpoints :
                entry.numpoints
            # time_interval: TOML config overrides catalogue entry default
            time_interval =
                config.time_interval !== nothing ? config.time_interval :
                entry.time_interval
            uneven_sampling_times = Float64[]
        end

        objective = TolerantObjective(
            model,
            outputs,
            entry.ic,
            p_true,
            time_interval,
            numpoints,
            entry.distance_function,
            entry.aggregate_distances;
            solver = solver,
            abstol = abstol,
            reltol = reltol,
            uneven_sampling_times = uneven_sampling_times,
        )

        # Domain center: config p_center → p_true (which may itself be overridden)
        p_center = config.p_center !== nothing ? config.p_center : p_true
        if config.p_center !== nothing
            println(
                io,
                "  p_center: $p_center (from config, domain centered away from p_true)",
            )
        end

        # Build bounds from config domain specification
        if config.radius !== nothing
            bounds = build_bounds(p_center, config.radius)
        elseif config.radii !== nothing
            bounds = build_bounds(p_center, config.radii)
        elseif config.bounds !== nothing
            bounds = config.bounds
        else
            bounds = entry.bounds
        end

        obj_name = config.entry_name
    else
        # ── 2b. Analytical mode: look up function from FUNCTION_REGISTRY ──
        func_name = config.analytical_function
        dim = config.dimension
        println(io, "  Analytical function: $func_name ($(dim)D)")

        bench = Globtim.get_benchmark_config_by_name(func_name, dim)
        objective = bench.objective
        obj_name = bench.name

        # Bounds: TOML config overrides registry defaults
        if config.bounds !== nothing
            bounds = config.bounds
        else
            bounds = bench.bounds
        end

        println(io, "  Bounds: $(join(["[$(lb), $(ub)]" for (lb, ub) in bounds], " × "))")
    end

    # ── 3. Build ExperimentParams ──
    experiment_params = Globtim.config_to_experiment_params(config)
    println(io, "  GN=$(config.GN), degrees=$(config.degree_range), basis=$(config.basis)")
    println(io, "  Dimension: $(length(bounds))D")

    # ── 4. Resolve output directory ──
    output_dir = if config.output_dir !== nothing
        config.output_dir
    else
        mktempdir()
    end

    # ── 5. Copy TOML config to output directory for reproducibility ──
    mkpath(output_dir)
    cp(
        realpath(path),
        joinpath(output_dir, "experiment_config.toml");
        force = true,
        follow_symlinks = true,
    )

    # ── 6. Run Globtim experiment ──
    println(io, "Running Globtim experiment...")
    exp_result = Globtim.run_standard_experiment(
        objective_function = objective,
        objective_name = obj_name,
        bounds = bounds,
        experiment_config = experiment_params,
        output_dir = output_dir,
        metadata = Dict{String,Any}(
            "config_file" => abspath(path),
            "config_name" => config.name,
            # Reconstruction metadata — enables standalone post-processing
            "analytical_function" => config.analytical_function,   # e.g. "Levy" or nothing
            "dimension" => config.dimension,                       # e.g. 3 or nothing
            "catalogue_path" => config.catalogue_path,             # e.g. "path/to/catalogue.jsonl" or nothing
            "entry_name" => config.entry_name,                     # e.g. "lv4d_basic" or nothing
        ),
        true_params = p_true,
    )

    degree_results = exp_result[:degree_results]
    println(io, "  Completed: $(length(degree_results)) degrees processed")

    # ── 7. Optional NelderMead/BFGS refinement ──
    refinement_results = nothing
    if config.refinement_enabled
        println(io, "Running refinement...")

        # Switch to tight tolerances for ODE refinement (Vern9/1e-10)
        if objective isa TolerantObjective
            set_tolerance!(objective, 1e-10)
            set_solver!(objective, Vern9())
        end

        # Build refinement config
        # Method dispatch (4xgu): NelderMead (default), BFGS, LBFGS supported in
        # the cluster path. Newton variants are not yet in the cluster runner —
        # use experiments/sandbox/run_refinement_shootout.jl for Newton comparison.
        ref_method = if config.refinement_method == "BFGS"
            BFGS()
        elseif config.refinement_method == "LBFGS"
            LBFGS()
        else
            NelderMead()
        end

        ref_max_time =
            config.refinement_max_time !== nothing ? config.refinement_max_time : 60.0
        ref_grad_method =
            config.refinement_gradient_method !== nothing ?
            Symbol(config.refinement_gradient_method) : :finitediff
        ref_grad_tol =
            config.refinement_gradient_tolerance !== nothing ?
            config.refinement_gradient_tolerance : 1e-4
        # Newton-step-norm tolerance (62qv): scale-invariant CP criterion. Default
        # 1e-4 matches "a Newton step would move us less than 1e-4."
        ref_step_tol =
            config.refinement_step_tolerance !== nothing ?
            config.refinement_step_tolerance : 1e-4

        ref_config = GlobtimPostProcessing.RefinementConfig(
            method = ref_method,
            max_time_per_point = ref_max_time,
            gradient_method = ref_grad_method,
            gradient_tolerance = ref_grad_tol,
            step_tolerance = ref_step_tol,
            bounds = bounds,
        )

        refinement_results = GlobtimPostProcessing.run_degree_analyses(
            degree_results,
            objective,
            output_dir,
            ref_config;
            gradient_method = ref_grad_method,
            io = io,
        )

        println(io, "  Refinement complete: $(length(refinement_results)) degrees refined")
    end

    # ── 8. Optional Newton-based CP validation ──
    known_cps = nothing
    if config.analysis_enabled
        # Find highest successful degree with critical points
        highest_dr = nothing
        for dr in reverse(degree_results)
            if dr.status == "success" && dr.n_critical_points > 0
                highest_dr = dr
                break
            end
        end

        if highest_dr === nothing
            println(io, "  WARNING: No successful degrees with CPs — skipping analysis")
        else
            # Resolve analysis parameters (config overrides, then defaults)
            ana_refinement_goal =
                config.analysis_refinement_goal !== nothing ?
                Symbol(config.analysis_refinement_goal) : :minimum
            ana_grad_method =
                config.analysis_gradient_method !== nothing ?
                Symbol(config.analysis_gradient_method) : :forwarddiff
            ana_newton_tol =
                config.analysis_newton_tol !== nothing ? config.analysis_newton_tol : 1e-8
            ana_newton_max_iter =
                config.analysis_newton_max_iterations !== nothing ?
                config.analysis_newton_max_iterations : 200
            ana_max_time_pt =
                config.analysis_max_time_per_point !== nothing ?
                config.analysis_max_time_per_point : 60.0
            ana_hessian_tol =
                config.analysis_hessian_tol !== nothing ? config.analysis_hessian_tol : 1e-6
            ana_dedup_fraction =
                config.analysis_dedup_fraction !== nothing ?
                config.analysis_dedup_fraction : 0.02
            ana_top_k = config.analysis_top_k
            ana_accept_tol =
                config.analysis_accept_tol !== nothing ? config.analysis_accept_tol : 1e-2
            ana_f_accept_tol = config.analysis_f_accept_tol  # nothing by default

            # Phase header
            goal_label =
                ana_refinement_goal == :minimum ? "NelderMead (f-minimization)" :
                "Newton (∇f = 0)"
            println(io)
            println(
                io,
                "══ CP Refinement: $goal_label ══════════════════════════════════════",
            )
            println(
                io,
                "  Source: degree $(highest_dr.degree) ($(highest_dr.n_critical_points) raw CPs)",
            )
            top_k_str = ana_top_k === nothing ? "none (refining all)" : string(ana_top_k)
            Printf.@printf(
                io,
                "  Params: goal=%s, tol=%.0e, accept_tol=%.0e, max_iter=%d, gradient=%s\n",
                ana_refinement_goal,
                ana_newton_tol,
                ana_accept_tol,
                ana_newton_max_iter,
                ana_grad_method
            )
            if ana_f_accept_tol !== nothing
                Printf.@printf(
                    io,
                    "  f_accept_tol=%.0e (accept CPs with f(x) below this value)\n",
                    ana_f_accept_tol
                )
            end
            println(io, "  Pre-filter: top_k=$top_k_str, dedup=$(ana_dedup_fraction)")
            println(io)

            ana_method = if ana_refinement_goal == :minimum
                GlobtimPostProcessing.OptimNelderMead(;
                    gradient_method = ana_grad_method,
                    hessian_tol = ana_hessian_tol,
                    max_iterations = ana_newton_max_iter,
                    max_time = ana_max_time_pt,
                )
            else
                GlobtimPostProcessing.NewtonCP(;
                    gradient_method = ana_grad_method,
                    tol = ana_newton_tol,
                    accept_tol = ana_accept_tol,
                    f_accept_tol = ana_f_accept_tol,
                    max_iterations = ana_newton_max_iter,
                    hessian_tol = ana_hessian_tol,
                )
            end
            discovery_result = GlobtimPostProcessing.build_known_cps_from_refinement(
                objective,
                highest_dr.critical_points,
                bounds;
                method = ana_method,
                dedup_fraction = ana_dedup_fraction,
                top_k = ana_top_k,
            )
            known_cps = discovery_result !== nothing ? discovery_result.known_cps : nothing
            # Note: build_known_cps_from_refinement already prints dedup summary with CP counts
        end
    end

    # ── 9. Assemble results ──
    return Dict{Symbol,Any}(
        :config => config,
        :experiment_result => exp_result,
        :degree_results => degree_results,
        :output_dir => output_dir,
        :entry => entry,
        :refinement_results => refinement_results,
        :known_cps => known_cps,
        :objective => objective,
        :bounds => bounds,
    )
end

end # module DynamicObjectivesGlobtimExt
