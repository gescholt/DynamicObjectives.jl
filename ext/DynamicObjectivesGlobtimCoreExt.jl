"""
DynamicObjectivesGlobtimCoreExt — config → (objective, bounds) resolution

Narrow package extension triggered by `Globtim` ALONE (unlike
`DynamicObjectivesGlobtimExt`, which needs Globtim + GlobtimPostProcessing +
Optim and carries the full HC/refinement pipeline). It provides
`DynamicObjectives.build_experiment_objective`, the single canonical
model-resolution step shared by:

  - DynamicObjectivesGlobtimExt.run_experiment_from_config (full pipeline)
  - the per-axis audit driver (audit, no HC)
  - the per-axis counterfactual driver (cf, no HC)

The audit/cf drivers deliberately depend only on DynamicObjectives + Globtim;
keeping this builder in its own extension preserves that property (rmow).

Triggered by: using DynamicObjectives; using Globtim
"""
module DynamicObjectivesGlobtimCoreExt

using DynamicObjectives
using Globtim

import DynamicObjectives:
    build_experiment_objective, _resolve_solver, load_catalogue, build_bounds,
    TolerantObjective, DISTANCE_REGISTRY, Tsit5

"""
    build_experiment_objective(config::Globtim.ExperimentPipelineConfig; io=stdout)
        -> (; objective, bounds, obj_name, entry, p_true)

Resolve a parsed TOML experiment config into an objective function and domain
bounds — catalogue mode (ODE `TolerantObjective`) or analytical mode
(`Globtim.FUNCTION_REGISTRY` benchmark). This is the single source of truth
for the resolution semantics:

  - `p_true`:  config overrides catalogue entry
  - solver:    config `[solver]` overrides (Tsit5 / 1e-4 defaults)
  - sampling:  `sample_times` overrides `time_interval` + `numpoints`
  - distance:  `[model].distance_function_override` resolves via
               `DISTANCE_REGISTRY` (opt-in swap, e.g. L2_squared)
  - domain:    centered on `p_center` (default `p_true`), precedence
               `radius` → `radii` → `bounds` → `entry.bounds`

In analytical mode `entry` and `p_true` are `nothing`.

Progress/diagnostic lines go to `io`; pass `io = devnull` to silence.
"""
function DynamicObjectives.build_experiment_objective(
    config::Globtim.ExperimentPipelineConfig;
    io::IO = stdout,
)
    if config.catalogue_path !== nothing
        # ── Catalogue mode: load entry, build TolerantObjective ──
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
            numpoints =
                config.solver_numpoints !== nothing ? config.solver_numpoints :
                entry.numpoints
            time_interval =
                config.time_interval !== nothing ? config.time_interval :
                entry.time_interval
            uneven_sampling_times = Float64[]
        end

        # Honour [model].distance_function_override — opt-in swap of the
        # catalogue entry's distance_function (e.g., L2_squared in place of
        # L2_norm) without editing the catalogue file.
        distance_fn = if config.distance_function_override !== nothing
            resolved = get(DISTANCE_REGISTRY, config.distance_function_override, nothing)
            resolved === nothing && error(
                "distance_function_override = \"$(config.distance_function_override)\" " *
                "not in DISTANCE_REGISTRY (have: " *
                join(sort(collect(keys(DISTANCE_REGISTRY))), ", ") *
                ")",
            )
            println(
                io,
                "  distance_function_override: \"$(config.distance_function_override)\" " *
                "(catalogue entry was using its default)",
            )
            resolved
        else
            entry.distance_function
        end

        objective = TolerantObjective(
            model,
            outputs,
            entry.ic,
            p_true,
            time_interval,
            numpoints,
            distance_fn,
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

        bounds = if config.radius !== nothing
            build_bounds(p_center, config.radius)
        elseif config.radii !== nothing
            build_bounds(p_center, config.radii)
        elseif config.bounds !== nothing
            config.bounds
        else
            entry.bounds
        end

        return (;
            objective,
            bounds,
            obj_name = config.entry_name,
            entry,
            p_true,
        )
    else
        # ── Analytical mode: look up function from FUNCTION_REGISTRY ──
        func_name = config.analytical_function
        dim = config.dimension
        println(io, "  Analytical function: $func_name ($(dim)D)")

        bench = Globtim.get_benchmark_config_by_name(func_name, dim)
        bounds = config.bounds !== nothing ? config.bounds : bench.bounds

        println(io, "  Bounds: $(join(["[$(lb), $(ub)]" for (lb, ub) in bounds], " × "))")

        return (;
            objective = bench.objective,
            bounds,
            obj_name = bench.name,
            entry = nothing,
            p_true = nothing,
        )
    end
end

end # module
