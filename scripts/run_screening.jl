#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════════════════════
# run_screening.jl — TOML-driven entry point for screen_and_probe (bead 20p7)
#
# Replaces the per-model `run_*_screening.jl` example scripts with a single
# config-driven driver. The TOML's [screening] section carries every parameter
# previously hardcoded in those scripts (model factory name, IC, bounds, tspan,
# n_candidates, n_probes, top_n, solver, catalogue output path).
#
# Usage:
#   julia --project=profiles/dev pkg/DynamicObjectives/scripts/run_screening.jl CONFIG.toml
#
# The config file must contain a [screening] section with at least:
#   - model_factory   string key into DynamicObjectives.MODEL_REGISTRY
#   - ic              [x0_1, x0_2, ...]
#   - bounds          [[lo, hi], ...]
#   - time_interval   [t_start, t_end]
#   - catalogue_path  output JSONL path (relative paths resolve via globtim_results/)
#
# All other [screening] fields are optional; defaults come from
# experiments/sandbox/sandbox_shared.jl::SCREENING_*_DEFAULT.
#
# See examples/configs/screening_lv4d.toml for a worked example.
# ═══════════════════════════════════════════════════════════════════════════════

using Printf
using DynamicObjectives
using Globtim: load_experiment_config

# The solvers used by _resolve_solver below (Tsit5, AutoTsit5, Vern7, Vern9,
# Rosenbrock23, Rodas5, TRBDF2) all come from `using DynamicObjectives`, which
# imports them from the granular OrdinaryDiffEq* sub-packages and re-exports
# them. This file used to `using OrdinaryDiffEq` — the umbrella package, which
# is NOT in profiles/cluster, so the script could not run on the cluster at all.
# Do not reinstate that import; add any new solver to DynamicObjectives instead.

# Pull the screening_options helper + SCREENING_*_DEFAULT constants.
include(joinpath(@__DIR__, "..", "..", "..", "experiments", "sandbox", "sandbox_shared.jl"))

# ── Solver lookup ────────────────────────────────────────────────────────────
# Mirror of Globtim.KNOWN_SCREENING_SOLVERS — string keys → solver objects.
# Adding a new key here also requires adding it to KNOWN_SCREENING_SOLVERS.
function _resolve_solver(name::String)
    name == "Tsit5" && return Tsit5()
    name == "AutoTsit5_Rosenbrock23" && return AutoTsit5(Rosenbrock23())
    name == "Vern7" && return Vern7()
    name == "Vern9" && return Vern9()
    name == "Rosenbrock23" && return Rosenbrock23()
    name == "Rodas5" && return Rodas5()
    name == "TRBDF2" && return TRBDF2()
    error("Unknown screening solver '$name'. Update _resolve_solver in run_screening.jl.")
end

function _resolve_model_factory(name::String)
    if !haskey(DynamicObjectives.MODEL_REGISTRY, name)
        known = join(sort(collect(keys(DynamicObjectives.MODEL_REGISTRY))), ", ")
        error(
            "model_factory '$name' not registered in DynamicObjectives.MODEL_REGISTRY. Known factories: $known",
        )
    end
    return DynamicObjectives.MODEL_REGISTRY[name]
end

function main(config_path::String)
    isfile(config_path) || error("Config file not found: $config_path")
    config = load_experiment_config(config_path)
    opts = screening_options(config)

    # Required fields
    isnothing(opts.model_factory) &&
        error("[screening] model_factory is required (TOML field model_factory)")
    isnothing(opts.ic) && error("[screening] ic is required (TOML field ic)")
    isnothing(opts.bounds) && error("[screening] bounds is required (TOML field bounds)")
    isnothing(opts.time_interval) &&
        error("[screening] time_interval is required (TOML field time_interval)")
    isnothing(opts.catalogue_path) &&
        error("[screening] catalogue_path is required (TOML field catalogue_path)")

    model_fn = _resolve_model_factory(opts.model_factory)
    solver = _resolve_solver(opts.solver)

    println("="^100)
    @printf("  TOML-driven screening: %s\n", config.name)
    println("="^100)
    println()
    @printf("  Model factory: %s\n", opts.model_factory)
    @printf("  IC: %s\n", string(opts.ic))
    @printf("  Bounds: %s\n", string(opts.bounds))
    @printf("  Tspan: %s\n", string(opts.time_interval))
    @printf("  Solver: %s\n", opts.solver)
    @printf(
        "  N: %d candidates × %d probes (top %d → catalogue)\n",
        opts.n_candidates,
        opts.n_probes,
        opts.top_n,
    )
    println()

    model, _params, _states, outputs = model_fn()

    screen_start = time()
    result = screen_and_probe(
        model,
        outputs,
        opts.ic,
        opts.bounds,
        opts.time_interval;
        n_candidates = opts.n_candidates,
        n_probes = opts.n_probes,
        numpoints_screen = opts.numpoints_screen,
        numpoints_probe = opts.numpoints_probe,
        min_finite_fraction = opts.min_finite_fraction,
        max_noise_ratio = opts.max_noise_ratio,
        solver = solver,
        verbose = true,
    )
    screen_time = time() - screen_start

    sweep = result.sweep
    ranked = result.ranking.ranked_indices

    println()
    println("─"^100)
    println("  Screening Results")
    println("─"^100)
    @printf(
        "  Candidates: %d generated, %d passed (%.1f%% pass rate)\n",
        sweep.diagnostics.n_total,
        sweep.diagnostics.n_valid,
        100.0 * sweep.pass_rate,
    )
    @printf(
        "  Ranked: %d candidates, %d filtered out\n",
        length(ranked),
        length(result.ranking.filtered_out),
    )
    @printf("  Screening time: %.1fs\n", screen_time)
    println()

    if sweep.diagnostics.n_valid == 0
        error("No valid candidates found — try wider bounds or a different time span.")
    end

    n_show = min(opts.top_n + 5, length(ranked))
    println("  Top $n_show candidates (by dynamic range):")
    @printf(
        "  %-4s  %-40s  %14s  %12s  %s\n",
        "Rank",
        "p_true",
        "dynamic_range",
        "variance",
        "Selected"
    )
    println("  " * "-"^90)
    for i in 1:n_show
        idx = ranked[i]
        p = result.sweep.valid[idx]
        probe = result.probes[idx]
        marker = i <= opts.top_n ? " *" : "  "
        @printf(
            "  %-4d  [%s]  %14.4f  %12.2e  %s\n",
            i,
            join([@sprintf("%7.4f", x) for x in p], ", "),
            probe.dynamic_range,
            probe.variance,
            marker,
        )
    end
    println("  (* = selected for catalogue)")
    println()

    name_prefix = something(opts.name_prefix, config.name)
    description = something(
        opts.description,
        "Auto-discovered via screen_and_probe (n_candidates=$(opts.n_candidates), n_probes=$(opts.n_probes))",
    )

    entries = screening_to_catalogue(
        result,
        model_fn,
        opts.ic,
        opts.bounds,
        opts.time_interval;
        top_n = opts.top_n,
        name_prefix = name_prefix,
        description = description,
    )

    mkpath(dirname(opts.catalogue_path))
    save_catalogue(opts.catalogue_path, entries)
    @printf("  Saved %d entries to %s\n", length(entries), opts.catalogue_path)
    println()

    println("="^100)
    println("  SCREENING COMPLETE")
    println("="^100)
    @printf(
        "  Pass rate: %.1f%% (%d/%d)\n",
        100.0 * sweep.pass_rate,
        sweep.diagnostics.n_valid,
        sweep.diagnostics.n_total
    )
    @printf("  Saved: %d catalogue entries → %s\n", length(entries), opts.catalogue_path)
    @printf("  Wall time: %.1fs\n", screen_time)
    println()
end

if length(ARGS) == 0
    println(
        stderr,
        "Usage: julia --project=profiles/dev pkg/DynamicObjectives/scripts/run_screening.jl CONFIG.toml",
    )
    exit(1)
end

main(ARGS[1])
