#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════════════════════
# run_3d_landscape_search.jl — Creative search for interesting 3D landscapes
#
# Screens all new 3D model candidates, scores them, and prints a comparison table.
# Three creative approaches to find 3D landscapes matching LV2D quality:
#
#   1. LV2D Reparam 3D: nonlinear polynomial map (θ₁,θ₂,θ₃) → (a,b) on LV2D v1
#   2. Rosenzweig-MacArthur 3D: Holling type II functional response (natural nonlinearity)
#   3. Goodwin 3D Product Obs: product observation y₁ = x₁·x₃ on Goodwin oscillator
#
# Pipeline per model:
#   1. screen_and_probe(500, 10) — discover viable p_true candidates
#   2. screening_to_catalogue — pick top 5, save to JSONL
#   3. score_landscape_grid — compute interestingness metrics
#   4. Print comparison table
#
# Usage:
#   julia --project=profiles/dev --threads=auto pkg/DynamicObjectives/examples/run_3d_landscape_search.jl
#
# Target: curvature > 0.15, plateau < 10%, local minima > 10
# Gold standard: LV2D_paper_1 — 12 minima, curvature 0.255, 1% plateau
# ═══════════════════════════════════════════════════════════════════════════════

using Printf
using DynamicObjectives

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))

# ── Model configurations ────────────────────────────────────────────────────

struct ModelConfig
    name::String
    define_fn::Function
    ic::Vector{Float64}
    bounds::Vector{Tuple{Float64,Float64}}
    tspan::Vector{Float64}
    solver::Any
    catalogue_path::String
    name_prefix::String
    n_candidates::Int
    n_probes::Int
    top_n::Int
end

const MODELS = [
    # Approach 1: LV2D reparameterized via polynomial map
    # a = θ₁ + θ₃², b = θ₂·θ₃ — centered on known-good (a*,b*) ≈ (0.2, 0.4)
    ModelConfig(
        "LV2D Reparam 3D (polynomial)",
        define_lv2d_reparam_3d_model,
        [0.3, 0.6],                                                  # LV2D IC
        [(-0.3, 0.7), (-0.5, 1.5), (-0.8, 0.8)],                  # θ₁, θ₂, θ₃
        [0.0, 1.0],                                                  # LV2D_paper_1 time span
        Tsit5(),
        joinpath(REPO_ROOT, "globtim_results", "lv2d_reparam_3d_catalogue.jsonl"),
        "LV2D_Reparam_3D",
        500,
        10,
        5,
    ),

    # Approach 2: Rosenzweig-MacArthur predator-prey
    # Holling type II: a·x/(1+a·h·x) creates natural nonlinearity
    ModelConfig(
        "Rosenzweig-MacArthur 3D",
        define_rosenzweig_macarthur_3d_model,
        [0.5, 0.3],                                                  # prey=0.5, predator=0.3
        [(0.5, 3.0), (0.5, 4.0), (0.5, 3.0)],                     # r, a, K
        [0.0, 10.0],                                                 # longer time for ecology
        Tsit5(),
        joinpath(REPO_ROOT, "globtim_results", "rosenzweig_macarthur_3d_catalogue.jsonl"),
        "RMA_3D",
        500,
        10,
        5,
    ),

    # Approach 3: Goodwin 3D with product observation y₁ = x₁·x₃
    ModelConfig(
        "Goodwin 3D Product Obs",
        define_goodwin_3d_product_obs_model,
        [0.3, 0.5, 0.6],                                            # Goodwin IC (from 3D screening)
        [(0.5, 5.0), (0.1, 2.0), (0.1, 2.0)],                     # k1, k2, k4
        [0.0, 10.0],                                                 # Goodwin time span
        Tsit5(),
        joinpath(REPO_ROOT, "globtim_results", "goodwin_3d_product_obs_catalogue.jsonl"),
        "Goodwin3D_ProdObs",
        500,
        10,
        5,
    ),
]

# ── Screen a single model ──────────────────────────────────────────────────

function screen_model(cfg::ModelConfig)
    println()
    println("═"^100)
    @printf("  %s\n", cfg.name)
    println("═"^100)

    model, params, states, outputs = cfg.define_fn()
    @printf("  IC: %s\n", string(cfg.ic))
    @printf("  Bounds: %s\n", string(cfg.bounds))
    @printf("  Tspan: %s, Solver: %s\n", string(cfg.tspan), string(typeof(cfg.solver)))
    @printf("  Screening: %d candidates, %d probes\n", cfg.n_candidates, cfg.n_probes)
    println()

    t0 = time()
    result = screen_and_probe(
        model,
        outputs,
        cfg.ic,
        cfg.bounds,
        cfg.tspan;
        n_candidates = cfg.n_candidates,
        n_probes = cfg.n_probes,
        solver = cfg.solver,
        verbose = true,
    )
    screen_time = time() - t0

    sweep = result.sweep
    ranked = result.ranking.ranked_indices

    @printf(
        "\n  Results: %d/%d passed (%.1f%%), %d ranked, %.1fs\n",
        sweep.diagnostics.n_valid,
        sweep.diagnostics.n_total,
        100.0 * sweep.pass_rate,
        length(ranked),
        screen_time
    )

    if sweep.diagnostics.n_valid == 0
        @printf("  ⚠ No valid candidates — all trajectories diverged\n")
        return nothing
    end

    # Show top candidates
    n_show = min(cfg.top_n + 3, length(ranked))
    println("\n  Top $n_show candidates:")
    @printf("  %-4s  %-40s  %14s  %s\n", "Rank", "p_true", "dynamic_range", "Sel")
    println("  " * "-"^70)
    for i in 1:n_show
        idx = ranked[i]
        p = result.sweep.valid[idx]
        probe = result.probes[idx]
        marker = i <= cfg.top_n ? " *" : "  "
        @printf(
            "  %-4d  [%s]  %14.4f  %s\n",
            i,
            join([@sprintf("%.4f", x) for x in p], ", "),
            probe.dynamic_range,
            marker
        )
    end

    # Save catalogue
    entries = screening_to_catalogue(
        result,
        cfg.define_fn,
        cfg.ic,
        cfg.bounds,
        cfg.tspan;
        top_n = cfg.top_n,
        name_prefix = cfg.name_prefix,
        description = "Auto-discovered via run_3d_landscape_search.jl ($(cfg.n_candidates) candidates). $(cfg.name).",
    )
    mkpath(dirname(cfg.catalogue_path))
    save_catalogue(cfg.catalogue_path, entries)
    @printf("\n  Saved %d entries → %s\n", length(entries), cfg.catalogue_path)

    return entries
end

# ── Score catalogue entries ─────────────────────────────────────────────────

struct ScoreResult
    name::String
    model_name::String
    n_local_minima::Int
    curvature_score::Float64
    plateau_fraction::Float64
    dynamic_range::Float64
    basin_fraction::Float64
    n_deceptive::Int
    interestingness::Float64
    eval_time_ms::Float64
end

function score_entries(catalogue_path::String, model_name::String)
    if !isfile(catalogue_path)
        @printf("  ⚠ Catalogue not found: %s\n", catalogue_path)
        return ScoreResult[]
    end

    entries = load_catalogue(catalogue_path)
    scores = ScoreResult[]

    for entry in entries
        @printf("  Scoring %s...", entry.name)
        try
            model, _, _, outputs = entry.model_fn()
            result = score_landscape_grid(
                model,
                outputs,
                entry.ic,
                entry.p_true,
                entry.bounds,
                entry.time_interval;
                numpoints = entry.numpoints,
                distance_function = entry.distance_function,
                aggregate_distances = entry.aggregate_distances,
            )
            score = interestingness_score(result)
            push!(
                scores,
                ScoreResult(
                    entry.name,
                    model_name,
                    result.n_local_minima,
                    result.curvature_score,
                    result.plateau_fraction,
                    result.dynamic_range,
                    result.basin_fraction,
                    result.n_deceptive,
                    score,
                    result.evaluation_time_ms,
                ),
            )
            @printf(
                " min=%d, curv=%.3f, plat=%.2f, score=%.3f (%.0fms)\n",
                result.n_local_minima,
                result.curvature_score,
                result.plateau_fraction,
                score,
                result.evaluation_time_ms
            )
        catch e
            @printf(" ERROR: %s\n", first(split(sprint(showerror, e), '\n')))
        end
    end
    return scores
end

# ── Comparison table ────────────────────────────────────────────────────────

function print_comparison(all_scores::Vector{ScoreResult})
    isempty(all_scores) && (println("No scores to compare."); return)

    sorted = sort(all_scores, by = s -> -s.interestingness)

    println()
    println("═"^120)
    println("  3D Landscape Search — Comparison Table")
    println("  Target: curvature > 0.15, plateau < 10%, minima > 10")
    println("  Reference: LV2D_paper_1 — 12 minima, curvature 0.255, 1% plateau")
    println("═"^120)
    @printf(
        "%-4s %-24s %-28s %7s %6s %6s %8s %7s %6s %9s\n",
        "Rank",
        "Entry",
        "Model",
        "Minima",
        "Curv",
        "Plat",
        "DynRange",
        "Basin",
        "Decep",
        "Score"
    )
    println("─"^120)

    for (i, s) in enumerate(sorted)
        marker =
            (
                s.curvature_score > 0.15 &&
                s.plateau_fraction < 0.10 &&
                s.n_local_minima > 10
            ) ? "★ " : "  "
        @printf(
            "%s%-2d %-24s %-28s %7d %6.3f %6.2f %8.2f %7.2f %6d %9.3f\n",
            marker,
            i,
            s.name,
            s.model_name,
            s.n_local_minima,
            s.curvature_score,
            s.plateau_fraction,
            s.dynamic_range,
            s.basin_fraction,
            s.n_deceptive,
            s.interestingness
        )
    end

    println("─"^120)
    n_target = count(
        s ->
            s.curvature_score > 0.15 && s.plateau_fraction < 0.10 && s.n_local_minima > 10,
        sorted,
    )
    @printf(
        "  %d/%d entries meet all three targets (marked with ★)\n",
        n_target,
        length(sorted)
    )
    println("═"^120)
end

# ── Main ────────────────────────────────────────────────────────────────────

function main()
    println("="^100)
    println("  Creative Search for Interesting 3D Landscapes")
    println("  3 approaches × 500 candidates × 5 entries = 15 landscapes to score")
    println("="^100)

    # Step 1: Screen all models
    all_entries = Dict{String,Any}()
    for cfg in MODELS
        entries = screen_model(cfg)
        all_entries[cfg.name] = entries
    end

    # Step 2: Score all catalogues
    println()
    println("═"^100)
    println("  Scoring all catalogues")
    println("═"^100)

    all_scores = ScoreResult[]
    for cfg in MODELS
        if all_entries[cfg.name] !== nothing
            println("\n  $(cfg.name):")
            scores = score_entries(cfg.catalogue_path, cfg.name)
            append!(all_scores, scores)
        end
    end

    # Step 3: Print comparison
    print_comparison(all_scores)

    # Summary
    println()
    println("  Next steps:")
    println(
        "    1. Run full grid scoring: julia --project=profiles/dev --threads=auto experiments/sandbox/run_grid_scoring.jl",
    )
    println(
        "    2. Visualize best candidates: julia --project=profiles/dev --threads=auto pkg/DynamicObjectives/examples/level_set_animation_3d.jl <catalogue>.jsonl",
    )
    println(
        "    3. Iterate on winners: adjust time span, IC, bounds, or try LV2D reparam variant B (trig)",
    )
    println()
end

main()
