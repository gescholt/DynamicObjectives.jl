# Run with: julia --project=profiles/dev --threads=auto pkg/DynamicObjectives/examples/run_coupled_lv2d_3d_screening.jl
#
# Coupled LV2D 3D: parameter landscape screening
#
# Two LV2D v1 subsystems coupled via shared parameter θ₃ + state feed-forward.
# Subsystem A: (a, b) = (θ₁, θ₃)
# Subsystem B: (a, b) = (0.5·θ₃, θ₂) + ε·(x₁, x₂) feed-forward
#
# Non-stiff LV dynamics — uses Tsit5().
#
# Pipeline:
#   1. screen_and_probe(200, 10) — discover viable p_true candidates
#   2. screening_to_catalogue — pick top 5, save to JSONL
#   3. Print screening summary + top candidates with dynamic_range scores
#
# Next step: copy catalogue to paper/catalogue/ for 3D explorer:
#   cp globtim_results/coupled_lv2d_3d_catalogue.jsonl pkg/DynamicObjectives/paper/catalogue/

using Printf
using DynamicObjectives

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

const N_CANDIDATES = 200
const N_PROBES = 10
const TOP_N = 5

# Parameters: θ₁ (a of A), θ₂ (b of B), θ₃ (b of A / determines a of B)
# Bounds match LV2D_paper_1: a∈(-0.3, 0.7), b∈(-0.1, 0.9)
const SCREEN_IC = [0.3, 0.6, 0.3, 0.6]                                # both subsystems at LV2D IC
const SCREEN_BOUNDS = [(-0.3, 0.7), (-0.1, 0.9), (-0.1, 0.9)]        # θ₁, θ₂, θ₃
const SCREEN_TSPAN = [0.0, 1.0]                                       # LV2D_paper_1 uses T=1; longer diverges

# Non-stiff LV dynamics
const SCREEN_SOLVER = Tsit5()

# Repo root (scripts live in pkg/DynamicObjectives/examples/)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))

# Output paths
const CATALOGUE_PATH =
    joinpath(REPO_ROOT, "globtim_results", "coupled_lv2d_3d_catalogue.jsonl")

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Screen the parameter landscape
# ═══════════════════════════════════════════════════════════════════════════════

println("="^100)
println("  Coupled LV2D 3D: Parameter Landscape Screening")
println("="^100)
println()

model, params, states, outputs = define_coupled_lv2d_3d_model()

println("─"^100)
println("  Step 1: screen_and_probe ($N_CANDIDATES candidates, $N_PROBES probes)")
@printf("  IC: %s  (both subsystems at LV2D IC)\n", string(SCREEN_IC))
@printf("  Bounds: %s  (θ₁, θ₂, θ₃)\n", string(SCREEN_BOUNDS))
@printf("  Tspan: %s\n", string(SCREEN_TSPAN))
@printf("  Solver: Tsit5()  (non-stiff LV)\n")
@printf("  Outputs: x₁ + x₃  (one state per subsystem)\n")
println("─"^100)
println()

screen_start = time()
result = screen_and_probe(
    model,
    outputs,
    SCREEN_IC,
    SCREEN_BOUNDS,
    SCREEN_TSPAN;
    n_candidates = N_CANDIDATES,
    n_probes = N_PROBES,
    solver = SCREEN_SOLVER,
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
    100.0 * sweep.pass_rate
)
@printf(
    "  Ranked: %d candidates, %d filtered out\n",
    length(ranked),
    length(result.ranking.filtered_out)
)
@printf("  Screening time: %.1fs\n", screen_time)
println()

if sweep.diagnostics.n_valid == 0
    error(
        "No valid candidates found — all $N_CANDIDATES trajectories diverged. " *
        "Try wider bounds or a different time span.",
    )
end

# Show ranked candidates
n_show = min(TOP_N + 5, length(ranked))  # show extras for context
println("  Top $n_show candidates (by dynamic range):")
@printf(
    "  %-4s  %-30s  %14s  %12s  %s\n",
    "Rank",
    "p_true (θ₁, θ₂, θ₃)",
    "dynamic_range",
    "variance",
    "Selected"
)
println("  " * "-"^80)
for i in 1:n_show
    idx = ranked[i]
    p = result.sweep.valid[idx]
    probe = result.probes[idx]
    marker = i <= TOP_N ? " *" : "  "
    @printf(
        "  %-4d  [%s]  %14.4f  %12.2e  %s\n",
        i,
        join([@sprintf("%6.4f", x) for x in p], ", "),
        probe.dynamic_range,
        probe.variance,
        marker,
    )
end
println("  (* = selected for catalogue)")
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Convert to catalogue entries + save
# ═══════════════════════════════════════════════════════════════════════════════

println("─"^100)
println("  Step 2: screening_to_catalogue (top $TOP_N)")
println("─"^100)

entries = screening_to_catalogue(
    result,
    define_coupled_lv2d_3d_model,
    SCREEN_IC,
    SCREEN_BOUNDS,
    SCREEN_TSPAN;
    top_n = TOP_N,
    name_prefix = "CoupledLV2D_3D",
    description = "Auto-discovered via screen_and_probe ($N_CANDIDATES candidates, $N_PROBES probes). " *
                  "Solver: Tsit5(). Coupled LV2D subsystems with shared θ₃ + state feed-forward (ε=0.1).",
)

mkpath(dirname(CATALOGUE_PATH))
save_catalogue(CATALOGUE_PATH, entries)
@printf("  Saved %d entries to %s\n", length(entries), CATALOGUE_PATH)
println()

# Print catalogue entries for reference
println("  Catalogue entries:")
for (i, entry) in enumerate(entries)
    @printf(
        "    %d. %s: p_true = [%s]  (θ₁=%.4f, θ₂=%.4f, θ₃=%.4f)\n",
        i,
        entry.name,
        join([@sprintf("%.4f", x) for x in entry.p_true], ", "),
        entry.p_true[1],
        entry.p_true[2],
        entry.p_true[3]
    )
end
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════════

println("="^100)
println("  SCREENING COMPLETE")
println("="^100)
@printf("  Model:       Coupled LV2D 3D (θ₁, θ₂, θ₃ — x₁+x₃ observed)\n")
@printf("  Solver:      Tsit5() — non-stiff\n")
@printf(
    "  Pass rate:   %.1f%% (%d/%d)\n",
    100.0 * sweep.pass_rate,
    sweep.diagnostics.n_valid,
    sweep.diagnostics.n_total
)
@printf("  Ranked:      %d candidates\n", length(ranked))
@printf("  Saved:       %d catalogue entries → %s\n", length(entries), CATALOGUE_PATH)
@printf("  Wall time:   %.1fs\n", screen_time)
println()
println("  Next: copy catalogue for 3D explorer:")
println(
    "    cp globtim_results/coupled_lv2d_3d_catalogue.jsonl pkg/DynamicObjectives/paper/catalogue/",
)
println()
