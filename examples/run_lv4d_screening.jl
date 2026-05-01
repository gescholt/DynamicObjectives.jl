# Run with: julia --project=. Dynamic_objectives/examples/run_lv4d_screening.jl
#
# Constrained LV4D: parameter landscape screening
#
# Screens the constrained 4-species Lotka-Volterra model for viable p_true
# candidates and saves the top results to a catalogue for downstream use.
#
# Pipeline:
#   1. screen_and_probe(200, 10) — discover viable p_true candidates
#   2. screening_to_catalogue — pick top 3, save to JSONL
#   3. Print screening summary + top candidates with dynamic_range scores

using Printf
using Dynamic_objectives

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

const N_CANDIDATES = 200
const N_PROBES = 10
const TOP_N = 3

# From test_configs.jl LV_4D_Constrained entry
const SCREEN_IC = [0.8, 1.2, 0.8, 1.2]
const SCREEN_BOUNDS = [(-0.1, 0.1), (-0.1, 0.1), (-0.1, 0.1), (-0.1, 0.1)]
const SCREEN_TSPAN = [0.0, 10.0]

# Repo root (scripts moved from root to Dynamic_objectives/examples/)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

# Output paths
const CATALOGUE_PATH = joinpath(REPO_ROOT, "globtim_results", "lv4d_catalogue.jsonl")

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Screen the parameter landscape
# ═══════════════════════════════════════════════════════════════════════════════

println("=" ^ 100)
println("  Constrained LV4D: Parameter Landscape Screening")
println("=" ^ 100)
println()

model, params, states, outputs = define_constrained_lotka_volterra_4D()

println("─" ^ 100)
println("  Step 1: screen_and_probe ($N_CANDIDATES candidates, $N_PROBES probes)")
@printf("  IC: %s\n", string(SCREEN_IC))
@printf("  Bounds: %s\n", string(SCREEN_BOUNDS))
@printf("  Tspan: %s\n", string(SCREEN_TSPAN))
println("─" ^ 100)
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
    verbose = true,
)
screen_time = time() - screen_start

sweep = result.sweep
ranked = result.ranking.ranked_indices

println()
println("─" ^ 100)
println("  Screening Results")
println("─" ^ 100)
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
    error("No valid candidates found — all $N_CANDIDATES trajectories diverged")
end

# Show all ranked candidates
n_show = min(TOP_N + 5, length(ranked))  # show extras for context
println("  Top $n_show candidates (by dynamic range):")
@printf(
    "  %-4s  %-40s  %14s  %12s  %s\n",
    "Rank",
    "p_true",
    "dynamic_range",
    "variance",
    "Selected"
)
println("  " * "-" ^ 90)
for i in 1:n_show
    idx = ranked[i]
    p = result.sweep.valid[idx]
    probe = result.probes[idx]
    marker = i <= TOP_N ? " *" : "  "
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

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Convert to catalogue entries + save
# ═══════════════════════════════════════════════════════════════════════════════

println("─" ^ 100)
println("  Step 2: screening_to_catalogue (top $TOP_N)")
println("─" ^ 100)

entries = screening_to_catalogue(
    result,
    define_constrained_lotka_volterra_4D,
    SCREEN_IC,
    SCREEN_BOUNDS,
    SCREEN_TSPAN;
    top_n = TOP_N,
    name_prefix = "LV4D",
    description = "Auto-discovered via screen_and_probe ($N_CANDIDATES candidates, $N_PROBES probes)",
)

mkpath(dirname(CATALOGUE_PATH))
save_catalogue(CATALOGUE_PATH, entries)
@printf("  Saved %d entries to %s\n", length(entries), CATALOGUE_PATH)
println()

# Print catalogue entries for reference
println("  Catalogue entries:")
for (i, entry) in enumerate(entries)
    @printf(
        "    %d. %s: p_true = [%s]\n",
        i,
        entry.name,
        join([@sprintf("%.4f", x) for x in entry.p_true], ", ")
    )
end
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════════

println("=" ^ 100)
println("  SCREENING COMPLETE")
println("=" ^ 100)
@printf("  Model:       Constrained LV4D (4 skew-symmetric eps parameters)\n")
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
