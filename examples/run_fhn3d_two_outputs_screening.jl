# Run with: julia --project=profiles/dev --threads=auto pkg/Dynamic_objectives/examples/run_fhn3d_two_outputs_screening.jl
#
# FitzHugh-Nagumo 3D dual-output: parameter landscape screening
#
# Same FHN dynamics as run_fhn3d_screening.jl, but measuring BOTH voltage V
# and recovery variable R. The dual observation eliminates the structural
# non-identifiability that collapses single-output level sets to hyperplanes.
#
# FHN is a stiff system — uses AutoTsit5(Rosenbrock23()) instead of plain Tsit5.
#
# Pipeline:
#   1. screen_and_probe(200, 10) — discover viable p_true candidates
#   2. screening_to_catalogue — pick top 5, save to JSONL
#   3. Print screening summary + top candidates with dynamic_range scores
#
# Next step: copy catalogue to paper/catalogue/ for 3D explorer:
#   cp globtim_results/fhn3d_two_outputs_catalogue.jsonl pkg/Dynamic_objectives/paper/catalogue/

using Printf
using Dynamic_objectives

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

const N_CANDIDATES = 200
const N_PROBES = 10
const TOP_N = 5

# FHN parameters: g (coupling), a (threshold), b (recovery)
# g > 0.1 required (appears in denominator 1/g)
const SCREEN_IC = [0.0, 0.0]                                        # resting state (V, R)
const SCREEN_BOUNDS = [(0.1, 2.0), (0.0, 2.0), (0.0, 2.0)]         # g, a, b
const SCREEN_TSPAN = [0.0, 20.0]                                    # moderate; 50.0 causes stiffness failures

# Stiff-aware solver — plain Tsit5 fails on FHN
const SCREEN_SOLVER = AutoTsit5(Rosenbrock23())

# Repo root (scripts live in pkg/Dynamic_objectives/examples/)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))

# Output paths
const CATALOGUE_PATH =
    joinpath(REPO_ROOT, "globtim_results", "fhn3d_two_outputs_catalogue.jsonl")

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Screen the parameter landscape
# ═══════════════════════════════════════════════════════════════════════════════

println("=" ^ 100)
println("  FitzHugh-Nagumo 3D Dual-Output (V+R): Parameter Landscape Screening")
println("=" ^ 100)
println()

model, params, states, outputs = define_fitzhugh_nagumo_3D_model_two_outputs()

println("─" ^ 100)
println("  Step 1: screen_and_probe ($N_CANDIDATES candidates, $N_PROBES probes)")
@printf("  IC: %s  (V=0, R=0 resting state)\n", string(SCREEN_IC))
@printf("  Bounds: %s  (g, a, b)\n", string(SCREEN_BOUNDS))
@printf("  Tspan: %s\n", string(SCREEN_TSPAN))
@printf("  Solver: AutoTsit5(Rosenbrock23())  (stiff-aware)\n")
@printf("  Outputs: V + R  (dual-output — both states observed)\n")
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
    solver = SCREEN_SOLVER,
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
    "p_true (g, a, b)",
    "dynamic_range",
    "variance",
    "Selected"
)
println("  " * "-" ^ 80)
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

println("─" ^ 100)
println("  Step 2: screening_to_catalogue (top $TOP_N)")
println("─" ^ 100)

entries = screening_to_catalogue(
    result,
    define_fitzhugh_nagumo_3D_model_two_outputs,
    SCREEN_IC,
    SCREEN_BOUNDS,
    SCREEN_TSPAN;
    top_n = TOP_N,
    name_prefix = "FHN3D_2out",
    description = "Auto-discovered via screen_and_probe ($N_CANDIDATES candidates, $N_PROBES probes). " *
                  "Solver: AutoTsit5(Rosenbrock23()). Dual-output (V+R) neuronal excitability model.",
)

mkpath(dirname(CATALOGUE_PATH))
save_catalogue(CATALOGUE_PATH, entries)
@printf("  Saved %d entries to %s\n", length(entries), CATALOGUE_PATH)
println()

# Print catalogue entries for reference
println("  Catalogue entries:")
for (i, entry) in enumerate(entries)
    @printf(
        "    %d. %s: p_true = [%s]  (g=%.4f, a=%.4f, b=%.4f)\n",
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

println("=" ^ 100)
println("  SCREENING COMPLETE")
println("=" ^ 100)
@printf("  Model:       FitzHugh-Nagumo 3D dual-output (g, a, b — V+R observed)\n")
@printf("  Solver:      AutoTsit5(Rosenbrock23()) — stiff-aware\n")
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
    "    cp globtim_results/fhn3d_two_outputs_catalogue.jsonl pkg/Dynamic_objectives/paper/catalogue/",
)
println()
