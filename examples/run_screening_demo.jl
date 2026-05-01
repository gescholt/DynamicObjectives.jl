#!/usr/bin/env julia
# End-to-end screening pipeline demo
# Usage: julia --project=Dynamic_objectives Dynamic_objectives/examples/run_screening_demo.jl
#
# Demonstrates the full p_true discovery workflow:
#   1. Define model + parameter space
#   2. screen_and_probe: sweep candidates, filter bounded, probe landscape, rank
#   3. screening_to_catalogue: convert top candidates to CatalogueEntry
#   4. Save/load catalogue via JSONL persistence
#   5. Verify: construct objective from loaded entry, sanity-check at p_true

using Dynamic_objectives
using Printf

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Define model and parameter space
# ═══════════════════════════════════════════════════════════════════════════════

display_section("Screening Pipeline Demo: DAISY 4D")

model, params, states, outputs = define_daisy_ex3_model_4D()
ic = [1.0, 0.5, 0.3, 0.0]
bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)]
tspan = [0.0, 20.0]

display_results(
    [
        "Model" => "DAISY Ex3 4D",
        "Parameters" => length(bounds),
        "Bounds" => join(["($(b[1]), $(b[2]))" for b in bounds], ", "),
        "Time span" => "$(tspan[1]) → $(tspan[2])",
        "IC" => join([@sprintf("%.1f", x) for x in ic], ", "),
    ],
    title = "Model Configuration",
)

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Run screening pipeline
# ═══════════════════════════════════════════════════════════════════════════════

display_section("Step 2: screen_and_probe")
println("  Generating candidates, filtering for bounded trajectories,")
println("  probing landscape structure, ranking by dynamic range...\n")

N_CANDIDATES = 200
N_PROBES = 10
TOP_N = 5

result = screen_and_probe(
    model,
    outputs,
    ic,
    bounds,
    tspan;
    n_candidates = N_CANDIDATES,
    n_probes = N_PROBES,
    verbose = true,
)

sweep = result.sweep
ranked = result.ranking.ranked_indices

display_results(
    [
        "Candidates generated" => sweep.diagnostics.n_total,
        "Passed boundedness" => sweep.diagnostics.n_valid,
        "Pass rate" => @sprintf("%.1f%%", 100.0 * sweep.pass_rate),
        "Rejected (solver)" => sweep.diagnostics.n_rejected_solver_failure,
        "Rejected (amplitude)" => sweep.diagnostics.n_rejected_amplitude,
        "Screening time" => @sprintf("%.1fs", sweep.diagnostics.elapsed_seconds),
        "Ranked candidates" => length(ranked),
        "Filtered out" => length(result.ranking.filtered_out),
    ],
    title = "Screening Results",
)

# Print top candidates
display_section("Top $TOP_N Candidates (by dynamic range)")
n_show = min(TOP_N, length(ranked))
for i in 1:n_show
    idx = ranked[i]
    p = result.sweep.valid[idx]
    probe = result.probes[idx]
    @printf(
        "  #%d: p_true = [%s]  dynamic_range=%.2f  variance=%.2e  finite=%d/%d\n",
        i,
        join([@sprintf("%.4f", x) for x in p], ", "),
        probe.dynamic_range,
        probe.variance,
        probe.n_finite,
        probe.n_finite + probe.n_inf,
    )
end
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Step 3: Convert to catalogue entries
# ═══════════════════════════════════════════════════════════════════════════════

display_section("Step 3: screening_to_catalogue")

entries = screening_to_catalogue(
    result,
    define_daisy_ex3_model_4D,
    ic,
    bounds,
    tspan;
    top_n = TOP_N,
    name_prefix = "DAISY_4D_screened",
    description = "Auto-discovered via screening demo ($(N_CANDIDATES) candidates)",
)

println("  Created $(length(entries)) CatalogueEntry objects:")
for e in entries
    @printf(
        "    %-25s  p_true=[%s]  dim=%d\n",
        e.name,
        join([@sprintf("%.4f", x) for x in e.p_true], ", "),
        dimension(e),
    )
end
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Save and reload via JSONL
# ═══════════════════════════════════════════════════════════════════════════════

display_section("Step 4: JSONL persistence round-trip")

catalogue_path = joinpath(tempdir(), "screening_demo_catalogue.jsonl")
save_catalogue(catalogue_path, entries)
println("  Saved $(length(entries)) entries to: $catalogue_path")

loaded = load_catalogue(catalogue_path)
println("  Loaded $(length(loaded)) entries back")

summary = catalogue_summary(catalogue_path)
display_results(
    [
        "Total entries" => summary.total,
        "Unique models" => length(summary.models),
        "Dimensions" =>
            join(["$(d)D: $(n)" for (d, n) in sort(collect(summary.dimensions))], ", "),
    ],
    title = "Catalogue Summary",
)

# ═══════════════════════════════════════════════════════════════════════════════
# Step 5: Verify — construct objective from loaded entry
# ═══════════════════════════════════════════════════════════════════════════════

display_section("Step 5: Verify loaded entries")
println("  Constructing objective from each loaded entry and evaluating at p_true.\n")

for (i, entry) in enumerate(loaded)
    objective = create_objective(entry)
    val_at_true = objective(entry.p_true)
    status = val_at_true < 1e-6 ? "PASS" : (val_at_true < 1e-3 ? "WARN" : "FAIL")
    @printf("  %s  %-25s  f(p_true) = %.6e\n", status, entry.name, val_at_true)
end

# Clean up temp file
rm(catalogue_path, force = true)

println("\nDone. Pipeline validated end-to-end.")
