# Run with: julia --project=. Dynamic_objectives/examples/run_daisy4d_multi_ptrue.jl
#
# DAISY Ex3 4D: multi-p_true validation via screening pipeline
#
# Tests whether Globtim CP recovery generalizes across the DAISY Ex3 parameter
# landscape, not just for one lucky p_true.
#
# Pipeline:
#   1. screen_and_probe(200, 10) — discover viable p_true candidates
#   2. screening_to_catalogue — pick top 5, save to JSONL
#   3. For each p_true: run Globtim at radius 0.005, degrees 4/6/8, GN=12
#      with TolerantObjective (Tsit5/1e-4 grid -> Vern9/1e-10 post-processing)
#   4. Print per-p_true recovery table + aggregate statistics

using Printf
using LinearAlgebra
using Dynamic_objectives
using Globtim
using GlobtimPostProcessing

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

const DOMAIN_RADIUS = 0.005    # good recovery at this radius
const GN = 12
const DEGREE_RANGE = 4:2:8
const N_CANDIDATES = 200
const N_PROBES = 10
const TOP_N = 5

# Screening uses wider bounds / longer tspan (proven in run_screening_demo.jl)
const SCREEN_IC = [1.0, 0.5, 0.3, 0.0]
const SCREEN_BOUNDS = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)]
const SCREEN_TSPAN = [0.0, 20.0]

# Globtim runs use tighter time interval
const GLOBTIM_TSPAN = [0.0, 10.0]
const GLOBTIM_IC = [1.0, 2.0, 1.0, 1.0]
const GLOBTIM_NUMPOINTS = 21

# Repo root (scripts moved from root to Dynamic_objectives/examples/)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

# Output paths
const CATALOGUE_PATH = joinpath(REPO_ROOT, "globtim_results", "daisy4d_catalogue.jsonl")
const RESULTS_DIR = joinpath(REPO_ROOT, "globtim_results", "daisy4d_multi_ptrue")

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Screen the parameter landscape
# ═══════════════════════════════════════════════════════════════════════════════

println("=" ^ 100)
println("  DAISY Ex3 4D: Multi-p_true Validation via Screening Pipeline")
println("  (Step 2 of Globtim 3D-4D validation)")
println("=" ^ 100)
println()

model, params, states, outputs = define_daisy_ex3_model_4D()

println("─" ^ 100)
println("  Step 1: screen_and_probe ($N_CANDIDATES candidates, $N_PROBES probes)")
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

# Show top candidates
n_show = min(TOP_N + 2, length(ranked))  # show a couple extra for context
println("  Top $n_show candidates (by dynamic range):")
for i in 1:n_show
    idx = ranked[i]
    p = result.sweep.valid[idx]
    probe = result.probes[idx]
    marker = i <= TOP_N ? "*" : " "
    @printf(
        "  %s #%d: p_true = [%s]  dynamic_range=%.2f  variance=%.2e\n",
        marker,
        i,
        join([@sprintf("%.4f", x) for x in p], ", "),
        probe.dynamic_range,
        probe.variance,
    )
end
println("  (* = selected for Globtim validation)")
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Convert to catalogue entries + save
# ═══════════════════════════════════════════════════════════════════════════════

println("─" ^ 100)
println("  Step 2: screening_to_catalogue (top $TOP_N)")
println("─" ^ 100)

entries = screening_to_catalogue(
    result,
    define_daisy_ex3_model_4D,
    SCREEN_IC,
    SCREEN_BOUNDS,
    SCREEN_TSPAN;
    top_n = TOP_N,
    name_prefix = "DAISY_4D",
    description = "Auto-discovered via cc0 screening ($N_CANDIDATES candidates, $N_PROBES probes)",
)

mkpath(dirname(CATALOGUE_PATH))
save_catalogue(CATALOGUE_PATH, entries)
@printf("  Saved %d entries to %s\n", length(entries), CATALOGUE_PATH)
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Step 3: Run Globtim for each p_true
# ═══════════════════════════════════════════════════════════════════════════════

println("─" ^ 100)
println("  Step 3: Globtim CP recovery for each p_true")
@printf(
    "  Radius: %.4f, Degrees: %s, GN: %d (%d^%d = %d grid points)\n",
    DOMAIN_RADIUS,
    string(collect(DEGREE_RANGE)),
    GN,
    GN,
    4,
    GN^4
)
println("─" ^ 100)
println()

all_results = CandidateResult[]

for (ci, entry) in enumerate(entries)
    p_true = entry.p_true
    println("── Candidate $ci/$(length(entries)): $(entry.name)")
    @printf("   p_true = [%s]\n", join([@sprintf("%.4f", x) for x in p_true], ", "))

    # Build domain around this p_true
    bounds = build_bounds(p_true, DOMAIN_RADIUS)
    config = Globtim.ExperimentParams(
        GN = GN,
        degree_range = DEGREE_RANGE,
        domain_size = DOMAIN_RADIUS,
    )
    output_dir = joinpath(RESULTS_DIR, @sprintf("candidate_%02d", ci))

    known_cps = KnownCriticalPoints([p_true], [0.0], [:min], bounds)

    # Build TolerantObjective from this entry's p_true
    # Note: we use GLOBTIM_IC and GLOBTIM_TSPAN, not the screening IC/tspan
    objective = TolerantObjective(
        model,
        outputs,
        GLOBTIM_IC,
        p_true,
        GLOBTIM_TSPAN,
        GLOBTIM_NUMPOINTS,
        L2_norm;
        solver = Tsit5(),
        abstol = 1e-4,
        reltol = 1e-4,
    )

    # Run Globtim
    cand_start = time()
    exp_result = run_standard_experiment(
        objective_function = objective,
        objective_name = @sprintf("daisy4d_c%02d", ci),
        bounds = bounds,
        experiment_config = config,
        output_dir = output_dir,
        metadata = Dict{String,Any}("candidate" => ci, "p_true" => p_true),
        true_params = p_true,
    )

    degree_results = exp_result[:degree_results]

    # Switch to tight tolerances for capture analysis
    set_tolerance!(objective, 1e-10)
    set_solver!(objective, Vern9())

    degree_capture_results = compute_degree_capture_results(degree_results, known_cps)
    capture_by_deg = Dict{Int,CaptureResult}()
    for (deg, cr) in degree_capture_results
        capture_by_deg[deg] = cr
    end

    cand_time = time() - cand_start

    # Extract per-degree results
    for dr in degree_results
        raw_obj = dr.best_objective !== nothing ? dr.best_objective : NaN
        raw_rec = dr.best_estimate !== nothing ? norm(dr.best_estimate - p_true) : NaN

        cap_1pct = NaN
        cap_5pct = NaN
        if haskey(capture_by_deg, dr.degree)
            cr = capture_by_deg[dr.degree]
            cap_1pct = capture_rate_at(cr, 0.01)
            cap_5pct = capture_rate_at(cr, 0.05)
        end

        push!(
            all_results,
            CandidateResult(
                entry.name,
                p_true,
                dr.degree,
                dr.n_critical_points,
                raw_obj,
                raw_rec,
                dr.l2_approx_error,
                dr.relative_l2_error,
                cap_1pct,
                cap_5pct,
            ),
        )
    end

    @printf("   Wall time: %.1fs\n", cand_time)
    for dr in degree_results
        rec = dr.best_estimate !== nothing ? norm(dr.best_estimate - p_true) : NaN
        rec_str = isnan(rec) ? "N/A" : @sprintf("%.3e", rec)
        @printf(
            "     deg %d: %d CPs, raw recovery = %s, L2 = %.3e\n",
            dr.degree,
            dr.n_critical_points,
            rec_str,
            dr.l2_approx_error
        )
    end
    println()
end

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Summary tables
# ═══════════════════════════════════════════════════════════════════════════════

print_candidate_summary_table(
    all_results,
    entries;
    title = "MULTI-P_TRUE SUMMARY: DAISY Ex3 4D",
    radius = DOMAIN_RADIUS,
    degree_range = DEGREE_RANGE,
    name_width = 25,
    p_true_precision = 3,
)
