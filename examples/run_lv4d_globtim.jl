# Run with: julia --project=. DynamicObjectives/examples/run_lv4d_globtim.jl
#
# Constrained LV4D: Globtim CP recovery on screening candidates
#
# Loads top candidates from the LV4D catalogue (produced by run_lv4d_screening.jl)
# and runs Globtim CP recovery at radius 0.005, degrees 4/6/8, GN=12.
#
# Pipeline per candidate:
#   1. Build TolerantObjective at Tsit5/1e-4 (coarse for grid eval)
#   2. run_standard_experiment — polynomial fit + HomotopyContinuation CP solving
#   3. Switch to Vern9/1e-10 → compute_degree_capture_results
#   4. Print per-candidate per-degree recovery table

using Printf
using LinearAlgebra
using DynamicObjectives
using Globtim
using GlobtimPostProcessing

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

const DOMAIN_RADIUS = 0.005    # good recovery at this radius
const GN = 12
const DEGREE_RANGE = 4:2:8

# IC and tspan for Globtim objective evaluation
# Use the same as screening (catalogue entries store tspan=[0,10], ic=[0.8,1.2,0.8,1.2])
const GLOBTIM_NUMPOINTS = 21

# Repo root (scripts moved from root to DynamicObjectives/examples/)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

# Input/output paths
const CATALOGUE_PATH = joinpath(REPO_ROOT, "globtim_results", "lv4d_catalogue.jsonl")
const RESULTS_DIR = joinpath(REPO_ROOT, "globtim_results", "lv4d_globtim")

# ═══════════════════════════════════════════════════════════════════════════════
# Load catalogue
# ═══════════════════════════════════════════════════════════════════════════════

println("="^120)
println("  Constrained LV4D: Globtim CP Recovery")
println("="^120)
println()

entries = load_catalogue(CATALOGUE_PATH)
@printf("  Loaded %d candidates from %s\n", length(entries), CATALOGUE_PATH)
for (i, entry) in enumerate(entries)
    @printf(
        "    %d. %s: p_true = [%s]\n",
        i,
        entry.name,
        join([@sprintf("%.4f", x) for x in entry.p_true], ", ")
    )
end
println()

if isempty(entries)
    error("No catalogue entries found at $CATALOGUE_PATH")
end

model, params, states, outputs = define_constrained_lotka_volterra_4D()

# ═══════════════════════════════════════════════════════════════════════════════
# Run Globtim for each candidate
# ═══════════════════════════════════════════════════════════════════════════════

println("─"^120)
@printf(
    "  Globtim CP recovery: radius=%.4f, degrees=%s, GN=%d (%d^%d = %d grid points)\n",
    DOMAIN_RADIUS,
    string(collect(DEGREE_RANGE)),
    GN,
    GN,
    4,
    GN^4
)
println("─"^120)
println()

all_results = CandidateResult[]

for (ci, entry) in enumerate(entries)
    p_true = entry.p_true
    ic = entry.ic
    tspan = entry.time_interval

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

    # Build TolerantObjective — coarse for grid eval
    objective = TolerantObjective(
        model,
        outputs,
        ic,
        p_true,
        tspan,
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
        objective_name = @sprintf("lv4d_c%02d", ci),
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
# Summary tables
# ═══════════════════════════════════════════════════════════════════════════════

print_candidate_summary_table(
    all_results,
    entries;
    title = "LV4D GLOBTIM RECOVERY SUMMARY",
    radius = DOMAIN_RADIUS,
    degree_range = DEGREE_RANGE,
)
