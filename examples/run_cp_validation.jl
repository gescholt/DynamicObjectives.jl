# Run with: julia --project=. Dynamic_objectives/examples/run_cp_validation.jl
#
# Validate polynomial CPs as true objective stationary points
#
# For each HC-discovered critical point on the polynomial approximation, evaluate
# the gradient of the TRUE ODE objective (Vern9/1e-10) and classify:
#   - Is ||∇f(cp)|| small? (true stationary point vs polynomial artifact)
#   - Newton refinement → Hessian eigenvalue classification (min/saddle/max/degenerate)
#
# Runs on all catalogue entries: 3 LV4D + 5 DAISY Ex3 4D.
#
# Pipeline per candidate:
#   1. Build TolerantObjective at Tsit5/1e-4 (coarse for grid eval)
#   2. run_standard_experiment — polynomial fit + HomotopyContinuation CP solving
#   3. Switch to Vern9/1e-10
#   4. validate_critical_points at thresholds 1e-3, 1e-6, 1e-9
#   5. refine_to_critical_points (Newton on ∇f=0, Hessian classification)
#   6. Print per-degree validation summary

using Printf
using LinearAlgebra
using Dynamic_objectives
using Globtim
using GlobtimPostProcessing

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

const DOMAIN_RADIUS = 0.005
const GN = 12
const DEGREE_RANGE = 4:2:8
const NUMPOINTS = 21

const GRADIENT_THRESHOLDS = [1e-3, 1e-6, 1e-9]

# Newton refinement parameters
const NEWTON_TOL = 1e-8
const NEWTON_MAX_ITER = 100
const HESSIAN_TOL = 1e-6

# Repo root (scripts moved from root to Dynamic_objectives/examples/)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

# Catalogue paths
const LV4D_CATALOGUE = joinpath(REPO_ROOT, "globtim_results", "lv4d_catalogue.jsonl")
const DAISY_CATALOGUE = joinpath(REPO_ROOT, "globtim_results", "daisy4d_catalogue.jsonl")
const RESULTS_DIR = joinpath(REPO_ROOT, "globtim_results", "cp_validation")

# DAISY Ex3 uses different IC/tspan for Globtim than for screening
const DAISY_GLOBTIM_IC = [1.0, 2.0, 1.0, 1.0]
const DAISY_GLOBTIM_TSPAN = [0.0, 10.0]

# ═══════════════════════════════════════════════════════════════════════════════
# Result types
# ═══════════════════════════════════════════════════════════════════════════════

struct DegreeValidation
    system::String
    candidate::String
    degree::Int
    n_total_cps::Int
    # Gradient validation at multiple thresholds
    n_valid_1e3::Int
    n_valid_1e6::Int
    n_valid_1e9::Int
    grad_min::Float64
    grad_median::Float64
    grad_max::Float64
    # Newton refinement classification
    n_converged::Int
    n_min::Int
    n_max::Int
    n_saddle::Int
    n_degenerate::Int
    # Raw recovery (for context)
    raw_recovery::Float64
    l2_error::Float64
end

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: run validation for one candidate
# ═══════════════════════════════════════════════════════════════════════════════

function validate_candidate(;
    system::String,
    model,
    outputs,
    entry::CatalogueEntry,
    ic::Vector{Float64},
    tspan::Vector{Float64},
    ci::Int,
    total::Int,
)::Vector{DegreeValidation}
    p_true = entry.p_true
    println("── $system candidate $ci/$total: $(entry.name)")
    @printf("   p_true = [%s]\n", join([@sprintf("%.6f", x) for x in p_true], ", "))

    bounds = [(p - DOMAIN_RADIUS, p + DOMAIN_RADIUS) for p in p_true]
    config = Globtim.ExperimentParams(
        GN = GN,
        degree_range = DEGREE_RANGE,
        domain_size = DOMAIN_RADIUS,
    )
    output_dir = joinpath(RESULTS_DIR, lowercase(system), @sprintf("candidate_%02d", ci))

    # Build TolerantObjective — coarse for grid eval
    objective = TolerantObjective(
        model,
        outputs,
        ic,
        p_true,
        tspan,
        NUMPOINTS,
        L2_norm;
        solver = Tsit5(),
        abstol = 1e-4,
        reltol = 1e-4,
    )

    # Run Globtim
    cand_start = time()
    exp_result = run_standard_experiment(
        objective_function = objective,
        objective_name = @sprintf("%s_c%02d", lowercase(system), ci),
        bounds = bounds,
        experiment_config = config,
        output_dir = output_dir,
        metadata = Dict{String,Any}(
            "system" => system,
            "candidate" => ci,
            "p_true" => p_true,
        ),
        true_params = p_true,
    )

    degree_results = exp_result[:degree_results]
    globtim_time = time() - cand_start
    @printf("   Globtim: %.1fs\n", globtim_time)

    # Switch to tight tolerances for validation
    set_tolerance!(objective, 1e-10)
    set_solver!(objective, Vern9())

    results = DegreeValidation[]

    for dr in degree_results
        if dr.status != "success" || isempty(dr.critical_points)
            @printf(
                "     deg %d: %s (%d CPs) — skipped\n",
                dr.degree,
                dr.status,
                dr.n_critical_points
            )
            continue
        end

        cps = dr.critical_points
        n_cps = length(cps)
        @printf("     deg %d: %d CPs — validating gradients...\n", dr.degree, n_cps)

        # 1. Gradient validation at multiple thresholds
        #    compute_gradient_norms once, then threshold
        grad_start = time()
        norms = compute_gradient_norms(cps, objective; gradient_method = :finitediff)
        grad_time = time() - grad_start

        finite_norms = filter(isfinite, norms)
        n_valid_1e3 = count(n -> n < 1e-3, norms)
        n_valid_1e6 = count(n -> n < 1e-6, norms)
        n_valid_1e9 = count(n -> n < 1e-9, norms)

        grad_min = isempty(finite_norms) ? Inf : minimum(finite_norms)
        grad_med =
            isempty(finite_norms) ? Inf :
            sort(finite_norms)[div(length(finite_norms) + 1, 2)]
        grad_max = isempty(finite_norms) ? Inf : maximum(finite_norms)

        @printf(
            "       gradients: %.1fs, valid@1e-3=%d/%d, @1e-6=%d/%d, @1e-9=%d/%d\n",
            grad_time,
            n_valid_1e3,
            n_cps,
            n_valid_1e6,
            n_cps,
            n_valid_1e9,
            n_cps
        )
        @printf(
            "       ||∇f||: min=%.2e, median=%.2e, max=%.2e\n",
            grad_min,
            grad_med,
            grad_max
        )

        # 2. Newton refinement + Hessian classification
        newton_start = time()
        newton_results = refine_to_critical_points(
            objective,
            cps;
            gradient_method = :finitediff,
            tol = NEWTON_TOL,
            max_iterations = NEWTON_MAX_ITER,
            bounds = bounds,
            hessian_tol = HESSIAN_TOL,
        )
        newton_time = time() - newton_start

        n_converged = count(r -> r.converged, newton_results)
        n_min = count(r -> r.converged && r.cp_type == :min, newton_results)
        n_max = count(r -> r.converged && r.cp_type == :max, newton_results)
        n_saddle = count(r -> r.converged && r.cp_type == :saddle, newton_results)
        n_degenerate = count(r -> r.converged && r.cp_type == :degenerate, newton_results)

        @printf(
            "       Newton: %.1fs, converged=%d/%d [min=%d, max=%d, saddle=%d, degen=%d]\n",
            newton_time,
            n_converged,
            n_cps,
            n_min,
            n_max,
            n_saddle,
            n_degenerate
        )

        # Raw recovery for context
        raw_rec = dr.best_estimate !== nothing ? norm(dr.best_estimate - p_true) : NaN

        push!(
            results,
            DegreeValidation(
                system,
                entry.name,
                dr.degree,
                n_cps,
                n_valid_1e3,
                n_valid_1e6,
                n_valid_1e9,
                grad_min,
                grad_med,
                grad_max,
                n_converged,
                n_min,
                n_max,
                n_saddle,
                n_degenerate,
                raw_rec,
                dr.l2_approx_error,
            ),
        )
    end

    println()
    return results
end

# ═══════════════════════════════════════════════════════════════════════════════
# Load catalogues
# ═══════════════════════════════════════════════════════════════════════════════

println("=" ^ 120)
println("  CP Validation: Polynomial CPs as True Objective Stationary Points")
println("=" ^ 120)
println()

lv4d_entries = load_catalogue(LV4D_CATALOGUE)
daisy_entries = load_catalogue(DAISY_CATALOGUE)
@printf("  Loaded %d LV4D candidates from %s\n", length(lv4d_entries), LV4D_CATALOGUE)
@printf(
    "  Loaded %d DAISY Ex3 candidates from %s\n",
    length(daisy_entries),
    DAISY_CATALOGUE
)
println()

if isempty(lv4d_entries)
    error("No LV4D catalogue entries found at $LV4D_CATALOGUE")
end
if isempty(daisy_entries)
    error("No DAISY catalogue entries found at $DAISY_CATALOGUE")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Run LV4D validation
# ═══════════════════════════════════════════════════════════════════════════════

println("─" ^ 120)
println(
    "  LV4D: $(length(lv4d_entries)) candidates, radius=$DOMAIN_RADIUS, degrees=$(collect(DEGREE_RANGE)), GN=$GN",
)
println("─" ^ 120)
println()

lv4d_model, _, _, lv4d_outputs = define_constrained_lotka_volterra_4D()

all_validations = DegreeValidation[]

for (ci, entry) in enumerate(lv4d_entries)
    dvs = validate_candidate(
        system = "LV4D",
        model = lv4d_model,
        outputs = lv4d_outputs,
        entry = entry,
        ic = entry.ic,
        tspan = entry.time_interval,
        ci = ci,
        total = length(lv4d_entries),
    )
    append!(all_validations, dvs)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Run DAISY Ex3 validation
# ═══════════════════════════════════════════════════════════════════════════════

println("─" ^ 120)
println(
    "  DAISY Ex3: $(length(daisy_entries)) candidates, radius=$DOMAIN_RADIUS, degrees=$(collect(DEGREE_RANGE)), GN=$GN",
)
println("─" ^ 120)
println()

daisy_model, _, _, daisy_outputs = define_daisy_ex3_model_4D()

for (ci, entry) in enumerate(daisy_entries)
    dvs = validate_candidate(
        system = "DAISY_Ex3",
        model = daisy_model,
        outputs = daisy_outputs,
        entry = entry,
        ic = DAISY_GLOBTIM_IC,
        tspan = DAISY_GLOBTIM_TSPAN,
        ci = ci,
        total = length(daisy_entries),
    )
    append!(all_validations, dvs)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 1: Per-candidate, per-degree validation
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("=" ^ 140)
println("  CP VALIDATION SUMMARY")
@printf(
    "  Radius: %.4f, Degrees: %s, GN: %d\n",
    DOMAIN_RADIUS,
    string(collect(DEGREE_RANGE)),
    GN
)
@printf(
    "  Gradient thresholds: %s\n",
    join([@sprintf("%.0e", t) for t in GRADIENT_THRESHOLDS], ", ")
)
println("=" ^ 140)
println()

# Header
@printf(
    "%-10s  %-22s  %4s  %5s  %12s  %10s  %-21s  %-21s  %-21s  %8s  %6s  %6s  %6s  %6s\n",
    "System",
    "Candidate",
    "Deg",
    "#CPs",
    "Raw Recov",
    "L2 Error",
    "Valid@1e-3",
    "Valid@1e-6",
    "Valid@1e-9",
    "Newton",
    "min",
    "max",
    "sadl",
    "degn"
)
println("-" ^ 140)

for v in all_validations
    rec_str = isnan(v.raw_recovery) ? "N/A" : @sprintf("%.3e", v.raw_recovery)

    v3 = @sprintf(
        "%3d/%3d (%3.0f%%)",
        v.n_valid_1e3,
        v.n_total_cps,
        100.0 * v.n_valid_1e3 / v.n_total_cps
    )
    v6 = @sprintf(
        "%3d/%3d (%3.0f%%)",
        v.n_valid_1e6,
        v.n_total_cps,
        100.0 * v.n_valid_1e6 / v.n_total_cps
    )
    v9 = @sprintf(
        "%3d/%3d (%3.0f%%)",
        v.n_valid_1e9,
        v.n_total_cps,
        100.0 * v.n_valid_1e9 / v.n_total_cps
    )

    newton_str = @sprintf("%3d/%3d", v.n_converged, v.n_total_cps)

    @printf(
        "%-10s  %-22s  %4d  %5d  %12s  %10.3e  %-21s  %-21s  %-21s  %8s  %6d  %6d  %6d  %6d\n",
        v.system,
        v.candidate,
        v.degree,
        v.n_total_cps,
        rec_str,
        v.l2_error,
        v3,
        v6,
        v9,
        newton_str,
        v.n_min,
        v.n_max,
        v.n_saddle,
        v.n_degenerate
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 2: Gradient norm distribution per system/degree
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("─" ^ 120)
println("  Gradient Norm Distribution (||∇f|| on true ODE objective)")
println("─" ^ 120)
println()

@printf(
    "%-10s  %-22s  %4s  %5s  %12s  %12s  %12s\n",
    "System",
    "Candidate",
    "Deg",
    "#CPs",
    "||∇f|| min",
    "||∇f|| med",
    "||∇f|| max"
)
println("-" ^ 85)

for v in all_validations
    @printf(
        "%-10s  %-22s  %4d  %5d  %12.3e  %12.3e  %12.3e\n",
        v.system,
        v.candidate,
        v.degree,
        v.n_total_cps,
        v.grad_min,
        v.grad_median,
        v.grad_max
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 3: Aggregate by system
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("─" ^ 120)
println("  Aggregate by System (across all candidates and degrees)")
println("─" ^ 120)
println()

for system in ["LV4D", "DAISY_Ex3"]
    rows = filter(v -> v.system == system, all_validations)
    if isempty(rows)
        continue
    end

    total_cps = sum(v -> v.n_total_cps, rows)
    total_v3 = sum(v -> v.n_valid_1e3, rows)
    total_v6 = sum(v -> v.n_valid_1e6, rows)
    total_v9 = sum(v -> v.n_valid_1e9, rows)
    total_converged = sum(v -> v.n_converged, rows)
    total_min = sum(v -> v.n_min, rows)
    total_max = sum(v -> v.n_max, rows)
    total_saddle = sum(v -> v.n_saddle, rows)
    total_degenerate = sum(v -> v.n_degenerate, rows)

    @printf("  %s:\n", system)
    @printf(
        "    Total polynomial CPs:    %d (across %d candidate-degree runs)\n",
        total_cps,
        length(rows)
    )
    @printf(
        "    True stationary @1e-3:   %d/%d (%.1f%%)\n",
        total_v3,
        total_cps,
        100.0 * total_v3 / total_cps
    )
    @printf(
        "    True stationary @1e-6:   %d/%d (%.1f%%)\n",
        total_v6,
        total_cps,
        100.0 * total_v6 / total_cps
    )
    @printf(
        "    True stationary @1e-9:   %d/%d (%.1f%%)\n",
        total_v9,
        total_cps,
        100.0 * total_v9 / total_cps
    )
    @printf(
        "    Newton converged:        %d/%d (%.1f%%)\n",
        total_converged,
        total_cps,
        100.0 * total_converged / total_cps
    )
    @printf(
        "    Classification:          min=%d, max=%d, saddle=%d, degenerate=%d\n",
        total_min,
        total_max,
        total_saddle,
        total_degenerate
    )

    artifact_rate = 1.0 - total_v3 / total_cps
    @printf("    Artifact rate (@1e-3):   %.1f%%\n", 100.0 * artifact_rate)
    println()
end

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 4: Aggregate by degree (across both systems)
# ═══════════════════════════════════════════════════════════════════════════════

println("─" ^ 120)
println("  Aggregate by Degree (across both systems)")
println("─" ^ 120)
println()

@printf(
    "  %4s  %6s  %12s  %12s  %12s  %8s  %6s  %6s  %6s  %6s\n",
    "Deg",
    "#CPs",
    "Valid@1e-3",
    "Valid@1e-6",
    "Valid@1e-9",
    "Newton",
    "min",
    "max",
    "sadl",
    "degn"
)
println("  " * "-" ^ 100)

for deg in collect(DEGREE_RANGE)
    rows = filter(v -> v.degree == deg, all_validations)
    if isempty(rows)
        continue
    end

    total_cps = sum(v -> v.n_total_cps, rows)
    total_v3 = sum(v -> v.n_valid_1e3, rows)
    total_v6 = sum(v -> v.n_valid_1e6, rows)
    total_v9 = sum(v -> v.n_valid_1e9, rows)
    total_converged = sum(v -> v.n_converged, rows)
    total_min = sum(v -> v.n_min, rows)
    total_max = sum(v -> v.n_max, rows)
    total_saddle = sum(v -> v.n_saddle, rows)
    total_degenerate = sum(v -> v.n_degenerate, rows)

    v3_str = @sprintf("%4d/%4d %3.0f%%", total_v3, total_cps, 100.0 * total_v3 / total_cps)
    v6_str = @sprintf("%4d/%4d %3.0f%%", total_v6, total_cps, 100.0 * total_v6 / total_cps)
    v9_str = @sprintf("%4d/%4d %3.0f%%", total_v9, total_cps, 100.0 * total_v9 / total_cps)
    nw_str = @sprintf("%4d/%4d", total_converged, total_cps)

    @printf(
        "  %4d  %6d  %12s  %12s  %12s  %8s  %6d  %6d  %6d  %6d\n",
        deg,
        total_cps,
        v3_str,
        v6_str,
        v9_str,
        nw_str,
        total_min,
        total_max,
        total_saddle,
        total_degenerate
    )
end

println()
