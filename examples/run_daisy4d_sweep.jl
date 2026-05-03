# Run with: julia --project=. DynamicObjectives/examples/run_daisy4d_sweep.jl
#
# DAISY Ex3 4D: domain radius x degree sweep with capture analysis
#
# Sweeps domain radii [0.01, 0.005, 0.001] x degrees [4, 6, 8] around the
# known minimum at p_true = [0.2, 0.3, 0.5, 0.6].
# For each (radius, degree) combination, runs the full 3-step pipeline:
#   Step 1: Setup — TolerantObjective with coarse Tsit5/1e-4
#   Step 2: Run Globtim — polynomial fit + HC critical point solving
#   Step 3: Post-process — switch to Vern9/1e-10, gradient validation,
#           Nelder-Mead refinement, capture analysis
#
# Output: summary table of radius x degree -> {n_CPs, best_raw, best_refined,
#         recovery_error, L2_error, capture_at_1%, capture_at_5%, wall_time}
#
# Expected runtime: ~8-12 minutes (3 radii x 3 degrees x ~40s grid eval each)

using Printf
using LinearAlgebra
using DynamicObjectives
using Globtim
using GlobtimPostProcessing

# ═══════════════════════════════════════════════════════════════════════════════
# Problem constants
# ═══════════════════════════════════════════════════════════════════════════════

model, params, states, outputs = define_daisy_ex3_model_4D()

const P_TRUE = [0.2, 0.3, 0.5, 0.6]
const IC = [1.0, 2.0, 1.0, 1.0]
const TIME_INTERVAL = [0.0, 10.0]
const NUM_POINTS = 21
const DIM = length(P_TRUE)
const GN = 12
const DEGREE_RANGE = 4:2:8

# Radii to sweep (large to small)
const RADII = [0.01, 0.005, 0.001]

# Repo root (scripts moved from root to DynamicObjectives/examples/)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

# ═══════════════════════════════════════════════════════════════════════════════
# Result storage
# ═══════════════════════════════════════════════════════════════════════════════

struct SweepRow
    radius::Float64
    degree::Int
    l2_error::Float64
    rel_l2::Float64
    n_cps::Int
    best_raw_obj::Float64       # best raw objective value (NaN if no CPs)
    raw_recovery_error::Float64 # ||best_raw_estimate - p_true|| (NaN if no CPs)
    best_refined_obj::Float64   # best refined objective value (NaN if no refinement)
    recovery_error::Float64     # ||best_refined - p_true|| (NaN if no estimate)
    capture_1pct::Float64       # capture rate at 1% tolerance (NaN if no capture data)
    capture_5pct::Float64       # capture rate at 5% tolerance (NaN if no capture data)
    wall_time::Float64          # total Globtim time for this degree (poly + HC + processing + IO)
    # Per-stage timing breakdown
    poly_time::Float64          # polynomial construction time
    solve_time::Float64         # HC critical point solving time
    process_time::Float64       # critical point processing time
    refinement_time::Float64    # Nelder-Mead refinement time
    capture_time::Float64       # capture analysis time
end

all_rows = SweepRow[]

# ═══════════════════════════════════════════════════════════════════════════════
# Banner
# ═══════════════════════════════════════════════════════════════════════════════

println("="^90)
println("  DAISY Ex3 4D: Domain Radius x Degree Sweep with Capture Analysis")
println("  (Step 1 of Globtim 3D-4D validation)")
println("="^90)
println()
@printf("  %-16s %s\n", "Model:", "DAISY Ex3 4D (4 params: p1, p3, p4, p6)")
@printf("  %-16s %s\n", "P_TRUE:", string(P_TRUE))
@printf("  %-16s %s\n", "IC:", string(IC))
@printf(
    "  %-16s [%.1f, %.1f], %d points\n",
    "Time:",
    TIME_INTERVAL[1],
    TIME_INTERVAL[2],
    NUM_POINTS
)
@printf("  %-16s %d^%d = %d points per run\n", "Grid:", GN, DIM, GN^DIM)
@printf("  %-16s %s\n", "Degrees:", string(collect(DEGREE_RANGE)))
@printf("  %-16s %s\n", "Radii:", string(RADII))
@printf("  %-16s Tsit5/1e-4 (grid) -> Vern9/1e-10 (post-processing)\n", "ODE solver:")
@printf("  %-16s FiniteDiff (numerical)\n", "Gradient:")
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Sweep: for each radius, run full pipeline
# ═══════════════════════════════════════════════════════════════════════════════

total_start = time()

for (ri, radius) in enumerate(RADII)
    bounds = [(p - radius, p + radius) for p in P_TRUE]
    config =
        Globtim.ExperimentParams(GN = GN, degree_range = DEGREE_RANGE, domain_size = radius)
    output_dir =
        joinpath(REPO_ROOT, "globtim_results", "daisy4d_sweep", @sprintf("r_%.4f", radius))

    known_cps = KnownCriticalPoints([P_TRUE], [0.0], [:min], bounds)

    # Refinement config: slightly inset bounds to avoid boundary issues
    const_eps = radius * 1e-6
    refinement_config = ode_refinement_config(
        max_time_per_point = 5.0,
        bounds = [(b[1] + const_eps, b[2] - const_eps) for b in bounds],
        show_progress = false,
    )

    # Step 1: Build TolerantObjective with coarse tolerances
    objective = TolerantObjective(
        model,
        outputs,
        IC,
        P_TRUE,
        TIME_INTERVAL,
        NUM_POINTS,
        L2_norm;
        solver = Tsit5(),
        abstol = 1e-4,
        reltol = 1e-4,
    )

    println("─"^90)
    @printf(
        "  Radius %.4f (domain width %.4f, %d/%d)\n",
        radius,
        2 * radius,
        ri,
        length(RADII)
    )
    println("─"^90)

    # Step 2: Run Globtim (at coarse tolerances)
    radius_start = time()
    result = run_standard_experiment(
        objective_function = objective,
        objective_name = @sprintf("daisy4d_r%.4f", radius),
        bounds = bounds,
        experiment_config = config,
        output_dir = output_dir,
        metadata = Dict{String,Any}("sweep_radius" => radius, "ic" => IC),
        true_params = P_TRUE,
    )

    degree_results = result[:degree_results]

    # Step 3: Switch to tight tolerances for post-processing
    set_tolerance!(objective, 1e-10)
    set_solver!(objective, Vern9())

    # Per-degree: refinement + capture analysis
    refinement_start = time()
    degree_analyses = run_degree_analyses(
        degree_results,
        objective,
        output_dir,
        refinement_config;
        gradient_method = :finitediff,
    )
    refinement_wall = time() - refinement_start

    capture_start = time()
    degree_capture_results = compute_degree_capture_results(degree_results, known_cps)
    capture_wall = time() - capture_start

    radius_wall = time() - radius_start

    # Build lookup dicts from vector-of-tuple results
    analyses_by_deg = Dict{Int,Tuple}()
    for da in degree_analyses
        analyses_by_deg[da[1]] = da
    end
    capture_by_deg = Dict{Int,CaptureResult}()
    for (deg, cr) in degree_capture_results
        capture_by_deg[deg] = cr
    end

    # Distribute capture wall time evenly across degrees (capture is fast, no per-degree timing)
    n_success =
        count(dr -> dr.status == "success" && dr.n_critical_points > 0, degree_results)
    capture_per_deg = n_success > 0 ? capture_wall / n_success : 0.0

    # Extract per-degree rows
    for dr in degree_results
        raw_obj = dr.best_objective !== nothing ? dr.best_objective : NaN
        raw_rec_err = dr.best_estimate !== nothing ? norm(dr.best_estimate - P_TRUE) : NaN

        # Find refined objective for this degree
        refined_obj = NaN
        recovery_err = NaN
        ref_time = 0.0
        if haskey(analyses_by_deg, dr.degree)
            ref = analyses_by_deg[dr.degree][3]  # RefinedExperimentResult
            ref_time = ref.total_time
            if ref.n_converged > 0
                refined_obj = ref.best_refined_value
                best_pt = ref.refined_points[ref.best_refined_idx]
                recovery_err = norm(best_pt - P_TRUE)
            end
        end

        # Capture rates at 1% and 5% tolerance fraction
        cap_1pct = NaN
        cap_5pct = NaN
        cap_time = 0.0
        if haskey(capture_by_deg, dr.degree)
            cr = capture_by_deg[dr.degree]
            cap_1pct = capture_rate_at(cr, 0.01)
            cap_5pct = capture_rate_at(cr, 0.05)
            cap_time = capture_per_deg
        end

        push!(
            all_rows,
            SweepRow(
                radius,
                dr.degree,
                dr.l2_approx_error,
                dr.relative_l2_error,
                dr.n_critical_points,
                raw_obj,
                raw_rec_err,
                refined_obj,
                recovery_err,
                cap_1pct,
                cap_5pct,
                dr.total_computation_time,
                dr.polynomial_construction_time,
                dr.critical_point_solving_time,
                dr.critical_point_processing_time,
                ref_time,
                cap_time,
            ),
        )
    end

    # Print per-radius summary
    @printf("  Wall time: %.1fs\n", radius_wall)
    for dr in degree_results
        rec_str =
            dr.recovery_error !== nothing ? @sprintf("%.3e", dr.recovery_error) : "N/A"
        @printf(
            "    deg %d: %d CPs, raw recovery = %s, L2 = %.3e\n",
            dr.degree,
            dr.n_critical_points,
            rec_str,
            dr.l2_approx_error
        )
    end
    println()
end

total_wall = time() - total_start

# ═══════════════════════════════════════════════════════════════════════════════
# Full summary table
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("="^120)
println("  SWEEP SUMMARY: DAISY Ex3 4D — Radius x Degree")
@printf("  Total wall time: %.1fs\n", total_wall)
println("="^120)
println()

# Table header
@printf(
    "%-8s  %4s  %10s  %10s  %5s  %11s  %12s  %11s  %12s  %8s  %8s  %7s\n",
    "Radius",
    "Deg",
    "L2 Error",
    "Rel L2",
    "#CPs",
    "Raw f(x)",
    "Raw Recov",
    "Ref f(x)",
    "Ref Recov",
    "Cap@1%",
    "Cap@5%",
    "Time(s)"
)
println("-"^130)

for row in all_rows
    rel_str = isnan(row.rel_l2) ? "N/A" : @sprintf("%.3e", row.rel_l2)
    raw_str = isnan(row.best_raw_obj) ? "N/A" : @sprintf("%.3e", row.best_raw_obj)
    raw_rc_str =
        isnan(row.raw_recovery_error) ? "N/A" : @sprintf("%.3e", row.raw_recovery_error)
    ref_str = isnan(row.best_refined_obj) ? "N/A" : @sprintf("%.3e", row.best_refined_obj)
    rec_str = isnan(row.recovery_error) ? "N/A" : @sprintf("%.3e", row.recovery_error)
    c1_str = isnan(row.capture_1pct) ? "---" : @sprintf("%.0f%%", 100 * row.capture_1pct)
    c5_str = isnan(row.capture_5pct) ? "---" : @sprintf("%.0f%%", 100 * row.capture_5pct)
    @printf(
        "%-8.4f  %4d  %10.3e  %10s  %5d  %11s  %12s  %11s  %12s  %8s  %8s  %7.1f\n",
        row.radius,
        row.degree,
        row.l2_error,
        rel_str,
        row.n_cps,
        raw_str,
        raw_rc_str,
        ref_str,
        rec_str,
        c1_str,
        c5_str,
        row.wall_time
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Condensed view: best per radius
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("─"^130)
println("  Best per radius (lowest raw recovery error across degrees):")
println("─"^130)
@printf(
    "%-8s  %4s  %10s  %10s  %5s  %11s  %12s  %12s  %8s  %8s\n",
    "Radius",
    "Deg",
    "L2 Error",
    "Rel L2",
    "#CPs",
    "Raw f(x)",
    "Raw Recov",
    "Ref Recov",
    "Cap@1%",
    "Cap@5%"
)
println("-"^105)

for radius in RADII
    rows = filter(r -> r.radius == radius && !isnan(r.raw_recovery_error), all_rows)
    if isempty(rows)
        @printf(
            "%-8.4f  %4s  %10s  %10s  %5s  %11s  %12s  %12s  %8s  %8s\n",
            radius,
            "---",
            "---",
            "---",
            "---",
            "---",
            "---",
            "---",
            "---",
            "---"
        )
        continue
    end
    best = rows[argmin([r.raw_recovery_error for r in rows])]
    rel_str = isnan(best.rel_l2) ? "N/A" : @sprintf("%.3e", best.rel_l2)
    raw_str = isnan(best.best_raw_obj) ? "N/A" : @sprintf("%.3e", best.best_raw_obj)
    raw_rc_str = @sprintf("%.3e", best.raw_recovery_error)
    ref_rc_str = isnan(best.recovery_error) ? "N/A" : @sprintf("%.3e", best.recovery_error)
    c1_str = isnan(best.capture_1pct) ? "---" : @sprintf("%.0f%%", 100 * best.capture_1pct)
    c5_str = isnan(best.capture_5pct) ? "---" : @sprintf("%.0f%%", 100 * best.capture_5pct)
    @printf(
        "%-8.4f  %4d  %10.3e  %10s  %5d  %11s  %12s  %12s  %8s  %8s\n",
        best.radius,
        best.degree,
        best.l2_error,
        rel_str,
        best.n_cps,
        raw_str,
        raw_rc_str,
        ref_rc_str,
        c1_str,
        c5_str
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Per-stage timing breakdown
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("="^120)
println("  TIMING BREAKDOWN: Per-Stage Wall Time (seconds)")
println("="^120)
println()

@printf(
    "%-8s  %4s  %8s  %8s  %8s  %8s  %8s  %8s  %6s\n",
    "Radius",
    "Deg",
    "Poly(s)",
    "HC(s)",
    "Proc(s)",
    "Ref(s)",
    "Cap(s)",
    "Total(s)",
    "%Ref"
)
println("-"^90)

for row in all_rows
    total =
        row.poly_time +
        row.solve_time +
        row.process_time +
        row.refinement_time +
        row.capture_time
    pct_ref = total > 0 ? 100 * row.refinement_time / total : 0.0
    @printf(
        "%-8.4f  %4d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %5.0f%%\n",
        row.radius,
        row.degree,
        row.poly_time,
        row.solve_time,
        row.process_time,
        row.refinement_time,
        row.capture_time,
        total,
        pct_ref
    )
end

println()
