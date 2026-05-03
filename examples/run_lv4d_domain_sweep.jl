# Run with: julia --project=. DynamicObjectives/examples/run_lv4d_domain_sweep.jl
#
# LV4D Domain × Degree Sweep
#
# Answers: "How large an LV4D domain can we reliably solve on the laptop at
# degree ≤ 10 in ≤ 2 min/combo?"  Also detects multiple critical points at
# larger radii.
#
# Sweeps radii [0.005, 0.01, 0.02, 0.05, 0.1] × degrees [4, 6, 8, 10]
# for one catalogue LV4D candidate (LV4D_1) at GN=12.
#
# For each (radius, degree) combination:
#   1. Setup — TolerantObjective with coarse Tsit5/1e-4
#   2. Run Globtim — polynomial fit + HC critical point solving
#   3. Post-process — Vern9/1e-10, refinement, capture analysis
#   4. Multi-CP analysis — count in-domain CPs, cluster to find distinct minima
#
# Output: summary tables for accuracy, timing, and multi-CP discovery.
#
# Expected runtime: ~30-50 min total (20 combos, budget ~2 min each;
#   small-radius combos finish in seconds, large-radius deg-10 may hit ~2 min)

using Printf
using LinearAlgebra
using DynamicObjectives
using Globtim
using GlobtimPostProcessing

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

# Load first catalogue candidate
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const CATALOGUE_PATH = joinpath(REPO_ROOT, "globtim_results", "lv4d_catalogue.jsonl")

entries = load_catalogue(CATALOGUE_PATH)
isempty(entries) && error("No catalogue entries at $CATALOGUE_PATH")
const ENTRY = entries[1]  # LV4D_1

model, params, states, outputs = define_constrained_lotka_volterra_4D()

const P_TRUE = ENTRY.p_true
const IC = ENTRY.ic
const TIME_INTERVAL = ENTRY.time_interval
const NUM_POINTS = 21        # grid density for Globtim objective
const DIM = length(P_TRUE)
const GN = 12

# Sweep axes
const RADII = [0.005, 0.01, 0.02, 0.05, 0.1]
const DEGREE_RANGE = 4:2:10

# Time budget per combo — skip refinement if Globtim alone exceeds this
const MAX_WALL_SECONDS = 180.0   # 3 min hard cutoff (generous beyond 2 min target)

# Multi-CP clustering: two CPs are "same" if closer than this fraction of domain diameter
const CLUSTER_FRACTION = 0.05

# ═══════════════════════════════════════════════════════════════════════════════
# Result storage
# ═══════════════════════════════════════════════════════════════════════════════

struct SweepRow
    radius::Float64
    degree::Int
    l2_error::Float64
    rel_l2::Float64
    n_cps_total::Int           # all real CPs from HC
    n_cps_in_domain::Int       # CPs inside the domain bounds
    n_distinct_minima::Int     # spatially distinct low-objective CPs in domain
    best_raw_obj::Float64
    raw_recovery_error::Float64
    best_refined_obj::Float64
    recovery_error::Float64
    capture_1pct::Float64
    capture_5pct::Float64
    wall_time::Float64
    poly_time::Float64
    solve_time::Float64
    process_time::Float64
    refinement_time::Float64
    capture_time::Float64
    status::String             # "ok", "timeout", "no_cps", "error"
end

# Distinct minimum info for the multi-CP table
struct MultiCPInfo
    radius::Float64
    degree::Int
    n_in_domain::Int
    n_distinct::Int
    # For each distinct minimum: (objective_value, distance_to_p_true, coordinates)
    minima::Vector{Tuple{Float64,Float64,Vector{Float64}}}
end

all_rows = SweepRow[]
all_multicp = MultiCPInfo[]

# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""Check if a point lies inside the domain bounds."""
function in_domain(pt::Vector{Float64}, bounds::Vector{Tuple{Float64,Float64}})
    all(bounds[i][1] <= pt[i] <= bounds[i][2] for i in eachindex(pt))
end

"""
Cluster in-domain CPs by spatial proximity and return distinct groups.
Each group is represented by the member with the lowest objective value.
"""
function cluster_critical_points(
    points::Vector{Vector{Float64}},
    values::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}};
    cluster_fraction::Float64 = 0.05,
)
    domain_diameter = norm([b[2] - b[1] for b in bounds])
    threshold = cluster_fraction * domain_diameter

    # Filter to in-domain points
    in_mask = [in_domain(pt, bounds) for pt in points]
    idx_in = findall(in_mask)

    if isempty(idx_in)
        return (
            n_in_domain = 0,
            n_distinct = 0,
            representatives = Vector{Float64}[],
            representative_values = Float64[],
            cluster_sizes = Int[],
        )
    end

    pts_in = points[idx_in]
    vals_in = values[idx_in]

    # Greedy clustering: sort by objective value, assign each to nearest existing
    # cluster or start a new one
    order = sortperm(vals_in)
    clusters = Vector{Int}[]        # each cluster = indices into pts_in
    centroids = Vector{Float64}[]   # representative point per cluster

    for i in order
        assigned = false
        for (ci, centroid) in enumerate(centroids)
            if norm(pts_in[i] - centroid) < threshold
                push!(clusters[ci], i)
                assigned = true
                break
            end
        end
        if !assigned
            push!(clusters, [i])
            push!(centroids, pts_in[i])
        end
    end

    # Representative = lowest objective in each cluster (first by sort order)
    rep_pts = Vector{Float64}[]
    rep_vals = Float64[]
    sizes = Int[]
    for cluster in clusters
        best = cluster[argmin([vals_in[j] for j in cluster])]
        push!(rep_pts, pts_in[best])
        push!(rep_vals, vals_in[best])
        push!(sizes, length(cluster))
    end

    return (
        n_in_domain = length(idx_in),
        n_distinct = length(clusters),
        representatives = rep_pts,
        representative_values = rep_vals,
        cluster_sizes = sizes,
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Banner
# ═══════════════════════════════════════════════════════════════════════════════

println("="^100)
println("  LV4D Domain × Degree Sweep")
println("  How large a domain can we solve reliably at degree ≤ 10?")
println("="^100)
println()
@printf("  %-18s %s\n", "Candidate:", ENTRY.name)
@printf("  %-18s [%s]\n", "P_TRUE:", join([@sprintf("%.4f", x) for x in P_TRUE], ", "))
@printf("  %-18s %s\n", "IC:", string(IC))
@printf(
    "  %-18s [%.1f, %.1f], %d Globtim points\n",
    "Time:",
    TIME_INTERVAL[1],
    TIME_INTERVAL[2],
    NUM_POINTS
)
@printf("  %-18s %d^%d = %d points per run\n", "Grid:", GN, DIM, GN^DIM)
@printf("  %-18s %s\n", "Degrees:", string(collect(DEGREE_RANGE)))
@printf("  %-18s %s\n", "Radii:", string(RADII))
@printf("  %-18s %.0fs\n", "Budget per combo:", MAX_WALL_SECONDS)
@printf("  %-18s Tsit5/1e-4 (grid) → Vern9/1e-10 (post-processing)\n", "ODE solver:")
println()

# ═══════════════════════════════════════════════════════════════════════════════
# Sweep
# ═══════════════════════════════════════════════════════════════════════════════

total_start = time()

for (ri, radius) in enumerate(RADII)
    bounds = build_bounds(P_TRUE, radius)
    domain_diameter = norm([2 * radius for _ in 1:DIM])

    println("─"^100)
    @printf(
        "  Radius %.4f  (domain width %.4f, diameter %.4f, %d/%d)\n",
        radius,
        2 * radius,
        domain_diameter,
        ri,
        length(RADII)
    )
    println("─"^100)

    for degree in DEGREE_RANGE
        combo_start = time()
        @printf("    deg %2d: ", degree)

        # Fresh objective per combo (avoids stale solver state between radii)
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

        known_cps = KnownCriticalPoints([P_TRUE], [0.0], [:min], bounds)

        config = Globtim.ExperimentParams(
            GN = GN,
            degree_range = degree:2:degree,
            domain_size = radius,
        )
        output_dir = mktempdir()

        status = "ok"
        l2_err = NaN
        rel_l2 = NaN
        n_cps_total = 0
        n_cps_in_domain = 0
        n_distinct = 0
        best_raw = NaN
        raw_rec = NaN
        best_ref = NaN
        ref_rec = NaN
        cap_1 = NaN
        cap_5 = NaN
        poly_t = 0.0
        solve_t = 0.0
        proc_t = 0.0
        ref_t = 0.0
        cap_t = 0.0

        cluster_result = nothing

        try
            # Run Globtim (coarse tolerances)
            exp_result = run_standard_experiment(
                objective_function = objective,
                objective_name = @sprintf("lv4d_r%.4f_d%d", radius, degree),
                bounds = bounds,
                experiment_config = config,
                output_dir = output_dir,
                metadata = Dict{String,Any}("radius" => radius, "degree" => degree),
                true_params = P_TRUE,
            )

            dr = exp_result[:degree_results][1]  # single degree
            combo_so_far = time() - combo_start

            l2_err = dr.l2_approx_error
            rel_l2 = dr.relative_l2_error
            n_cps_total = dr.n_critical_points
            poly_t = dr.polynomial_construction_time
            solve_t = dr.critical_point_solving_time
            proc_t = dr.critical_point_processing_time

            if dr.status != "success" || n_cps_total == 0
                status = "no_cps"
                @printf("FAILED (no CPs), L2=%.2e, %.1fs\n", l2_err, combo_so_far)
            else
                best_raw = dr.best_objective !== nothing ? dr.best_objective : NaN
                raw_rec =
                    dr.best_estimate !== nothing ? norm(dr.best_estimate - P_TRUE) : NaN

                # Multi-CP clustering
                cluster_result = cluster_critical_points(
                    dr.critical_points,
                    dr.objective_values,
                    bounds;
                    cluster_fraction = CLUSTER_FRACTION,
                )
                n_cps_in_domain = cluster_result.n_in_domain
                n_distinct = cluster_result.n_distinct

                # Check time budget before refinement
                if combo_so_far > MAX_WALL_SECONDS
                    status = "timeout"
                    @printf(
                        "TIMEOUT after Globtim (%.1fs), %d CPs (%d in-domain, %d distinct), L2=%.2e\n",
                        combo_so_far,
                        n_cps_total,
                        n_cps_in_domain,
                        n_distinct,
                        l2_err
                    )
                else
                    # Switch to tight tolerances for post-processing
                    set_tolerance!(objective, 1e-10)
                    set_solver!(objective, Vern9())

                    # Refinement
                    const_eps = radius * 1e-6
                    refinement_config = ode_refinement_config(
                        max_time_per_point = 5.0,
                        bounds = [(b[1] + const_eps, b[2] - const_eps) for b in bounds],
                        show_progress = false,
                    )

                    ref_start = time()
                    degree_analyses = run_degree_analyses(
                        [dr],
                        objective,
                        output_dir,
                        refinement_config;
                        gradient_method = :finitediff,
                    )
                    ref_t = time() - ref_start

                    if !isempty(degree_analyses)
                        ref = degree_analyses[1][3]  # RefinedExperimentResult
                        if ref.n_converged > 0
                            best_ref = ref.best_refined_value
                            best_pt = ref.refined_points[ref.best_refined_idx]
                            ref_rec = norm(best_pt - P_TRUE)
                        end
                    end

                    # Capture analysis
                    cap_start = time()
                    degree_capture_results = compute_degree_capture_results([dr], known_cps)
                    cap_t = time() - cap_start

                    for (deg, cr) in degree_capture_results
                        cap_1 = capture_rate_at(cr, 0.01)
                        cap_5 = capture_rate_at(cr, 0.05)
                    end

                    wall = time() - combo_start
                    rec_str =
                        isnan(ref_rec) ?
                        (isnan(raw_rec) ? "N/A" : @sprintf("%.2e(raw)", raw_rec)) :
                        @sprintf("%.2e", ref_rec)
                    @printf(
                        "%d CPs (%d in-domain, %d distinct), recovery=%s, L2=%.2e, %.1fs\n",
                        n_cps_total,
                        n_cps_in_domain,
                        n_distinct,
                        rec_str,
                        l2_err,
                        wall
                    )
                end
            end
        catch e
            status = "error"
            wall = time() - combo_start
            @printf("ERROR: %s (%.1fs)\n", sprint(showerror, e), wall)
        end

        wall = time() - combo_start
        push!(
            all_rows,
            SweepRow(
                radius,
                degree,
                l2_err,
                rel_l2,
                n_cps_total,
                n_cps_in_domain,
                n_distinct,
                best_raw,
                raw_rec,
                best_ref,
                ref_rec,
                cap_1,
                cap_5,
                wall,
                poly_t,
                solve_t,
                proc_t,
                ref_t,
                cap_t,
                status,
            ),
        )

        # Store multi-CP info if we found distinct minima
        if cluster_result !== nothing && cluster_result.n_distinct > 0
            minima_info = Tuple{Float64,Float64,Vector{Float64}}[]
            for (val, pt) in
                zip(cluster_result.representative_values, cluster_result.representatives)
                push!(minima_info, (val, norm(pt - P_TRUE), pt))
            end
            push!(
                all_multicp,
                MultiCPInfo(
                    radius,
                    degree,
                    cluster_result.n_in_domain,
                    cluster_result.n_distinct,
                    minima_info,
                ),
            )
        end
    end
    println()
end

total_wall = time() - total_start

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 1: Accuracy × Timing
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("="^140)
println("  TABLE 1: ACCURACY & TIMING — Radius × Degree")
@printf("  Total wall time: %.1fs (%.1f min)\n", total_wall, total_wall / 60)
println("="^140)
println()

@printf(
    "%-8s  %4s  %6s  %10s  %10s  %5s  %5s  %5s  %11s  %12s  %12s  %8s  %8s  %7s\n",
    "Radius",
    "Deg",
    "Status",
    "L2 Error",
    "Rel L2",
    "#CPs",
    "InDom",
    "Dist",
    "Raw f(x)",
    "Raw Recov",
    "Ref Recov",
    "Cap@1%",
    "Cap@5%",
    "Time(s)"
)
println("-"^140)

for row in all_rows
    rel_str = isnan(row.rel_l2) ? "N/A" : @sprintf("%.3e", row.rel_l2)
    raw_str = isnan(row.best_raw_obj) ? "N/A" : @sprintf("%.3e", row.best_raw_obj)
    raw_rc_str =
        isnan(row.raw_recovery_error) ? "N/A" : @sprintf("%.3e", row.raw_recovery_error)
    ref_rc_str = isnan(row.recovery_error) ? "N/A" : @sprintf("%.3e", row.recovery_error)
    c1_str = isnan(row.capture_1pct) ? "---" : @sprintf("%.0f%%", 100 * row.capture_1pct)
    c5_str = isnan(row.capture_5pct) ? "---" : @sprintf("%.0f%%", 100 * row.capture_5pct)
    @printf(
        "%-8.4f  %4d  %6s  %10.3e  %10s  %5d  %5d  %5d  %11s  %12s  %12s  %8s  %8s  %7.1f\n",
        row.radius,
        row.degree,
        row.status,
        row.l2_error,
        rel_str,
        row.n_cps_total,
        row.n_cps_in_domain,
        row.n_distinct_minima,
        raw_str,
        raw_rc_str,
        ref_rc_str,
        c1_str,
        c5_str,
        row.wall_time
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 2: Timing Breakdown
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("="^120)
println("  TABLE 2: TIMING BREAKDOWN (seconds)")
println("="^120)
println()

@printf(
    "%-8s  %4s  %6s  %8s  %8s  %8s  %8s  %8s  %8s  %6s\n",
    "Radius",
    "Deg",
    "Status",
    "Poly(s)",
    "HC(s)",
    "Proc(s)",
    "Ref(s)",
    "Cap(s)",
    "Total(s)",
    "%HC"
)
println("-"^100)

for row in all_rows
    total =
        row.poly_time +
        row.solve_time +
        row.process_time +
        row.refinement_time +
        row.capture_time
    pct_hc = total > 0 ? 100 * row.solve_time / total : 0.0
    @printf(
        "%-8.4f  %4d  %6s  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %5.0f%%\n",
        row.radius,
        row.degree,
        row.status,
        row.poly_time,
        row.solve_time,
        row.process_time,
        row.refinement_time,
        row.capture_time,
        total,
        pct_hc
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 3: Multi-CP Discovery
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("="^120)
println("  TABLE 3: MULTI-CRITICAL-POINT DISCOVERY")
println("  (rows with >1 distinct in-domain minimum)")
println("="^120)
println()

multi_rows = filter(m -> m.n_distinct > 1, all_multicp)

if isempty(multi_rows)
    println(
        "  No multi-CP cases found. All (radius, degree) combos yielded 0 or 1 in-domain minimum.",
    )
else
    @printf(
        "%-8s  %4s  %5s  %5s  %-60s\n",
        "Radius",
        "Deg",
        "InDom",
        "Dist",
        "Distinct Minima: (obj, dist_to_p_true)"
    )
    println("-"^120)

    for m in multi_rows
        min_strs = [@sprintf("(%.3e, d=%.4f)", val, dist) for (val, dist, _) in m.minima]
        @printf(
            "%-8.4f  %4d  %5d  %5d  %s\n",
            m.radius,
            m.degree,
            m.n_in_domain,
            m.n_distinct,
            join(min_strs, "  ")
        )
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Summary Table 4: Best per radius
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("="^120)
println("  TABLE 4: BEST PER RADIUS (lowest recovery error, status=ok only)")
println("="^120)
println()

@printf(
    "%-8s  %4s  %10s  %5s  %5s  %12s  %12s  %8s  %7s\n",
    "Radius",
    "Deg",
    "L2 Error",
    "#CPs",
    "Dist",
    "Raw Recov",
    "Ref Recov",
    "Cap@5%",
    "Time(s)"
)
println("-"^100)

for radius in RADII
    rows = filter(
        r -> r.radius == radius && r.status == "ok" && !isnan(r.raw_recovery_error),
        all_rows,
    )
    if isempty(rows)
        @printf(
            "%-8.4f  %4s  %10s  %5s  %5s  %12s  %12s  %8s  %7s\n",
            radius,
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
    # Prefer refined recovery if available, else raw
    best = rows[argmin([
        isnan(r.recovery_error) ? r.raw_recovery_error : r.recovery_error for r in rows
    ])]
    raw_rc_str =
        isnan(best.raw_recovery_error) ? "N/A" : @sprintf("%.3e", best.raw_recovery_error)
    ref_rc_str = isnan(best.recovery_error) ? "N/A" : @sprintf("%.3e", best.recovery_error)
    c5_str = isnan(best.capture_5pct) ? "---" : @sprintf("%.0f%%", 100 * best.capture_5pct)
    @printf(
        "%-8.4f  %4d  %10.3e  %5d  %5d  %12s  %12s  %8s  %7.1f\n",
        best.radius,
        best.degree,
        best.l2_error,
        best.n_cps_in_domain,
        best.n_distinct_minima,
        raw_rc_str,
        ref_rc_str,
        c5_str,
        best.wall_time
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Conclusion
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("="^100)
println("  NOTES")
println("="^100)
println()
println(
    "  Grid: GN=$GN -> $(GN^DIM) eval points per degree (~2ms/eval -> ~$(round(Int, GN^DIM * 0.002))s grid eval)",
)
println(
    "  HC solve cost scales as (d-1)^$DIM total-degree paths (deg 10 -> 9^4 = 6561 paths)",
)
println("  Larger radii -> worse polynomial fit (higher L2) but more CPs in domain")
println("  The 'frontier' is where L2 error starts degrading AND timing hits budget")
println("  Multi-CP rows show where the landscape has structure beyond a single basin")
println()
