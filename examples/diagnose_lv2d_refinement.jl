#!/usr/bin/env julia
#
# Diagnostic script: Refinement methods on LV2D landscapes
#
# Compares log_L2_norm vs L2_norm objectives across:
#   A. Landscape characterization (objective values, gradient, Hessian, grid scan)
#   B. Newton refinement (FiniteDiff vs ForwardDiff) vs NelderMead (gradient-free)
#   C. Valley detection (Hessian eigenanalysis, FiniteDiff + ForwardDiff)
#
# Usage:
#   julia --project=profiles/dev --threads=auto \
#       pkg/DynamicObjectives/examples/diagnose_lv2d_refinement.jl
#

using DynamicObjectives
using GlobtimPostProcessing:
    refine_to_critical_point, refine_point, refinement_method, method_name
using ForwardDiff, FiniteDiff, LinearAlgebra, Printf

const CATALOGUE_PATH =
    joinpath(@__DIR__, "..", "paper", "catalogue", "lv2d_all_catalogue.jsonl")

# ═══════════════════════════════════════════════════════════════════════════════
# Load both catalogue entries
# ═══════════════════════════════════════════════════════════════════════════════

println("Loading catalogue...")
entries = load_catalogue(CATALOGUE_PATH)

entry_log = entries[findfirst(e -> e.name == "LV2D_SciML_d10_log", entries)]
entry_l2 = entries[findfirst(e -> e.name == "LV2D_SciML_d10", entries)]
@assert entry_log.name == "LV2D_SciML_d10_log"
@assert entry_l2.name == "LV2D_SciML_d10"
@assert entry_l2.distance_function === L2_norm "Expected L2_norm, got $(entry_l2.distance_function)"
@assert entry_l2.p_true == entry_log.p_true "p_true mismatch between entries"

obj_log = create_objective(entry_log)
obj_l2 = create_objective(entry_l2)
p_true = collect(Float64, entry_log.p_true)  # [1.0, 1.5]
bounds_tuples = [Tuple(b) for b in entry_log.bounds]
lb = [b[1] for b in bounds_tuples]
ub = [b[2] for b in bounds_tuples]

println("="^72)
println("  Refinement Diagnostics: log_L2_norm vs L2_norm")
println("="^72)
println("  p_true       = $p_true")
println("  bounds       = $bounds_tuples")
println("  f_log(p_true) = $(obj_log(p_true))")
println("  f_L2(p_true)  = $(obj_l2(p_true))")
println("="^72)

# ═══════════════════════════════════════════════════════════════════════════════
# Shared test infrastructure
# ═══════════════════════════════════════════════════════════════════════════════

struct PointDiag
    p::Vector{Float64}
    f::Float64
    grad::Vector{Float64}
    grad_norm::Float64
    hessian::Matrix{Float64}
    eigenvalues::Vector{Float64}
    eigenvectors::Matrix{Float64}
    condition_number::Float64
    is_positive_definite::Bool
    valley_dim::Int          # number of near-zero eigenvalues
end

function diagnose_point(
    obj,
    p::Vector{Float64};
    eig_threshold::Float64 = 1e-3,
    method::Symbol = :finitediff,
)
    f = obj(p)
    g = if method == :forwarddiff
        ForwardDiff.gradient(obj, p)
    elseif method == :finitediff
        FiniteDiff.finite_difference_gradient(obj, p)
    else
        error("Unknown method: $method. Use :forwarddiff or :finitediff")
    end
    H = if method == :forwarddiff
        ForwardDiff.hessian(obj, p)
    elseif method == :finitediff
        FiniteDiff.finite_difference_hessian(obj, p)
    else
        error("Unknown method: $method. Use :forwarddiff or :finitediff")
    end
    eig = eigen(Symmetric(H))
    λ = eig.values
    V = eig.vectors
    cond_num = maximum(abs, λ) / max(minimum(abs, λ), 1e-30)
    pd = all(λ .> 0)
    vdim = count(abs.(λ) .< eig_threshold)
    PointDiag(copy(p), f, g, norm(g), H, λ, V, cond_num, pd, vdim)
end

# Starting points shared across all refinement tests
start_points = [
    ([1.0, 1.5], "p_true exact"),
    ([1.001, 1.501], "p_true + 0.001"),
    ([1.01, 1.51], "p_true + 0.01"),
    ([1.05, 1.55], "p_true + 0.05"),
    ([0.95, 1.45], "p_true - 0.05"),
    ([1.1, 1.6], "p_true + 0.1"),
    ([0.9, 1.4], "p_true - 0.1"),
    ([1.2, 1.7], "p_true + 0.2"),
    ([0.8, 1.3], "p_true - 0.2"),
    ([1.5, 2.0], "p_true + 0.5"),
    ([0.5, 1.0], "p_true - 0.5"),
    ([2.0, 2.0], "near p_true"),
    ([5.0, 5.0], "center"),
    ([0.5, 0.5], "corner low"),
    ([9.0, 9.0], "corner high"),
]

function in_bounds(p)
    all(p[i] >= bounds_tuples[i][1] && p[i] <= bounds_tuples[i][2] for i in 1:length(p))
end

# ═════════════════════════════════════════════════════════════════════════════
#  PART A: Landscape Characterization (both objectives, side by side)
# ═════════════════════════════════════════════════════════════════════════════

println("\n" * "═"^72)
println("  PART A: Landscape Characterization")
println("═"^72)

# ─── Test 1: Radial offsets ─────────────────────────────────────────────────

println("\n" * "─"^72)
println("  TEST 1: Objective along radial offsets from p_true")
println("─"^72)
@printf("  %-12s  %16s  %16s\n", "δ", "f_log", "f_L2")
println("  " * "─"^48)

for δ in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
    p = p_true .+ δ
    @printf("  %-12.3f  %+16.6f  %+16.6f\n", δ, obj_log(p), obj_l2(p))
end

# ─── Test 2: Gradient & Hessian at p_true ───────────────────────────────────

println("\n" * "─"^72)
println("  TEST 2: Gradient & Hessian at p_true = $p_true")
println("─"^72)

# FiniteDiff diagnostics
d_log = diagnose_point(obj_log, p_true; method = :finitediff)
d_l2 = diagnose_point(obj_l2, p_true; method = :finitediff)

# ForwardDiff diagnostics (only L2_norm is ForwardDiff-compatible via SciMLBase.remake)
d_l2_fd = diagnose_point(obj_l2, p_true; method = :forwarddiff)

for (label, d, d_fd) in [("log_L2_norm", d_log, nothing), ("L2_norm", d_l2, d_l2_fd)]
    println("\n  [$label] (FiniteDiff)")
    println("    f(p_true)   = $(d.f)")
    println("    grad        = $(d.grad)")
    println("    ||grad||    = $(d.grad_norm)")
    println("    Hessian:")
    @printf("      [%+.6e  %+.6e]\n", d.hessian[1, 1], d.hessian[1, 2])
    @printf("      [%+.6e  %+.6e]\n", d.hessian[2, 1], d.hessian[2, 2])
    println("    eigenvalues = $(d.eigenvalues)")
    @printf("    condition   = %.2e\n", d.condition_number)
    println("    pos. def.?  = $(d.is_positive_definite)")
    println("    valley_dim  = $(d.valley_dim)")

    if d_fd !== nothing
        println()
        println("  [$label] (ForwardDiff)")
        println("    f(p_true)   = $(d_fd.f)")
        println("    grad        = $(d_fd.grad)")
        println("    ||grad||    = $(d_fd.grad_norm)")
        println("    Hessian:")
        @printf("      [%+.6e  %+.6e]\n", d_fd.hessian[1, 1], d_fd.hessian[1, 2])
        @printf("      [%+.6e  %+.6e]\n", d_fd.hessian[2, 1], d_fd.hessian[2, 2])
        println("    eigenvalues = $(d_fd.eigenvalues)")
        @printf("    condition   = %.2e\n", d_fd.condition_number)
        println("    pos. def.?  = $(d_fd.is_positive_definite)")
        println("    valley_dim  = $(d_fd.valley_dim)")
        println()
        # Agreement check (away from cusp this should be ~1e-9)
        grad_relerr = norm(d.grad - d_fd.grad) / max(norm(d.grad), norm(d_fd.grad), 1e-30)
        @printf("    FD vs FiniteDiff gradient rel. error = %.2e\n", grad_relerr)
        @printf("    NOTE: At p_true, L2_norm = sqrt(v'v) → 0. The gradient of sqrt(x)\n")
        @printf("          at x=0 is undefined. ForwardDiff propagates this faithfully;\n")
        @printf("          FiniteDiff averages over the cusp → much smaller ||grad||.\n")
    end
end

# ─── Test 3: Grid definiteness scan ────────────────────────────────────────

function print_grid(obj, label)
    println("\n" * "─"^72)
    println("  TEST 3: Hessian definiteness grid [$label]")
    println("─"^72)
    println("  Legend: + = PD, ~ = indefinite, - = ND, V = valley (near-zero eig)")
    println()

    grid_n = 10
    p1_grid = range(0.1, 9.9, length = grid_n)
    p2_grid = range(0.1, 9.9, length = grid_n)

    n_pd = 0
    n_indef = 0
    n_nd = 0
    n_valley = 0

    print("  p₂\\p₁ ")
    for p1 in p1_grid
        @printf(" %4.1f", p1)
    end
    println()
    print("        ")
    println("─"^(5 * grid_n))

    for p2 in reverse(p2_grid)
        @printf("  %4.1f |", p2)
        for p1 in p1_grid
            p = [p1, p2]
            try
                H = FiniteDiff.finite_difference_hessian(obj, p)
                λ = eigvals(Symmetric(H))
                has_valley = any(abs.(λ) .< 1e-3)
                if has_valley
                    print("    V")
                    n_valley += 1
                elseif all(λ .> 0)
                    print("    +")
                    n_pd += 1
                elseif all(λ .< 0)
                    print("    -")
                    n_nd += 1
                else
                    print("    ~")
                    n_indef += 1
                end
            catch e
                print("    !")
                @warn "Hessian evaluation failed at p=[$p1, $p2]" exception =
                    (e, catch_backtrace())
            end
        end
        println()
    end
    println()
    println("  PD: $n_pd  Indef: $n_indef  ND: $n_nd  Valley: $n_valley  / $(grid_n^2)")
end

print_grid(obj_log, "log_L2_norm")
print_grid(obj_l2, "L2_norm")

# ═════════════════════════════════════════════════════════════════════════════
#  PART B: Newton vs NelderMead Refinement (L2_norm only)
# ═════════════════════════════════════════════════════════════════════════════
#
# log_L2_norm is proven unsuitable for gradient-based refinement (Test 2 shows
# ||grad||=2.4 and cond(H)~1e10 at p_true). Focus refinement comparison on L2_norm.

println("\n" * "═"^72)
println("  PART B: Refinement Comparison on L2_norm")
println("═"^72)

# ─── Test 4: Newton refinement ──────────────────────────────────────────────

println("\n" * "─"^72)
println("  TEST 4: Newton refinement (trust-region, FiniteDiff) [L2_norm]")
println("─"^72)
@printf(
    "  %-20s  %-8s  %10s  %10s  %5s  %-22s  %-10s  %s\n",
    "Start",
    "CP Type",
    "f(final)",
    "||grad||",
    "Iters",
    "Final point",
    "Status",
    "d(p_true)"
)
println("  " * "─"^110)

for (start, slabel) in start_points
    r = refine_to_critical_point(
        obj_l2,
        start;
        gradient_method = :finitediff,
        trace = false,
        bounds = bounds_tuples,
        max_iterations = 200,
    )
    status = r.converged ? "CONVERGED" : "NOT CONV"
    dist = norm(r.point .- p_true)
    @printf(
        "  %-20s  %-8s  %+.3e  %.2e  %5d  [%+.5f, %+.5f]  %-10s  %.2e\n",
        slabel,
        r.cp_type,
        r.objective_value,
        r.gradient_norm,
        r.iterations,
        r.point[1],
        r.point[2],
        status,
        dist
    )
end

# ─── Test 4b: Newton refinement with ForwardDiff ────────────────────────────

println("\n" * "─"^72)
println("  TEST 4b: Newton refinement (trust-region, ForwardDiff) [L2_norm]")
println("─"^72)
@printf(
    "  %-20s  %-8s  %10s  %10s  %5s  %-22s  %-10s  %s\n",
    "Start",
    "CP Type",
    "f(final)",
    "||grad||",
    "Iters",
    "Final point",
    "Status",
    "d(p_true)"
)
println("  " * "─"^110)

for (start, slabel) in start_points
    r = refine_to_critical_point(
        obj_l2,
        start;
        gradient_method = :forwarddiff,
        trace = false,
        bounds = bounds_tuples,
        max_iterations = 200,
    )
    status = r.converged ? "CONVERGED" : "NOT CONV"
    dist = norm(r.point .- p_true)
    @printf(
        "  %-20s  %-8s  %+.3e  %.2e  %5d  [%+.5f, %+.5f]  %-10s  %.2e\n",
        slabel,
        r.cp_type,
        r.objective_value,
        r.gradient_norm,
        r.iterations,
        r.point[1],
        r.point[2],
        status,
        dist
    )
end

# ─── Test 4c: Newton refinement with mode=:minimize ────────────────────────

println("\n" * "─"^72)
println(
    "  TEST 4c: Newton refinement (modified Newton, mode=:minimize, FiniteDiff) [L2_norm]",
)
println("─"^72)
@printf(
    "  %-20s  %-8s  %10s  %10s  %5s  %-22s  %-10s  %s\n",
    "Start",
    "CP Type",
    "f(final)",
    "||grad||",
    "Iters",
    "Final point",
    "Status",
    "d(p_true)"
)
println("  " * "─"^110)

for (start, slabel) in start_points
    r = refine_to_critical_point(
        obj_l2,
        start;
        gradient_method = :finitediff,
        trace = false,
        bounds = bounds_tuples,
        max_iterations = 200,
        mode = :minimize,
    )
    status = r.converged ? "CONVERGED" : "NOT CONV"
    dist = norm(r.point .- p_true)
    @printf(
        "  %-20s  %-8s  %+.3e  %.2e  %5d  [%+.5f, %+.5f]  %-10s  %.2e\n",
        slabel,
        r.cp_type,
        r.objective_value,
        r.gradient_norm,
        r.iterations,
        r.point[1],
        r.point[2],
        status,
        dist
    )
end

# ─── Test 5: NelderMead refinement ──────────────────────────────────────────

println("\n" * "─"^72)
println("  TEST 5: NelderMead refinement (gradient-free) [L2_norm]")
println("─"^72)
@printf(
    "  %-20s  %-8s  %10s  %10s  %5s  %-22s  %-10s  %s\n",
    "Start",
    "CP Type",
    "f(final)",
    "||grad||",
    "Iters",
    "Final point",
    "Status",
    "d(p_true)"
)
println("  " * "─"^110)

nm_method = refinement_method(:neldermead; gradient_method = :finitediff)
for (start, slabel) in start_points
    r = refine_point(nm_method, obj_l2, start; bounds = bounds_tuples)
    status = r.converged ? "CONVERGED" : "NOT CONV"
    dist = norm(r.point .- p_true)
    @printf(
        "  %-20s  %-8s  %+.3e  %.2e  %5d  [%+.5f, %+.5f]  %-10s  %.2e\n",
        slabel,
        r.cp_type,
        r.objective_value,
        r.gradient_norm,
        r.iterations,
        r.point[1],
        r.point[2],
        status,
        dist
    )
end

# ═════════════════════════════════════════════════════════════════════════════
#  PART C: Valley Detection (FiniteDiff + ForwardDiff eigenanalysis)
# ═════════════════════════════════════════════════════════════════════════════

println("\n" * "═"^72)
println("  PART C: Valley Detection")
println("═"^72)
println()
println("  Valley dimension = number of Hessian eigenvalues with |λ| < 1e-3")
println("  (threshold matches ValleyWalkConfig default)")
println()

# Collect test points: p_true + NelderMead endpoints on L2_norm
valley_test_points = [
    (p_true, "p_true"),
    ([0.0 + 1e-8, 2.91], "Newton endpoint (log)"),
    ([5.54, 5.58], "Newton saddle (log)"),
]

nm_valley = refinement_method(:neldermead; gradient_method = :finitediff)
for (start, slabel) in
    [([1.05, 1.55], "near+0.05"), ([0.8, 1.3], "near-0.2"), ([5.0, 5.0], "center")]
    r = refine_point(nm_valley, obj_l2, start; bounds = bounds_tuples)
    push!(valley_test_points, (r.point, "NM_L2 from $slabel"))
end

for (label, obj) in [("log_L2_norm", obj_log), ("L2_norm", obj_l2)]
    println("\n" * "─"^72)
    println("  TEST 6: Valley detection [$label] (FiniteDiff)")
    println("─"^72)
    @printf(
        "  %-28s  %10s  %10s  %6s  %s\n",
        "Point",
        "f",
        "||grad||",
        "V-dim",
        "Eigenvalues"
    )
    println("  " * "─"^90)

    for (pt, plabel) in valley_test_points
        if !in_bounds(pt)
            @printf("  %-28s  OUT OF BOUNDS\n", plabel)
            continue
        end
        d = diagnose_point(obj, pt; method = :finitediff)
        @printf(
            "  %-28s  %+.3e  %.2e  %4d   [%s]\n",
            plabel,
            d.f,
            d.grad_norm,
            d.valley_dim,
            join([@sprintf("%+.3e", e) for e in d.eigenvalues], ", ")
        )
    end
end

# ForwardDiff valley detection on L2_norm only (log_L2_norm is not ForwardDiff-compatible)
println("\n" * "─"^72)
println("  TEST 6b: Valley detection [L2_norm] (ForwardDiff)")
println("─"^72)
@printf("  %-28s  %10s  %10s  %6s  %s\n", "Point", "f", "||grad||", "V-dim", "Eigenvalues")
println("  " * "─"^90)

for (pt, plabel) in valley_test_points
    if !in_bounds(pt)
        @printf("  %-28s  OUT OF BOUNDS\n", plabel)
        continue
    end
    d = diagnose_point(obj_l2, pt; method = :forwarddiff)
    @printf(
        "  %-28s  %+.3e  %.2e  %4d   [%s]\n",
        plabel,
        d.f,
        d.grad_norm,
        d.valley_dim,
        join([@sprintf("%+.3e", e) for e in d.eigenvalues], ", ")
    )
end

# ═════════════════════════════════════════════════════════════════════════════
#  SUMMARY
# ═════════════════════════════════════════════════════════════════════════════

println("\n" * "═"^72)
println("  SUMMARY")
println("═"^72)

# ─── Landscape at p_true ────────────────────────────────────────────────────

println("\n  Landscape at p_true = $p_true:")
println("                        log_L2_norm    L2_norm (FD)   L2_norm (AD)")
println("  ──────────────────────────────────────────────────────────────────────")
@printf("  f(p_true)       %+16.6e  %+16.6e  %+16.6e\n", d_log.f, d_l2.f, d_l2_fd.f)
@printf(
    "  ||grad(p_true)||%16.2e  %16.2e  %16.2e\n",
    d_log.grad_norm,
    d_l2.grad_norm,
    d_l2_fd.grad_norm
)
@printf(
    "  λ_min           %+16.3e  %+16.3e  %+16.3e\n",
    d_log.eigenvalues[1],
    d_l2.eigenvalues[1],
    d_l2_fd.eigenvalues[1]
)
@printf(
    "  λ_max           %+16.3e  %+16.3e  %+16.3e\n",
    d_log.eigenvalues[end],
    d_l2.eigenvalues[end],
    d_l2_fd.eigenvalues[end]
)
@printf(
    "  cond(H)         %16.1e  %16.1e  %16.1e\n",
    d_log.condition_number,
    d_l2.condition_number,
    d_l2_fd.condition_number
)
@printf(
    "  PD?             %16s  %16s  %16s\n",
    d_log.is_positive_definite ? "yes" : "NO",
    d_l2.is_positive_definite ? "yes" : "NO",
    d_l2_fd.is_positive_definite ? "yes" : "NO"
)
@printf(
    "  valley_dim      %16d  %16d  %16d\n",
    d_log.valley_dim,
    d_l2.valley_dim,
    d_l2_fd.valley_dim
)
println()
println("  FD = FiniteDiff,  AD = ForwardDiff")
println("  NOTE: At p_true, L2_norm = sqrt(v'v) → 0. ForwardDiff propagates the")
println("        sqrt(0) singularity faithfully (large ||grad||); FiniteDiff averages")
println("        over the cusp → much smaller ||grad||. Both are correct for their method.")

# ─── Refinement comparison (all methods) ──────────────────────────────────

println("\n" * "─"^72)
println(
    "  COMPARATIVE REFINEMENT (L2_norm, $(length(start_points)) starting points, d < 0.01 = found)",
)
println("─"^72)

# Define all methods to compare using unified refinement_method factory
methods_to_test = [
    refinement_method(:newton_cp; gradient_method = :finitediff, max_iterations = 200),
    refinement_method(
        :newton_minimize;
        gradient_method = :finitediff,
        max_iterations = 200,
    ),
    refinement_method(:neldermead; gradient_method = :finitediff),
    refinement_method(:lbfgs; gradient_method = :finitediff),
    refinement_method(:bfgs; gradient_method = :finitediff),
    refinement_method(:conjugategradient; gradient_method = :finitediff),
]

# Collect results: (name, n_found, mean_f, mean_iters, mean_time)
method_stats = []

for m in methods_to_test
    mname = String(method_name(m))
    n_found = 0
    sum_f = 0.0
    sum_iters = 0
    sum_time = 0.0
    n_ok = 0

    for (start, _) in start_points
        t0 = time()
        r = refine_point(m, obj_l2, start; bounds = bounds_tuples)
        elapsed = time() - t0
        n_ok += 1
        sum_f += r.objective_value
        sum_iters += r.iterations
        sum_time += elapsed
        if norm(r.point .- p_true) < 0.01
            n_found += 1
        end
    end

    mean_f = n_ok > 0 ? sum_f / n_ok : NaN
    mean_iters = n_ok > 0 ? sum_iters / n_ok : NaN
    mean_time = n_ok > 0 ? sum_time / n_ok : NaN
    push!(method_stats, (mname, n_found, mean_f, mean_iters, mean_time, n_ok))
end

# Print summary table
n_pts = length(start_points)
@printf(
    "  %-24s  %12s  %12s  %10s  %10s\n",
    "Method",
    "Found p_true",
    "Mean f(final)",
    "Mean iters",
    "Mean time"
)
println("  " * "─"^72)
for (mname, n_found, mean_f, mean_iters, mean_time, n_ok) in method_stats
    err_str = n_ok < n_pts ? " ($(n_pts - n_ok) err)" : ""
    @printf(
        "  %-24s  %5d / %-5d  %12.3e  %10.1f  %9.3fs%s\n",
        mname,
        n_found,
        n_pts,
        mean_f,
        mean_iters,
        mean_time,
        err_str
    )
end

println("\n" * "="^72)
println("  Diagnostics complete.")
println("="^72)
