# ============================================================================
# Grid-Based Interestingness Scoring
# ============================================================================
#
# Evaluates ODE objective functions on structured grids to compute landscape
# complexity metrics: local minima count, dynamic range, basin fraction,
# deceptive minima. Complements the cheap random-probe approach in screening.jl
# with high-fidelity grid-based analysis on top-N candidates.
#
# Typical usage:
#   result = score_landscape_grid(model, outputs, ic, p_true, bounds, tspan)
#   score  = interestingness_score(result)
#
# Or batch over screening results:
#   scores = score_top_candidates(screening_result, model, outputs, ic, bounds, tspan)

# ============================================================================
# Result type
# ============================================================================

"""
    GridScoreResult

Results from grid-based landscape evaluation. Contains both raw grid data
and derived interestingness metrics.

# Fields
- `grid_values::Array{Float64}`: objective values on structured grid (N-dimensional array)
- `grid_axes::Vector{Vector{Float64}}`: axis coordinates per dimension
- `grid_size::Vector{Int}`: number of points per dimension
- `n_local_minima::Int`: count of grid-based local minima
- `local_minima_indices::Vector{CartesianIndex}`: grid indices of local minima
- `local_minima_values::Vector{Float64}`: objective values at local minima
- `local_minima_points::Vector{Vector{Float64}}`: parameter-space coordinates of local minima
- `dynamic_range::Float64`: log10(max/min) of finite nonzero grid values
- `basin_fraction::Float64`: fraction of grid points whose nearest local minimum is closest to p_true
- `n_deceptive::Int`: local minima with low objective but far from p_true
- `deceptive_points::Vector{Vector{Float64}}`: locations of deceptive minima
- `n_finite::Int`: grid points with finite objective values
- `n_total::Int`: total grid points evaluated
- `evaluation_time_ms::Float64`: wall time for grid evaluation
- `directional_variation::Vector{Float64}`: per-dimension variation ratios (mean |∂f/∂d| / value range)
- `curvature_score::Float64`: min directional variation normalized to [0,1] (1=curved, 0=flat)
- `basin_depth_ratio::Float64`: best competing minimum value / p_true minimum value (1.0 = equally deep)
- `conditioning::Float64`: log10 condition number of approximate Hessian at p_true
- `ruggedness::Float64`: second-order / first-order variation ratio (high = bumpy)
- `plateau_fraction::Float64`: fraction of grid with near-zero gradient that aren't minima
"""
struct GridScoreResult{N}
    grid_values::Array{Float64}
    grid_axes::Vector{Vector{Float64}}
    grid_size::Vector{Int}

    n_local_minima::Int
    local_minima_indices::Vector{CartesianIndex{N}}
    local_minima_values::Vector{Float64}
    local_minima_points::Vector{Vector{Float64}}

    dynamic_range::Float64
    basin_fraction::Float64

    n_deceptive::Int
    deceptive_points::Vector{Vector{Float64}}

    n_finite::Int
    n_total::Int
    evaluation_time_ms::Float64

    directional_variation::Vector{Float64}
    curvature_score::Float64

    basin_depth_ratio::Float64
    conditioning::Float64
    ruggedness::Float64
    plateau_fraction::Float64
end

# ============================================================================
# Core grid evaluation
# ============================================================================

"""
    score_landscape_grid(model, outputs, ic, p_true, bounds, time_interval; kwargs...) -> GridScoreResult

Evaluate an ODE-based objective function on a structured Cartesian grid and
compute landscape interestingness metrics.

# Arguments
- `model`: ModelingToolkit system
- `outputs`: observation equations
- `ic`: initial conditions
- `p_true`: true parameter vector
- `bounds`: per-dimension bounds `[(lo1,hi1), (lo2,hi2), ...]`
- `time_interval`: ODE integration interval `[t0, tf]`

# Keyword arguments
- `points_per_dim::Union{Int,Nothing}=nothing`: grid resolution per dimension.
   Auto-selects based on dimension: 20 (2D), 12 (3D), 8 (4D), 6 (5D+).
- `numpoints::Int=30`: time samples for ODE trajectory comparison
- `distance_function=L2_norm`: trajectory distance metric
- `aggregate_distances=sum`: how to aggregate per-output distances
- `solver=Tsit5()`: ODE solver (moderate tolerance for speed)
- `abstol::Float64=1e-6`: absolute ODE tolerance
- `reltol::Float64=1e-6`: relative ODE tolerance
- `deceptive_threshold::Float64=2.0`: a local minimum is "deceptive" if
   its objective value < `deceptive_threshold * global_min_value` ...
- `deceptive_distance::Float64=0.1`: ... AND its normalized distance from
   p_true exceeds this fraction of the domain diagonal
"""
function score_landscape_grid(
    model,
    outputs,
    ic::Vector{Float64},
    p_true::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}},
    time_interval;
    points_per_dim::Union{Int,Nothing} = nothing,
    numpoints::Int = 30,
    distance_function = L2_norm,
    aggregate_distances = sum,
    solver = Tsit5(),
    abstol::Float64 = 1e-6,
    reltol::Float64 = 1e-6,
    deceptive_threshold::Float64 = 2.0,
    deceptive_distance::Float64 = 0.1,
)
    dim = length(bounds)

    # Auto-select resolution
    ppd = if points_per_dim !== nothing
        points_per_dim
    else
        _default_points_per_dim(dim)
    end

    # Build axis coordinates
    grid_axes = [collect(range(b[1], b[2]; length = ppd)) for b in bounds]
    grid_size = fill(ppd, dim)

    # Build objective function
    obj = make_error_distance(
        model,
        outputs,
        ic,
        p_true,
        time_interval,
        numpoints,
        distance_function,
        aggregate_distances;
        return_inf_on_error = true,
        solver = solver,
        abstol = abstol,
        reltol = reltol,
    )

    # Evaluate on grid
    t0 = time_ns()
    grid_values = _evaluate_grid(obj, grid_axes)
    elapsed_ms = (time_ns() - t0) / 1e6

    n_total = length(grid_values)
    finite_mask = isfinite.(grid_values)
    n_finite = count(finite_mask)

    # Detect local minima
    minima_indices = _detect_local_minima(grid_values)
    minima_values = [grid_values[idx] for idx in minima_indices]
    minima_points = [_index_to_point(idx, grid_axes) for idx in minima_indices]

    # Dynamic range
    dyn_range = _compute_dynamic_range(grid_values)

    # Basin fraction (what fraction of grid is nearest to p_true's basin)
    basin_frac = _compute_basin_fraction(
        grid_values,
        grid_axes,
        minima_indices,
        minima_points,
        p_true,
    )

    # Deceptive minima
    deceptive_pts = _detect_deceptive_minima(
        minima_points,
        minima_values,
        p_true,
        bounds,
        deceptive_threshold,
        deceptive_distance,
    )

    # Directional variation (curvature/degeneracy)
    dir_var, curv_score = _compute_directional_variation(grid_values)

    # New metrics
    depth_ratio = _compute_basin_depth_ratio(minima_values, minima_points, p_true)
    cond = _compute_conditioning(grid_values, grid_axes, p_true)
    rugged = _compute_ruggedness(grid_values)
    plateau = _compute_plateau_fraction(grid_values, minima_indices)

    return GridScoreResult(
        grid_values,
        grid_axes,
        grid_size,
        length(minima_indices),
        minima_indices,
        minima_values,
        minima_points,
        dyn_range,
        basin_frac,
        length(deceptive_pts),
        deceptive_pts,
        n_finite,
        n_total,
        elapsed_ms,
        dir_var,
        curv_score,
        depth_ratio,
        cond,
        rugged,
        plateau,
    )
end

"""
    _default_points_per_dim(dim::Int) -> Int

Auto-select grid resolution based on dimensionality.
Balances coverage against cost (total points = ppd^dim).
"""
function _default_points_per_dim(dim::Int)::Int
    if dim <= 2
        20
    elseif dim == 3
        12
    elseif dim == 4
        8
    else
        6
    end
end

# ============================================================================
# Grid evaluation
# ============================================================================

"""
    _evaluate_grid(obj, grid_axes) -> Array{Float64}

Evaluate `obj(p)` on the Cartesian product of `grid_axes`.
Returns an N-dimensional array of objective values.
"""
function _evaluate_grid(obj, grid_axes::Vector{Vector{Float64}})::Array{Float64}
    dim = length(grid_axes)
    sizes = Tuple(length(ax) for ax in grid_axes)
    values = Array{Float64}(undef, sizes)

    for idx in CartesianIndices(values)
        p = [grid_axes[d][idx[d]] for d in 1:dim]
        values[idx] = obj(p)
    end

    return values
end

# ============================================================================
# Local minima detection
# ============================================================================

"""
    _detect_local_minima(grid_values::Array{Float64}) -> Vector{CartesianIndex}

Find grid points that are strictly lower than all face-adjacent neighbors
(2d neighbors in d dimensions). Infinite values are never minima.
"""
function _detect_local_minima(grid_values::Array{Float64})
    dims = size(grid_values)
    ndim = ndims(grid_values)
    minima = Vector{CartesianIndex{ndim}}()

    for idx in CartesianIndices(grid_values)
        val = grid_values[idx]
        !isfinite(val) && continue

        is_min = true
        for d in 1:ndim
            for offset in (-1, 1)
                neighbor_i = idx[d] + offset
                (neighbor_i < 1 || neighbor_i > dims[d]) && continue

                # Build neighbor index
                neighbor_idx =
                    CartesianIndex(ntuple(k -> k == d ? neighbor_i : idx[k], ndim))

                neighbor_val = grid_values[neighbor_idx]
                # Non-finite neighbors don't block minimum status
                if isfinite(neighbor_val) && neighbor_val <= val
                    is_min = false
                    break
                end
            end
            !is_min && break
        end

        is_min && push!(minima, idx)
    end

    return minima
end

# ============================================================================
# Metric computations
# ============================================================================

"""
    _index_to_point(idx::CartesianIndex, grid_axes) -> Vector{Float64}

Convert a CartesianIndex into parameter-space coordinates.
"""
function _index_to_point(
    idx::CartesianIndex,
    grid_axes::Vector{Vector{Float64}},
)::Vector{Float64}
    return [grid_axes[d][idx[d]] for d in 1:length(grid_axes)]
end

"""
    _compute_dynamic_range(grid_values::Array{Float64}) -> Float64

Compute log10(max/min) of positive finite grid values.
Returns 0.0 if fewer than 2 distinct positive finite values exist.
"""
function _compute_dynamic_range(grid_values::Array{Float64})::Float64
    finite_positive = filter(v -> isfinite(v) && v > 0, vec(grid_values))
    length(finite_positive) < 2 && return 0.0

    vmin = minimum(finite_positive)
    vmax = maximum(finite_positive)
    vmin == vmax && return 0.0

    return log10(vmax / vmin)
end

"""
    _compute_basin_fraction(grid_values, grid_axes, minima_indices, minima_points, p_true) -> Float64

Fraction of finite grid points whose nearest local minimum (by Euclidean distance
in parameter space) is the one closest to `p_true`.

Returns 0.0 if no local minima exist.
"""
function _compute_basin_fraction(
    grid_values::Array{Float64},
    grid_axes::Vector{Vector{Float64}},
    minima_indices::Vector{<:CartesianIndex},
    minima_points::Vector{Vector{Float64}},
    p_true::Vector{Float64},
)::Float64
    isempty(minima_points) && return 0.0

    # Find which minimum is closest to p_true
    dists_to_ptrue = [norm(mp .- p_true) for mp in minima_points]
    best_idx = argmin(dists_to_ptrue)

    # Count grid points whose nearest minimum is the p_true-closest one
    n_in_basin = 0
    n_finite = 0

    for idx in CartesianIndices(grid_values)
        !isfinite(grid_values[idx]) && continue
        n_finite += 1

        pt = _index_to_point(idx, grid_axes)
        dists = [norm(pt .- mp) for mp in minima_points]
        nearest = argmin(dists)

        if nearest == best_idx
            n_in_basin += 1
        end
    end

    n_finite == 0 && return 0.0
    return n_in_basin / n_finite
end

"""
    _detect_deceptive_minima(minima_points, minima_values, p_true, bounds,
                            threshold, distance_frac) -> Vector{Vector{Float64}}

Identify local minima that are "deceptive": they have low objective value
(within `threshold` multiplicative factor of the global grid minimum) but are
far from `p_true` (normalized distance > `distance_frac` of the domain diagonal).
"""
function _detect_deceptive_minima(
    minima_points::Vector{Vector{Float64}},
    minima_values::Vector{Float64},
    p_true::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}},
    threshold::Float64,
    distance_frac::Float64,
)::Vector{Vector{Float64}}
    isempty(minima_values) && return Vector{Float64}[]

    # Domain diagonal for normalization
    diag = norm([b[2] - b[1] for b in bounds])
    diag == 0 && return Vector{Float64}[]

    # Global minimum value on grid
    global_min = minimum(filter(isfinite, minima_values))

    deceptive = Vector{Float64}[]
    for (pt, val) in zip(minima_points, minima_values)
        !isfinite(val) && continue

        # Low objective value (close to global min)?
        is_low = if global_min > 0
            val < threshold * global_min
        elseif global_min < 0
            # For negative objectives (e.g. log_L2_norm), "low" means more negative
            val < global_min / threshold
        else
            val == 0.0
        end

        # Far from p_true?
        is_far = norm(pt .- p_true) / diag > distance_frac

        if is_low && is_far
            push!(deceptive, pt)
        end
    end

    return deceptive
end

# ============================================================================
# Directional variation (curvature/degeneracy detection)
# ============================================================================

"""
    _compute_directional_variation(grid_values::Array{Float64}) -> (Vector{Float64}, Float64)

Compute per-dimension variation ratios to detect degenerate (flat) landscapes.

For each dimension d, computes the mean absolute finite difference along that axis,
normalized by the total value range. A landscape that is flat along dimension d
will have near-zero variation ratio for that dimension.

Returns `(variations, curvature_score)` where:
- `variations`: per-dimension variation ratios
- `curvature_score`: `clamp(min(variations) / 0.1, 0, 1)` — 0 means flat, 1 means curved
"""
function _compute_directional_variation(
    grid_values::Array{Float64},
)::Tuple{Vector{Float64},Float64}
    ndim = ndims(grid_values)
    dims = size(grid_values)

    # Compute value range over finite values
    finite_vals = filter(isfinite, vec(grid_values))
    if length(finite_vals) < 2
        return (zeros(ndim), 0.0)
    end
    val_range = maximum(finite_vals) - minimum(finite_vals)
    if val_range == 0.0
        return (zeros(ndim), 0.0)
    end

    variations = zeros(ndim)

    for d in 1:ndim
        diffs_sum = 0.0
        n_diffs = 0

        for idx in CartesianIndices(grid_values)
            idx[d] >= dims[d] && continue  # no next neighbor along d

            # Build neighbor index (next along dimension d)
            neighbor_idx = CartesianIndex(ntuple(k -> k == d ? idx[k] + 1 : idx[k], ndim))

            v1 = grid_values[idx]
            v2 = grid_values[neighbor_idx]

            if isfinite(v1) && isfinite(v2)
                diffs_sum += abs(v2 - v1)
                n_diffs += 1
            end
        end

        if n_diffs > 0
            variations[d] = (diffs_sum / n_diffs) / val_range
        end
    end

    curvature_score = clamp(minimum(variations) / 0.1, 0.0, 1.0)
    return (variations, curvature_score)
end

# ============================================================================
# Basin depth ratio
# ============================================================================

"""
    _compute_basin_depth_ratio(minima_values, minima_points, p_true) -> Float64

Ratio of best competing minimum value to the p_true minimum value.
Close to 1.0 means a competing minimum is equally deep — harder optimization.
Returns 0.0 if ≤1 minimum exists.
"""
function _compute_basin_depth_ratio(
    minima_values::Vector{Float64},
    minima_points::Vector{Vector{Float64}},
    p_true::Vector{Float64},
)::Float64
    length(minima_values) <= 1 && return 0.0

    # Find which minimum is closest to p_true
    dists = [norm(mp .- p_true) for mp in minima_points]
    true_idx = argmin(dists)
    true_val = minima_values[true_idx]

    # Collect other minima values
    other_values = [minima_values[i] for i in eachindex(minima_values) if i != true_idx]
    filter!(isfinite, other_values)
    isempty(other_values) && return 0.0

    !isfinite(true_val) && return 0.0

    if true_val > 0
        # Positive objectives: ratio of best competitor to true minimum
        best_other = minimum(other_values)
        best_other <= 0 && return 0.0
        return clamp(true_val / best_other, 0.0, 1.0)
    elseif true_val < 0
        # Negative objectives (e.g. log_L2_norm): more-negative competitor → higher score
        best_other = minimum(other_values)  # most negative competitor
        best_other >= 0 && return 0.0
        return clamp(best_other / true_val, 0.0, 1.0)
    else
        # true_val == 0: any finite other minimum gives some depth ratio
        return any(v != 0 for v in other_values) ? 0.5 : 0.0
    end
end

# ============================================================================
# Conditioning (Hessian condition number at p_true)
# ============================================================================

"""
    _compute_conditioning(grid_values, grid_axes, p_true) -> Float64

Approximate log10 condition number of the Hessian at the grid point nearest to p_true.
Uses 2nd-order central finite differences on the grid.
Returns 0.0 if the point is degenerate or computation fails.
"""
function _compute_conditioning(
    grid_values::Array{Float64},
    grid_axes::Vector{Vector{Float64}},
    p_true::Vector{Float64},
)::Float64
    ndim = length(grid_axes)
    dims = size(grid_values)

    # Find nearest grid point to p_true
    best_idx = _nearest_grid_index(grid_axes, p_true)

    # Shift to interior if on boundary (need ±1 neighbors for central differences)
    idx_arr = [best_idx[d] for d in 1:ndim]
    for d in 1:ndim
        idx_arr[d] = clamp(idx_arr[d], 2, dims[d] - 1)
    end
    ci = CartesianIndex(Tuple(idx_arr))

    !isfinite(grid_values[ci]) && return 0.0

    # Build approximate Hessian via central finite differences
    H = zeros(ndim, ndim)

    for d in 1:ndim
        h_d = grid_axes[d][2] - grid_axes[d][1]
        h_d == 0 && return 0.0

        # Diagonal: H[d,d] = (f[i+1] - 2f[i] + f[i-1]) / h_d^2
        idx_plus = CartesianIndex(ntuple(k -> k == d ? ci[k] + 1 : ci[k], ndim))
        idx_minus = CartesianIndex(ntuple(k -> k == d ? ci[k] - 1 : ci[k], ndim))

        fp = grid_values[idx_plus]
        fm = grid_values[idx_minus]
        fc = grid_values[ci]
        (!isfinite(fp) || !isfinite(fm) || !isfinite(fc)) && return 0.0

        H[d, d] = (fp - 2fc + fm) / (h_d^2)
    end

    for d1 in 1:ndim
        h1 = grid_axes[d1][2] - grid_axes[d1][1]
        for d2 in (d1+1):ndim
            h2 = grid_axes[d2][2] - grid_axes[d2][1]

            # Off-diagonal: (f[+1,+1] - f[+1,-1] - f[-1,+1] + f[-1,-1]) / (4·h1·h2)
            idx_pp = CartesianIndex(
                ntuple(k -> k == d1 ? ci[k] + 1 : k == d2 ? ci[k] + 1 : ci[k], ndim),
            )
            idx_pm = CartesianIndex(
                ntuple(k -> k == d1 ? ci[k] + 1 : k == d2 ? ci[k] - 1 : ci[k], ndim),
            )
            idx_mp = CartesianIndex(
                ntuple(k -> k == d1 ? ci[k] - 1 : k == d2 ? ci[k] + 1 : ci[k], ndim),
            )
            idx_mm = CartesianIndex(
                ntuple(k -> k == d1 ? ci[k] - 1 : k == d2 ? ci[k] - 1 : ci[k], ndim),
            )

            fpp = grid_values[idx_pp]
            fpm = grid_values[idx_pm]
            fmp = grid_values[idx_mp]
            fmm = grid_values[idx_mm]

            (!isfinite(fpp) || !isfinite(fpm) || !isfinite(fmp) || !isfinite(fmm)) &&
                return 0.0

            H[d1, d2] = (fpp - fpm - fmp + fmm) / (4 * h1 * h2)
            H[d2, d1] = H[d1, d2]
        end
    end

    # Compute condition number from eigenvalues
    evals = try
        eigvals(H)
    catch
        return 0.0
    end

    abs_evals = abs.(evals)
    any(e -> !isfinite(e) || e <= 0, abs_evals) && return 0.0

    return log10(maximum(abs_evals) / minimum(abs_evals))
end

"""
    _nearest_grid_index(grid_axes, point) -> CartesianIndex

Find the grid index nearest to `point` in parameter space.
"""
function _nearest_grid_index(grid_axes::Vector{Vector{Float64}}, point::Vector{Float64})
    ndim = length(grid_axes)
    idx_arr = Vector{Int}(undef, ndim)
    for d in 1:ndim
        _, idx_arr[d] = findmin(abs.(grid_axes[d] .- point[d]))
    end
    return CartesianIndex(Tuple(idx_arr))
end

# ============================================================================
# Ruggedness (second-order / first-order variation ratio)
# ============================================================================

"""
    _compute_ruggedness(grid_values::Array{Float64}) -> Float64

Ratio of second-order to first-order finite differences, averaged across dimensions.
High ratio = bumpy/rough surface; low ratio = smooth.
"""
function _compute_ruggedness(grid_values::Array{Float64})::Float64
    ndim = ndims(grid_values)
    dims = size(grid_values)
    ratios = Float64[]

    for d in 1:ndim
        first_sum = 0.0
        first_n = 0
        second_sum = 0.0
        second_n = 0

        for idx in CartesianIndices(grid_values)
            v0 = grid_values[idx]
            !isfinite(v0) && continue

            # First-order: |f[i+1] - f[i]|
            if idx[d] < dims[d]
                next_idx = CartesianIndex(ntuple(k -> k == d ? idx[k] + 1 : idx[k], ndim))
                v1 = grid_values[next_idx]
                if isfinite(v1)
                    first_sum += abs(v1 - v0)
                    first_n += 1
                end
            end

            # Second-order: |f[i+2] - 2f[i+1] + f[i]|
            if idx[d] + 2 <= dims[d]
                next1_idx = CartesianIndex(ntuple(k -> k == d ? idx[k] + 1 : idx[k], ndim))
                next2_idx = CartesianIndex(ntuple(k -> k == d ? idx[k] + 2 : idx[k], ndim))
                v1 = grid_values[next1_idx]
                v2 = grid_values[next2_idx]
                if isfinite(v1) && isfinite(v2)
                    second_sum += abs(v2 - 2v1 + v0)
                    second_n += 1
                end
            end
        end

        if first_n > 0 && first_sum > 0
            first_mean = first_sum / first_n
            second_mean = second_n > 0 ? second_sum / second_n : 0.0
            push!(ratios, second_mean / first_mean)
        end
    end

    isempty(ratios) && return 0.0
    return sum(ratios) / length(ratios)
end

# ============================================================================
# Plateau fraction
# ============================================================================

"""
    _compute_plateau_fraction(grid_values, minima_indices) -> Float64

Fraction of finite grid points with near-zero gradient that are NOT local minima.
High values indicate flat plateaus in the landscape.
"""
function _compute_plateau_fraction(
    grid_values::Array{Float64},
    minima_indices::Vector{<:CartesianIndex},
)::Float64
    ndim = ndims(grid_values)
    dims = size(grid_values)

    # Compute value range for threshold
    finite_vals = filter(isfinite, vec(grid_values))
    length(finite_vals) < 2 && return 0.0
    val_range = maximum(finite_vals) - minimum(finite_vals)
    val_range == 0 && return 0.0
    ε = 0.01 * val_range

    minima_set = Set(minima_indices)

    n_plateau = 0
    n_finite = 0

    for idx in CartesianIndices(grid_values)
        v = grid_values[idx]
        !isfinite(v) && continue
        n_finite += 1

        # Compute gradient magnitude via forward differences (use available neighbors)
        grad_sq = 0.0
        has_any_neighbor = false
        for d in 1:ndim
            if idx[d] < dims[d]
                next_idx = CartesianIndex(ntuple(k -> k == d ? idx[k] + 1 : idx[k], ndim))
                vn = grid_values[next_idx]
                if isfinite(vn)
                    grad_sq += (vn - v)^2
                    has_any_neighbor = true
                end
            elseif idx[d] > 1
                prev_idx = CartesianIndex(ntuple(k -> k == d ? idx[k] - 1 : idx[k], ndim))
                vp = grid_values[prev_idx]
                if isfinite(vp)
                    grad_sq += (v - vp)^2
                    has_any_neighbor = true
                end
            end
        end

        !has_any_neighbor && continue

        grad_mag = sqrt(grad_sq)
        if grad_mag < ε && !(idx in minima_set)
            n_plateau += 1
        end
    end

    n_finite == 0 && return 0.0
    return n_plateau / n_finite
end

# ============================================================================
# Composite interestingness score
# ============================================================================

"""
    interestingness_score(gs::GridScoreResult) -> Float64

Compute a composite interestingness score from grid-based metrics.
Higher values indicate more interesting/challenging landscapes for
optimization benchmarking.

!!! warning "What this score actually ranks (bead cbyn.1)"
    It ranks ONE of the two difficulty modes in the corpus — *deep wells in a
    flat sea* — and is close to blind to the other, *curved multimodal*. Do not
    read it as a general "hard vs boring" axis, and do not select a benchmark
    set by top-N of it. Use [`difficulty_profile`](@ref) to stratify.

    Measured over the 73-entry corpus (`experiments/sandbox/results/grid_scores.json`):

    - `curvature_score` carries the LARGEST weight (0.13) but supplies the
      second-SMALLEST share of the score's spread (6.3%), because its realized
      distribution is crushed: `clamp(min_variation / 0.1, 0, 1)` divides by 0.1
      when the corpus median `min_variation` is 0.0089, so **72 of 73 entries
      can never saturate it** and the median contribution is 0.012 of a 0.13
      budget.
    - The `Curvature ... penalizes flat/degenerate landscapes` line below is
      wrong as stated: a term with a POSITIVE weight can only fail to reward,
      never penalize. Flatness is in fact rewarded outright by
      `plateau_fraction` (+0.08).
    - Empirically `corr(curvature_score, interestingness) = -0.534`. This is NOT
      a weighting bug and cannot be reweighted away — curvature is negatively
      correlated with *every* other term (−0.15 to −0.66) and with their
      aggregate at −0.605. Dropping `plateau_norm` and recalibrating curvature
      together still leave it at −0.364. The nine terms span two opposed
      factors, so no all-positive linear scalar can rank both.
    - `dynamic_range` and `plateau_fraction` correlate at +0.91 with each other,
      i.e. 0.18 of the weight budget is spent twice on one property.

    Left numerically unchanged on purpose: 73 scores are persisted and callers
    compare against them. The defect is in what it CLAIMS to measure, not in its
    arithmetic.

# Scoring components (all normalized to [0, 1] range):
- **Multimodality** (weight 0.20): more local minima = more interesting
- **Dynamic range** (weight 0.10): wider value spread = more structure
- **Basin narrowness** (weight 0.12): smaller p_true basin = harder.
  Uses slope 4 so credit drops to zero once the p_true basin covers ≥25% of
  the domain — a wide basin is not a "narrow" basin.
- **Deceptiveness** (weight 0.10): more deceptive minima = harder
- **Curvature** (weight 0.13): penalizes flat/degenerate landscapes
- **Basin depth ratio** (weight 0.12): competing minima equally deep = harder.
  Multiplied by a multimodal-context factor `min(1, (n_local_minima-1)/4)` so
  that a near-unity ratio against 0–1 competitors (degenerate case) earns no
  credit — the metric is only meaningful when there are ≥5 competing minima.
- **Conditioning** (weight 0.10): ill-conditioned valley = harder
- **Ruggedness** (weight 0.05): bumpy surface = harder
- **Plateau fraction** (weight 0.08): flat plateaus = harder

Returns a value in [0, 1] where 1 is maximally interesting.
"""
function interestingness_score(gs::GridScoreResult)::Float64
    # Multimodality: saturates around 15 minima
    multi = min(1.0, (gs.n_local_minima - 1) / 14.0)
    multi = max(0.0, multi)

    # Dynamic range: saturates around 6 orders of magnitude
    dyn = min(1.0, gs.dynamic_range / 6.0)

    # Basin narrowness: slope 4 so basin_fraction ≥ 0.25 scores 0.
    # Previous formula `1 - basin_fraction` had a ~0.9 floor on typical
    # ODE landscapes (median basin_fraction ≈ 0.07), inflating scores by
    # ~0.11 across the board.
    basin = max(0.0, 1.0 - 4.0 * gs.basin_fraction)

    # Deceptiveness: saturates around 5 deceptive minima
    decept = min(1.0, gs.n_deceptive / 5.0)

    # Curvature: penalizes flat/degenerate landscapes (already [0,1])
    curv = gs.curvature_score

    # Basin depth ratio: weighted by multimodal context so ratios against
    # 0–1 competitors (degenerate near-unity values) earn no credit. The
    # raw ratio saturates at 1.0 on most ODE landscapes and was the second
    # structural floor in the pre-rebalance formula.
    depth_context = min(1.0, max(0, gs.n_local_minima - 1) / 4.0)
    depth = gs.basin_depth_ratio * depth_context

    # Conditioning: saturates at condition number 10^4
    cond_norm = min(1.0, gs.conditioning / 4.0)

    # Ruggedness: saturates at ratio 1.0
    rugged_norm = min(1.0, gs.ruggedness / 1.0)

    # Plateau fraction: saturates at 70% plateau
    plateau_norm = min(1.0, gs.plateau_fraction / 0.7)

    return 0.20 * multi +
           0.10 * dyn +
           0.12 * basin +
           0.10 * decept +
           0.13 * curv +
           0.12 * depth +
           0.10 * cond_norm +
           0.05 * rugged_norm +
           0.08 * plateau_norm
end

"Corpus 75th percentile of min(directional_variation); see [`structure_score`](@ref)."
const CURVATURE_SCALE = 0.02

"""
    structure_score(gs::GridScoreResult) -> Float64

Score the difficulty axis that [`interestingness_score`](@ref) is blind to:
**curved multimodality** — many local minima on a genuinely curved surface,
rather than narrow wells punched through a flat plateau.

Deliberately built from the terms that oppose the interestingness composite
(`corr(curvature_score, interestingness) = -0.534`), so the two scores are
near-independent by construction and a benchmark set can be stratified on both.

Curvature is recalibrated against the corpus rather than the arbitrary 0.1 of
`curvature_score`: dividing by 0.1 leaves 72 of 73 entries unable to saturate.
`CURVATURE_SCALE = 0.02` sits near the corpus 75th percentile of
`min(directional_variation)`, so the term spans its range instead of hugging 0.

Returns a value in `[0, 1]`. Higher means more curved-multimodal.
"""
function structure_score(gs::GridScoreResult)::Float64
    # Curvature on a corpus-realistic scale (see CURVATURE_SCALE).
    curv = clamp(gs.curvature_score * 0.1 / CURVATURE_SCALE, 0.0, 1.0)

    # Multimodality, same saturation as the interestingness composite.
    multi = clamp((gs.n_local_minima - 1) / 14.0, 0.0, 1.0)

    # Flatness counts AGAINST this axis — the term interestingness_score cannot
    # express, because there every weight is positive.
    not_flat = 1.0 - clamp(gs.plateau_fraction / 0.7, 0.0, 1.0)

    # Deceptive minima on a curved surface are what makes this mode hard.
    decept = clamp(gs.n_deceptive / 5.0, 0.0, 1.0)

    return 0.40 * curv + 0.30 * multi + 0.20 * not_flat + 0.10 * decept
end

"""
    difficulty_profile(gs::GridScoreResult) -> NamedTuple

Both difficulty axes plus a coarse mode label, for stratifying a benchmark set
instead of cutting top-N from a single scalar (bead cbyn.1, BENCH-4/t25x).

Returns `(; resolution, structure, mode)`:

- `resolution` — [`interestingness_score`](@ref). Deep narrow wells in a flat
  sea. Stresses polynomial RESOLUTION.
- `structure` — [`structure_score`](@ref). Dense multimodality on a curved
  surface. Stresses ROOT COUNT and separation.
- `mode` — `:deep_well`, `:curved`, `:both` or `:easy`, by comparing each axis
  against 0.5.

Both modes are genuinely hard; they are hard in different ways, and a set drawn
from one axis alone will systematically miss the other.
"""
function difficulty_profile(gs::GridScoreResult)
    r = interestingness_score(gs)
    s = structure_score(gs)
    mode = if r >= 0.5 && s >= 0.5
        :both
    elseif r >= 0.5
        :deep_well
    elseif s >= 0.5
        :curved
    else
        :easy
    end
    return (; resolution = r, structure = s, mode = mode)
end

# ============================================================================
# Batch scoring over screening results
# ============================================================================

"""
    score_top_candidates(screening_result, model, outputs, ic, bounds, time_interval;
                         top_n=20, kwargs...) -> Vector{Union{GridScoreResult, Nothing}}

Run `score_landscape_grid` on the top-N candidates from a `ScreeningResult`.
Returns one `GridScoreResult` per candidate (or `nothing` if evaluation fails).

# Arguments
- `screening_result::ScreeningResult`: output from `screen_and_probe`
- `model`, `outputs`, `ic`: model specification
- `bounds`: parameter domain
- `time_interval`: ODE time span
- `top_n::Int=20`: how many top-ranked candidates to score

All other keyword arguments are forwarded to `score_landscape_grid`.
"""
function score_top_candidates(
    screening_result::ScreeningResult,
    model,
    outputs,
    ic::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}},
    time_interval;
    top_n::Int = 20,
    kwargs...,
)::Vector{Union{GridScoreResult,Nothing}}
    ranked = screening_result.ranking.ranked_indices
    valid = screening_result.sweep.valid
    n = min(top_n, length(ranked))

    results = Vector{Union{GridScoreResult,Nothing}}(nothing, n)

    for i in 1:n
        candidate_idx = ranked[i]
        p_true = valid[candidate_idx]

        try
            results[i] = score_landscape_grid(
                model,
                outputs,
                ic,
                p_true,
                bounds,
                time_interval;
                kwargs...,
            )
        catch e
            @warn "Grid scoring failed for candidate $i (index $candidate_idx)" exception =
                (e, catch_backtrace())
            results[i] = nothing
        end
    end

    return results
end

# ============================================================================
# Display
# ============================================================================

"""
    print_grid_score(gs::GridScoreResult; io::IO=stdout)

Print a human-readable summary of grid-based landscape scoring results.
"""
function print_grid_score(gs::GridScoreResult; io::IO = stdout)
    println(io, "Grid-Based Landscape Score")
    println(io, "="^50)
    println(io, "  Grid size:       $(join(gs.grid_size, " x ")) = $(gs.n_total) points")
    println(
        io,
        "  Finite evals:    $(gs.n_finite) / $(gs.n_total) ($(round(100*gs.n_finite/gs.n_total; digits=1))%)",
    )
    println(io, "  Eval time:       $(round(gs.evaluation_time_ms; digits=1)) ms")
    println(io)
    println(io, "  Local minima:    $(gs.n_local_minima)")
    if !isempty(gs.local_minima_values)
        best = minimum(gs.local_minima_values)
        worst = maximum(gs.local_minima_values)
        println(io, "    Best:          $(round(best; sigdigits=4))")
        println(io, "    Worst:         $(round(worst; sigdigits=4))")
    end
    println(
        io,
        "  Dynamic range:   $(round(gs.dynamic_range; digits=2)) orders of magnitude",
    )
    println(
        io,
        "  Basin fraction:  $(round(gs.basin_fraction * 100; digits=1))% (p_true basin)",
    )
    println(io, "  Deceptive min:   $(gs.n_deceptive)")
    println(
        io,
        "  Dir. variation:  [$(join([round(v; digits=3) for v in gs.directional_variation], ", "))]",
    )
    println(
        io,
        "  Curvature score: $(round(gs.curvature_score; digits=3)) (1=curved, 0=flat)",
    )
    println(
        io,
        "  Depth ratio:     $(round(gs.basin_depth_ratio; digits=3)) (1=equally deep competitor)",
    )
    println(
        io,
        "  Conditioning:    $(round(gs.conditioning; digits=2)) (log10 condition number)",
    )
    println(
        io,
        "  Ruggedness:      $(round(gs.ruggedness; digits=3)) (2nd/1st order ratio)",
    )
    println(io, "  Plateau frac:    $(round(gs.plateau_fraction * 100; digits=1))%")
    println(io)
    println(io, "  Interestingness: $(round(interestingness_score(gs); digits=3))")
    println(io, "="^50)
end
