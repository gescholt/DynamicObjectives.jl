# Parallel Grid Evaluation
# Thread-safe parallel evaluation of objectives on N-dimensional Cartesian grids.
# Uses a factory pattern: each thread gets its own objective instance.
# Note: since sample_data now uses SciMLBase.remake (no shared mutable state),
# a single objective instance is inherently thread-safe. The factory pattern
# is retained for performance (avoids remake contention on shared problem).

"""
    evaluate_grid_threaded(make_obj, grid_axes; show_progress=true) -> Array{Float64}

Evaluate an objective on an N-dimensional Cartesian grid using multiple threads.

# Thread safety

`make_obj` is a **zero-argument factory function** that returns a fresh callable
objective. It is called once per thread to produce independent instances.
While `sample_data` now uses `SciMLBase.remake` (thread-safe), the factory
pattern is retained for performance to avoid contention on shared problem state.

# Arguments
- `make_obj::Function`: Factory `() -> obj` where `obj(p::Vector{Float64}) -> Float64`.
- `grid_axes::Vector{Vector{Float64}}`: One vector of sample coordinates per dimension.
  The Cartesian product defines the full grid (size = `prod(length.(grid_axes))`).

# Keyword Arguments
- `show_progress::Bool = true`: Print periodic progress updates to stdout.

# Returns
- `Array{Float64}` of shape `(length(grid_axes[1]), ..., length(grid_axes[N]))`.
  Failed evaluations are stored as `Inf` (consistent with `return_inf_on_error`).

# Example
```julia
entry = load_catalogue("fhn3d_catalogue.jsonl")[1]
make_obj = () -> create_objective(entry)
axes = [collect(range(lo, hi, length=20)) for (lo, hi) in entry.bounds]
values = evaluate_grid_threaded(make_obj, axes)
```
"""
function evaluate_grid_threaded(
    make_obj::Function,
    grid_axes::Vector{Vector{Float64}};
    show_progress::Bool = true
)
    dim = length(grid_axes)
    dim >= 1 || throw(ArgumentError("grid_axes must have at least 1 dimension"))
    sizes = Tuple(length(ax) for ax in grid_axes)
    n_total = prod(sizes)

    # Build one objective per outer-loop iteration. Each iteration runs on one
    # thread at a time, so there are no data races. We cannot use threadid()
    # because Julia 1.12's task migration means threadid() can exceed nthreads().
    n_outer = sizes[1]
    objectives = Vector{Any}(undef, n_outer)
    for i in 1:n_outer
        objectives[i] = make_obj()
    end

    values = fill(Inf, sizes)
    completed = Threads.Atomic{Int}(0)
    t_start = time()

    # Progress reporter (runs on a separate task)
    progress_task = if show_progress
        @async begin
            while true
                n_done = completed[]
                if n_done >= n_total
                    break
                end
                elapsed = time() - t_start
                rate = n_done / max(elapsed, 1e-6)
                remaining = (n_total - n_done) / max(rate, 1e-6)
                pct = 100.0 * n_done / n_total
                print("\r  Grid eval: $n_done/$n_total ($(@sprintf("%.1f", pct))%) — " *
                      "$(@sprintf("%.0f", rate)) evals/s — ETA $(@sprintf("%.0f", remaining))s   ")
                sleep(2.0)
            end
        end
    else
        nothing
    end

    # Parallelize over the first axis (outer loop), sequential inner axes.
    # This matches the proven pattern from landscape_explorer.jl:compute_distance_grid.
    Threads.@threads for i_outer in 1:n_outer
        obj = objectives[i_outer]
        for idx in CartesianIndices(sizes[2:end])
            # Build the full index tuple: (i_outer, idx[1], idx[2], ...)
            full_idx = CartesianIndex(i_outer, Tuple(idx)...)
            p = [grid_axes[d][full_idx[d]] for d in 1:dim]
            val = try
                result = obj(p)
                isfinite(result) ? result : Inf
            catch e
                if e isa InterruptException
                    rethrow(e)
                end
                Inf
            end
            values[full_idx] = val
            Threads.atomic_add!(completed, 1)
        end
    end

    # Stop progress reporter
    if show_progress
        completed[] = n_total  # signal completion
        sleep(0.1)  # let reporter print final state
        elapsed = time() - t_start
        rate = n_total / max(elapsed, 1e-6)
        n_finite = count(isfinite, values)
        println("\r  Grid eval: $n_total/$n_total (100.0%) — done in $(@sprintf("%.1f", elapsed))s " *
                "($(@sprintf("%.0f", rate)) evals/s) — $n_finite/$n_total finite values   ")
    end

    return values
end
