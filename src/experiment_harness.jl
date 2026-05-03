# Candidate-level parallelism for catalogue experiments (bead 1yt).
# Fan-out across entries; complements the inner threading from iyj
# (thread_evals for per-entry grid evaluation).

"""
    run_catalogue_experiments(entries, runner; parallel=false, max_workers=Threads.nthreads())

Apply `runner(entry)` to every entry in `entries` and return the collected
results in the same order as `entries`.

# Thread safety

When `parallel=true`, entries are distributed across chunked `Threads.@spawn`
tasks (same pattern as `thread_evals` in globtim's `MainGenerate` /
`adaptive_subdivision`). The caller-supplied `runner` **must be thread-safe**
— typically this means either (a) it builds any mutable per-entry state
inside the closure, or (b) if it reuses an `ODEProblem` / solver across
calls, the underlying APIs are re-entrant (e.g. `SciMLBase.remake` rather
than in-place parameter mutation).

# Arguments

- `entries::AbstractVector`: any vector of catalogue-like items. `runner`
  is called once per element.
- `runner`: unary callable. What it returns is opaque to this function;
  all returned values are collected into a vector.

# Keyword arguments

- `parallel::Bool = false`: serial (`map`) by default; set true to fan out.
- `max_workers::Int = Threads.nthreads()`: controls chunk granularity
  (same role as `Threads.nthreads()` in `thread_evals`). When ≤ 1 or when
  `entries` has ≤ 1 element, the serial path is used regardless of
  `parallel`.

# Example

```julia
entries = load_catalogue("data/catalogue.jsonl")
results = run_catalogue_experiments(
    entries,
    entry -> run_standard_experiment(create_objective(entry), build_bounds(entry));
    parallel = true,
)
```
"""
function run_catalogue_experiments(
    entries::AbstractVector,
    runner;
    parallel::Bool = false,
    max_workers::Int = Threads.nthreads(),
)
    n = length(entries)
    if !parallel || max_workers <= 1 || n <= 1
        return map(runner, entries)
    end

    first_result = runner(entries[1])
    results = Vector{typeof(first_result)}(undef, n)
    results[1] = first_result

    chunk_size = max(1, cld(n - 1, 4 * max_workers))
    @sync for chunk_start in 2:chunk_size:n
        chunk_end = min(chunk_start + chunk_size - 1, n)
        Threads.@spawn begin
            for i in chunk_start:chunk_end
                results[i] = runner(entries[i])
            end
        end
    end

    return results
end
