# Error Metrics
# Distance functions and error metric construction for parameter estimation

"""
L1 norm distance function (scaled by 100)

Arguments:
- Y_true: True values
- Y_test: Test values

Returns:
    100 * ||Y_true - Y_test||₁

Note: The 100x scaling converts the L1 norm to a percentage-like scale for
readability in experiment summaries. Not a mathematical normalization.
"""
const L1_NORM_DISPLAY_SCALE = 100
L1_norm(Y_true, Y_test) = L1_NORM_DISPLAY_SCALE * norm(Y_true - Y_test, 1)

"""
L2 norm distance function

Arguments:
- Y_true: True values
- Y_test: Test values

Returns:
    ||Y_true - Y_test||₂
"""
L2_norm(Y_true, Y_test) = norm(Y_true - Y_test, 2)

"""
Squared Euclidean distance — `Σᵢ (Y_true[i] - Y_test[i])²`.

Same minimizer as `L2_norm` (squaring is monotone on ℝ≥₀), but **C² at zero**
where `L2_norm` is Lipschitz/conical. Use this for parameter-estimation error
functions when the polynomial-approximation framework (Safey/Scholten/Trélat,
HAL hal-05160251) is intended to apply: Thm 1's Morse-class smoothness
assumption `c ≥ max(3, βn+1)` requires C² landscapes, which `L2_norm` violates
at the global minimum.

Empirically (radius_sweep_squared_l2_fhn3d_tight.md, 2026-05-27): switching
fhn3d_tight from `L2_norm` to `L2_squared` reduced ‖w_d - f‖_∞ by 60×–14000×
at the same degrees and restored `f_max ∝ r²` scaling (quadratic well)
from the Lipschitz `f_max ∝ r` of `L2_norm`.

Arguments:
- Y_true: True values
- Y_test: Test values

Returns:
    Σᵢ (Y_true[i] - Y_test[i])²
"""
L2_squared(Y_true, Y_test) = sum(abs2, Y_true .- Y_test)

"""
Log L2 norm distance function

Arguments:
- Y_true: True values
- Y_test: Test values

Returns:
    log₂(||Y_true - Y_test||₂ + ε)
"""
log_L2_norm(Y_true, Y_test) = log2(norm(Y_true - Y_test, 2) + eps(eltype(Y_true)))

"""
    make_error_distance(model, outputs, ic, p_true, time_interval,
                       numpoints, distance_function, aggregate_distances, add_noise_in_time_series;
                       return_inf_on_error, eval_timeout, solver, abstol, reltol)

Construct an error function comparing model predictions against reference data.

The returned closure uses `SciMLBase.remake` to set parameters, which:
- Supports ForwardDiff.Dual types for automatic differentiation
- Is thread-safe (each call creates a new problem instance)

Arguments:
- `model`: ModelingToolkit System
- `outputs`: Vector of measurement equations
- `ic`: Vector of initial conditions for the ODE system
- `p_true`: Vector of true parameter values
- `time_interval`: Time interval [start_time, end_time] for simulation
- `numpoints`: Number of time points to sample (default: 5)
- `distance_function`: Function to compute distance (default: L2_norm).
    The function should take two vectors (true and predicted) and return a scalar distance value.
- `aggregate_distances`: Function to aggregate multiple distance values (default: sum)
- `add_noise_in_time_series`: Function to add noise to time series (default: nothing = no noise)
- `return_inf_on_error`: If true (default), return Inf when ODE solving fails (for optimization compatibility).
    If false, throw detailed error messages (for debugging).
- `eval_timeout`: Maximum time in seconds for a single objective function evaluation (default: nothing = no timeout).
    If evaluation exceeds this time, returns Inf (when return_inf_on_error=true) or throws TimeoutError.
- `solver`: ODE solver algorithm (default: Vern9()). Options include Tsit5(), Vern7(), Vern9(),
    AutoTsit5(Rosenbrock23()), etc. Lower-order solvers like Tsit5() are faster at coarser tolerances.
- `abstol`: Absolute tolerance for ODE solver (default: 1e-10). Loosening to 1e-6 can give 2-5x speedup.
- `reltol`: Relative tolerance for ODE solver (default: 1e-10). Loosening to 1e-6 can give 2-5x speedup.

Returns:
    Function that computes error between predictions and reference data
"""
function make_error_distance(
    model::ModelingToolkit.AbstractSystem,
    outputs::Vector{ModelingToolkit.Equation},
    ic::AbstractVector{<:Real},
    p_true::Vector{T},
    time_interval,
    numpoints::Int = 5,
    distance_function = L2_norm,
    aggregate_distances = sum,
    add_noise_in_time_series = nothing;
    return_inf_on_error::Bool = true,  # For optimization compatibility
    eval_timeout::Union{Float64,Nothing} = nothing,  # Timeout per evaluation
    solver = Vern9(),
    abstol::Real = 1e-10,
    reltol::Real = 1e-10,
    uneven_sampling_times::AbstractVector{<:Real} = Float64[],  # Explicit sample times (overrides uniform spacing)
) where {T}
    @assert length(p_true) == length(ModelingToolkit.parameters(model)) "Parameter vector length mismatch"
    @assert length(ic) == length(ModelingToolkit.unknowns(model)) "Initial conditions length mismatch"
    @assert length(outputs) > 0 "At least one output variable must be specified"
    @assert numpoints > 0 "Number of points must be greater than zero"
    @assert time_interval[2] > time_interval[1] "End time must be greater than start time"

    # Determine sampling mode: explicit times vs uniform spacing
    uneven_sampling = !isempty(uneven_sampling_times)
    if uneven_sampling
        @assert length(uneven_sampling_times) == numpoints "uneven_sampling_times length ($(length(uneven_sampling_times))) must match numpoints ($numpoints)"
        @assert issorted(uneven_sampling_times) "uneven_sampling_times must be in ascending order"
    end

    # Generate reference solution once during function creation
    # Model is already completed by @mtkcompile, so pass directly to ODEProblem
    problem = ODEProblem(
        model,
        merge(
            Dict(ModelingToolkit.unknowns(model) .=> ic),
            Dict(ModelingToolkit.parameters(model) .=> p_true),
        ),
        time_interval,
    )

    data_sample_true = sample_data(
        problem,
        model,
        outputs,
        time_interval,
        p_true,
        ic,
        numpoints;
        uneven_sampling = uneven_sampling,
        uneven_sampling_times = uneven_sampling_times,
    )

    if add_noise_in_time_series === nothing
        add_noise_in_time_series = x -> x  # No noise function, return input as is
    end
    for (key, values) in data_sample_true
        if key == "t"
            continue  # Skip time array
        end
        data_sample_true[key] = add_noise_in_time_series(values)
    end

    # ── Build the integrator freelist for the fast path (vemc) ─────────────
    # A Channel of `init`'d (integrator, out_buffer) pairs + a single `setp`
    # closure replace the per-call `remake(problem; p=...)+solve(...)`
    # allocation. take!/put! per call (~100 ns) is noise against the ~15 μs
    # solve, and unlike `pool[Threads.threadid()]` it cannot alias two
    # migrating tasks onto one integrator. The closure branches on T2:
    # Float64 → pool (~15 μs/call), Dual/other → legacy path.
    # See experiments/sandbox/spike_integrator_reuse.jl for validation.
    sampling_times_vec = uneven_sampling ?
        collect(uneven_sampling_times) :
        collect(range(time_interval[1], time_interval[2], length = numpoints))
    integrators_pool, setp_fn = _build_integrator_pool(
        problem,
        model,
        solver,
        abstol,
        reltol,
        sampling_times_vec,
        outputs,
        numpoints,
    )

    function error_distance(p_test::AbstractVector{T2}) where {T2}

        use_pool = (integrators_pool !== nothing) && (T2 === Float64)

        # Timeout wrapper: use @async + timedwait if eval_timeout is set
        if eval_timeout !== nothing
            # Run computation in async task
            task = @async begin
                use_pool ? _compute_error_distance_pool(
                    p_test,
                    outputs,
                    data_sample_true,
                    distance_function,
                    aggregate_distances,
                    integrators_pool,
                    setp_fn,
                    numpoints,
                ) : _compute_error_distance(
                    p_test,
                    outputs,
                    time_interval,
                    numpoints,
                    data_sample_true,
                    problem,
                    model,
                    ic,
                    distance_function,
                    aggregate_distances,
                    solver,
                    abstol,
                    reltol,
                )
            end

            # Wait for task with timeout
            timed_result = timedwait(() -> istaskdone(task), eval_timeout)

            if timed_result == :timed_out
                # Kill the task (best effort - may not actually interrupt ODE solver)
                Base.throwto(task, InterruptException())

                if return_inf_on_error
                    # Use T2(Inf) to preserve ForwardDiff Dual types
                    return T2(Inf)
                else
                    error(
                        "Objective evaluation timed out after $(eval_timeout)s for p_test=$p_test",
                    )
                end
            else
                # Task completed - check for errors
                try
                    return fetch(task)
                catch e
                    if isa(e, InterruptException)
                        rethrow(e)
                    end
                    if return_inf_on_error
                        return T2(Inf)
                    else
                        rethrow(e)
                    end
                end
            end
        end

        # No timeout - original behavior
        # Wrap in try-catch when return_inf_on_error is true (for optimization)
        if return_inf_on_error
            try
                return use_pool ? _compute_error_distance_pool(
                    p_test,
                    outputs,
                    data_sample_true,
                    distance_function,
                    aggregate_distances,
                    integrators_pool,
                    setp_fn,
                    numpoints,
                ) : _compute_error_distance(
                    p_test,
                    outputs,
                    time_interval,
                    numpoints,
                    data_sample_true,
                    problem,
                    model,
                    ic,
                    distance_function,
                    aggregate_distances,
                    solver,
                    abstol,
                    reltol,
                )
            catch e
                if isa(e, InterruptException)
                    rethrow(e)
                end
                # Return typed Inf to preserve ForwardDiff Dual types
                return T2(Inf)
            end
        else
            # Let errors propagate for debugging
            return use_pool ? _compute_error_distance_pool(
                p_test,
                outputs,
                data_sample_true,
                distance_function,
                aggregate_distances,
                integrators_pool,
                setp_fn,
                numpoints,
            ) : _compute_error_distance(
                p_test,
                outputs,
                time_interval,
                numpoints,
                data_sample_true,
                problem,
                model,
                ic,
                distance_function,
                aggregate_distances,
                solver,
                abstol,
                reltol,
            )
        end
    end

    # Inner function with actual computation and error checking
    function _compute_error_distance(
        p_test,
        measured_data,
        time_interval,
        datasize,
        data_sample_true,
        problem,
        model,
        ic,
        distance_function,
        aggregate_distances,
        solver,
        abstol,
        reltol,
    )
        # Check datasize consistency
        if datasize != length(data_sample_true["t"])
            error(
                "Datasize mismatch: requested $datasize but reference data has $(length(data_sample_true["t"])) points",
            )
        end

        # Sample data with test parameters using specified solver and tolerances
        # Uses same sampling times (uniform or explicit) as the reference data
        data_sample_test = try
            sample_data(
                problem,
                model,
                measured_data,
                time_interval,
                p_test,
                ic,
                datasize;
                uneven_sampling = uneven_sampling,
                uneven_sampling_times = uneven_sampling_times,
                solver = solver,
                abstol = abstol,
                reltol = reltol,
            )
        catch e
            if isa(e, InterruptException)
                rethrow(e)
            end
            # Provide context about which parameters caused the sampling failure
            error(
                "ODE sampling failed for parameters p_test=$p_test: $(sprint(showerror, e))",
            )
        end

        if isempty(data_sample_test)
            error("sample_data returned empty dataset for parameters p_test=$p_test")
        end

        # Validate data quality
        for (key, values) in data_sample_test
            if key == "t"
                continue  # Skip time array
            end
            if any(isnan.(values))
                error(
                    "NaN detected in sampled data for variable '$key' with parameters p_test=$p_test",
                )
            end
            if any(isinf.(values))
                error(
                    "Inf detected in sampled data for variable '$key' with parameters p_test=$p_test",
                )
            end
            if any(isnan.(data_sample_true[key]))
                error("NaN detected in reference data for variable '$key'")
            end
            if any(isinf.(data_sample_true[key]))
                error("Inf detected in reference data for variable '$key'")
            end
        end

        return aggregate_distances([
            distance_function(data_sample_true[key], data_sample_test[key]) for
            key in keys(data_sample_true) if key != "t"
        ])
    end

    # Fast path: borrow an (integrator, out_buffer) pair from the freelist.
    # Mirrors _compute_error_distance's validation + aggregation, but the inner
    # `sample_data!` does setp/reinit!/solve! against a pre-built integrator
    # instead of remake+solve. Float64-only — gated by the outer closure.
    #
    # take!/put! rather than `pool[Threads.threadid()]`: under task migration
    # (Julia ≥1.7) a task can move threads between the threadid() read and the
    # solve, aliasing two tasks onto one integrator — which corrupts the solve
    # silently (the return_inf_on_error closure absorbs it as Inf). The
    # freelist hands each in-flight call exclusive ownership of a pair.
    function _compute_error_distance_pool(
        p_test::AbstractVector{Float64},
        measured_data,
        data_sample_true,
        distance_function,
        aggregate_distances,
        integrators_pool,
        setp_fn,
        datasize,
    )
        if datasize != length(data_sample_true["t"])
            error(
                "Datasize mismatch: requested $datasize but reference data has $(length(data_sample_true["t"])) points",
            )
        end

        integrator, out = take!(integrators_pool)
        try
            try
                sample_data!(out, integrator, p_test, setp_fn, measured_data)
            catch e
                if isa(e, InterruptException)
                    rethrow(e)
                end
                error(
                    "ODE sampling (integrator-pool path) failed for parameters p_test=$p_test: $(sprint(showerror, e))",
                )
            end

            for (key, values) in out
                if key == "t"
                    continue
                end
                if any(isnan, values)
                    error(
                        "NaN detected in sampled data for variable '$key' with parameters p_test=$p_test",
                    )
                end
                if any(isinf, values)
                    error(
                        "Inf detected in sampled data for variable '$key' with parameters p_test=$p_test",
                    )
                end
                if any(isnan, data_sample_true[key])
                    error("NaN detected in reference data for variable '$key'")
                end
                if any(isinf, data_sample_true[key])
                    error("Inf detected in reference data for variable '$key'")
                end
            end

            return aggregate_distances([
                distance_function(data_sample_true[key], out[key]) for
                key in keys(data_sample_true) if key != "t"
            ])
        finally
            # A pair borrowed by a failed call goes back too: sample_data!'s
            # reinit! resets the integrator state on next use.
            put!(integrators_pool, (integrator, out))
        end
    end

    return error_distance
end

# ─────────────────────────────────────────────────────────────────────────────
# Integrator-pool builder (vemc — per-thread reuse of init/reinit!/solve!)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _build_integrator_pool(problem, model, solver, abstol, reltol, saveat, outputs, numpoints)
        -> (pool::Channel, setp_fn)

Build the freelist used by the `make_error_distance` Float64 fast path: a
`Channel` of `(integrator, out_buffer)` pairs — each an `init`'d integrator
sharing the same `problem` template plus a pre-allocated output buffer
(`OrderedDict{Any,Vector{Float64}}`) — and a single `setp` closure for
in-place parameter mutation. Callers `take!` a pair for exclusive use and
`put!` it back; a Channel rather than `pool[Threads.threadid()]` indexing
because task migration can move a task between the threadid() read and the
solve, aliasing two tasks onto one integrator (silent-Inf corruption under
`thread_evals`).

If construction fails (e.g., the model is incompatible with `setp` or `init`),
raises — callers should treat that as a real configuration error, not a
silent fallback. Per-call AD compatibility (Dual eltype) is handled separately
at the call site.
"""
_solver_pool_safe(solver) = nameof(typeof(solver)) !== :CompositeAlgorithm

function _build_integrator_pool(
    problem,
    model,
    solver,
    abstol,
    reltol,
    saveat::AbstractVector{<:Real},
    outputs::Vector{ModelingToolkit.Equation},
    numpoints::Int,
)
    # Bail out for composite / auto-switching solvers (e.g. AutoTsit5(Rosenbrock23())).
    # CompositeAlgorithm carries stiffness-detector state in AutoSwitch that
    # `reinit!` does not reset, causing ~1e-5 drift across sequential solves
    # (diagnostic: experiments/sandbox/diagnose_vemc_drift.jl). Returning
    # (nothing, nothing) signals the caller to use the legacy remake+solve
    # path, which is allocation-heavy but correct.
    if !_solver_pool_safe(solver)
        return nothing, nothing
    end

    # Freelist capacity: Threads.maxthreadid() (default + interactive pools)
    # bounds how many tasks the schedulers can run simultaneously, so take!
    # never blocks under the chunked-@spawn fan-out; if it ever did (more
    # in-flight callers than pairs), the call waits for a free pair instead
    # of aliasing state.
    n = Threads.maxthreadid()
    proto = SciMLBase.init(
        problem,
        solver;
        saveat = saveat,
        abstol = abstol,
        reltol = reltol,
        verbose = false,
        maxiters = 1000000,
    )

    setp_fn = SciMLBase.setp(problem, ModelingToolkit.parameters(model))

    saveat_vec = collect(Float64, saveat)
    make_buffer() =
        let buf = DataStructures.OrderedDict{Any,Vector{Float64}}()
            for v in outputs
                buf[Num(v.lhs)] = Vector{Float64}(undef, numpoints)
            end
            buf["t"] = copy(saveat_vec)
            buf
        end

    pool = Channel{Tuple{typeof(proto),DataStructures.OrderedDict{Any,Vector{Float64}}}}(n)
    put!(pool, (proto, make_buffer()))
    for _ in 2:n
        integ = SciMLBase.init(
            problem,
            solver;
            saveat = saveat,
            abstol = abstol,
            reltol = reltol,
            verbose = false,
            maxiters = 1000000,
        )
        put!(pool, (integ, make_buffer()))
    end

    return pool, setp_fn
end

#==============================================================================#
#                      TOLERANT OBJECTIVE WRAPPER                               #
#==============================================================================#

"""
    TolerantObjective

A mutable wrapper around an ODE objective function that allows runtime adjustment
of solver tolerances and solver algorithm. This enables adaptive subdivision to use
coarse tolerances during early phases and tight tolerances for final accuracy.

The wrapper is callable, so it satisfies the `f::Function` interface expected by
globtim's `adaptive_refine`, `two_phase_refine`, and `MainGenerate`.

# Usage

```julia
# Create base objective
base_obj = make_error_distance(model, outputs, ic, p_true, time_interval, numpoints,
                                L2_norm; return_inf_on_error=true)

# Wrap with mutable tolerances
tol_obj = TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
                            abstol=1e-10, reltol=1e-10)

# Use like a regular function
value = tol_obj([1.0, 0.5])

# Switch to coarse mode for subdivision Phase 1
set_tolerance!(tol_obj, 1e-6)

# Switch to tight mode for Phase 2
set_tolerance!(tol_obj, 1e-10)

# Change solver
set_solver!(tol_obj, Tsit5())
```
"""
mutable struct TolerantObjective
    _make_args::NamedTuple  # Frozen args for rebuilding: model, outputs, ic, p_true, etc.
    _error_func::Function   # Current callable (rebuilt when tolerances change)
    solver::Any             # Current ODE solver
    abstol::Float64         # Current absolute tolerance
    reltol::Float64         # Current relative tolerance
end

"""
    TolerantObjective(model, outputs, ic, p_true, time_interval, numpoints;
                       distance_function=L2_norm, aggregate_distances=sum,
                      add_noise_in_time_series=nothing,
                      return_inf_on_error=true, eval_timeout=nothing,
                      solver=Vern9(), abstol=1e-10, reltol=1e-10)

Construct a TolerantObjective with mutable solver settings.

All positional and keyword arguments match `make_error_distance` — the difference is that
`solver`, `abstol`, and `reltol` can be changed after construction via `set_tolerance!`
and `set_solver!`.
"""
function TolerantObjective(
    model::ModelingToolkit.AbstractSystem,
    outputs::Vector{ModelingToolkit.Equation},
    ic::AbstractVector{<:Real},
    p_true::Vector{T},
    time_interval,
    numpoints::Int = 5,
    distance_function = L2_norm,
    aggregate_distances = sum,
    add_noise_in_time_series = nothing;
    return_inf_on_error::Bool = true,
    eval_timeout::Union{Float64,Nothing} = nothing,
    solver = Vern9(),
    abstol::Real = 1e-10,
    reltol::Real = 1e-10,
    uneven_sampling_times::AbstractVector{<:Real} = Float64[],
) where {T}
    make_args = (
        model = model,
        outputs = outputs,
        ic = ic,
        p_true = p_true,
        time_interval = time_interval,
        numpoints = numpoints,
        distance_function = distance_function,
        aggregate_distances = aggregate_distances,
        add_noise_in_time_series = add_noise_in_time_series,
        return_inf_on_error = return_inf_on_error,
        eval_timeout = eval_timeout,
        uneven_sampling_times = uneven_sampling_times,
    )

    error_func = make_error_distance(
        model,
        outputs,
        ic,
        p_true,
        time_interval,
        numpoints,
        distance_function,
        aggregate_distances,
        add_noise_in_time_series;
        return_inf_on_error = return_inf_on_error,
        eval_timeout = eval_timeout,
        solver = solver,
        abstol = abstol,
        reltol = reltol,
        uneven_sampling_times = uneven_sampling_times,
    )

    return TolerantObjective(make_args, error_func, solver, abstol, reltol)
end

# Make TolerantObjective callable — this is the key interface
# For 1D ODE models, globtim passes a scalar; the inner objective expects a vector.
# Wrap scalars into a static 1-element vector (allocation-free via StaticArrays).
function (obj::TolerantObjective)(x)
    x_vec = x isa Number ? SVector{1}(x) : x
    return obj._error_func(x_vec)
end

"""
    set_tolerance!(obj::TolerantObjective, tol::Float64)

Set both abstol and reltol to the same value and rebuild the inner objective.
"""
function set_tolerance!(obj::TolerantObjective, tol::Float64)
    set_tolerance!(obj, tol, tol)
end

"""
    set_tolerance!(obj::TolerantObjective, abstol::Float64, reltol::Float64)

Set abstol and reltol independently and rebuild the inner objective.
"""
function set_tolerance!(obj::TolerantObjective, abstol::Float64, reltol::Float64)
    if obj.abstol == abstol && obj.reltol == reltol
        return obj  # No change needed
    end
    obj.abstol = abstol
    obj.reltol = reltol
    _rebuild_error_func!(obj)
    return obj
end

"""
    set_solver!(obj::TolerantObjective, solver)

Change the ODE solver algorithm and rebuild the inner objective.
"""
function set_solver!(obj::TolerantObjective, solver)
    obj.solver = solver
    _rebuild_error_func!(obj)
    return obj
end

"""
    get_tolerance(obj::TolerantObjective) -> (abstol, reltol)

Return current tolerance settings.
"""
get_tolerance(obj::TolerantObjective) = (obj.abstol, obj.reltol)

"""
    get_solver(obj::TolerantObjective) -> solver

Return current solver.
"""
get_solver(obj::TolerantObjective) = obj.solver

# Internal: rebuild the error function closure with current settings
function _rebuild_error_func!(obj::TolerantObjective)
    a = obj._make_args
    obj._error_func = make_error_distance(
        a.model,
        a.outputs,
        a.ic,
        a.p_true,
        a.time_interval,
        a.numpoints,
        a.distance_function,
        a.aggregate_distances,
        a.add_noise_in_time_series;
        return_inf_on_error = a.return_inf_on_error,
        eval_timeout = a.eval_timeout,
        solver = obj.solver,
        abstol = obj.abstol,
        reltol = obj.reltol,
        uneven_sampling_times = a.uneven_sampling_times,
    )
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Named aggregation strategies (bead 0iq)
# ─────────────────────────────────────────────────────────────────────────────

"""
    AGGREGATION_STRATEGIES :: Dict{Symbol, Function}

Registry of named strategies for combining per-output distance vectors into a
scalar objective. Pass `AGGREGATION_STRATEGIES[name]` (or use `resolve_aggregation(name)`)
as the `aggregate_distances` positional argument to `make_error_distance`.

Current strategies:
- `:sum`     — sum of per-output distances (default of `make_error_distance`)
- `:mean`    — arithmetic mean
- `:maximum` — worst-case output error
- `:minimum` — best-case output error
- `:first`   — first output only (legacy behaviour prior to bead pam)
- `:rms`     — root-mean-square of per-output distances
"""
const AGGREGATION_STRATEGIES = Dict{Symbol,Function}(
    :sum => sum,
    :mean => mean,
    :maximum => maximum,
    :minimum => minimum,
    :first => first,
    :rms => ds -> sqrt(mean(abs2, ds)),
)

"""
    resolve_aggregation(name::Symbol) -> Function

Look up a named aggregation strategy in `AGGREGATION_STRATEGIES`. Throws
`ArgumentError` listing valid names if the symbol is unknown. Useful when
catalogue entries or TOML configs persist strategy choices as symbols rather
than function references.
"""
function resolve_aggregation(name::Symbol)
    haskey(AGGREGATION_STRATEGIES, name) || throw(
        ArgumentError(
            "Unknown aggregation strategy :$name. " *
            "Valid options: $(sort!(collect(keys(AGGREGATION_STRATEGIES))))",
        ),
    )
    return AGGREGATION_STRATEGIES[name]
end

# ─────────────────────────────────────────────────────────────────────────────
# Partial observability helper (bead dds)
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_partial_observability_distance(model, outputs, output_indices, ic, p_true,
                                         time_interval,
                                         numpoints=5, distance_function=L2_norm,
                                         aggregate_distances=sum,
                                         add_noise_in_time_series=nothing;
                                         kwargs...)

Thin wrapper around `make_error_distance` that selects a subset of the model's
measurement equations via `output_indices`. Simulates partial observability —
cases where only some state variables can be measured in practice — so the
caller can study how identifiability changes under reduced output sets.

The same effect is achieved by passing `outputs[output_indices]` directly to
`make_error_distance`; this helper exists mainly to document the pattern and
validate the index range.

Arguments match `make_error_distance` except for `output_indices`, a
collection of integer indices into `outputs`. Throws `ArgumentError` if the
indices are empty or out of range.
"""
function make_partial_observability_distance(
    model::ModelingToolkit.AbstractSystem,
    outputs::Vector{<:ModelingToolkit.Equation},
    output_indices::AbstractVector{<:Integer},
    ic::AbstractVector{<:Real},
    p_true::Vector{T},
    time_interval,
    args...;
    kwargs...,
) where {T}
    isempty(output_indices) && throw(ArgumentError("output_indices must not be empty"))
    (minimum(output_indices) >= 1 && maximum(output_indices) <= length(outputs)) || throw(
        ArgumentError(
            "output_indices $(collect(output_indices)) out of range 1:$(length(outputs))",
        ),
    )
    subset = outputs[output_indices]
    return make_error_distance(model, subset, ic, p_true, time_interval, args...; kwargs...)
end
