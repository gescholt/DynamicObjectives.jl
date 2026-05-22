# Data Generation
# Functions for generating synthetic time series data from ODE systems

"""
    sample_data(problem, model, measured_data, time_interval, p_true, u0, num_points; kwargs...)

Generate synthetic time series data from an ODE system with specified parameters.

Uses `SciMLBase.remake` to set parameters, which:
- Supports ForwardDiff.Dual types for automatic differentiation
- Is thread-safe (creates a new problem instead of mutating shared state)

Arguments:
- `problem`: ODEProblem template (not mutated — remade with new parameters)
- `model`: ModelingToolkit System representing the differential equations
- `measured_data`: Vector of measurement equations
- `time_interval`: [start_time, end_time] for simulation
- `p_true`: Parameter values (Vector or SVector — may contain ForwardDiff.Dual types)
- `u0`: Vector of initial conditions
- `num_points`: Number of time points to sample

Optional kwargs:
- `uneven_sampling`: Boolean for non-uniform time sampling
- `uneven_sampling_times`: Vector of specific sampling times
- `solver`: ODE solver (default: Vern9())
- `inject_noise`: Boolean for adding measurement noise
- `mean_noise`: Mean of Gaussian noise
- `stddev_noise`: Standard deviation of Gaussian noise
- `abstol`: Absolute tolerance for solver
- `reltol`: Relative tolerance for solver

Returns:
    OrderedDict containing time series data for each measured variable.
    Value element type matches the ODE solution (Float64 normally, Dual under ForwardDiff).
"""
function sample_data(
    problem,
    model,
    measured_data::Vector{ModelingToolkit.Equation},
    time_interval::Vector{T},
    p_true,
    u0,
    num_points::Int;
    uneven_sampling = false,
    uneven_sampling_times = Vector{T}(),
    solver = Vern9(),
    inject_noise = false,
    mean_noise = zero(T),
    stddev_noise = one(T),
    abstol = convert(T, 1e-10),
    reltol = convert(T, 1e-10),
) where {T<:Number}
    @assert length(time_interval) == 2 "Time interval must be [start_time, end_time]"

    if uneven_sampling
        if length(uneven_sampling_times) == 0
            error("No uneven sampling times provided")
        end
        if length(uneven_sampling_times) != num_points
            error("Uneven sampling times must be of length num_points")
        end
        sampling_times = uneven_sampling_times
    else
        sampling_times = range(time_interval[1], time_interval[2], length = num_points)
    end

    # Use remake to set parameters — supports ForwardDiff Dual types and is thread-safe.
    # Dict-based parameter mapping is required for MTK v10 @mtkcompile models.
    param_syms = ModelingToolkit.parameters(model)
    remade_problem = SciMLBase.remake(problem; p = Dict(param_syms .=> p_true))

    local solution_true
    # Disable all logging below Error level to suppress ODE solver warnings
    solution_true = Logging.with_logger(Logging.NullLogger()) do
        ModelingToolkit.solve(
            remade_problem,
            solver,
            saveat = sampling_times;
            abstol,
            reltol,
            verbose = false,
            maxiters = 1000000,
        )
    end

    # Determine element type from solution (Float64 normally, Dual under ForwardDiff)
    first_values = solution_true[Num(measured_data[1].rhs)]
    V = eltype(first_values)

    data_sample = DataStructures.OrderedDict{Any,Vector{V}}(
        Num(v.lhs) => solution_true[Num(v.rhs)] for v in measured_data
    )

    if inject_noise
        for (key, sample) in data_sample
            data_sample[key] = sample + randn(num_points) .* stddev_noise .+ mean_noise
        end
    end

    data_sample["t"] = sampling_times
    return data_sample
end

"""
    sample_data!(out, integrator, p_test, setp_fn, measured_data)

In-place variant of `sample_data` for the audit / Tolerant-objective hot path.

Reuses a pre-allocated integrator (built via `SciMLBase.init` at construction time)
and a pre-allocated output buffer `out::OrderedDict`. The parameter cache is mutated
via `setp_fn` (no allocation), then `reinit!` + `solve!` reuse the integrator's
working memory.

Per spike (`experiments/sandbox/spike_integrator_reuse.jl`): 210× faster, 168× less
allocation than the `remake(problem; p=...)+solve(...)` path, at bit-identical output.

Arguments:
- `out::OrderedDict{Any,Vector{Float64}}`: pre-allocated; each `Num(v.lhs)` key
  must already hold a `Vector{Float64}` of the right length. `out["t"]` is preserved.
- `integrator`: built once via `SciMLBase.init(problem, solver; saveat, abstol, reltol)`.
- `p_test::AbstractVector{Float64}`: current parameter vector.
- `setp_fn`: closure returned by `SciMLBase.setp(problem, parameters(model))`.
- `measured_data::Vector{ModelingToolkit.Equation}`: same as for `sample_data`.

Returns `out` (mutated). The integrator's internal `sol` buffer is overwritten on the
next call — do not retain references across calls.
"""
function sample_data!(
    out::DataStructures.OrderedDict,
    integrator,
    p_test::AbstractVector{Float64},
    setp_fn,
    measured_data::Vector{ModelingToolkit.Equation},
)
    setp_fn(integrator, p_test)
    Logging.with_logger(Logging.NullLogger()) do
        SciMLBase.reinit!(integrator, integrator.sol.prob.u0; tstops = Float64[])
        SciMLBase.solve!(integrator)
    end
    sol = integrator.sol
    for v in measured_data
        key = Num(v.lhs)
        buf = out[key]
        rhs_vals = sol[Num(v.rhs)]
        copyto!(buf, rhs_vals)
    end
    return out
end
