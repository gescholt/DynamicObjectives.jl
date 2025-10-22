# Data Generation
# Functions for generating synthetic time series data from ODE systems

"""
    sample_data(problem, model, measured_data, time_interval, p_true, u0, num_points; kwargs...)

Generate synthetic time series data from an ODE system with specified parameters.

Arguments:
- `problem`: ODEProblem instance
- `model`: ModelingToolkit ODESystem representing the differential equations
- `measured_data`: Vector of measurement equations
- `time_interval`: [start_time, end_time] for simulation
- `p_true`: Parameter values (Vector or SVector)
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
    OrderedDict containing time series data for each measured variable
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
    reltol = convert(T, 1e-10)
) where {T <: Number}

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

    problem.p.tunable .= p_true
    local solution_true
    # Disable all logging below Error level to suppress ODE solver warnings
    solution_true = Logging.with_logger(Logging.NullLogger()) do
        ModelingToolkit.solve(problem, solver, saveat = sampling_times; abstol, reltol, verbose=false, maxiters=1000000)
    end

    data_sample = DataStructures.OrderedDict{Any, Vector{T}}(
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
