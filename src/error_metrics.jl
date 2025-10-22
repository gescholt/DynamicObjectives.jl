# Error Metrics
# Distance functions and error metric construction for parameter estimation

"""
L1 norm distance function (scaled by 100)

Arguments:
- Y_true: True values
- Y_test: Test values

Returns:
    100 * ||Y_true - Y_test||₁
"""
L1_norm(Y_true, Y_test) = 100 * norm(Y_true - Y_test, 1)

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
Log L2 norm distance function

Arguments:
- Y_true: True values
- Y_test: Test values

Returns:
    log₂(||Y_true - Y_test||₂ + ε)
"""
log_L2_norm(Y_true, Y_test) = log2(norm(Y_true - Y_test, 2) + eps(eltype(Y_true)))

"""
    make_error_distance(model, outputs, initial_conditions, p_true, time_interval,
                       numpoints, distance_function, aggregate_distances, add_noise_in_time_series;
                       return_inf_on_error, eval_timeout)

Construct an error function comparing model predictions against reference data.

Arguments:
- `model`: ModelingToolkit ODESystem
- `outputs`: Vector of measurement equations
- `initial_conditions`: Vector of initial conditions for the ODE system
- `p_true`: Vector of true parameter values
- `time_interval`: Time interval [start_time, end_time] for simulation
- `numpoints`: Number of time points to sample (default: 5)
- `distance_function`: Function to compute distance (default: L2_norm).
    The function should take two vectors (true and predicted) and return a scalar distance value.
- `aggregate_distances`: Function to aggregate multiple distance values (default: first)
- `add_noise_in_time_series`: Function to add noise to time series (default: nothing = no noise)
- `return_inf_on_error`: If true (default), return Inf when ODE solving fails (for optimization compatibility).
    If false, throw detailed error messages (for debugging).
- `eval_timeout`: Maximum time in seconds for a single objective function evaluation (default: nothing = no timeout).
    If evaluation exceeds this time, returns Inf (when return_inf_on_error=true) or throws TimeoutError.

Returns:
    Function that computes error between predictions and reference data
"""
function make_error_distance(
    model::ModelingToolkit.ODESystem,
    outputs::Vector{ModelingToolkit.Equation},
    initial_conditions::Vector{Float64},
    p_true::Vector{T},
    time_interval,
    numpoints::Int = 5,
    distance_function = L2_norm,
    aggregate_distances = first,
    add_noise_in_time_series = nothing;
    return_inf_on_error::Bool = true,  # For optimization compatibility
    eval_timeout::Union{Float64, Nothing} = nothing  # Timeout per evaluation
) where {T}
    @assert length(p_true) == length(ModelingToolkit.parameters(model)) "Parameter vector length mismatch"
    @assert length(initial_conditions) == length(ModelingToolkit.unknowns(model)) "Initial conditions length mismatch"
    @assert length(outputs) > 0 "At least one output variable must be specified"
    @assert numpoints > 0 "Number of points must be greater than zero"
    @assert time_interval[2] > time_interval[1] "End time must be greater than start time"

    # Generate reference solution once during function creation
    problem = ODEProblem(
        ModelingToolkit.complete(model),
        merge(
            Dict(ModelingToolkit.unknowns(model) .=> initial_conditions),
            Dict(ModelingToolkit.parameters(model) .=> p_true)
        ),
        time_interval
    )

    data_sample_true = sample_data(
        problem,
        model,
        outputs,
        time_interval,
        p_true,
        initial_conditions,
        numpoints
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

    function Error_distance(
        p_test::Union{SVector{N, T2}, Vector{T2}}
    ) where {T2, N}

        # Timeout wrapper: use @async + timedwait if eval_timeout is set
        if eval_timeout !== nothing
            # Run computation in async task
            task = @async begin
                _compute_error_distance(
                    p_test, outputs, time_interval, numpoints,
                    data_sample_true, problem, model, initial_conditions,
                    distance_function, aggregate_distances
                )
            end

            # Wait for task with timeout
            timed_result = timedwait(() -> istaskdone(task), eval_timeout)

            if timed_result == :timed_out
                # Kill the task (best effort - may not actually interrupt ODE solver)
                Base.throwto(task, InterruptException())

                if return_inf_on_error
                    return Inf  # Guide optimizer away from slow regions
                else
                    error("Objective evaluation timed out after $(eval_timeout)s for p_test=$p_test")
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
                        return Inf
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
                return _compute_error_distance(
                    p_test, outputs, time_interval, numpoints,
                    data_sample_true, problem, model, initial_conditions,
                    distance_function, aggregate_distances
                )
            catch e
                if isa(e, InterruptException)
                    rethrow(e)
                end
                # Return Inf to guide optimizer away from this region
                return Inf
            end
        else
            # Let errors propagate for debugging
            return _compute_error_distance(
                p_test, outputs, time_interval, numpoints,
                data_sample_true, problem, model, initial_conditions,
                distance_function, aggregate_distances
            )
        end
    end

    # Inner function with actual computation and error checking
    function _compute_error_distance(
        p_test, measured_data, time_interval, datasize,
        data_sample_true, problem, model, initial_conditions,
        distance_function, aggregate_distances
    )
        # Check datasize consistency
        if datasize != length(data_sample_true["t"])
            error("Datasize mismatch: requested $datasize but reference data has $(length(data_sample_true["t"])) points")
        end

        # Sample data with test parameters
        data_sample_test = try
            sample_data(
                problem,
                model,
                measured_data,
                time_interval,
                p_test,
                initial_conditions,
                datasize
            )
        catch e
            if isa(e, InterruptException)
                rethrow(e)
            end
            # Provide context about which parameters caused the sampling failure
            error("ODE sampling failed for parameters p_test=$p_test: $(sprint(showerror, e))")
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
                error("NaN detected in sampled data for variable '$key' with parameters p_test=$p_test")
            end
            if any(isinf.(values))
                error("Inf detected in sampled data for variable '$key' with parameters p_test=$p_test")
            end
            if any(isnan.(data_sample_true[key]))
                error("NaN detected in reference data for variable '$key'")
            end
            if any(isinf.(data_sample_true[key]))
                error("Inf detected in reference data for variable '$key'")
            end
        end

        return aggregate_distances(
            [
            distance_function(data_sample_true[key], data_sample_test[key])
            for key in keys(data_sample_true) if key != "t"
        ]
        )
    end

    return Error_distance
end
