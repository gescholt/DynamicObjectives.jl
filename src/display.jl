# Display Infrastructure
# Pretty and readable output for Dynamic_objectives data structures

using PrettyTables
using Term
using UnicodePlots
using DataStructures: OrderedDict
using Printf

"""
    DisplayConfig

Configuration for display output customization.

Fields:
- `use_color`: Enable colored output (default: true)
- `table_backend`: Backend for tables (:unicode, :ascii, :markdown) (default: :unicode)
- `plot_width`: Width of plots in characters (default: 60)
- `plot_height`: Height of plots in characters (default: 20)
- `precision`: Decimal precision for numerical output (default: 4)
"""
Base.@kwdef mutable struct DisplayConfig
    use_color::Bool = true
    table_backend::Symbol = :unicode
    plot_width::Int = 60
    plot_height::Int = 20
    precision::Int = 4
end

# Global default configuration
const DEFAULT_CONFIG = DisplayConfig()

"""
    set_display_config!(; kwargs...)

Update the global display configuration.

# Examples
```julia
set_display_config!(use_color=false, precision=6)
```
"""
function set_display_config!(; kwargs...)
    for (key, value) in kwargs
        setfield!(DEFAULT_CONFIG, key, value)
    end
end

"""
    display_model_summary(model; config=DEFAULT_CONFIG)

Display a formatted summary of an ODE model.

# Arguments
- `model`: ModelingToolkit ODESystem
- `config`: DisplayConfig instance (optional)

# Examples
```julia
model, params, states, outputs = define_daisy_ex3_model_4D()
display_model_summary(model)
```
"""
function display_model_summary(model::ModelingToolkit.ODESystem; config::DisplayConfig=DEFAULT_CONFIG)
    params = ModelingToolkit.parameters(model)
    states = ModelingToolkit.unknowns(model)
    eqs = ModelingToolkit.equations(model)

    # Create summary panel
    panel_content = ""
    panel_content *= "Model: $(nameof(model))\n"
    panel_content *= "━"^50 * "\n"
    panel_content *= "States ($(length(states))): $(join(string.(states), ", "))\n"
    panel_content *= "Parameters ($(length(params))): $(join(string.(params), ", "))\n"
    panel_content *= "Equations: $(length(eqs))\n"

    if config.use_color
        panel = Term.Panel(
            panel_content,
            title="ODE Model Summary",
            title_style="bold cyan",
            style="cyan",
            fit=true
        )
        println(panel)
    else
        println("\n=== ODE Model Summary ===")
        println(panel_content)
        println("="^50)
    end
end

"""
    display_parameters(param_names, param_values; labels=["Parameters"], config=DEFAULT_CONFIG)

Display parameter values in a formatted table.

# Arguments
- `param_names`: Vector of parameter names (as symbols or strings)
- `param_values`: Vector or matrix of parameter values
- `labels`: Column labels for parameter sets (default: ["Parameters"])
- `config`: DisplayConfig instance (optional)

# Examples
```julia
param_names = [:α, :β, :γ, :δ]
p_true = [0.1, 0.2, 0.3, 0.4]
p_test = [0.15, 0.25, 0.35, 0.45]
display_parameters(param_names, [p_true p_test], labels=["True", "Test"])
```
"""
function display_parameters(
    param_names::Vector,
    param_values;
    labels::Vector{String}=["Value"],
    config::DisplayConfig=DEFAULT_CONFIG
)
    # Convert to matrix if vector
    if param_values isa AbstractVector
        param_values = reshape(param_values, :, 1)
    end

    # Create data matrix with parameter names
    data = hcat(string.(param_names), param_values)
    header = vcat(["Parameter"], labels)

    # Format numbers with specified precision
    formatters = (v, i, j) -> begin
        if j > 1 && v isa Number
            return string(round(v, digits=config.precision))
        end
        return string(v)
    end

    # Select backend
    backend_map = Dict(
        :unicode => Val(:unicode),
        :ascii => Val(:ascii),
        :markdown => Val(:markdown)
    )
    backend = get(backend_map, config.table_backend, Val(:unicode))

    if config.use_color
        println()
        pretty_table(
            data,
            header=header,
            backend=backend,
            formatters=formatters,
            header_crayon=crayon"bold cyan",
            border_crayon=crayon"cyan"
        )
        println()
    else
        println()
        pretty_table(
            data,
            header=header,
            backend=backend,
            formatters=formatters
        )
        println()
    end
end

"""
    display_time_series(data_dict; variables=nothing, show_plot=true, config=DEFAULT_CONFIG)

Display time series data with optional plots.

# Arguments
- `data_dict`: OrderedDict with "t" key for time and variable names as keys
- `variables`: Vector of variables to display (nothing = all) (optional)
- `show_plot`: Whether to show plots (default: true)
- `config`: DisplayConfig instance (optional)

# Examples
```julia
model, params, states, outputs = define_daisy_ex3_model_4D()
p_true = [0.1, 0.2, 0.3, 0.4]
ic = [1.0, 2.0, 1.0, 1.0]
problem = ODEProblem(...)
data = sample_data(problem, model, outputs, [0.0, 10.0], p_true, ic, 25)
display_time_series(data)
```
"""
function display_time_series(
    data_dict::OrderedDict;
    variables=nothing,
    show_plot::Bool=true,
    config::DisplayConfig=DEFAULT_CONFIG
)
    if !haskey(data_dict, "t")
        error("Time series data must contain 't' key for time values")
    end

    t = data_dict["t"]

    # Select variables to display
    var_keys = if variables === nothing
        [k for k in keys(data_dict) if k != "t"]
    else
        variables
    end

    if isempty(var_keys)
        @warn "No variables to display"
        return
    end

    # Display header
    if config.use_color
        header_panel = Term.Panel(
            "Time Series Data ($(length(t)) points)",
            title="Data Summary",
            title_style="bold green",
            style="green",
            fit=true
        )
        println(header_panel)
    else
        println("\n=== Time Series Data ($(length(t)) points) ===\n")
    end

    # Display table with first/last few points
    n_preview = min(5, length(t))
    preview_indices = vcat(1:n_preview, (length(t)-n_preview+1):length(t))
    preview_indices = unique(sort(preview_indices))

    # Build data matrix
    data_matrix = zeros(length(preview_indices), length(var_keys) + 1)
    data_matrix[:, 1] = t[preview_indices]

    for (i, var) in enumerate(var_keys)
        data_matrix[:, i+1] = data_dict[var][preview_indices]
    end

    # Format with ellipsis if truncated
    row_labels = string.(preview_indices)
    if length(t) > 2 * n_preview
        mid_idx = n_preview
        row_labels[mid_idx] = "⋮"
        data_matrix[mid_idx, :] .= NaN
    end

    header_row = vcat(["Time"], string.(var_keys))

    formatters = (v, i, j) -> begin
        if isnan(v)
            return "⋮"
        elseif v isa Number
            return string(round(v, digits=config.precision))
        end
        return string(v)
    end

    backend_map = Dict(
        :unicode => Val(:unicode),
        :ascii => Val(:ascii),
        :markdown => Val(:markdown)
    )
    backend = get(backend_map, config.table_backend, Val(:unicode))

    println()
    if config.use_color
        pretty_table(
            data_matrix,
            header=header_row,
            row_labels=row_labels,
            backend=backend,
            formatters=formatters,
            header_crayon=crayon"bold green",
            border_crayon=crayon"green"
        )
    else
        pretty_table(
            data_matrix,
            header=header_row,
            row_labels=row_labels,
            backend=backend,
            formatters=formatters
        )
    end

    # Display plots
    if show_plot
        println()
        for var in var_keys
            y = data_dict[var]

            # Check if data is valid for plotting
            if all(isfinite.(y))
                plot_title = "$(var) vs Time"

                if config.use_color
                    p = lineplot(
                        t, y,
                        title=plot_title,
                        xlabel="Time",
                        ylabel=string(var),
                        width=config.plot_width,
                        height=config.plot_height ÷ 2,
                        border=:solid,
                        color=:cyan
                    )
                else
                    p = lineplot(
                        t, y,
                        title=plot_title,
                        xlabel="Time",
                        ylabel=string(var),
                        width=config.plot_width,
                        height=config.plot_height ÷ 2,
                        border=:solid
                    )
                end

                println(p)
                println()
            else
                @warn "Skipping plot for $var: contains non-finite values"
            end
        end
    end
end

"""
    display_comparison(data_true, data_test; labels=["Reference", "Test"], config=DEFAULT_CONFIG)

Display side-by-side comparison of two time series datasets.

# Arguments
- `data_true`: OrderedDict with reference time series
- `data_test`: OrderedDict with test time series
- `labels`: Labels for the two datasets (default: ["Reference", "Test"])
- `config`: DisplayConfig instance (optional)

# Examples
```julia
display_comparison(data_true, data_test, labels=["True Parameters", "Test Parameters"])
```
"""
function display_comparison(
    data_true::OrderedDict,
    data_test::OrderedDict;
    labels::Vector{String}=["Reference", "Test"],
    config::DisplayConfig=DEFAULT_CONFIG
)
    if !haskey(data_true, "t") || !haskey(data_test, "t")
        error("Both datasets must contain 't' key")
    end

    t = data_true["t"]
    var_keys = [k for k in keys(data_true) if k != "t"]

    # Display header
    if config.use_color
        header = Term.Panel(
            "Comparing $(labels[1]) vs $(labels[2])",
            title="Time Series Comparison",
            title_style="bold magenta",
            style="magenta",
            fit=true
        )
        println(header)
    else
        println("\n=== Time Series Comparison: $(labels[1]) vs $(labels[2]) ===\n")
    end

    # Display plots comparing the two datasets
    for var in var_keys
        y_true = data_true[var]
        y_test = data_test[var]

        if all(isfinite.(y_true)) && all(isfinite.(y_test))
            plot_title = "$(var): $(labels[1]) vs $(labels[2])"

            p = lineplot(
                t, y_true,
                name=labels[1],
                title=plot_title,
                xlabel="Time",
                ylabel=string(var),
                width=config.plot_width,
                height=config.plot_height ÷ 2,
                border=:solid
            )

            lineplot!(p, t, y_test, name=labels[2])

            println(p)
            println()
        else
            @warn "Skipping comparison plot for $var: contains non-finite values"
        end
    end
end

"""
    display_error_metrics(errors; labels=nothing, config=DEFAULT_CONFIG)

Display error/distance metrics in a formatted table.

# Arguments
- `errors`: Dictionary or vector of error values
- `labels`: Labels for error metrics (optional)
- `config`: DisplayConfig instance (optional)

# Examples
```julia
errors = Dict("L1" => 0.0123, "L2" => 0.0045, "log_L2" => -7.89)
display_error_metrics(errors)
```
"""
function display_error_metrics(
    errors;
    labels=nothing,
    config::DisplayConfig=DEFAULT_CONFIG
)
    # Convert to dictionary if needed
    error_dict = if errors isa AbstractDict
        errors
    elseif labels !== nothing
        Dict(string(labels[i]) => errors[i] for i in 1:length(errors))
    else
        Dict("Error_$i" => errors[i] for i in 1:length(errors))
    end

    # Create table
    metric_names = collect(keys(error_dict))
    metric_values = [error_dict[k] for k in metric_names]

    data = hcat(metric_names, metric_values)
    header = ["Metric", "Value"]

    formatters = (v, i, j) -> begin
        if j == 2 && v isa Number
            if abs(v) < 1e-3 || abs(v) > 1e4
                return @sprintf("%.4e", v)
            else
                return string(round(v, digits=config.precision))
            end
        end
        return string(v)
    end

    backend_map = Dict(
        :unicode => Val(:unicode),
        :ascii => Val(:ascii),
        :markdown => Val(:markdown)
    )
    backend = get(backend_map, config.table_backend, Val(:unicode))

    if config.use_color
        println()
        panel = Term.Panel(
            "Error Metrics",
            title="Metrics",
            title_style="bold yellow",
            style="yellow",
            fit=true
        )
        println(panel)

        pretty_table(
            data,
            header=header,
            backend=backend,
            formatters=formatters,
            header_crayon=crayon"bold yellow",
            border_crayon=crayon"yellow"
        )
        println()
    else
        println("\n=== Error Metrics ===\n")
        pretty_table(
            data,
            header=header,
            backend=backend,
            formatters=formatters
        )
        println()
    end
end

"""
    display_optimization_result(result; config=DEFAULT_CONFIG)

Display optimization results in a formatted summary.

# Arguments
- `result`: Dictionary or named tuple with optimization results
  Expected keys: :params, :error, :iterations, :converged, :time_elapsed
- `config`: DisplayConfig instance (optional)

# Examples
```julia
result = (
    params=[0.1, 0.2, 0.3, 0.4],
    error=0.0012,
    iterations=150,
    converged=true,
    time_elapsed=2.5
)
display_optimization_result(result)
```
"""
function display_optimization_result(
    result;
    config::DisplayConfig=DEFAULT_CONFIG
)
    # Extract fields
    params = get(result, :params, nothing)
    error_val = get(result, :error, nothing)
    iterations = get(result, :iterations, nothing)
    converged = get(result, :converged, nothing)
    time_elapsed = get(result, :time_elapsed, nothing)

    # Build summary string
    summary = ""

    if converged !== nothing
        status = converged ? "✓ Converged" : "✗ Not Converged"
        summary *= "Status: $status\n"
    end

    if error_val !== nothing
        if abs(error_val) < 1e-3 || abs(error_val) > 1e4
            summary *= @sprintf("Final Error: %.4e\n", error_val)
        else
            summary *= "Final Error: $(round(error_val, digits=config.precision))\n"
        end
    end

    if iterations !== nothing
        summary *= "Iterations: $iterations\n"
    end

    if time_elapsed !== nothing
        summary *= "Time Elapsed: $(round(time_elapsed, digits=2))s\n"
    end

    if params !== nothing
        summary *= "\nOptimal Parameters:\n"
        for (i, p) in enumerate(params)
            summary *= "  p[$i] = $(round(p, digits=config.precision))\n"
        end
    end

    # Display
    if config.use_color
        panel = Term.Panel(
            summary,
            title="Optimization Result",
            title_style="bold green",
            style="green",
            fit=true
        )
        println(panel)
    else
        println("\n=== Optimization Result ===")
        println(summary)
        println("="^50)
    end
end

"""
    display_optimization_progress(iteration, params, error; config=DEFAULT_CONFIG)

Display real-time optimization progress (for interactive use).

# Arguments
- `iteration`: Current iteration number
- `params`: Current parameter values
- `error`: Current error value
- `config`: DisplayConfig instance (optional)

# Examples
```julia
# In optimization callback
display_optimization_progress(iter, current_params, current_error)
```
"""
function display_optimization_progress(
    iteration::Int,
    params::AbstractVector,
    error::Real;
    config::DisplayConfig=DEFAULT_CONFIG
)
    param_str = join([round(p, digits=config.precision) for p in params], ", ")

    if abs(error) < 1e-3 || abs(error) > 1e4
        error_str = @sprintf("%.4e", error)
    else
        error_str = string(round(error, digits=config.precision))
    end

    if config.use_color
        print("\r")  # Carriage return for same-line update
        print("Iter: ", Term.apply_style("$iteration", "bold cyan"), " | ")
        print("Error: ", Term.apply_style(error_str, "bold yellow"), " | ")
        print("Params: ", Term.apply_style("[$(param_str)]", "bold white"))
    else
        print("\rIter: $iteration | Error: $error_str | Params: [$(param_str)]")
    end
    flush(stdout)
end
