"""
DynamicObjectivesPrettyTablesExt — terminal tables

Loaded when PrettyTables is loaded next to DynamicObjectives. It supplies the
methods of the table-printing display functions (`display_parameters`,
`display_error_metrics`, `display_results`, `display_gradient_analysis`,
`display_quality_summary`, `display_degree_comparison`,
`display_top_critical_points`) and the preview table of `display_time_series`.
The core package holds the docstrings and stubs that raise an error naming
PrettyTables when this extension is absent.

Written against the PrettyTables 2.x API (`backend = Val(...)`, `header`,
`*_crayon` keywords); compat is pinned to "2".

Triggered by: using DynamicObjectives; using PrettyTables
"""
module DynamicObjectivesPrettyTablesExt

using DynamicObjectives
using DynamicObjectives: DisplayConfig, DEFAULT_CONFIG
using LinearAlgebra: norm
using PrettyTables
using Printf

const _BACKENDS =
    Dict(:text => Val(:text), :ascii => Val(:ascii), :markdown => Val(:markdown))

_backend(config::DisplayConfig) = get(_BACKENDS, config.table_backend, Val(:text))

function DynamicObjectives.display_parameters(
    param_names::Vector,
    param_values;
    labels::Vector{String} = ["Value"],
    config::DisplayConfig = DEFAULT_CONFIG,
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
            return string(round(v, digits = config.precision))
        end
        return string(v)
    end

    backend = _backend(config)

    println()
    if config.use_color
        pretty_table(
            data,
            header = header,
            backend = backend,
            formatters = formatters,
            header_crayon = crayon"bold cyan",
            border_crayon = crayon"cyan",
        )
    else
        pretty_table(data, header = header, backend = backend, formatters = formatters)
    end
    println()
end

function DynamicObjectives._display_time_series_table(
    t::AbstractVector,
    data_dict::AbstractDict,
    var_keys,
    config::DisplayConfig,
)
    # Display table with first/last few points (or all points if dataset is small)
    n_preview = min(5, length(t))
    show_ellipsis = length(t) > 2 * n_preview

    if show_ellipsis
        preview_indices = vcat(1:n_preview, (length(t)-n_preview+1):length(t))
        preview_indices = unique(sort(preview_indices))
    else
        preview_indices = collect(1:length(t))
    end

    # Build data matrix
    data_matrix = zeros(length(preview_indices), length(var_keys) + 1)
    data_matrix[:, 1] = t[preview_indices]

    for (i, var) in enumerate(var_keys)
        data_matrix[:, i+1] = data_dict[var][preview_indices]
    end

    # Format with ellipsis if truncated
    row_labels = string.(preview_indices)
    if show_ellipsis
        mid_idx = n_preview
        row_labels[mid_idx] = "⋮"
        data_matrix[mid_idx, :] .= NaN
    end

    header_row = vcat(["Time"], string.(var_keys))

    formatters = (v, i, j) -> begin
        if isnan(v)
            return "⋮"
        elseif v isa Number
            return string(round(v, digits = config.precision))
        end
        return string(v)
    end

    title = "Time Series Data ($(length(t)) points)"
    backend = _backend(config)

    println()
    if config.use_color
        pretty_table(
            data_matrix,
            header = header_row,
            row_labels = row_labels,
            backend = backend,
            formatters = formatters,
            title = title,
            title_crayon = crayon"bold green",
            header_crayon = crayon"bold green",
            border_crayon = crayon"green",
        )
    else
        pretty_table(
            data_matrix,
            header = header_row,
            row_labels = row_labels,
            backend = backend,
            formatters = formatters,
            title = title,
        )
    end
end

function DynamicObjectives.display_error_metrics(
    errors;
    labels = nothing,
    config::DisplayConfig = DEFAULT_CONFIG,
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
                return string(round(v, digits = config.precision))
            end
        end
        return string(v)
    end

    backend = _backend(config)

    println()
    if config.use_color
        pretty_table(
            data,
            header = header,
            backend = backend,
            formatters = formatters,
            title = "Error Metrics",
            title_crayon = crayon"bold yellow",
            header_crayon = crayon"bold yellow",
            border_crayon = crayon"yellow",
        )
    else
        pretty_table(
            data,
            header = header,
            backend = backend,
            formatters = formatters,
            title = "Error Metrics",
        )
    end
    println()
end

function DynamicObjectives.display_results(
    results::Union{Vector{<:Pair},AbstractDict};
    title::String = "Results",
    config::DisplayConfig = DEFAULT_CONFIG,
)
    # Convert to vector of pairs if dict
    pairs = results isa AbstractDict ? collect(results) : results

    # Build data matrix
    data = Matrix{Any}(undef, length(pairs), 2)
    for (i, (key, value)) in enumerate(pairs)
        data[i, 1] = string(key)
        if value isa AbstractFloat
            if abs(value) < 1e-3 || abs(value) > 1e4
                data[i, 2] = @sprintf("%.4e", value)
            else
                data[i, 2] = round(value, digits = config.precision)
            end
        else
            data[i, 2] = value
        end
    end

    backend = _backend(config)

    println()
    if config.use_color
        pretty_table(
            data,
            header = ["Metric", "Value"],
            header_crayon = crayon"bold cyan",
            border_crayon = crayon"cyan",
            backend = backend,
            alignment = [:l, :r],
            title = title,
            title_crayon = crayon"bold white",
        )
    else
        pretty_table(
            data,
            header = ["Metric", "Value"],
            backend = backend,
            alignment = [:l, :r],
            title = title,
        )
    end
    println()
end

function DynamicObjectives.display_gradient_analysis(
    norms::Vector{Float64};
    tolerance::Float64 = 1e-6,
    title::String = "Gradient Norm Analysis",
    config::DisplayConfig = DEFAULT_CONFIG,
)
    # Count valid and invalid points
    valid_count = count(n -> n <= tolerance, norms)
    invalid_count = length(norms) - valid_count
    total_count = length(norms)

    # Compute statistics
    min_norm = minimum(norms)
    mean_norm = sum(norms) / length(norms)
    max_norm = maximum(norms)

    results = [
        "Valid points" => "$valid_count/$total_count",
        "Invalid points" => "$invalid_count/$total_count",
        "Min ||∇f||" => min_norm,
        "Mean ||∇f||" => mean_norm,
        "Max ||∇f||" => max_norm,
        "Tolerance" => tolerance,
    ]

    DynamicObjectives.display_results(results; title = title, config = config)
end

function DynamicObjectives.display_quality_summary(
    refined_results;
    title::String = "Critical Point Quality",
    config::DisplayConfig = DEFAULT_CONFIG,
)
    n_raw = refined_results.n_raw
    n_converged = refined_results.n_converged
    success_rate = 100.0 * n_converged / n_raw
    best_value = refined_results.best_refined_value
    mean_improvement = refined_results.mean_improvement

    results = [
        "Raw critical points" => n_raw,
        "Converged points" => n_converged,
        "Success rate (%)" => success_rate,
        "Mean improvement" => mean_improvement,
        "Best objective value" => best_value,
    ]

    DynamicObjectives.display_results(results; title = title, config = config)
end

function DynamicObjectives.display_degree_comparison(
    degree_results::Vector;
    title::String = "Degree Comparison",
    config::DisplayConfig = DEFAULT_CONFIG,
)
    n_degrees = length(degree_results)
    data = Matrix{Any}(undef, n_degrees, 3)

    for (i, result) in enumerate(degree_results)
        data[i, 1] = result.degree
        data[i, 2] = result.n_critical_points

        best_obj = result.best_objective
        if abs(best_obj) < 1e-3 || abs(best_obj) > 1e4
            data[i, 3] = @sprintf("%.4e", best_obj)
        else
            data[i, 3] = round(best_obj, digits = config.precision)
        end
    end

    backend = _backend(config)

    println()
    if config.use_color
        pretty_table(
            data,
            header = ["Degree", "Critical Points", "Best Objective"],
            header_crayon = crayon"bold cyan",
            border_crayon = crayon"cyan",
            backend = backend,
            alignment = [:c, :c, :r],
            title = title,
            title_crayon = crayon"bold white",
        )
    else
        pretty_table(
            data,
            header = ["Degree", "Critical Points", "Best Objective"],
            backend = backend,
            alignment = [:c, :c, :r],
            title = title,
        )
    end
    println()
end

function DynamicObjectives.display_top_critical_points(
    points::Vector,
    objective_values::Vector{<:Real},
    param_names::Vector;
    n::Int = 5,
    true_params::Union{Vector{<:Real},Nothing} = nothing,
    title::String = "Top Critical Points",
    config::DisplayConfig = DEFAULT_CONFIG,
)
    sorted_indices = sortperm(objective_values)
    n_display = min(n, length(sorted_indices))

    if n_display == 0
        println("No critical points to display")
        return
    end

    n_params = length(param_names)

    header = ["Rank", "Objective"]
    for name in param_names
        push!(header, string(name))
    end
    if true_params !== nothing
        push!(header, "Rel. Error")
    end

    data = Matrix{Any}(undef, n_display, length(header))

    for (row, idx) in enumerate(sorted_indices[1:n_display])
        point = points[idx]
        obj_val = objective_values[idx]

        data[row, 1] = row  # Rank

        if abs(obj_val) < 1e-3 || abs(obj_val) > 1e4
            data[row, 2] = @sprintf("%.4e", obj_val)
        else
            data[row, 2] = round(obj_val, digits = config.precision)
        end

        for (j, p) in enumerate(point)
            data[row, 2+j] = round(p, digits = config.precision)
        end

        if true_params !== nothing
            rel_error = norm(point .- true_params) / norm(true_params)
            data[row, end] = @sprintf("%.2f%%", rel_error * 100)
        end
    end

    backend = _backend(config)
    alignment =
        vcat([:c, :r], fill(:r, n_params), true_params !== nothing ? [:r] : Symbol[])

    println()
    if config.use_color
        pretty_table(
            data,
            header = header,
            header_crayon = crayon"bold cyan",
            border_crayon = crayon"cyan",
            backend = backend,
            alignment = alignment,
            title = title,
            title_crayon = crayon"bold white",
        )
    else
        pretty_table(
            data,
            header = header,
            backend = backend,
            alignment = alignment,
            title = title,
        )
    end
    println()

    top_values = objective_values[sorted_indices[1:n_display]]
    println("  Best: $(@sprintf("%.4e", minimum(top_values)))")
    println("  Worst (of top $n_display): $(@sprintf("%.4e", maximum(top_values)))")

    if true_params !== nothing
        best_point = points[sorted_indices[1]]
        best_error = norm(best_point .- true_params) / norm(true_params)
        println("  Best point relative error: $(@sprintf("%.2f%%", best_error * 100))")
    end
end

end # module DynamicObjectivesPrettyTablesExt
