"""
DynamicObjectivesUnicodePlotsExt — terminal line plots

Loaded when UnicodePlots is loaded next to DynamicObjectives. It supplies
`display_comparison` and the plots of `display_time_series` (the preview table
of `display_time_series` comes from DynamicObjectivesPrettyTablesExt, so
`display_time_series(data; show_plot = false)` needs PrettyTables only).

Triggered by: using DynamicObjectives; using UnicodePlots
"""
module DynamicObjectivesUnicodePlotsExt

using DynamicObjectives
using DynamicObjectives: DisplayConfig, DEFAULT_CONFIG, OrderedDict
using UnicodePlots: lineplot, lineplot!

function DynamicObjectives._display_time_series_plots(
    t::AbstractVector,
    data_dict::AbstractDict,
    var_keys,
    config::DisplayConfig,
)
    println()
    for var in var_keys
        y = data_dict[var]

        if all(isfinite.(y))
            plot_title = "$(var) vs Time"
            color_kw = config.use_color ? (color = :cyan,) : NamedTuple()
            p = lineplot(
                t,
                y;
                title = plot_title,
                xlabel = "Time",
                ylabel = string(var),
                width = config.plot_width,
                height = config.plot_height ÷ 2,
                border = :solid,
                color_kw...,
            )
            println(p)
            println()
        else
            @warn "Skipping plot for $var: contains non-finite values"
        end
    end
end

function DynamicObjectives.display_comparison(
    data_true::OrderedDict,
    data_test::OrderedDict;
    labels::Vector{String} = ["Reference", "Test"],
    config::DisplayConfig = DEFAULT_CONFIG,
)
    if !haskey(data_true, "t") || !haskey(data_test, "t")
        error("Both datasets must contain 't' key")
    end

    t = data_true["t"]
    var_keys = [k for k in keys(data_true) if k != "t"]

    # Plain header: a styled panel would make this function depend on Term too.
    println("\n=== Time Series Comparison: $(labels[1]) vs $(labels[2]) ===\n")

    for var in var_keys
        y_true = data_true[var]
        y_test = data_test[var]

        if all(isfinite.(y_true)) && all(isfinite.(y_test))
            plot_title = "$(var): $(labels[1]) vs $(labels[2])"

            p = lineplot(
                t,
                y_true,
                name = labels[1],
                title = plot_title,
                xlabel = "Time",
                ylabel = string(var),
                width = config.plot_width,
                height = config.plot_height ÷ 2,
                border = :solid,
            )

            lineplot!(p, t, y_test, name = labels[2])

            println(p)
            println()
        else
            @warn "Skipping comparison plot for $var: contains non-finite values"
        end
    end
end

end # module DynamicObjectivesUnicodePlotsExt
