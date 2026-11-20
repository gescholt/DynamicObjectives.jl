"""
DynamicObjectivesTermExt — styled terminal panels, headers and progress output

Loaded when Term is loaded next to DynamicObjectives. It supplies the methods of
`display_model_summary`, `display_optimization_result`,
`display_optimization_progress`, `display_section`, `display_subsection`,
`display_step`, `stop_spinner!`, and the internal renderers `render_progress`
(behind `progress_bar`, `update_progress!`, `finish_progress!`, `with_progress`)
and `render_spinner` (behind `spinner`, `spin!`, `with_spinner`). The `ProgressBar`
and `Spinner` structs stay in the core package, since an extension cannot define
exported types.

Triggered by: using DynamicObjectives; using Term
"""
module DynamicObjectivesTermExt

using DynamicObjectives
using DynamicObjectives:
    DisplayConfig, DEFAULT_CONFIG, ProgressBar, Spinner, SPINNER_FRAMES, format_duration
using DynamicObjectives: ModelingToolkit
using Printf
using Term: Panel, apply_style

function DynamicObjectives.display_model_summary(
    model::ModelingToolkit.AbstractSystem;
    config::DisplayConfig = DEFAULT_CONFIG,
)
    params = ModelingToolkit.parameters(model)
    states = ModelingToolkit.unknowns(model)
    eqs = ModelingToolkit.equations(model)

    panel_content = ""
    panel_content *= "Model: $(nameof(model))\n"
    panel_content *= "━"^50 * "\n"
    panel_content *= "States ($(length(states))): $(join(string.(states), ", "))\n"
    panel_content *= "Parameters ($(length(params))): $(join(string.(params), ", "))\n"
    panel_content *= "Equations: $(length(eqs))\n"

    if config.use_color
        panel = Panel(
            panel_content,
            title = "ODE Model Summary",
            title_style = "bold cyan",
            style = "cyan",
            fit = true,
        )
        println(panel)
    else
        println("\n=== ODE Model Summary ===")
        println(panel_content)
        println("="^50)
    end
end

function DynamicObjectives.display_optimization_result(
    result;
    config::DisplayConfig = DEFAULT_CONFIG,
)
    params = get(result, :params, nothing)
    error_val = get(result, :error, nothing)
    iterations = get(result, :iterations, nothing)
    converged = get(result, :converged, nothing)
    time_elapsed = get(result, :time_elapsed, nothing)

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

    if config.use_color
        panel = Panel(
            summary,
            title = "Optimization Result",
            title_style = "bold green",
            style = "green",
            fit = true,
        )
        println(panel)
    else
        println("\n=== Optimization Result ===")
        println(summary)
        println("="^50)
    end
end

function DynamicObjectives.display_optimization_progress(
    iteration::Int,
    params::AbstractVector,
    error::Real;
    config::DisplayConfig = DEFAULT_CONFIG,
)
    param_str = join([round(p, digits = config.precision) for p in params], ", ")

    if abs(error) < 1e-3 || abs(error) > 1e4
        error_str = @sprintf("%.4e", error)
    else
        error_str = string(round(error, digits = config.precision))
    end

    if config.use_color
        print("\r")  # Carriage return for same-line update
        print("Iter: ", apply_style("$iteration", "bold cyan"), " | ")
        print("Error: ", apply_style(error_str, "bold yellow"), " | ")
        print("Params: ", apply_style("[$(param_str)]", "bold white"))
    else
        print("\rIter: $iteration | Error: $error_str | Params: [$(param_str)]")
    end
    flush(stdout)
end

function DynamicObjectives.display_section(
    title::String;
    subtitle::String = "",
    config::DisplayConfig = DEFAULT_CONFIG,
)
    content = subtitle == "" ? title : "$title\n$subtitle"

    if config.use_color
        panel = Panel(
            content,
            title_style = "bold white",
            style = "bold blue",
            fit = true,
            width = 80,
        )
        println()
        println(panel)
    else
        println()
        println("="^80)
        println(title)
        if subtitle != ""
            println(subtitle)
        end
        println("="^80)
    end
end

function DynamicObjectives.display_subsection(
    title::String;
    config::DisplayConfig = DEFAULT_CONFIG,
)
    if config.use_color
        println()
        println(apply_style("▶ $title", "bold cyan"))
        println(apply_style("─"^80, "dim"))
    else
        println()
        println(title)
        println("-"^80)
    end
end

function DynamicObjectives.render_progress(pb::ProgressBar)
    pct = pb.current / pb.total
    filled = round(Int, pct * pb.width)
    empty = pb.width - filled

    elapsed = time() - pb.start_time
    if pb.current > 0 && !pb.finished
        eta = elapsed / pb.current * (pb.total - pb.current)
        eta_str = format_duration(eta)
    else
        eta_str = "--:--"
    end

    elapsed_str = format_duration(elapsed)

    if pb.config.use_color
        bar_filled = apply_style("█"^filled, "bold green")
        bar_empty = apply_style("░"^empty, "dim")
        pct_str = apply_style(@sprintf("%5.1f%%", pct * 100), "bold cyan")
        count_str = apply_style("[$(pb.current)/$(pb.total)]", "bold white")
        time_str = apply_style("$(elapsed_str)<$(eta_str)", "dim")
        desc_str = apply_style(pb.description, "bold")
    else
        bar_filled = "█"^filled
        bar_empty = "░"^empty
        pct_str = @sprintf("%5.1f%%", pct * 100)
        count_str = "[$(pb.current)/$(pb.total)]"
        time_str = "$(elapsed_str)<$(eta_str)"
        desc_str = pb.description
    end

    print("\r")
    print("$desc_str $count_str |$bar_filled$bar_empty| $pct_str $time_str")
    print(" "^10)  # Clear any trailing characters
    flush(stdout)
end

function DynamicObjectives.stop_spinner!(
    sp::Spinner;
    message::Union{String,Nothing} = nothing,
    success::Bool = true,
)
    sp.active = false
    final_msg = message !== nothing ? message : sp.message
    elapsed = time() - sp.start_time

    icon = success ? "✓" : "✗"
    icon_style = success ? "bold green" : "bold red"

    print("\r")
    if sp.config.use_color
        print(apply_style(icon, icon_style), " ")
        print(apply_style(final_msg, "bold"), " ")
        print(apply_style("($(format_duration(elapsed)))", "dim"))
    else
        print("$icon $final_msg ($(format_duration(elapsed)))")
    end
    print(" "^20)
    println()
    flush(stdout)
end

function DynamicObjectives.render_spinner(sp::Spinner)
    elapsed = time() - sp.start_time
    frame_char = SPINNER_FRAMES[sp.frame]

    print("\r")
    if sp.config.use_color
        print(apply_style(frame_char, "bold cyan"), " ")
        print(sp.message, " ")
        print(apply_style("($(format_duration(elapsed)))", "dim"))
    else
        print("$frame_char $(sp.message) ($(format_duration(elapsed)))")
    end
    print(" "^10)
    flush(stdout)
end

function DynamicObjectives.display_step(
    step_num::Int,
    total_steps::Int,
    title::String;
    status::Symbol = :pending,
    config::DisplayConfig = DEFAULT_CONFIG,
)
    icons = Dict(
        :pending => ("○", "dim"),
        :running => ("◐", "bold yellow"),
        :completed => ("●", "bold green"),
        :failed => ("✗", "bold red"),
    )

    icon, style = get(icons, status, ("?", "dim"))

    step_str = "[$step_num/$total_steps]"

    if config.use_color
        println(
            apply_style(icon, style),
            " ",
            apply_style(step_str, "bold cyan"),
            " ",
            apply_style(title, status == :running ? "bold" : "default"),
        )
    else
        println("$icon $step_str $title")
    end
end

end # module DynamicObjectivesTermExt
