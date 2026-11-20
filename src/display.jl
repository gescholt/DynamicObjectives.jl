# Display Infrastructure
# Terminal output for DynamicObjectives data structures.
#
# The renderers live in three package extensions, one per display package, so a
# user loads only the package behind the feature they use:
#
#   DynamicObjectivesPrettyTablesExt  (using PrettyTables)  tables
#   DynamicObjectivesTermExt          (using Term)          panels, headers, progress
#   DynamicObjectivesUnicodePlotsExt  (using UnicodePlots)  line plots
#
# This file keeps the configuration, the progress/spinner state types, the
# pure-string helpers (`format_duration`, `format_critical_points_markdown`), the
# docstrings, and a catch-all stub per moved function. A stub raises an error
# naming the package to load; the extension adds a more specific method, so
# dispatch prefers it and nothing is overwritten. No display function has a
# plain-text fallback when its package is missing.

"""
    DisplayConfig

Configuration for display output customization.

Fields:
- `use_color`: Enable colored output (default: true)
- `table_backend`: Backend for tables (:text, :markdown, :latex) (default: :text)
- `plot_width`: Width of plots in characters (default: 60)
- `plot_height`: Height of plots in characters (default: 20)
- `precision`: Decimal precision for numerical output (default: 4)
"""
Base.@kwdef mutable struct DisplayConfig
    use_color::Bool = true
    table_backend::Symbol = :text
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

# ==============================================================================
# Extension stubs
# ==============================================================================

const _DISPLAY_EXTS = Dict(
    "PrettyTables" => (:DynamicObjectivesPrettyTablesExt, "the terminal tables"),
    "Term" => (:DynamicObjectivesTermExt, "the styled panels and progress output"),
    "UnicodePlots" => (:DynamicObjectivesUnicodePlotsExt, "the terminal line plots"),
)

# Called by every stub. If the extension is loaded, the stub was reached only
# because no extension method matches the argument types: report a MethodError
# rather than a misleading "requires X.jl".
function _display_requires(fn::Function, pkg::String, args::Tuple)
    ext, what = _DISPLAY_EXTS[pkg]
    if Base.get_extension(@__MODULE__, ext) !== nothing
        throw(MethodError(fn, args))
    end
    error(
        "$(nameof(fn)) requires $pkg.jl. Add `using $pkg` to load $what " *
        "(package extension $ext).",
    )
end

# ------------------------------------------------------------------------------
# PrettyTables: tables
# ------------------------------------------------------------------------------

"""
    display_parameters(param_names, param_values; labels=["Value"], config=DEFAULT_CONFIG)

Display parameter values in a formatted table. Requires `using PrettyTables`.

# Arguments
- `param_names`: Vector of parameter names (as symbols or strings)
- `param_values`: Vector or matrix of parameter values
- `labels`: Column labels for parameter sets (default: ["Value"])
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using PrettyTables
param_names = [:α, :β, :γ, :δ]
p_true = [0.1, 0.2, 0.3, 0.4]
p_test = [0.15, 0.25, 0.35, 0.45]
display_parameters(param_names, [p_true p_test], labels=["True", "Test"])
```
"""
display_parameters(args...; kwargs...) =
    _display_requires(display_parameters, "PrettyTables", args)

"""
    display_error_metrics(errors; labels=nothing, config=DEFAULT_CONFIG)

Display error/distance metrics in a table titled "Error Metrics". Requires
`using PrettyTables`.

# Arguments
- `errors`: Dictionary or vector of error values
- `labels`: Labels for error metrics (optional)
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using PrettyTables
errors = Dict("L1" => 0.0123, "L2" => 0.0045, "log_L2" => -7.89)
display_error_metrics(errors)
```
"""
display_error_metrics(args...; kwargs...) =
    _display_requires(display_error_metrics, "PrettyTables", args)

"""
    display_results(results; title="Results", config=DEFAULT_CONFIG)

Display stage/workflow results as a formatted key-value table. Requires
`using PrettyTables`.

# Arguments
- `results`: Vector of Pairs (key => value) or Dict
- `title`: Table title (default: "Results")
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using PrettyTables
display_results([
    "Critical points found" => 42,
    "Best objective value" => 1.23e-6,
    "Time elapsed (s)" => 12.5
], title="Stage 1 Results")
```
"""
display_results(args...; kwargs...) = _display_requires(display_results, "PrettyTables", args)

"""
    display_gradient_analysis(norms; tolerance=1e-6, title="Gradient Norm Analysis", config=DEFAULT_CONFIG)

Display gradient norm analysis table showing validation metrics. Requires
`using PrettyTables`.

# Arguments
- `norms`: Vector of gradient norms (Float64)
- `tolerance`: Threshold for valid gradient (default: 1e-6)
- `title`: Table title (default: "Gradient Norm Analysis")
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using PrettyTables
grad_norms = [1.2e-12, 3.4e-09, 2.1e-05, 8.7e-08]
display_gradient_analysis(grad_norms, tolerance=1e-6)
```
"""
display_gradient_analysis(args...; kwargs...) =
    _display_requires(display_gradient_analysis, "PrettyTables", args)

"""
    display_quality_summary(refined_results; title="Critical Point Quality", config=DEFAULT_CONFIG)

Display critical point quality summary with combined metrics. Requires
`using PrettyTables`.

# Arguments
- `refined_results`: NamedTuple with refinement results (must have fields: n_raw, n_converged, best_refined_value, mean_improvement)
- `title`: Table title (default: "Critical Point Quality")
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using PrettyTables
results = (
    n_raw=81,
    n_converged=78,
    best_refined_value=1.23e-8,
    mean_improvement=2.5,
)
display_quality_summary(results)
```
"""
display_quality_summary(args...; kwargs...) =
    _display_requires(display_quality_summary, "PrettyTables", args)

"""
    display_degree_comparison(degree_results; title="Degree Comparison", config=DEFAULT_CONFIG)

Display degree-by-degree comparison table showing critical points and quality per
degree. Requires `using PrettyTables`.

# Arguments
- `degree_results`: Vector of results per degree (each must have: degree, n_critical_points, best_objective)
- `title`: Table title (default: "Degree Comparison")
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using PrettyTables
degree_results = [
    (degree=4, n_critical_points=15, best_objective=2.3e-5),
    (degree=5, n_critical_points=28, best_objective=1.1e-6),
    (degree=6, n_critical_points=38, best_objective=3.4e-8)
]
display_degree_comparison(degree_results)
```
"""
display_degree_comparison(args...; kwargs...) =
    _display_requires(display_degree_comparison, "PrettyTables", args)

"""
    display_top_critical_points(points, objective_values, param_names;
                                n=5, true_params=nothing, title="Top Critical Points",
                                config=DEFAULT_CONFIG)

Display a table of the top N critical points ranked by objective value. Requires
`using PrettyTables`. For a dependency-free Markdown string see
[`format_critical_points_markdown`](@ref).

# Arguments
- `points`: Vector of parameter vectors (each element is a parameter vector)
- `objective_values`: Vector of objective function values
- `param_names`: Vector of parameter names (symbols or strings)
- `n`: Number of top points to display (default: 5)
- `true_params`: Optional true parameter values for comparison
- `title`: Table title
- `config`: DisplayConfig instance

# Examples
```julia
using PrettyTables
display_top_critical_points(
    refined.refined_points,
    refined.refined_values,
    [:p1, :p2, :p3, :p4],
    n=5,
    true_params=p_true
)
```
"""
display_top_critical_points(args...; kwargs...) =
    _display_requires(display_top_critical_points, "PrettyTables", args)

# Internal: preview table of `display_time_series` (PrettyTables extension).
_display_time_series_table(args...) =
    _display_requires(_display_time_series_table, "PrettyTables", args)

# ------------------------------------------------------------------------------
# UnicodePlots: line plots
# ------------------------------------------------------------------------------

# Internal: per-variable plots of `display_time_series` (UnicodePlots extension).
_display_time_series_plots(args...) =
    _display_requires(_display_time_series_plots, "UnicodePlots", args)

"""
    display_comparison(data_true, data_test; labels=["Reference", "Test"], config=DEFAULT_CONFIG)

Plot two time series datasets against each other, one plot per variable, under a
plain-text header. Requires `using UnicodePlots`.

# Arguments
- `data_true`: OrderedDict with reference time series
- `data_test`: OrderedDict with test time series
- `labels`: Labels for the two datasets (default: ["Reference", "Test"])
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using UnicodePlots
display_comparison(data_true, data_test, labels=["True Parameters", "Test Parameters"])
```
"""
display_comparison(args...; kwargs...) =
    _display_requires(display_comparison, "UnicodePlots", args)

"""
    display_time_series(data_dict; variables=nothing, show_plot=true, config=DEFAULT_CONFIG)

Display time series data: a preview table (first and last points), then one line
plot per variable. The table needs `using PrettyTables`; the plots need
`using UnicodePlots`, so `show_plot = false` needs PrettyTables only.

# Arguments
- `data_dict`: OrderedDict with "t" key for time and variable names as keys
- `variables`: Vector of variables to display (nothing = all) (optional)
- `show_plot`: Whether to show plots (default: true)
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using PrettyTables, UnicodePlots
data = sample_data(problem, model, outputs, [0.0, 10.0], p_true, ic, 25)
display_time_series(data)
```
"""
function display_time_series(
    data_dict::OrderedDict;
    variables = nothing,
    show_plot::Bool = true,
    config::DisplayConfig = DEFAULT_CONFIG,
)
    if !haskey(data_dict, "t")
        error("Time series data must contain 't' key for time values")
    end

    t = data_dict["t"]

    var_keys = if variables === nothing
        [k for k in keys(data_dict) if k != "t"]
    else
        variables
    end

    if isempty(var_keys)
        @warn "No variables to display"
        return
    end

    # Fail before printing anything if the plots were asked for but cannot be drawn.
    if show_plot && Base.get_extension(@__MODULE__, :DynamicObjectivesUnicodePlotsExt) === nothing
        _display_time_series_plots(t, data_dict, var_keys, config)
    end

    _display_time_series_table(t, data_dict, var_keys, config)
    show_plot && _display_time_series_plots(t, data_dict, var_keys, config)
    return nothing
end

# ------------------------------------------------------------------------------
# Term: panels, headers, progress rendering
# ------------------------------------------------------------------------------

"""
    display_model_summary(model; config=DEFAULT_CONFIG)

Display a formatted summary of an ODE model in a panel. Requires `using Term`.

# Arguments
- `model`: ModelingToolkit System
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using Term
model, params, states, outputs = define_daisy_ex3_model_4D()
display_model_summary(model)
```
"""
display_model_summary(args...; kwargs...) =
    _display_requires(display_model_summary, "Term", args)

"""
    display_optimization_result(result; config=DEFAULT_CONFIG)

Display optimization results in a summary panel. Requires `using Term`.

# Arguments
- `result`: Dictionary or named tuple with optimization results
  Expected keys: :params, :error, :iterations, :converged, :time_elapsed
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using Term
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
display_optimization_result(args...; kwargs...) =
    _display_requires(display_optimization_result, "Term", args)

"""
    display_optimization_progress(iteration, params, error; config=DEFAULT_CONFIG)

Display real-time optimization progress on one updating line (for interactive use).
Requires `using Term`.

# Arguments
- `iteration`: Current iteration number
- `params`: Current parameter values
- `error`: Current error value
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using Term
# In optimization callback
display_optimization_progress(iter, current_params, current_error)
```
"""
display_optimization_progress(args...; kwargs...) =
    _display_requires(display_optimization_progress, "Term", args)

"""
    display_section(title; subtitle="", config=DEFAULT_CONFIG)

Display a workflow section header with consistent styling. Requires `using Term`.

# Arguments
- `title`: Main section title
- `subtitle`: Optional subtitle or description
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using Term
display_section("Step 1: Model Setup")
display_section("RESULTS", subtitle="Parameter Recovery Analysis")
```
"""
display_section(args...; kwargs...) = _display_requires(display_section, "Term", args)

"""
    display_subsection(title; config=DEFAULT_CONFIG)

Display a workflow subsection header with consistent styling. Requires `using Term`.

# Arguments
- `title`: Subsection title
- `config`: DisplayConfig instance (optional)

# Examples
```julia
using Term
display_subsection("Configuring experiment parameters")
```
"""
display_subsection(args...; kwargs...) = _display_requires(display_subsection, "Term", args)

"""
    display_step(step_num, total_steps, title; status=:pending, config=DEFAULT_CONFIG)

Display a workflow step with status indicator. Requires `using Term`.

# Arguments
- `step_num`: Current step number
- `total_steps`: Total number of steps
- `title`: Step title/description
- `status`: One of :pending, :running, :completed, :failed
- `config`: DisplayConfig instance

# Examples
```julia
using Term
display_step(1, 3, "Model Setup", status=:completed)
display_step(2, 3, "Running Optimization", status=:running)
display_step(3, 3, "Analyzing Results", status=:pending)
```
"""
display_step(args...; kwargs...) = _display_requires(display_step, "Term", args)

# Internal renderers (Term extension).
render_progress(args...) = _display_requires(render_progress, "Term", args)
render_spinner(args...) = _display_requires(render_spinner, "Term", args)

# ==============================================================================
# Progress Bar Infrastructure (state in core, rendering in the Term extension)
# ==============================================================================

"""
    ProgressBar

Mutable struct for tracking and displaying progress of iterative operations.
Rendering requires `using Term`.

# Fields
- `total`: Total number of steps
- `current`: Current step (0-based internally)
- `description`: Description shown before the bar
- `width`: Width of the progress bar in characters
- `start_time`: Time when progress started
- `config`: DisplayConfig for styling
"""
mutable struct ProgressBar
    total::Int
    current::Int
    description::String
    width::Int
    start_time::Float64
    config::DisplayConfig
    finished::Bool
end

"""
    progress_bar(total; description="Progress", width=40, config=DEFAULT_CONFIG)

Create a new progress bar and draw it. Requires `using Term`.

# Arguments
- `total`: Total number of steps
- `description`: Text shown before the bar (default: "Progress")
- `width`: Width of the bar in characters (default: 40)
- `config`: DisplayConfig instance (optional)

# Returns
- `ProgressBar` instance

# Examples
```julia
using Term
pb = progress_bar(100, description="Processing models")
for i in 1:100
    # do work
    update_progress!(pb)
end
finish_progress!(pb)
```
"""
function progress_bar(
    total::Int;
    description::String = "Progress",
    width::Int = 40,
    config::DisplayConfig = DEFAULT_CONFIG,
)
    pb = ProgressBar(total, 0, description, width, time(), config, false)
    render_progress(pb)
    return pb
end

"""
    update_progress!(pb::ProgressBar; increment=1, description=nothing)

Update the progress bar by incrementing the current count. Requires `using Term`.

# Arguments
- `pb`: ProgressBar instance
- `increment`: Amount to increment (default: 1)
- `description`: Optionally update the description
"""
function update_progress!(
    pb::ProgressBar;
    increment::Int = 1,
    description::Union{String,Nothing} = nothing,
)
    pb.current = min(pb.current + increment, pb.total)
    if description !== nothing
        pb.description = description
    end
    render_progress(pb)
end

"""
    finish_progress!(pb::ProgressBar; description=nothing)

Complete the progress bar and move to next line. Requires `using Term`.

# Arguments
- `pb`: ProgressBar instance
- `description`: Optional final description
"""
function finish_progress!(pb::ProgressBar; description::Union{String,Nothing} = nothing)
    pb.current = pb.total
    pb.finished = true
    if description !== nothing
        pb.description = description
    end
    render_progress(pb)
    println()  # Move to next line
end

"""
    format_duration(seconds)

Format seconds into MM:SS or HH:MM:SS string.
"""
function format_duration(seconds::Real)
    if !isfinite(seconds) || seconds < 0
        return "--:--"
    end

    total_secs = round(Int, seconds)
    hours = total_secs ÷ 3600
    mins = (total_secs % 3600) ÷ 60
    secs = total_secs % 60

    if hours > 0
        return @sprintf("%02d:%02d:%02d", hours, mins, secs)
    else
        return @sprintf("%02d:%02d", mins, secs)
    end
end

"""
    with_progress(f, collection; description="Progress", config=DEFAULT_CONFIG)

Execute a function over a collection with automatic progress tracking. Requires
`using Term`.

# Arguments
- `f`: Function to apply to each element (receives element and index)
- `collection`: Iterable collection
- `description`: Progress bar description
- `config`: DisplayConfig instance

# Returns
- Vector of results from applying f to each element

# Examples
```julia
using Term
results = with_progress(models, description="Testing models") do model, i
    run_test(model)
end
```
"""
function with_progress(
    f::Function,
    collection;
    description::String = "Progress",
    config::DisplayConfig = DEFAULT_CONFIG,
)
    n = length(collection)
    pb = progress_bar(n; description = description, config = config)
    results = []

    for (i, item) in enumerate(collection)
        result = f(item, i)
        push!(results, result)
        update_progress!(pb)
    end

    finish_progress!(pb)
    return results
end

# ==============================================================================
# Spinner for Indeterminate Progress (state in core, rendering in the Term extension)
# ==============================================================================

const SPINNER_FRAMES = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]

"""
    Spinner

Mutable struct for displaying an animated spinner for indeterminate progress.
Rendering requires `using Term`.

# Fields
- `message`: Message to display
- `frame`: Current animation frame
- `start_time`: When spinner started
- `config`: DisplayConfig for styling
- `active`: Whether spinner is active
"""
mutable struct Spinner
    message::String
    frame::Int
    start_time::Float64
    config::DisplayConfig
    active::Bool
end

"""
    spinner(message; config=DEFAULT_CONFIG)

Create and start a new spinner for indeterminate progress. Requires `using Term`.

# Arguments
- `message`: Message to display next to spinner
- `config`: DisplayConfig instance

# Returns
- `Spinner` instance

# Examples
```julia
using Term
sp = spinner("Loading data...")
# do work
stop_spinner!(sp; message = "Done!")
```
"""
function spinner(message::String; config::DisplayConfig = DEFAULT_CONFIG)
    sp = Spinner(message, 1, time(), config, true)
    render_spinner(sp)
    return sp
end

"""
    spin!(sp::Spinner; message=nothing)

Advance the spinner animation by one frame. Requires `using Term`.

# Arguments
- `sp`: Spinner instance
- `message`: Optionally update the message
"""
function spin!(sp::Spinner; message::Union{String,Nothing} = nothing)
    if !sp.active
        return
    end
    sp.frame = (sp.frame % length(SPINNER_FRAMES)) + 1
    if message !== nothing
        sp.message = message
    end
    render_spinner(sp)
end

"""
    stop_spinner!(sp::Spinner; message=nothing, success=true)

Stop the spinner and display final status. Requires `using Term`.

# Arguments
- `sp`: Spinner instance
- `message`: Final message (optional, uses current if not provided)
- `success`: Whether operation succeeded (affects icon)
"""
stop_spinner!(args...; kwargs...) = _display_requires(stop_spinner!, "Term", args)

"""
    with_spinner(f, message; config=DEFAULT_CONFIG)

Execute a function while displaying a spinner. Requires `using Term`.

# Arguments
- `f`: Function to execute (no arguments)
- `message`: Message to display
- `config`: DisplayConfig instance

# Returns
- Result of function f

# Examples
```julia
using Term
result = with_spinner("Computing gradients...") do
    expensive_computation()
end
```
"""
function with_spinner(f::Function, message::String; config::DisplayConfig = DEFAULT_CONFIG)
    sp = spinner(message; config = config)

    # Start a background task to animate the spinner
    spinner_task = @async begin
        while sp.active
            sleep(0.1)
            spin!(sp)
        end
    end

    try
        result = f()
        stop_spinner!(sp; success = true)
        return result
    catch e
        stop_spinner!(sp; message = "Failed: $message", success = false)
        rethrow(e)
    end
end

# ==============================================================================
# Critical Points as Markdown (no display package needed)
# ==============================================================================

"""
    format_critical_points_markdown(points, objective_values, param_names;
                                     n=5, true_params=nothing)

Format top N critical points as a Markdown table string.

# Returns
- String containing Markdown-formatted table
"""
function format_critical_points_markdown(
    points::Vector,
    objective_values::Vector{<:Real},
    param_names::Vector;
    n::Int = 5,
    true_params::Union{Vector{<:Real},Nothing} = nothing,
    precision::Int = 4,
)
    io = IOBuffer()

    sorted_indices = sortperm(objective_values)
    n_display = min(n, length(sorted_indices))

    if n_display == 0
        return "No critical points found."
    end

    header_parts = ["Rank", "Objective"]
    for name in param_names
        push!(header_parts, string(name))
    end
    if true_params !== nothing
        push!(header_parts, "Rel. Error")
    end

    println(io, "| ", join(header_parts, " | "), " |")
    println(io, "| ", join(fill("---:", length(header_parts)), " | "), " |")

    for (row, idx) in enumerate(sorted_indices[1:n_display])
        point = points[idx]
        obj_val = objective_values[idx]

        row_parts = String[]
        push!(row_parts, string(row))

        if abs(obj_val) < 1e-3 || abs(obj_val) > 1e4
            push!(row_parts, @sprintf("%.4e", obj_val))
        else
            push!(row_parts, string(round(obj_val, digits = precision)))
        end

        for p in point
            push!(row_parts, string(round(p, digits = precision)))
        end

        if true_params !== nothing
            rel_error = norm(point .- true_params) / norm(true_params)
            push!(row_parts, @sprintf("%.2f%%", rel_error * 100))
        end

        println(io, "| ", join(row_parts, " | "), " |")
    end

    return String(take!(io))
end
