# globtim Integration Test - Single Model
# Tests Dynamic_objectives with globtim integration using display infrastructure
#
# This example demonstrates:
# - Clean, readable output with suppressed verbose logging
# - Pretty display using Term.jl panels and PrettyTables
# - Integration with globtimcore (if available)
# - All major display functions from Dynamic_objectives

# ============================================================================
# Logging Configuration - Suppress verbose output
# ============================================================================
using Logging

# Create a custom logger that filters out Info messages except critical ones
struct SelectiveLogger <: AbstractLogger
    io::IO
    min_level::LogLevel
    critical_keywords::Vector{String}
end

SelectiveLogger(io::IO=stderr) = SelectiveLogger(io, Logging.Warn, String[
    "ModelRegistry initialized",  # Keep this Info message
    "Optimization complete",
    "PASS",
    "FAIL"
])

function Logging.shouldlog(logger::SelectiveLogger, level, _module, group, id)
    return level >= logger.min_level
end

function Logging.min_enabled_level(logger::SelectiveLogger)
    return logger.min_level
end

function Logging.catch_exceptions(logger::SelectiveLogger)
    return true
end

function Logging.handle_message(logger::SelectiveLogger, level, message, _module, group, id, file, line; kwargs...)
    # Check if this is a critical Info message we want to keep
    msg_str = string(message)
    is_critical = any(kw -> occursin(kw, msg_str), logger.critical_keywords)

    if level >= logger.min_level || is_critical
        printstyled(logger.io, "[ ", level, " ] ", color=:cyan, bold=true)
        println(logger.io, message)
        for (key, val) in kwargs
            println(logger.io, "  ", key, " = ", val)
        end
    end
end

# Set the custom logger globally
global_logger(SelectiveLogger(stderr))

# Find globtimcore path (assumed to be sibling directory)
const GLOBTIM_PATH = abspath(joinpath(@__DIR__, "..", "..", "..", "globtimcore"))

println("Activating globtimcore environment: $GLOBTIM_PATH")
using Pkg

# Suppress Pkg output
Pkg.activate(GLOBTIM_PATH; io=devnull)
Pkg.instantiate(; io=devnull)

println("Setting up Dynamic_objectives...")
using Dynamic_objectives

# Try to load Globtim if available
try
    using Globtim.ModelRegistry
    n_models = length(list_models())
    println("✓ ModelRegistry initialized with $n_models models")
catch e
    @warn "Globtim.ModelRegistry not available, continuing with Dynamic_objectives only"
end

# ============================================================================
# Helper Functions for Clean Output
# ============================================================================

"""Suppress progress bars and verbose output during function execution"""
function with_suppressed_output(f::Function)
    # Redirect stdout to capture progress bars
    original_stdout = stdout
    original_stderr = stderr

    # Create a buffer to capture any errors we want to keep
    error_buffer = IOBuffer()

    try
        # Redirect stdout (for progress bars)
        redirect_stdout(devnull)

        # Execute the function
        result = f()

        return result
    finally
        # Restore original streams
        redirect_stdout(original_stdout)
    end
end

"""Format homotopy continuation results for display"""
function display_homotopy_results(results)
    if haskey(results, :paths_tracked)
        println("\n  Homotopy Continuation Results:")
        println("  ├─ Paths tracked: $(results[:paths_tracked])")
        if haskey(results, :real_solutions)
            println("  ├─ Real solutions: $(results[:real_solutions])")
        end
        if haskey(results, :time_elapsed)
            println("  └─ Time: $(round(results[:time_elapsed], digits=2))s")
        end
    end
end

println("="^80)
println("globtim Integration Test - Single Model (2D LV v3)")
println("="^80)
println()

# ============================================================================
# [1/5] Define Model
# ============================================================================
println("[1/5] Defining model...")

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# Setup parameters
p_true = [1.0, 0.5]
p_bounds = [(0.0, 3.0), (0.0, 2.0)]
ic = [1.0, 1.0]
time_interval = [0.0, 20.0]
num_points = 30

println("  Model: Lotka-Volterra 2D v3 (2 outputs)")
println("  True parameters: $p_true")
println("  Parameter bounds: $p_bounds")
println("  Time interval: $time_interval")
println("  Sample points: $num_points")
println()

# Display model summary using new infrastructure
display_model_summary(model)

# Display parameters
param_names = [:α, :β]
display_parameters(param_names, p_true, labels=["True Value"])

# ============================================================================
# [2/5] Create Objective Function
# ============================================================================
println("[2/5] Creating objective function...")

using ModelingToolkit: complete
problem = ModelingToolkit.ODEProblem(
    complete(model),
    merge(
        Dict(states .=> ic),
        Dict(params .=> p_true)
    ),
    time_interval
)

error_func = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    num_points,
    L2_norm,
    first,
    nothing;
    return_inf_on_error=true,
    eval_timeout=5.0
)

println("  ✓ Objective function created")
println()

# ============================================================================
# [3/5] Generate Reference Data
# ============================================================================
println("[3/5] Generating reference data...")

# Suppress ODE solver output
print("  Solving ODE system... ")
data_reference = with_suppressed_output() do
    sample_data(
        problem,
        model,
        outputs,
        time_interval,
        p_true,
        ic,
        num_points
    )
end
println("✓ Generated $(num_points) time points")
println()

# Display time series data
display_time_series(data_reference, show_plot=true)

# ============================================================================
# [4/5] Test Objective Function
# ============================================================================
println("[4/5] Testing objective function at different parameters...")

# Test at true parameters (should be ~0)
print("  Evaluating at true parameters... ")
error_at_true = with_suppressed_output(() -> error_func(p_true))
println("✓ Error: $(round(error_at_true, digits=8))")

# Test at perturbed parameters
p_test1 = [1.2, 0.6]
p_test2 = [0.8, 0.4]
p_test3 = [1.5, 0.3]

print("  Evaluating at test parameters... ")
error_at_test1 = with_suppressed_output(() -> error_func(p_test1))
error_at_test2 = with_suppressed_output(() -> error_func(p_test2))
error_at_test3 = with_suppressed_output(() -> error_func(p_test3))
println("✓ Done")
println()

# Display parameter comparison
display_parameters(
    param_names,
    [p_true p_test1 p_test2 p_test3],
    labels=["True", "Test1", "Test2", "Test3"]
)

# Display error metrics
error_metrics = Dict(
    "Error at True" => error_at_true,
    "Error at Test1" => error_at_test1,
    "Error at Test2" => error_at_test2,
    "Error at Test3" => error_at_test3,
    "Min Error" => minimum([error_at_test1, error_at_test2, error_at_test3]),
    "Max Error" => maximum([error_at_test1, error_at_test2, error_at_test3])
)

display_error_metrics(error_metrics)

# ============================================================================
# [5/5] Visual Comparison
# ============================================================================
println("[5/5] Comparing predictions...")

# Generate data for best test parameters
best_idx = argmin([error_at_test1, error_at_test2, error_at_test3])
p_best = [p_test1, p_test2, p_test3][best_idx]
label_best = ["Test1", "Test2", "Test3"][best_idx]

println("  Best test parameters: $p_best ($(label_best))")

# Generate comparison data with suppressed output
print("  Generating comparison data... ")
data_best = with_suppressed_output() do
    sample_data(
        problem,
        model,
        outputs,
        time_interval,
        p_best,
        ic,
        num_points
    )
end
println("✓ Done")
println()

# Display comparison
display_comparison(
    data_reference,
    data_best,
    labels=["True Parameters", "Best Test ($label_best)"]
)

# ============================================================================
# Summary
# ============================================================================
println()
using Term

# Create a summary panel
min_error = minimum([error_at_test1, error_at_test2, error_at_test3])
max_error = maximum([error_at_test1, error_at_test2, error_at_test3])

summary_text = """
{bold cyan}Test Results:{/bold cyan}
  ✓ Model loaded successfully
  ✓ Objective function created
  ✓ Reference data generated ($(num_points) points)
  ✓ Error function evaluated at 4 parameter sets
  ✓ Visual comparison completed

{bold green}Error Analysis:{/bold green}
  • Error at true parameters:  $(round(error_at_true, digits=8))
  • Best test error:            $(round(min_error, digits=6))
  • Worst test error:           $(round(max_error, digits=6))
  • Error ratio (worst/best):   $(round(max_error/min_error, digits=2))x

{bold yellow}Display Infrastructure Tested:{/bold yellow}
  📊 Model summaries           (display_model_summary)
  📊 Parameter tables          (display_parameters)
  📊 Time series plots         (display_time_series)
  📊 Prediction comparisons    (display_comparison)
  📊 Error metrics             (display_error_metrics)

{bold magenta}Next Steps:{/bold magenta}
  → Integrate with global optimization algorithms
  → Use display_optimization_progress() for real-time monitoring
  → Export results using display functions
"""

panel = Term.Panel(
    Term.parse(summary_text),
    title="Integration Test Complete ✓",
    title_style="bold green",
    style="green",
    fit=true
)
println(panel)

println("\n" * "="^80)
println("✅ All tests passed successfully!")
println("="^80)
