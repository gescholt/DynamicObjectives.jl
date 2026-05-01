# Display Infrastructure Demo
# This example demonstrates the pretty display capabilities of Dynamic_objectives

using Dynamic_objectives

println("\n" * "="^70)
println("Dynamic_objectives Display Infrastructure Demo")
println("="^70)

# ============================================================================
# 1. Model Summary Display
# ============================================================================
println("\n📊 1. MODEL SUMMARY DISPLAY\n")

model, params, states, outputs = define_daisy_ex3_model_4D()
display_model_summary(model)

# ============================================================================
# 2. Parameter Display
# ============================================================================
println("\n📊 2. PARAMETER DISPLAY\n")

# Define parameter values
param_names = [:α, :β, :γ, :δ]
p_true = [0.1, 0.2, 0.3, 0.4]
p_test = [0.15, 0.25, 0.35, 0.45]
p_optimal = [0.102, 0.198, 0.305, 0.397]

# Display single parameter set
display_parameters(param_names, p_true, labels = ["True Value"])

# Display multiple parameter sets for comparison
display_parameters(
    param_names,
    [p_true p_test p_optimal],
    labels = ["True", "Initial Guess", "Optimized"],
)

# ============================================================================
# 3. Time Series Data Display
# ============================================================================
println("\n📊 3. TIME SERIES DATA DISPLAY\n")

# Generate synthetic time series data
ic = [1.0, 2.0, 1.0, 1.0]
time_interval = [0.0, 10.0]
num_points = 25

# Create ODE problem and generate data
using Dynamic_objectives: ODEProblem
problem =
    ODEProblem(model, merge(Dict(states .=> ic), Dict(params .=> p_true)), time_interval)

data_true = sample_data(problem, model, outputs, time_interval, p_true, ic, num_points)

# Display time series with plots
display_time_series(data_true, show_plot = true)

# ============================================================================
# 4. Time Series Comparison
# ============================================================================
println("\n📊 4. TIME SERIES COMPARISON\n")

# Generate data with test parameters
data_test = sample_data(problem, model, outputs, time_interval, p_test, ic, num_points)

# Compare two datasets
display_comparison(data_true, data_test, labels = ["True Parameters", "Test Parameters"])

# ============================================================================
# 5. Error Metrics Display
# ============================================================================
println("\n📊 5. ERROR METRICS DISPLAY\n")

# Compute various error metrics
error_L1 = L1_norm(data_true[outputs[1].lhs], data_test[outputs[1].lhs])
error_L2 = L2_norm(data_true[outputs[1].lhs], data_test[outputs[1].lhs])
error_log_L2 = log_L2_norm(data_true[outputs[1].lhs], data_test[outputs[1].lhs])

error_dict = Dict(
    "L1 Norm" => error_L1,
    "L2 Norm" => error_L2,
    "Log-L2 Norm" => error_log_L2,
    "Relative Error" => error_L2 / maximum(abs.(data_true[outputs[1].lhs])),
)

display_error_metrics(error_dict)

# ============================================================================
# 6. Optimization Result Display
# ============================================================================
println("\n📊 6. OPTIMIZATION RESULT DISPLAY\n")

# Simulate an optimization result
optimization_result = (
    params = p_optimal,
    error = 0.0012,
    iterations = 150,
    converged = true,
    time_elapsed = 2.5,
)

display_optimization_result(optimization_result)

# ============================================================================
# 7. Optimization Progress Display (Interactive)
# ============================================================================
println("\n📊 7. OPTIMIZATION PROGRESS DISPLAY (Simulated)\n")

println("Simulating optimization progress updates...\n")

# Simulate 10 iterations
for iter in 1:10
    # Simulate parameter updates
    p_current = p_true .+ (p_test .- p_true) .* (1.0 - iter/10.0)
    error_current = 0.1 * exp(-iter/2)

    display_optimization_progress(iter, p_current, error_current)
    sleep(0.2)  # Simulate computation time
end

println("\n")  # New line after progress updates

# ============================================================================
# 8. Custom Display Configuration
# ============================================================================
println("\n📊 8. CUSTOM DISPLAY CONFIGURATION\n")

println("Example with custom configuration (ASCII, no color, higher precision):\n")

# Create custom config
custom_config = DisplayConfig(use_color = false, table_backend = :ascii, precision = 6)

# Display with custom config
display_parameters(
    param_names,
    [p_true p_optimal],
    labels = ["True", "Optimized"],
    config = custom_config,
)

# ============================================================================
# 9. Real-World Example: Parameter Estimation
# ============================================================================
println("\n📊 9. COMPLETE PARAMETER ESTIMATION WORKFLOW\n")

println("Demonstrating a complete parameter estimation workflow...\n")

# Create error function
error_func = make_error_distance(model, outputs, ic, p_true, time_interval, num_points)

# Test at true parameters (should be ~0)
error_at_true = error_func(p_true)
println("Error at true parameters: $(round(error_at_true, digits=8))")

# Test at perturbed parameters
error_at_test = error_func(p_test)
println("Error at test parameters: $(round(error_at_test, digits=6))\n")

# Display comprehensive results
display_parameters(param_names, [p_true p_test], labels = ["True", "Test"])

error_comparison = Dict(
    "Error at True Params" => error_at_true,
    "Error at Test Params" => error_at_test,
    "Ratio" => error_at_test / (error_at_true + 1e-10),
)

display_error_metrics(error_comparison)

# ============================================================================
# Summary
# ============================================================================
println("\n" * "="^70)
println("Demo Complete!")
println("="^70)
println("""
The display infrastructure provides:
  ✓ Model summaries with styled panels
  ✓ Parameter tables with flexible formatting
  ✓ Time series plots and data tables
  ✓ Side-by-side comparisons
  ✓ Error metric displays
  ✓ Optimization result summaries
  ✓ Real-time progress updates
  ✓ Customizable configurations

Future enhancements:
  → Interactive refinement controls
  → Live optimization monitoring
  → Export to various formats (HTML, LaTeX, etc.)
  → Integration with optimization libraries
""")
println("="^70)
