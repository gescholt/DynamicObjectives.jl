# Examples

This directory contains example scripts demonstrating various features of Dynamic_objectives.

## Basic Workflow

**File:** `basic_workflow.jl`

A simple, practical example showing a typical parameter estimation workflow:

1. Define an ODE model
2. Set up problem parameters
3. Generate reference data
4. Create objective function
5. Evaluate at different parameters
6. Visualize comparisons
7. Next steps for optimization

This is the best starting point for new users. Run with:

```bash
julia --project=. examples/basic_workflow.jl
```

## globtim Integration Test

**File:** `globtim_integration/run_single_model.jl`

Integration test demonstrating Dynamic_objectives working with globtimcore:

1. Activates globtimcore environment
2. Loads a 2D Lotka-Volterra model
3. Creates objective function
4. Generates reference data
5. Tests at multiple parameter sets
6. Compares predictions visually

Uses the display infrastructure throughout to show model information, parameters, time series, and comparisons. Run with:

```bash
julia examples/globtim_integration/run_single_model.jl
```

**Note:** Requires globtimcore to be installed in a sibling directory (`../../../globtimcore` relative to the script).

## Display Infrastructure Demo

**File:** `display_demo.jl`

A comprehensive demonstration of the display infrastructure, showing:

1. **Model Summary Display** - Pretty formatted ODE model information
2. **Parameter Display** - Tables showing parameter values and comparisons
3. **Time Series Data Display** - Data tables with plots
4. **Time Series Comparison** - Side-by-side visualization of different datasets
5. **Error Metrics Display** - Formatted error/distance metric tables
6. **Optimization Result Display** - Summary of optimization outcomes
7. **Optimization Progress Display** - Real-time iteration updates
8. **Custom Display Configuration** - Customizable output styles
9. **Complete Workflow Example** - End-to-end parameter estimation demo

### Running the Demo

```bash
julia --project=. examples/display_demo.jl
```

### Display Configuration

The display system can be customized:

```julia
using Dynamic_objectives

# Use default configuration
display_model_summary(model)

# Create custom configuration
config = DisplayConfig(
    use_color = false,        # Disable colors
    table_backend = :ascii,   # Use ASCII tables
    plot_width = 80,          # Wider plots
    plot_height = 25,         # Taller plots
    precision = 6             # More decimal places
)

# Use custom config
display_parameters(params, values, config=config)

# Or update global defaults
set_display_config!(use_color=false, precision=6)
```

## Display Functions Reference

### Model Information

```julia
display_model_summary(model)
```

Shows ODE system details: states, parameters, and equations.

### Parameters

```julia
# Single parameter set
display_parameters([:α, :β, :γ], [0.1, 0.2, 0.3])

# Multiple parameter sets for comparison
display_parameters(
    [:α, :β, :γ],
    [p_true p_test p_optimal],
    labels=["True", "Test", "Optimal"]
)
```

### Time Series

```julia
# Display data with plots
display_time_series(data_dict, show_plot=true)

# Display specific variables only
display_time_series(data_dict, variables=[:x, :y])

# Compare two datasets
display_comparison(data_true, data_test, labels=["Reference", "Prediction"])
```

### Error Metrics

```julia
# From dictionary
errors = Dict("L1" => 0.123, "L2" => 0.045)
display_error_metrics(errors)

# From vector with labels
display_error_metrics([0.123, 0.045], labels=["L1", "L2"])
```

### Optimization Results

```julia
# Complete optimization summary
result = (
    params = [0.1, 0.2, 0.3],
    error = 0.0012,
    iterations = 150,
    converged = true,
    time_elapsed = 2.5
)
display_optimization_result(result)

# Real-time progress updates (in optimization callback)
display_optimization_progress(iteration, current_params, current_error)
```

## Integration with Optimization

The display functions are designed to work seamlessly with optimization workflows:

```julia
using Dynamic_objectives
using Optimization

# Define your objective
model, params, states, outputs = define_daisy_ex3_model_4D()
error_func = make_error_distance(model, outputs, ic, p_true, time_interval, 25)

# Display initial setup
display_model_summary(model)
display_parameters(param_names, initial_guess, labels=["Initial"])

# Run optimization with progress display
function callback(p, error_val, iteration)
    if iteration % 10 == 0
        display_optimization_progress(iteration, p, error_val)
    end
    return false  # Continue optimization
end

# ... run optimization ...

# Display final results
display_optimization_result((
    params = optimal_params,
    error = final_error,
    iterations = total_iterations,
    converged = success,
    time_elapsed = elapsed_time
))

# Compare solutions
display_comparison(data_true, data_optimal, labels=["Target", "Fitted"])
```

## Future Enhancements

The display infrastructure is designed to be extended with:

- **Interactive refinement** - User prompts to refine search domains
- **Live dashboards** - Real-time monitoring of multiple optimization runs
- **Export capabilities** - Generate HTML/LaTeX reports
- **Plotting integrations** - Makie.jl for publication-quality figures
- **Web interfaces** - Pluto.jl notebooks with interactive controls
