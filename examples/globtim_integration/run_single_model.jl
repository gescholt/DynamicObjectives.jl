# globtim Integration Test - Single Model
# Tests Dynamic_objectives with globtim integration using display infrastructure

# Find globtimcore path (assumed to be sibling directory)
const GLOBTIM_PATH = abspath(joinpath(@__DIR__, "..", "..", "..", "globtimcore"))

println("Activating globtimcore environment: $GLOBTIM_PATH")
using Pkg
Pkg.activate(GLOBTIM_PATH)
Pkg.instantiate()

println("Setting up Dynamic_objectives...")
using Dynamic_objectives

# Try to load Globtim if available
try
    using Globtim.ModelRegistry
    @info "ModelRegistry initialized with $(length(list_models())) models"
catch e
    @warn "Globtim.ModelRegistry not available, continuing with Dynamic_objectives only"
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

data_reference = sample_data(
    problem,
    model,
    outputs,
    time_interval,
    p_true,
    ic,
    num_points
)

println("  ✓ Generated $(num_points) time points")
println()

# Display time series data
display_time_series(data_reference, show_plot=true)

# ============================================================================
# [4/5] Test Objective Function
# ============================================================================
println("[4/5] Testing objective function at different parameters...")

# Test at true parameters (should be ~0)
error_at_true = error_func(p_true)
println("  Error at true params: $(round(error_at_true, digits=8))")

# Test at perturbed parameters
p_test1 = [1.2, 0.6]
p_test2 = [0.8, 0.4]
p_test3 = [1.5, 0.3]

error_at_test1 = error_func(p_test1)
error_at_test2 = error_func(p_test2)
error_at_test3 = error_func(p_test3)

println("  Error at test1 params: $(round(error_at_test1, digits=6))")
println("  Error at test2 params: $(round(error_at_test2, digits=6))")
println("  Error at test3 params: $(round(error_at_test3, digits=6))")
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
println()

data_best = sample_data(
    problem,
    model,
    outputs,
    time_interval,
    p_best,
    ic,
    num_points
)

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
println("="^80)
println("Integration Test Complete!")
println("="^80)
println("""
Summary:
  ✓ Model loaded successfully
  ✓ Objective function created
  ✓ Reference data generated
  ✓ Error function evaluated at multiple parameter sets
  ✓ Visual comparison completed

Results:
  - Error at true parameters: $(round(error_at_true, digits=8))
  - Best test error: $(round(minimum([error_at_test1, error_at_test2, error_at_test3]), digits=6))
  - Worst test error: $(round(maximum([error_at_test1, error_at_test2, error_at_test3]), digits=6))

The display infrastructure successfully visualized:
  📊 Model summaries
  📊 Parameter tables
  📊 Time series plots
  📊 Data comparisons
  📊 Error metrics

Next steps:
  - Integrate with global optimization algorithms
  - Use display_optimization_progress() for real-time monitoring
  - Export results to reports using display functions
""")
println("="^80)
