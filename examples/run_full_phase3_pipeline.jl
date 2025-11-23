"""
Full Phase 3 Pipeline Example: globtimcore + globtimpostprocessing

This script demonstrates the complete 2-stage optimization pipeline:
1. Stage 1 (globtimcore): Find raw critical points using polynomial approximation
2. Stage 2 (globtimpostprocessing): Refine critical points using local optimization

# Requirements

This example requires two external packages in separate environments:
- globtimcore: For polynomial approximation and critical point finding
- globtimpostprocessing: For local refinement (implemented in Dynamic_objectives)

# Usage

```bash
# Activate Dynamic_objectives environment
julia --project=. examples/run_full_phase3_pipeline.jl
```

# Expected Results

For Lotka-Volterra 2D (EASY model):
- Stage 1: ~20% parameter recovery error (grid approximation)
- Stage 2: <1% parameter recovery error (after refinement)
- Improvement: 10-100x in objective value
- Total time: <1 minute

# Integration Patterns

This example uses **Pattern 1: Full Pipeline**
See REFINEMENT_INTEGRATION_GUIDE.md for alternative patterns.
"""

using Pkg
using Dynamic_objectives

println("="^80)
println("Phase 3 Full Pipeline Example")
println("="^80)
println()

# ============================================================================
# Configuration
# ============================================================================

# Model selection
MODEL_NAME = "LotkaVolterra_2D"
OUTPUT_DIR = joinpath(@__DIR__, "..", "results", "phase3_pipeline_$(MODEL_NAME)")

# Stage 1: globtimcore configuration
MAX_DEGREE = 18
GRID_SIZE = 8  # GN parameter

# Stage 2: Refinement configuration
REFINEMENT_TIMEOUT = 60.0  # seconds per point

println("Configuration:")
println("  Model: $MODEL_NAME")
println("  Output directory: $OUTPUT_DIR")
println("  Max degree: $MAX_DEGREE")
println("  Grid size: $GRID_SIZE")
println("  Refinement timeout: $(REFINEMENT_TIMEOUT)s")
println()

# ============================================================================
# Step 1: Create Dynamic_objectives model
# ============================================================================

println("Step 1: Creating model and objective function")
println("-"^80)

# Define Lotka-Volterra 2D model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# True parameters and initial conditions
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

# Time configuration
time_interval = [0.0, 20.0]
numpoints = 30

println("  Model: Lotka-Volterra 2D")
println("  True parameters: $p_true")
println("  Parameter bounds: $bounds")
println("  Time points: $numpoints")
println()

# Create objective function (1-argument signature)
objective = make_error_distance(
    model, outputs, ic, p_true,
    time_interval, numpoints,
    L2_norm,     # distance function
    first,       # aggregate function
    nothing;     # no noise
    return_inf_on_error = true,
    eval_timeout = 10.0  # timeout for stiff ODEs
)

println("✓ Objective function created")
println()

# ============================================================================
# Step 2: Stage 1 - Run globtimcore (find raw critical points)
# ============================================================================

println("Step 2: Stage 1 - Finding raw critical points with globtimcore")
println("-"^80)

# Check if globtimcore is available
globtimcore_available = false
try
    # Try to detect globtimcore installation
    globtimcore_path = get(ENV, "GLOBTIMCORE_PATH", "")
    if !isempty(globtimcore_path) && isdir(globtimcore_path)
        globtimcore_available = true
        println("  Found globtimcore at: $globtimcore_path")
    end
catch
end

if globtimcore_available
    println("  Running globtimcore polynomial approximation...")
    println()

    # OPTION A: Use external globtimcore package
    # This requires switching environments - see REFINEMENT_INTEGRATION_GUIDE.md
    println("  NOTE: This example uses standalone globtimcore integration.")
    println("  For full polynomial approximation, use the standalone script:")
    println("  examples/globtim_integration/run_single_model.jl")
    println()

    # For demonstration, use the simplified grid search
    mkpath(OUTPUT_DIR)

    result_stage1 = run_globtim_optimization(
        create_globtim_objective(model, outputs, ic, p_true, time_interval, numpoints),
        bounds,
        p_true;
        GN = GRID_SIZE,
        degree_range = 4:8,
        output_dir = OUTPUT_DIR,
        model_name = MODEL_NAME,
        show_progress = true,
        save_results = true
    )

    println()
    println("✓ Stage 1 complete")
    println("  Critical points found: $(result_stage1[:total_critical_points])")
    println("  Best objective (raw): $(round(result_stage1[:best_objective], digits=6))")
    println("  Recovery error (raw): $(round(100*result_stage1[:recovery_error], digits=2))%")
    println()

else
    # OPTION B: Use mock data for demonstration
    println("  ⚠ globtimcore not detected - using mock critical points for demo")
    println()

    # Create mock raw critical points
    mkpath(OUTPUT_DIR)

    # Generate grid points
    n_points = GRID_SIZE^2
    grid_x1 = range(bounds[1][1], bounds[1][2], length=GRID_SIZE)
    grid_x2 = range(bounds[2][1], bounds[2][2], length=GRID_SIZE)

    mock_points = []
    for x1 in grid_x1
        for x2 in grid_x2
            point = [x1, x2]
            value = objective(point)
            push!(mock_points, (point, value))
        end
    end

    # Sort by objective value and take best candidates
    sort!(mock_points, by = x -> x[2])
    best_points = mock_points[1:min(10, length(mock_points))]

    # Write to CSV (mock degree 8)
    using CSV, DataFrames
    mock_df = DataFrame(
        point_id = 1:length(best_points),
        p1 = [p[1][1] for p in best_points],
        p2 = [p[1][2] for p in best_points],
        objective_value = [p[2] for p in best_points]
    )

    CSV.write(joinpath(OUTPUT_DIR, "critical_points_deg_8.csv"), mock_df)

    println("✓ Mock critical points generated")
    println("  Points: $(length(best_points))")
    println("  Best objective (raw): $(round(minimum(p[2] for p in best_points), digits=6))")
    println()
end

# ============================================================================
# Step 3: Stage 2 - Refine critical points with globtimpostprocessing
# ============================================================================

println("Step 3: Stage 2 - Refining critical points with local optimization")
println("-"^80)

# Configure refinement
config = ode_refinement_config(
    max_time_per_point = REFINEMENT_TIMEOUT,
    optimization_method = :nelder_mead,
    max_iterations = 1000,
    x_tol = 1e-6,
    f_tol = 1e-6,
    verbose = true
)

println("  Refinement configuration:")
println("    Method: Nelder-Mead (gradient-free)")
println("    Timeout: $(config.max_time_per_point)s per point")
println("    Max iterations: $(config.max_iterations)")
println()

# Run refinement
result_stage2 = refine_experiment_results(
    OUTPUT_DIR,
    objective,
    config
)

println()
println("✓ Stage 2 complete")
println()

# ============================================================================
# Step 4: Display results
# ============================================================================

println("="^80)
println("Final Results: 2-Stage Pipeline")
println("="^80)
println()

println("Stage 1 (globtimcore - raw critical points):")
println("  Critical points found: $(result_stage2[:n_raw])")
println("  Best objective (raw): $(round(result_stage2[:best_raw_value], digits=6))")
println()

println("Stage 2 (globtimpostprocessing - refinement):")
println("  Converged: $(result_stage2[:n_converged])/$(result_stage2[:n_raw]) ($(round(100*result_stage2[:n_converged]/result_stage2[:n_raw], digits=1))%)")
println("  Best objective (refined): $(round(result_stage2[:best_refined_value], digits=6))")
println("  Mean improvement: $(round(result_stage2[:mean_improvement], digits=2))x")
println()

println("Overall improvement:")
println("  Objective improvement: $(round(result_stage2[:best_raw_value]/result_stage2[:best_refined_value], digits=2))x")
println()

# Display best refined point
best_idx = result_stage2[:best_refined_idx]
best_params = result_stage2[:refined_points][best_idx]
recovery_error = norm(best_params .- p_true) / norm(p_true)

println("Best refined parameters:")
println("  p_true:    $p_true")
println("  p_refined: $(round.(best_params, digits=4))")
println("  Recovery error: $(round(100*recovery_error, digits=2))%")
println()

println("Output files:")
println("  Directory: $OUTPUT_DIR")
for file in readdir(OUTPUT_DIR)
    println("    - $file")
end
println()

println("="^80)
println("Pipeline Complete!")
println("="^80)

# Success criteria
if recovery_error < 0.01
    println("✓ SUCCESS: Recovery error < 1%")
elseif recovery_error < 0.05
    println("⚠ PARTIAL SUCCESS: Recovery error < 5%")
else
    println("✗ NEEDS IMPROVEMENT: Recovery error > 5%")
    println("  Consider:")
    println("    - Increasing grid size (GN)")
    println("    - Increasing polynomial degree")
    println("    - Adjusting refinement tolerances")
end
