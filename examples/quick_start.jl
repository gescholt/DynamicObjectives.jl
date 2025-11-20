#!/usr/bin/env julia
# Quick Start: Your First Test with globtim
# Usage: julia --project=. examples/quick_start.jl

using Dynamic_objectives
using Printf

println("="^80)
println("Quick Start: Testing Dynamic_objectives with globtim")
println("="^80)

# ==============================================================================
# STEP 1: Choose a Model (starting with easiest)
# ==============================================================================

println("\nSTEP 1: Setting up model")
println("-"^80)

model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]

println("Model: LV 2D v3 with 2 outputs (EASIEST)")
println("  Parameters: 2")
println("  Outputs: 2 (full observability)")
println("  True params: $p_true")
println("  Bounds: $bounds")

# ==============================================================================
# STEP 2: Create Objective Function
# ==============================================================================

println("\nSTEP 2: Creating objective function")
println("-"^80)

error_func = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    [0.0, 20.0],  # time interval
    30,           # number of time points
    L2_norm,      # distance metric
    first,        # aggregation function
    nothing;      # no noise
    return_inf_on_error = true,
    eval_timeout = nothing
)

println("✓ Objective function created")

# ==============================================================================
# STEP 3: Verify Setup
# ==============================================================================

println("\nSTEP 3: Verifying setup")
println("-"^80)

error_at_true = error_func(p_true)
@printf "Error at true parameters: %.6e\n" error_at_true

if error_at_true < 1e-6
    println("✓ Setup verified (error ≈ 0)")
else
    error("Setup failed - error should be near 0 at true parameters")
end

# Test a few points
println("\nTesting objective function at different points:")
test_points = [
    [1.0, 0.5],    # true
    [1.1, 0.5],    # perturb p1
    [1.0, 0.6],    # perturb p2
    [0.5, 0.5],    # far from true
    [2.0, 1.0],    # far from true
]

for p in test_points
    err = error_func(p)
    @printf "  p = [%.1f, %.1f]  →  error = %.6e\n" p[1] p[2] err
end

println("\n✓ Objective function responds correctly to parameter changes")

# ==============================================================================
# STEP 4: INTEGRATE WITH YOUR OPTIMIZER
# ==============================================================================

println("\n" * "="^80)
println("STEP 4: Integrate with globtim")
println("="^80)

println("""
The objective function is ready to use with your globtim optimizer.

Replace the code below with your actual globtim call:

```julia
using Globtim  # Your optimizer package

# Run optimization
result = your_globtim_function(
    objective = error_func,
    bounds = bounds,
    # ... your specific options ...
)

# Extract results
p_best = result.minimizer      # Best parameters found
error_best = result.minimum    # Best error value
n_evals = result.f_calls       # Number of function evaluations

# Compare with true solution
println("Results:")
println("  Best params: \$p_best")
println("  True params: \$p_true")
println("  Best error: \$error_best")
println("  Parameter error: \$(norm(p_best - p_true))")
```

Expected Results:
  - p_best should be close to [1.0, 0.5]
  - error_best should be < 1e-6
  - parameter error should be < 0.01
""")

# ==============================================================================
# DEMO: Simple Random Search (for testing without globtim)
# ==============================================================================

println("\n" * "="^80)
println("DEMO: Simple Random Search (replace with globtim)")
println("="^80)

function simple_random_search(error_func, bounds, n_samples=100)
    best_p = nothing
    best_error = Inf

    for i in 1:n_samples
        # Random sample within bounds
        p = [bounds[j][1] + rand() * (bounds[j][2] - bounds[j][1]) for j in 1:length(bounds)]

        # Evaluate
        err = error_func(p)

        if err < best_error
            best_error = err
            best_p = p
        end
    end

    return best_p, best_error
end

println("\nRunning simple random search (100 samples)...")
p_best, error_best = simple_random_search(error_func, bounds, 100)

println("\nResults:")
@printf "  Best params: [%.4f, %.4f]\n" p_best[1] p_best[2]
@printf "  True params: [%.4f, %.4f]\n" p_true[1] p_true[2]
@printf "  Best error: %.6e\n" error_best
@printf "  Parameter error: %.6e\n" norm(p_best - p_true)

if error_best < 1e-3
    println("\n✓ Random search found a good solution!")
    println("  (Your globtim optimizer should do much better)")
else
    println("\n⚠ Random search did not find optimal (expected)")
    println("  Use proper optimizer to find true parameters")
end

# ==============================================================================
# NEXT STEPS
# ==============================================================================

println("\n" * "="^80)
println("Next Steps")
println("="^80)

println("""
1. Replace the demo above with your globtim optimizer

2. Once that works, try more models:
   - examples/single_test_template.jl - Test any single model
   - examples/batch_test_campaign.jl - Run full test suite
   - examples/test_configs.jl - Pre-configured test cases

3. See documentation:
   - TESTING_GUIDE.md - Complete guide for all 13 models
   - MODEL_CATALOG.md - Quick reference table
   - README.md - Package overview

4. Progressive difficulty:
   Easy (2D)   → LV_2D models
   Medium (3-4D) → LV_3D, DAISY models
   Hard → FitzHugh-Nagumo, identifiability tests
   Very Hard → 20D Generalized LV

Good luck with your testing campaign! 🚀
""")
