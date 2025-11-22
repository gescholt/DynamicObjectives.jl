# Examples - Testing Campaign Scripts

Executable scripts for running testing campaigns with globtim optimizer.

## Quick Start

### 1. First Time Setup

```bash
# Activate the package environment
cd /path/to/Dynamic_objectives
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Run Your First Test

```bash
julia --project=. examples/quick_start.jl
```

This will:
- Set up the easiest model (2D LV with 2 outputs)
- Create an objective function
- Verify everything works
- Show you where to integrate your globtim optimizer

## Available Scripts

### `quick_start.jl` - Your First Test ⭐
**Start here!**

Complete walkthrough for testing a single model:
- Sets up the easiest 2-parameter model
- Shows how to create and verify objective function
- Demonstrates integration points for globtim
- Includes simple random search demo

```bash
julia --project=. examples/quick_start.jl
```

### `single_test_template.jl` - Test Any Model

Template for testing individual models:
- Easy configuration section (just change model name and parameters)
- Automatic setup and verification
- Grid search demo (replace with globtim)
- Copy and modify for your specific tests

```bash
julia --project=. examples/single_test_template.jl
```

**To customize:**
1. Open the file in an editor
2. Modify the CONFIGURATION section (lines 15-30)
3. Run it

### `batch_test_campaign.jl` - Run Multiple Models

Batch testing script for systematic campaigns:
- Pre-configured with 6 representative models (Easy → Hard)
- Runs all tests automatically
- Generates summary report
- Saves results to file

```bash
julia --project=. examples/batch_test_campaign.jl
```

**Outputs:**
- Console summary with pass/fail for each model
- Detailed results table
- `test_results_YYYYMMDD_HHMMSS.txt` file

### `test_configs.jl` - Pre-configured Test Cases

Configuration library with all 13 models:
- Complete configs from TESTING_GUIDE.md
- Organized by difficulty phase
- Helper functions for loading configs
- Import into your own scripts

**Usage in your scripts:**
```julia
include("examples/test_configs.jl")

# List all available configs
list_configs()

# Get a specific config
config = get_config("LV_2D_v3_2outputs")

# Get all easy configs
easy = get_configs_by_phase("EASY")

# Create objective from config
error_func = create_objective(config)
```

## Integration with globtim

All scripts have placeholder functions where you integrate your optimizer. Look for:

```julia
# ==============================================================================
# GLOBTIM INTEGRATION PLACEHOLDER
# ==============================================================================
```

### Example Integration

Replace the placeholder with your actual globtim call:

```julia
using Globtim  # Your optimizer

function run_optimizer(error_func, bounds, config)
    result = globtim_optimize(
        objective = error_func,
        bounds = bounds,
        max_evals = 1000,
        # ... your globtim-specific options ...
    )

    return (
        minimizer = result.best_params,
        minimum = result.best_error,
        converged = result.success,
        iterations = result.iterations,
        f_calls = result.function_evaluations
    )
end
```

## Testing Workflow

### Phase 1: Single Model Validation

1. Start with `quick_start.jl`
2. Integrate your globtim optimizer
3. Verify it finds the true parameters
4. Try 2-3 more models using `single_test_template.jl`

### Phase 2: Systematic Campaign

1. Modify `batch_test_campaign.jl` with your optimizer
2. Run on Easy models first
3. Progress to Medium, then Hard
4. Analyze results file

### Phase 3: Full Evaluation

1. Load all configs from `test_configs.jl`
2. Run comprehensive campaign on all 13 models
3. Generate performance report
4. Identify strengths and weaknesses

## File Overview

```
examples/
├── README.md                    # This file
├── quick_start.jl              # ⭐ START HERE - First test
├── single_test_template.jl     # Template for individual tests
├── batch_test_campaign.jl      # Batch testing script
└── test_configs.jl             # Configuration library (all 13 models)
```

## Expected Results

### For Identifiable Models (most models)
- `error_best < 1e-6` (very close to zero)
- `norm(p_best - p_true) < 0.01` (parameters within 1%)
- Consistent convergence across runs

### For Non-Identifiable Models
These models have multiple parameter sets that give identical outputs:
- `Simple_1D_square`: ±a give same output
- `Simple_2D_product`: Only a*b and a+b are observable
- `Simple_2D_square`: ±b give same output

**Expected:** Optimizer may not find exact `p_true`, but should find equivalent parameters with same invariants.

## Troubleshooting

### "Package not found"
Make sure you've activated the project:
```bash
julia --project=. examples/script.jl
```

### "undefined variable error_func"
Make sure you're using the provided scripts as-is first. If you copied code, ensure all variables are in scope.

### Models timeout or fail
Some models need longer timeout:
- FitzHugh-Nagumo: oscillatory, may need `eval_timeout = 10.0`
- 20D Generalized LV: complex, use `eval_timeout = 15.0`

### Optimizer doesn't converge
1. Check bounds include true parameters
2. Try more evaluations
3. Consider different distance metrics (L2, L1, log_L2)
4. Some models are just hard (identifiability issues)

## Performance Tips

### Start Simple
- Begin with 2D models (fastest, easiest)
- Verify optimizer works correctly
- Then scale to harder problems

### Use Timeout
- Set `eval_timeout = 10.0` for hard models
- Prevents hanging on bad parameter regions
- Returns `Inf` to guide optimizer away

### Choose Distance Metric
- `L2_norm`: Default, works for most models
- `L1_norm`: Robust to outliers
- `log_L2_norm`: Good for wide error ranges (e.g., 20D model)

### Adjust Time Series
- More time points = better accuracy, slower evaluation
- Longer time interval = more information, may be unstable
- Balance speed vs. quality based on your needs

## Next Steps

1. ✅ Run `quick_start.jl` to verify setup
2. ✅ Integrate your globtim optimizer
3. ✅ Test on easy models (2D)
4. ✅ Progress to medium models (3-4D)
5. ✅ Challenge with hard models (FitzHugh-Nagumo, identifiability)
6. ✅ Stress test with 20D model

## See Also

- **[TESTING_GUIDE.md](../TESTING_GUIDE.md)** - Complete testing strategy
- **[MODEL_CATALOG.md](../MODEL_CATALOG.md)** - All 13 models reference
- **[README.md](../README.md)** - Package overview
- **[SETUP.md](../SETUP.md)** - Installation and setup

Good luck with your testing campaign! 🚀
