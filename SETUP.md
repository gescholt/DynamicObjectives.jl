# Quick Setup and Run Guide

## First-Time Setup

### 1. Install Dependencies

```julia
cd /path/to/Dynamic_objectives

# Activate the project environment
using Pkg
Pkg.activate(".")

# Install all dependencies
Pkg.instantiate()
```

This will install all required packages:
- DataStructures
- LinearAlgebra
- Logging
- ModelingToolkit
- OrdinaryDiffEq
- StaticArrays
- Test (for testing only)

### 2. Verify Installation

```julia
# Load the package
using Dynamic_objectives

# Quick smoke test
model, params, states, outputs = define_lotka_volterra_2D_model()
println("✓ Package loaded successfully!")
println("  Model has $(length(params)) parameters")
```

## Running Tests

### Option 1: Shell Script (Easiest)

```bash
cd /path/to/Dynamic_objectives

# Run all tests
./run_tests.sh

# Or specific test suite
./run_tests.sh basic      # Basic tests (fast)
./run_tests.sh validate   # All 13 models (slower)
```

### Option 2: Julia REPL

```julia
using Pkg
Pkg.activate("/path/to/Dynamic_objectives")

# Basic test suite
Pkg.test()

# Or run validation directly
include("test/validate_all_models.jl")
```

### Option 3: Command Line

```bash
cd /path/to/Dynamic_objectives

# Basic tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Validation suite
julia --project=. test/validate_all_models.jl
```

## Expected Test Results

### Basic Tests (test/runtests.jl)
```
Test Summary:           | Pass  Total  Time
Dynamic_objectives.jl   |    5      5  X.Xs
  DAISY LV4D Definition |    3      3  X.Xs
  Generalized LV4D ...  |    3      3  X.Xs
  Constrained LV4D ...  |    3      3  X.Xs
  Error Function ...    |    2      2  X.Xs
  Distance Functions    |    3      3  X.Xs
```

### Validation Suite (test/validate_all_models.jl)
```
================================================================================
Dynamic_objectives - Model Validation Suite
================================================================================
[1/13] Testing LV 2D v1... ✓ PASS (error @ true: 0.0, @ perturbed: 1.23)
[2/13] Testing LV 2D v2... ✓ PASS (error @ true: 0.0, @ perturbed: 1.45)
...
[13/13] Testing Simple 2D (b²)... ✓ PASS (error @ true: 0.0, @ perturbed: 2.34)
================================================================================
[14/14] Testing LV 4D Generalized (20 params)... ✓ PASS (error @ true: 0.0)
================================================================================
Validation Summary:
  PASSED: 14/14
  FAILED: 0/14
================================================================================

✅ All 14 models validated successfully!

Ready for globtim testing campaign!
See TESTING_GUIDE.md for recommended test configurations.
```

## Troubleshooting

### Error: "Package Test not found"
**Fixed!** This was resolved by adding Test to `[extras]` in Project.toml.

If you still see this, run:
```julia
using Pkg
Pkg.activate(".")
Pkg.add("Test")  # Shouldn't be needed, but forces Test installation
```

### Error: "Package X not found"
Make sure you've run `Pkg.instantiate()`:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Tests fail with ODE solver errors
Some models may be sensitive to initial conditions or parameters. Check:
1. Are you using the recommended configurations from TESTING_GUIDE.md?
2. Try increasing `eval_timeout` for slow models
3. Check if parameters are within reasonable bounds

### Validation takes too long
The 20D model can be slow. To skip it:
```julia
# Edit test/validate_all_models.jl and comment out the 20D test
# Or just run basic tests: ./run_tests.sh basic
```

## Next Steps After Setup

1. ✅ Run `./run_tests.sh` to verify everything works
2. 📖 Read [TESTING_GUIDE.md](TESTING_GUIDE.md) for complete model configurations
3. 🚀 Start with easy models from [MODEL_CATALOG.md](MODEL_CATALOG.md)
4. 🎯 Test with your globtim optimizer

## Quick Example to Test Setup

```julia
using Dynamic_objectives

# Simple 2-parameter model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# Create objective function
p_true = [1.0, 0.5]
ic = [1.0, 0.5]

error_func = make_error_distance(
    model, outputs, ic, p_true,
    [0.0, 20.0], 30,  # time interval, num points
    L2_norm, first, nothing;
    return_inf_on_error = true
)

# Test it works
@assert error_func(p_true) < 1e-6  "Error at true params should be ~0"
@assert error_func([1.1, 0.6]) > 0  "Error at other params should be > 0"

println("✅ Everything works! Ready to use with globtim.")
```

## Environment Info

- **Julia version required:** 1.10+
- **Key dependencies:**
  - ModelingToolkit v9 or v10
  - OrdinaryDiffEq v6
  - StaticArrays v1
  - DataStructures v0.18
- **Platform:** Cross-platform (Linux, macOS, Windows)

## Getting Help

If you encounter issues:
1. Check this guide's Troubleshooting section
2. Review [CLEANUP_PLAN.md](CLEANUP_PLAN.md) for repository structure
3. Check model-specific notes in [TESTING_GUIDE.md](TESTING_GUIDE.md)
4. Look at example usage in test files
