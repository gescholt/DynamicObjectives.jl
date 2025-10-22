# Migration Guide: DynamicalSystems.jl → Dynamic_objectives

This guide explains how to migrate experiments from the old `include("DynamicalSystems.jl")` approach to the new `Dynamic_objectives` package.

## Summary of Changes

- **Old**: `DynamicalSystems.jl` was loaded via `include()` in each experiment
- **New**: `Dynamic_objectives` is a proper Julia package with precompilation
- **Benefits**: Faster startup, no world age warnings, cleaner code

## Migration Steps

### 1. Update Experiment File

**Before (OLD):**
```julia
# CRITICAL: Load DynamicalSystems module FIRST
const PROJECT_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "Examples", "systems", "DynamicalSystems.jl"))

using Pkg
Pkg.activate(PROJECT_ROOT)
```

**After (NEW):**
```julia
const PROJECT_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

using Pkg
Pkg.activate(PROJECT_ROOT)
Pkg.instantiate()

using Dynamic_objectives
```

### 2. Remove invokelatest() Wrappers

**Before (OLD):**
```julia
model, params, states, outputs = Base.invokelatest() do
    DynamicalSystems.define_daisy_ex3_model_4D()
end

error_func = Base.invokelatest() do
    DynamicalSystems.make_error_distance(
        model, outputs, ic, p_true, time_interval, num_points,
        DynamicalSystems.L2_norm, first, nothing;
        return_inf_on_error = true, eval_timeout = 5.0
    )
end
```

**After (NEW):**
```julia
model, params, states, outputs = define_daisy_ex3_model_4D()

error_func = make_error_distance(
    model, outputs, ic, p_true, time_interval, num_points,
    L2_norm, first, nothing;
    return_inf_on_error = true, eval_timeout = 5.0
)
```

### 3. Update Function References

All functions are now imported directly from `Dynamic_objectives`:

| Old Reference | New Reference |
|--------------|---------------|
| `DynamicalSystems.define_daisy_ex3_model_4D()` | `define_daisy_ex3_model_4D()` |
| `DynamicalSystems.L2_norm` | `L2_norm` |
| `DynamicalSystems.make_error_distance` | `make_error_distance` |
| `DynamicalSystems.sample_data` | `sample_data` |

## Available Functions

All functions from the original `DynamicalSystems.jl` are now exported from `Dynamic_objectives`:

### ODE System Definitions

**DAISY Models:**
- `define_daisy_ex3_model_4D()`
- `define_daisy_ex3_model_4D_no_input()`

**Lotka-Volterra Models:**
- `define_generalized_lotka_volterra_4D()`
- `define_constrained_lotka_volterra_4D()`
- `define_lotka_volterra_3D_model()`
- `define_lotka_volterra_3D_model_v2()`
- `define_lotka_volterra_2D_model()`
- `define_lotka_volterra_2D_model_v2()`
- `define_lotka_volterra_2D_model_v3()`
- `define_lotka_volterra_2D_model_v3_two_outputs()`

**Other Models:**
- `define_fitzhugh_nagumo_3D_model()`
- `define_simple_1D_model_locally_identifiable()`
- `define_simple_2D_model_locally_identifiable()`
- `define_simple_2D_model_locally_identifiable_square()`

### Data Generation & Error Metrics

- `sample_data()` - Generate synthetic time series from ODE systems
- `make_error_distance()` - Create objective functions for parameter estimation
- `L1_norm()` - L1 distance metric
- `L2_norm()` - L2 distance metric
- `log_L2_norm()` - Log-L2 distance metric

## Example: Complete Migration

**File: `experiments/lv4d_2025/lv4d_experiment.jl`**

**Changes:**
1. ✅ Removed `include(joinpath(PROJECT_ROOT, "Examples", "systems", "DynamicalSystems.jl"))`
2. ✅ Added `using Dynamic_objectives` after `Pkg.instantiate()`
3. ✅ Removed all `Base.invokelatest()` wrappers
4. ✅ Changed `DynamicalSystems.define_daisy_ex3_model_4D()` → `define_daisy_ex3_model_4D()`
5. ✅ Changed `DynamicalSystems.L2_norm` → `L2_norm`

## ModelRegistry Integration

The `ModelRegistry` now detects if `Dynamic_objectives` is loaded:

```julia
using Dynamic_objectives  # Load in your experiment

# ModelRegistry will automatically register ODE models
using Globtim.ModelRegistry

# Query available models
models = list_models(category = :ode)
```

If you don't load `Dynamic_objectives`, you'll see a warning:
```
┌ Warning: Dynamic_objectives not loaded. ODE models will not be registered.
└ Load it with: using Dynamic_objectives
```

## Troubleshooting

### Issue: "Dynamic_objectives not found"

**Solution:** Ensure Dynamic_objectives is in your Project.toml:
```bash
cd /path/to/globtimcore
julia --project=. -e 'using Pkg; Pkg.develop(path="../Dynamic_objectives")'
```

### Issue: Still seeing world age warnings

**Solution:** Make sure you've:
1. Removed all `include("DynamicalSystems.jl")` calls
2. Removed all `Base.invokelatest()` wrappers
3. Added `using Dynamic_objectives` after `Pkg.activate()`

### Issue: Function not found (e.g., `define_daisy_ex3_model_4D`)

**Solution:** The function is exported from `Dynamic_objectives`, so just use it directly:
```julia
using Dynamic_objectives

# Direct call - no module prefix needed
model, params, states, outputs = define_daisy_ex3_model_4D()
```

## Benefits of Migration

✅ **Faster startup**: Module is precompiled once, not on every experiment run
✅ **No world age warnings**: Proper module boundaries eliminate Julia 1.12 warnings
✅ **Cleaner code**: No `include()` hacks or `invokelatest()` wrappers
✅ **Better organization**: Logical file structure with separate concerns
✅ **Independent testing**: Dedicated test suite for ODE models
✅ **Reusability**: Other projects can depend on `Dynamic_objectives`

## Questions?

See:
- [Dynamic_objectives README](/Users/ghscholt/GlobalOptim/Dynamic_objectives/README.md)
- [GitLab Repository](https://git.mpi-cbg.de/globaloptim/dynamic_objectives)
- [Example: lv4d_experiment.jl](/Users/ghscholt/GlobalOptim/globtimcore/experiments/lv4d_2025/lv4d_experiment.jl)
