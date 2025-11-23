# Dynamic_objectives Package Architecture

## Design Philosophy: Standalone with Optional Integration

Dynamic_objectives is designed to be **intentionally standalone** with **zero package dependencies** on globtimcore or globtimpostprocessing. This design has important benefits:

### Benefits of Standalone Design

1. **Fast Compilation**: No heavy dependencies means quick loading
2. **Independent Development**: Can be developed and tested independently
3. **No Circular Dependencies**: Avoids the circular dependency issues that plagued earlier versions
4. **Reusable**: Can be used with ANY optimizer, not just globtim
5. **Clean Architecture**: Clear separation of concerns

## Package Dependency Graph

```
┌─────────────────────────────────────────────────────────┐
│               Dynamic_objectives                        │
│     (ODE models, parameter estimation problems)         │
│                                                          │
│     Dependencies: ONLY standard packages                │
│     - DifferentialEquations, ModelingToolkit            │
│     - CSV, DataFrames, PrettyTables                     │
│     - NO globtim packages!                              │
└─────────────────────────────────────────────────────────┘
         │                            │
         │ calls (optional)           │ calls (optional)
         │ via dev setup              │ via dev setup
         ▼                            ▼
    ┌──────────┐              ┌──────────────────┐
    │globtim   │              │globtimpost       │
    │  core    │─────────────▶│  processing      │
    └──────────┘   produces   └──────────────────┘
                   CSV files
```

**Key Point**: Dynamic_objectives has NO formal dependency on globtim packages. The integration is achieved through local dev setup for testing purposes only.

## Integration Pattern: Local Dev Versions

When you want to use Dynamic_objectives with globtim for optimization, you use **local development versions**:

### What is `Pkg.develop()`?

```julia
Pkg.develop(path="../globtimcore")
```

This tells Julia's package manager:
- "Add globtimcore to this environment"
- "But don't download it from a registry"
- "Instead, use the local development version at this path"
- "Track changes to that local code"

### How This Enables Integration

After running setup:

```bash
./setup_dev_packages.jl
```

The following happens:

1. **Dynamic_objectives/Project.toml** is updated with dev entries:
   ```toml
   [deps]
   # ... existing deps ...
   Globtim = "..."  # Added by Pkg.develop()
   GlobtimPostProcessing = "..."  # Added by Pkg.develop()
   ```

2. **Dynamic_objectives/Manifest.toml** points to local paths:
   ```toml
   [[Globtim]]
   path = "../globtimcore"
   uuid = "..."

   [[GlobtimPostProcessing]]
   path = "../globtimpostprocessing"
   uuid = "..."
   ```

3. **Integration tests can now load the packages**:
   ```julia
   using Dynamic_objectives
   using Globtim  # Now available via dev
   using GlobtimPostProcessing  # Now available via dev
   ```

### Why Not Add Them as Formal Dependencies?

You might ask: "Why not just add them to Project.toml permanently?"

**Answer**: We want to keep the distinction:

- **Core package** (Dynamic_objectives): Standalone, no globtim dependencies
- **Testing/Integration**: Optional, requires explicit dev setup

This makes it clear that Dynamic_objectives is reusable with ANY optimizer, not just globtim.

## Testing Architecture

### Test Structure

```
Dynamic_objectives/test/
├── runtests.jl                    # Main test runner
├── test_error_metrics.jl          # Core functionality (no globtim)
├── test_model_definitions.jl      # Core functionality (no globtim)
├── test_globtim_integration.jl    # REQUIRES globtim dev setup ⚠️
└── validate_all_models.jl         # Core functionality (no globtim)
```

### Test Dependencies

**Basic Tests** (always work):
- `test_error_metrics.jl` - Tests error functions, distance metrics
- `test_model_definitions.jl` - Tests ODE models
- `validate_all_models.jl` - Validates all 13 models work

**Integration Tests** (require dev setup):
- `test_globtim_integration.jl` - Tests using Dynamic_objectives with globtim optimizer

### Running Tests

**Without dev setup** (basic tests only):
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
# Will skip integration tests, pass basic tests
```

**With dev setup** (all tests):
```bash
./setup_dev_packages.jl  # One-time setup
./run_tests.sh           # Now all tests pass
```

## Usage Patterns

### Pattern 1: Standalone Usage (No globtim)

```julia
using Dynamic_objectives

# Use with ANY optimizer (Optim.jl, NLopt.jl, etc.)
model, params, states, outputs = define_lotka_volterra_2D_model_v3()
error_func = make_error_distance(...)

using Optim
result = optimize(error_func, lower, upper, ParticleSwarm())
```

**No globtim required!**

### Pattern 2: Integration with globtim (Requires dev setup)

```julia
using Dynamic_objectives
using Globtim  # Available after ./setup_dev_packages.jl

# Create objective
error_func = make_error_distance(...)

# Use globtim optimizer
result = run_standard_experiment(
    error_func,
    bounds,
    StandardExperimentConfig(max_degree=18)
)
```

**Requires one-time dev setup.**

### Pattern 3: Full 2-Stage Pipeline (Requires dev setup)

```julia
using Dynamic_objectives
using Globtim
using GlobtimPostProcessing

# Stage 1: Find raw critical points
raw_result = run_standard_experiment(error_func, bounds, config)

# Stage 2: Refine critical points
refined_result = refine_experiment_results(
    raw_result.output_dir,
    error_func,
    ode_refinement_config()
)
```

**Requires one-time dev setup.**

## File Organization

### Core Package Files (No globtim dependencies)

```
src/
├── Dynamic_objectives.jl      # Main module
├── error_metrics.jl            # Distance functions
├── data_generation.jl          # Time series generation
├── display.jl                  # Pretty printing
└── systems/                    # ODE models
    ├── lotka_volterra.jl
    ├── daisy_benchmarks.jl
    ├── fitzhugh_nagumo.jl
    └── ...
```

### Integration Documentation (References globtim)

```
├── REFINEMENT_INTEGRATION_GUIDE.md   # How to use with globtim
├── TESTING_GUIDE.md                  # Testing with globtim
├── ARCHITECTURE.md                   # This file
└── setup_dev_packages.jl             # Dev setup script
```

### Test Files

```
test/
├── test_error_metrics.jl        # No globtim
├── test_model_definitions.jl    # No globtim
├── validate_all_models.jl       # No globtim
└── test_globtim_integration.jl  # Requires globtim (via dev setup)
```

## Summary: Key Takeaways

1. **Dynamic_objectives is standalone** - No formal dependencies on globtim
2. **Integration is optional** - Requires explicit dev setup
3. **One-time setup** - Run `./setup_dev_packages.jl` once
4. **Dev versions track local changes** - Always uses your latest local code
5. **Clean architecture** - Clear separation between core and integration

## Troubleshooting

### Error: "Package Globtim not found"

**Cause**: Dev setup not run

**Solution**:
```bash
./setup_dev_packages.jl
```

### Error: "Package at path '../globtimcore' not found"

**Cause**: Directory structure incorrect

**Solution**: Ensure directory structure is:
```
GlobalOptim/
├── Dynamic_objectives/
├── globtimcore/
└── globtimpostprocessing/
```

### Integration tests fail but basic tests pass

**Cause**: Dev setup not run or incomplete

**Solution**:
```bash
./setup_dev_packages.jl
```

## See Also

- **[README.md](README.md)** - Package overview and quick start
- **[TESTING_GUIDE.md](TESTING_GUIDE.md)** - Comprehensive testing strategy
- **[REFINEMENT_INTEGRATION_GUIDE.md](REFINEMENT_INTEGRATION_GUIDE.md)** - Integration patterns
- **[../docs/API_DESIGN_REFINEMENT.md](../docs/API_DESIGN_REFINEMENT.md)** - Overall architecture
