# Model Catalog

Complete list of all 13 time parameter estimation benchmark problems in Dynamic\_objectives.

## Summary Table

| # | Model Name | Function | Params | States | Outputs | Difficulty | Type | Notes |
|---|------------|----------|--------|--------|---------|------------|------|-------|
| 1 | LV 2D v1 | `define_lotka_volterra_2D_model()` | 2 | 2 | 1 | EASY | Predator-prey | c=1 |
| 2 | LV 2D v2 | `define_lotka_volterra_2D_model_v2()` | 2 | 2 | 1 | EASY | Predator-prey | c=0.1 |
| 3 | LV 2D v3 | `define_lotka_volterra_2D_model_v3()` | 2 | 2 | 1 | EASY | Predator-prey | c=0.5 |
| 4 | LV 2D v3 (2 out) | `define_lotka_volterra_2D_model_v3_two_outputs()` | 2 | 2 | 2 | EASY | Predator-prey | Full observability |
| 5 | LV 3D v1 | `define_lotka_volterra_3D_model()` | 3 | 2 | 1 | MEDIUM | Predator-prey | Classic formulation |
| 6 | LV 3D v2 | `define_lotka_volterra_3D_model_v2()` | 3 | 2 | 1 | MEDIUM | Predator-prey | Standard form |
| 7 | DAISY Ex3 (input) | `define_daisy_ex3_model_4D()` | 4 | 4 | 2 | MEDIUM | Linear | Time-dependent input |
| 8 | DAISY Ex3 (no input) | `define_daisy_ex3_model_4D_no_input()` | 4 | 3 | 2 | MEDIUM | Linear | Benchmark |
| 9 | LV 4D Constrained | `define_constrained_lotka_volterra_4D()` | 4 | 4 | 4 | MEDIUM-HARD | Competition | Skew-symmetric |
| 10 | FitzHugh-Nagumo | `define_fitzhugh_nagumo_3D_model()` | 3 | 2 | 1 | HARD | Neuronal | Oscillatory, may need timeout |
| 11 | Simple 1D (a²) | `define_simple_1D_model_locally_identifiable()` | 1 | 1 | 1 | MEDIUM-HARD | Identifiability test | Non-identifiable |
| 12 | Simple 2D (a*b) | `define_simple_2D_model_locally_identifiable()` | 2 | 1 | 1 | HARD | Identifiability test | Non-identifiable |
| 13 | Simple 2D (b²) | `define_simple_2D_model_locally_identifiable_square()` | 2 | 1 | 1 | HARD | Identifiability test | Non-identifiable |
| 14 | LV 4D Generalized | `define_generalized_lotka_volterra_4D()` | 20 | 4 | 4 | VERY HARD | Competition | Full interaction matrix |

## Grouped by Difficulty

### EASY (4 models)

Good for initial validation of optimizer setup:
- `define_lotka_volterra_2D_model()`
- `define_lotka_volterra_2D_model_v2()`
- `define_lotka_volterra_2D_model_v3()`
- `define_lotka_volterra_2D_model_v3_two_outputs()`

### MEDIUM (4 models)

Standard benchmarks for parameter estimation:
- `define_lotka_volterra_3D_model()`
- `define_lotka_volterra_3D_model_v2()`
- `define_daisy_ex3_model_4D()`
- `define_daisy_ex3_model_4D_no_input()`

### MEDIUM-HARD (2 models)

Challenging problems with structure:
- `define_constrained_lotka_volterra_4D()`
- `define_simple_1D_model_locally_identifiable()`

### HARD (3 models)

Advanced problems requiring robust optimization:
- `define_fitzhugh_nagumo_3D_model()` (oscillatory dynamics)
- `define_simple_2D_model_locally_identifiable()` (non-identifiable)
- `define_simple_2D_model_locally_identifiable_square()` (non-identifiable)

### VERY HARD (1 model)

Ultimate stress test:
- `define_generalized_lotka_volterra_4D()` (20 dimensions!)

## Grouped by Parameter Count

### 1 Parameter
- `define_simple_1D_model_locally_identifiable()`

### 2 Parameters
- `define_lotka_volterra_2D_model()`
- `define_lotka_volterra_2D_model_v2()`
- `define_lotka_volterra_2D_model_v3()`
- `define_lotka_volterra_2D_model_v3_two_outputs()`
- `define_simple_2D_model_locally_identifiable()`
- `define_simple_2D_model_locally_identifiable_square()`

### 3 Parameters
- `define_lotka_volterra_3D_model()`
- `define_lotka_volterra_3D_model_v2()`
- `define_fitzhugh_nagumo_3D_model()`

### 4 Parameters
- `define_daisy_ex3_model_4D()`
- `define_daisy_ex3_model_4D_no_input()`
- `define_constrained_lotka_volterra_4D()`

### 20 Parameters
- `define_generalized_lotka_volterra_4D()`

## Grouped by System Type

### Predator-Prey Dynamics (6 models)

Lotka-Volterra models with oscillatory or equilibrium behavior:
- All LV 2D and 3D models (6 total)

### Competition Models (2 models)

Multi-species competition with interaction matrices:
- `define_generalized_lotka_volterra_4D()` (full 4x4 interaction)
- `define_constrained_lotka_volterra_4D()` (perturbed structure)

### Linear Benchmark (2 models)

DAISY suite linear systems:
- `define_daisy_ex3_model_4D()`
- `define_daisy_ex3_model_4D_no_input()`

### Neuronal Models (1 model)

Excitable systems:
- `define_fitzhugh_nagumo_3D_model()`

### Identifiability Tests (3 models)

Non-identifiable systems (multiple parameter sets give same output):
- `define_simple_1D_model_locally_identifiable()`
- `define_simple_2D_model_locally_identifiable()`
- `define_simple_2D_model_locally_identifiable_square()`

## Special Characteristics

### Requires Timeout

Models that may need `eval_timeout` parameter:
- `define_generalized_lotka_volterra_4D()` (high-dimensional, stiff)
- `define_fitzhugh_nagumo_3D_model()` (oscillatory, limit cycles)

### Partial Observability

Models where not all states are measured (harder estimation):
- All LV 2D v1/v2/v3 single-output variants (7 total)
- All LV 3D models (2 total)
- FitzHugh-Nagumo (1 total)
- Simple 1D/2D identifiability models (3 total)

### Full Observability

All states are measured:
- `define_lotka_volterra_2D_model_v3_two_outputs()`
- `define_generalized_lotka_volterra_4D()`
- `define_constrained_lotka_volterra_4D()`
- Both DAISY models (2 total)

### Non-Identifiable

Multiple parameter sets produce identical outputs (test optimizer behavior):
- `define_simple_1D_model_locally_identifiable()` - ±a gives same output
- `define_simple_2D_model_locally_identifiable()` - (a,b) and (a',b') same if a*b and a+b match
- `define_simple_2D_model_locally_identifiable_square()` - ±b gives same output

## Recommended Testing Sequence

### Sequence 1: Progressive Difficulty

1. `define_lotka_volterra_2D_model_v3_two_outputs()` (2P, EASY, full obs)
2. `define_lotka_volterra_3D_model_v2()` (3P, MEDIUM)
3. `define_constrained_lotka_volterra_4D()` (4P, MEDIUM-HARD)
4. `define_fitzhugh_nagumo_3D_model()` (3P, HARD, oscillatory)
5. `define_generalized_lotka_volterra_4D()` (20P, VERY HARD)

### Sequence 2: Observability Study

1. `define_lotka_volterra_2D_model_v3()` (1 output, partial)
2. `define_lotka_volterra_2D_model_v3_two_outputs()` (2 outputs, full)

Compare convergence speed and accuracy.

### Sequence 3: Identifiability Study

1. `define_simple_1D_model_locally_identifiable()` (1P, ±a symmetry)
2. `define_simple_2D_model_locally_identifiable_square()` (2P, ±b symmetry)
3. `define_simple_2D_model_locally_identifiable()` (2P, product/sum only)

Check if optimizer finds correct invariants.

### Sequence 4: Scalability Study

1. `define_lotka_volterra_2D_model()` (2P)
2. `define_lotka_volterra_3D_model()` (3P)
3. `define_daisy_ex3_model_4D_no_input()` (4P)
4. `define_generalized_lotka_volterra_4D()` (20P)

Track how performance degrades with dimension.

## Quick Reference: Function Signatures

All model definition functions return:

```julia
(model, parameters, states, measured_quantities)
```

Where:
- `model`: ModelingToolkit ODESystem
- `parameters`: Vector of symbolic parameters
- `states`: Vector of symbolic state variables
- `measured_quantities`: Vector of measurement equations

## Usage Example

```julia
using Dynamic_objectives

# Select a model (e.g., #4 from table)
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# Set up problem (see testing_guide.md for recommended values)
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
bounds = [(0.0, 3.0), (0.0, 2.0)]
time_interval = [0.0, 20.0]
numpoints = 30

# Create objective
error_func = make_error_distance(
    model, outputs, ic, p_true, time_interval, numpoints,
    L2_norm, first, nothing;
    return_inf_on_error = true,
    eval_timeout = nothing  # Not needed for easy problems
)

# Verify setup
@assert error_func(p_true) < 1e-6

# Run optimizer on error_func with bounds
```

## See Also

- [Testing Guide](testing_guide.md) - Complete testing protocol with recommended parameters
- [Getting Started](getting_started.md) - Package overview and basic usage
- [Integration](integration.md) - Integration with globtim
