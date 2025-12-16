# Testing Guide

## Overview

This guide provides a comprehensive testing strategy for time parameter estimation problems using the Dynamic\_objectives package with globtim. All 13 models are suitable for testing parameter estimation algorithms.

## Setup: Loading Local Dev Versions (REQUIRED)

**Dynamic\_objectives is standalone** (no package dependencies), but to use it with globtim for testing, you need to set up local dev versions:

### One-Time Setup

```bash
cd /Users/ghscholt/GlobalOptim/Dynamic_objectives

julia --project=. -e '
using Pkg

# Add local dev versions of globtimcore and globtimpostprocessing
Pkg.develop(path="../globtimcore")
Pkg.develop(path="../globtimpostprocessing")

# Verify
using Globtim
using GlobtimPostProcessing
println("Setup complete!")
'
```

**This only needs to be done ONCE.** After this, the packages will always be available in the Dynamic\_objectives environment.

### Verifying Setup

```bash
# Run integration tests to verify
./run_tests.sh

# If you see errors like "Package Globtim not found", re-run the setup above
```

---

## Model Configurations

### 1. DAISY Benchmark Models

#### 1.1 DAISY Example 3 (4D with input)

**Function:** `define_daisy_ex3_model_4D()`

**Parameters:** 4 (p1, p3, p4, p6)

**Recommended Test Configuration:**

```julia
model, params, states, outputs = define_daisy_ex3_model_4D()

# Suggested true parameters
p_true = [0.5, 0.3, 0.4, 0.2]

# Suggested initial conditions
ic = [1.0, 0.5, 0.3, 0.0]

# Suggested parameter bounds
bounds = [
    (0.0, 2.0),   # p1
    (0.0, 1.0),   # p3
    (0.0, 1.0),   # p4
    (-0.5, 0.5)   # p6
]

# Suggested time interval
time_interval = [0.0, 20.0]
numpoints = 50

# Difficulty: MEDIUM
# Characteristics: Linear system, 2 outputs, time-dependent input
```

#### 1.2 DAISY Example 3 (4D no input)

**Function:** `define_daisy_ex3_model_4D_no_input()`

**Parameters:** 4 (p1, p3, p4, p6)

**Recommended Test Configuration:**

```julia
model, params, states, outputs = define_daisy_ex3_model_4D_no_input()

p_true = [0.5, 0.3, 0.4, 0.2]
ic = [1.0, 0.5, 0.3]

bounds = [
    (0.0, 2.0),   # p1
    (0.0, 1.0),   # p3
    (0.0, 1.0),   # p4
    (-0.5, 0.5)   # p6
]

time_interval = [0.0, 20.0]
numpoints = 50

# Difficulty: MEDIUM
# Characteristics: Linear system, 2 outputs, no input
```

---

### 2. Lotka-Volterra Systems

#### 2.1 Generalized LV 4D (20 parameters)

**Function:** `define_generalized_lotka_volterra_4D()`

**Parameters:** 20 (4 growth rates + 16 interaction matrix elements)

**Recommended Test Configuration:**

```julia
model, params, states, outputs = define_generalized_lotka_volterra_4D()

# Suggested true parameters (example: competitive system)
p_true = [
    # Growth rates
    1.0, 1.0, 1.0, 1.0,
    # Interaction matrix (row-major)
    -0.5, -0.1, -0.1, -0.1,  # B11-B14
    -0.1, -0.5, -0.1, -0.1,  # B21-B24
    -0.1, -0.1, -0.5, -0.1,  # B31-B34
    -0.1, -0.1, -0.1, -0.5   # B41-B44
]

ic = [0.5, 0.5, 0.5, 0.5]

# Suggested bounds
bounds = [
    # Growth rates
    (0.0, 3.0), (0.0, 3.0), (0.0, 3.0), (0.0, 3.0),
    # Interaction matrix elements
    (-2.0, 0.0), (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0),
    (-1.0, 1.0), (-2.0, 0.0), (-1.0, 1.0), (-1.0, 1.0),
    (-1.0, 1.0), (-1.0, 1.0), (-2.0, 0.0), (-1.0, 1.0),
    (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0), (-2.0, 0.0)
]

time_interval = [0.0, 15.0]
numpoints = 40

# Difficulty: VERY HARD
# Characteristics: High-dimensional (20D), nonlinear, 4 outputs
# Warning: May require eval_timeout due to stiffness
```

#### 2.2 Constrained LV 4D (4 parameters)

**Function:** `define_constrained_lotka_volterra_4D()`

**Parameters:** 4 (eps1, eps2, eps3, eps4 - skew-symmetric perturbations)

**Recommended Test Configuration:**

```julia
model, params, states, outputs = define_constrained_lotka_volterra_4D()

p_true = [0.01, 0.02, -0.01, 0.03]
ic = [0.8, 1.2, 0.8, 1.2]

bounds = [
    (-0.1, 0.1),   # eps1
    (-0.1, 0.1),   # eps2
    (-0.1, 0.1),   # eps3
    (-0.1, 0.1)    # eps4
]

time_interval = [0.0, 30.0]
numpoints = 60

# Difficulty: MEDIUM-HARD
# Characteristics: Constrained structure, coupled subsystems, 4 outputs
```

#### 2.3-2.8 LV 2D/3D Models

See [Model Catalog](model_catalog.md) for complete configurations.

---

### 3. Other Systems

#### 3.1 FitzHugh-Nagumo 3D

**Function:** `define_fitzhugh_nagumo_3D_model()`

**Parameters:** 3 (g, a, b)

**Recommended Test Configuration:**

```julia
model, params, states, outputs = define_fitzhugh_nagumo_3D_model()

p_true = [0.8, 0.7, 0.8]
ic = [0.0, 0.0]

bounds = [
    (0.1, 2.0),    # g (coupling)
    (0.0, 2.0),    # a (threshold)
    (0.0, 2.0)     # b (recovery)
]

time_interval = [0.0, 50.0]
numpoints = 100

# Difficulty: HARD
# Characteristics: Neuronal model, oscillatory behavior, 1 output
# Warning: May exhibit limit cycles, use eval_timeout
```

#### 3.2-3.4 Identifiability Test Models

**Function:** `define_simple_2D_model_locally_identifiable()`

**Characteristics:** NOT globally identifiable (a*b and a+b only)

> **Note:** Multiple parameter sets give same output - ideal for testing optimizer behavior. Expected: Optimizer should find parameter combinations where a*b and a+b match true values.

---

## Testing Protocol

### Basic Testing Setup

For each model, create an objective function using:

```julia
using Dynamic_objectives

# 1. Define model
model, params, states, outputs = define_XXX_model()

# 2. Set up problem
p_true = [...]  # From recommendations above
ic = [...]      # From recommendations above
time_interval = [...]
numpoints = ...

# 3. Create objective function
error_func = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    numpoints,
    L2_norm,              # Or L1_norm, log_L2_norm
    first,                # Aggregate function (for multi-output)
    nothing;              # No noise
    return_inf_on_error = true,
    eval_timeout = 10.0   # Recommended: 5-10 seconds
)

# 4. Test at true parameters (should be ≈ 0)
@assert error_func(p_true) < 1e-6

# 5. Test at perturbed parameters (should be > 0)
p_test = p_true .+ 0.1
@assert error_func(p_test) > 0
```

### Recommended Testing Campaign

#### Phase 1: Easy Problems (2 parameters)

Start with these to validate globtim setup:
1. `define_lotka_volterra_2D_model()`
2. `define_lotka_volterra_2D_model_v2()`
3. `define_lotka_volterra_2D_model_v3()`
4. `define_lotka_volterra_2D_model_v3_two_outputs()`

**Expected:** Quick convergence, low function evaluation counts

#### Phase 2: Medium Problems (3-4 parameters)

Test scalability and robustness:
1. `define_lotka_volterra_3D_model()`
2. `define_lotka_volterra_3D_model_v2()`
3. `define_daisy_ex3_model_4D()`
4. `define_daisy_ex3_model_4D_no_input()`
5. `define_constrained_lotka_volterra_4D()`
6. `define_fitzhugh_nagumo_3D_model()`

**Expected:** Moderate convergence time, may require more evaluations

#### Phase 3: Hard Problems (identifiability issues)

Test optimizer behavior on non-identifiable systems:
1. `define_simple_1D_model_locally_identifiable()`
2. `define_simple_2D_model_locally_identifiable()`
3. `define_simple_2D_model_locally_identifiable_square()`

**Expected:** Convergence to equivalent parameter sets, not necessarily p\_true

#### Phase 4: Very Hard Problems (high-dimensional)

Ultimate stress test:
1. `define_generalized_lotka_volterra_4D()` (20 parameters!)

**Expected:** Long computation time, may not fully converge

---

## Distance Function Comparison

Test each problem with different distance functions:

```julia
# L2 norm (standard Euclidean distance)
error_L2 = make_error_distance(..., L2_norm, ...)

# L1 norm (Manhattan distance, scaled by 100)
error_L1 = make_error_distance(..., L1_norm, ...)

# log-L2 norm (log-scale, good for wide-range errors)
error_log = make_error_distance(..., log_L2_norm, ...)
```

**Recommendation:** Start with L2\_norm, then compare with log\_L2\_norm for problems with wide error ranges.

---

## Timeout Configuration

Models that may require timeout (due to stiffness or oscillations):
- `define_generalized_lotka_volterra_4D()` (20D, nonlinear)
- `define_fitzhugh_nagumo_3D_model()` (oscillatory)
- Any model with large parameter values or long time intervals

**Recommended timeout:** 5-10 seconds per evaluation

```julia
error_func = make_error_distance(
    ...,
    eval_timeout = 10.0  # seconds
)
```

---

## Noise Injection Testing

Test robustness to noisy data:

```julia
# Add Gaussian noise to reference data
add_noise = y -> y .+ randn(length(y)) .* 0.05  # 5% noise

error_func = make_error_distance(
    model, outputs, ic, p_true, time_interval, numpoints,
    L2_norm, first, add_noise;
    return_inf_on_error = true
)
```

---

## Performance Metrics

For each test, record:
1. **Convergence success:** Did it find parameters within tolerance?
2. **Function evaluations:** How many evaluations to convergence?
3. **Wall-clock time:** Total optimization time
4. **Final error:** error\_func(p\_best)
5. **Parameter error:** norm(p\_best - p\_true) for identifiable models
6. **Timeout events:** How many evaluations timed out?

---

## Recommended Test Matrix

| Model | Difficulty | Params | Distance | Timeout | Priority |
|-------|-----------|--------|----------|---------|----------|
| LV 2D v3 (2 outputs) | EASY | 2 | L2 | No | HIGH |
| LV 2D v2 | EASY | 2 | L2 | No | HIGH |
| LV 3D v2 | MEDIUM | 3 | L2 | No | HIGH |
| DAISY Ex3 (no input) | MEDIUM | 4 | L2 | No | HIGH |
| FitzHugh-Nagumo | HARD | 3 | L2 | Yes (10s) | MEDIUM |
| Constrained LV 4D | MEDIUM-HARD | 4 | L2 | Optional | MEDIUM |
| Simple 2D (square) | HARD* | 2 | L2 | No | LOW |
| Generalized LV 4D | VERY HARD | 20 | log-L2 | Yes (10s) | LOW |

*Hard due to identifiability, not dimensionality

---

## Expected Results

**For globally identifiable models:**
- Optimum should be close to `p_true`
- `norm(p_best - p_true)` should be < 1% of parameter magnitudes
- `error_func(p_best)` should be < 1e-6 (no noise) or ≈ noise level

**For locally identifiable models:**
- Optimum may differ from `p_true`
- Check invariant quantities match:
  - Simple 1D: `a²` should match
  - Simple 2D (product): `a*b` and `a+b` should match
  - Simple 2D (square): `a` and `b²` should match

---

## Debugging Tips

If optimization fails:
1. **Verify setup:** Check `error_func(p_true)` is near zero
2. **Check bounds:** Ensure `p_true` is within bounds
3. **Reduce timeout:** May indicate stiff system
4. **Try different distance:** log\_L2\_norm may help with scaling
5. **Increase numpoints:** More data points may improve identifiability
6. **Check for NaN/Inf:** Review logs for failed ODE solves

---

## Advanced Testing

### Multi-output Aggregation

For models with multiple outputs, test different aggregation strategies:

```julia
# Take first output only
error_first = make_error_distance(..., aggregate_distances = first)

# Maximum error across outputs
error_max = make_error_distance(..., aggregate_distances = maximum)

# Sum of errors
error_sum = make_error_distance(..., aggregate_distances = sum)

# Mean error
error_mean = make_error_distance(..., aggregate_distances = mean)
```

### Custom Distance Functions

Create domain-specific distance metrics:

```julia
# Weighted L2 norm (emphasize early times)
function weighted_L2(y_true, y_test)
    weights = exp.(-0.1 * (0:length(y_true)-1))  # Decay weights
    return norm((y_true - y_test) .* weights, 2)
end

error_func = make_error_distance(..., weighted_L2, ...)
```

---

## Questions?

For issues or questions:
- Review code documentation in `src/error_metrics.jl`
- Check model definitions in `src/systems/*.jl`
- See [Getting Started](getting_started.md) for basic usage
- See [Integration](integration.md) for integration with globtimcore
