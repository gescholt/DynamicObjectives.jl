# Testing Guide for Dynamic_objectives

## Overview

This guide provides a comprehensive testing strategy for time parameter estimation problems using the Dynamic_objectives package with globtim. All 13 models are suitable for testing parameter estimation algorithms.

## Model Catalog

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
    # Interaction matrix elements (diagonal: negative, off-diagonal: negative or positive)
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

#### 2.3 LV 3D Model (variant 1)
**Function:** `define_lotka_volterra_3D_model()`

**Parameters:** 3 (a, b, c)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_lotka_volterra_3D_model()

p_true = [0.5, -0.3, 0.2]
ic = [1.0, 0.5]

bounds = [
    (-1.0, 2.0),   # a
    (-1.0, 0.0),   # b (negative for predator-prey)
    (-1.0, 1.0)    # c
]

time_interval = [0.0, 20.0]
numpoints = 40

# Difficulty: MEDIUM
# Characteristics: Classic predator-prey, 1 output (partial observability)
```

#### 2.4 LV 3D Model (variant 2)
**Function:** `define_lotka_volterra_3D_model_v2()`

**Parameters:** 3 (a, b, c)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_lotka_volterra_3D_model_v2()

p_true = [1.0, 0.5, 0.3]
ic = [1.0, 0.5]

bounds = [
    (0.0, 3.0),    # a (prey growth)
    (0.0, 2.0),    # b (interaction rate)
    (0.0, 1.0)     # c (predator conversion)
]

time_interval = [0.0, 20.0]
numpoints = 40

# Difficulty: MEDIUM
# Characteristics: Standard formulation, 1 output
```

#### 2.5 LV 2D Model (c=1)
**Function:** `define_lotka_volterra_2D_model()`

**Parameters:** 2 (a, b)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_lotka_volterra_2D_model()

p_true = [0.5, -0.3]
ic = [1.0, 0.5]

bounds = [
    (-1.0, 2.0),   # a
    (-1.0, 0.0)    # b
]

time_interval = [0.0, 20.0]
numpoints = 30

# Difficulty: EASY
# Characteristics: Low-dimensional, 1 output
```

#### 2.6 LV 2D Model v2 (c=0.1)
**Function:** `define_lotka_volterra_2D_model_v2()`

**Parameters:** 2 (a, b)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_lotka_volterra_2D_model_v2()

p_true = [0.5, -0.3]
ic = [1.0, 0.5]

bounds = [
    (-1.0, 2.0),   # a
    (-1.0, 0.0)    # b
]

time_interval = [0.0, 20.0]
numpoints = 30

# Difficulty: EASY
# Characteristics: Low-dimensional, 1 output, slow prey growth
```

#### 2.7 LV 2D Model v3 (c=0.5)
**Function:** `define_lotka_volterra_2D_model_v3()`

**Parameters:** 2 (a, b)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_lotka_volterra_2D_model_v3()

p_true = [1.0, 0.5]
ic = [1.0, 0.5]

bounds = [
    (0.0, 3.0),    # a
    (0.0, 2.0)     # b
]

time_interval = [0.0, 20.0]
numpoints = 30

# Difficulty: EASY
# Characteristics: Low-dimensional, 1 output
```

#### 2.8 LV 2D Model v3 (two outputs)
**Function:** `define_lotka_volterra_2D_model_v3_two_outputs()`

**Parameters:** 2 (a, b)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

p_true = [1.0, 0.5]
ic = [1.0, 0.5]

bounds = [
    (0.0, 3.0),    # a
    (0.0, 2.0)     # b
]

time_interval = [0.0, 20.0]
numpoints = 30

# Difficulty: EASY
# Characteristics: Low-dimensional, 2 outputs (full observability)
```

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

#### 3.2 Simple 2D Locally Identifiable (product form)
**Function:** `define_simple_2D_model_locally_identifiable()`

**Parameters:** 2 (a, b)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_simple_2D_model_locally_identifiable()

p_true = [2.0, 3.0]
ic = [1.0]

bounds = [
    (0.1, 5.0),    # a
    (0.1, 5.0)     # b
]

time_interval = [0.0, 5.0]
numpoints = 20

# Difficulty: HARD (identifiability issue)
# Characteristics: NOT globally identifiable (a*b and a+b only)
# Note: Multiple parameter sets give same output - ideal for testing optimizer behavior
# Expected: Optimizer should find parameter combinations where a*b and a+b match true values
```

#### 3.3 Simple 2D Locally Identifiable (square form)
**Function:** `define_simple_2D_model_locally_identifiable_square()`

**Parameters:** 2 (a, b)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_simple_2D_model_locally_identifiable_square()

p_true = [0.5, 2.0]
ic = [1.0]

bounds = [
    (-1.0, 2.0),   # a
    (-3.0, 3.0)    # b
]

time_interval = [0.0, 5.0]
numpoints = 20

# Difficulty: HARD (identifiability issue)
# Characteristics: NOT globally identifiable (b² only)
# Note: Both b and -b give same output - ideal for testing optimizer behavior
# Expected: Optimizer may find either b or -b
```

#### 3.4 Simple 1D Locally Identifiable
**Function:** `define_simple_1D_model_locally_identifiable()`

**Parameters:** 1 (a)

**Recommended Test Configuration:**
```julia
model, params, states, outputs = define_simple_1D_model_locally_identifiable()

p_true = [1.5]
ic = [1.0]

bounds = [
    (-3.0, 3.0)    # a
]

time_interval = [0.0, 5.0]
numpoints = 20

# Difficulty: MEDIUM-HARD (identifiability issue)
# Characteristics: NOT globally identifiable (a² only)
# Note: Both a and -a give same output
# Expected: Optimizer may find either a or -a
```

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

**Expected:** Convergence to equivalent parameter sets, not necessarily p_true

#### Phase 4: Very Hard Problems (high-dimensional)
Ultimate stress test:
1. `define_generalized_lotka_volterra_4D()` (20 parameters!)

**Expected:** Long computation time, may not fully converge

### Distance Function Comparison

Test each problem with different distance functions:

```julia
# L2 norm (standard Euclidean distance)
error_L2 = make_error_distance(..., L2_norm, ...)

# L1 norm (Manhattan distance, scaled by 100)
error_L1 = make_error_distance(..., L1_norm, ...)

# log-L2 norm (log-scale, good for wide-range errors)
error_log = make_error_distance(..., log_L2_norm, ...)
```

**Recommendation:** Start with L2_norm, then compare with log_L2_norm for problems with wide error ranges.

### Timeout Configuration

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

### Noise Injection Testing

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

### Performance Metrics

For each test, record:
1. **Convergence success:** Did it find parameters within tolerance?
2. **Function evaluations:** How many evaluations to convergence?
3. **Wall-clock time:** Total optimization time
4. **Final error:** error_func(p_best)
5. **Parameter error:** norm(p_best - p_true) for identifiable models
6. **Timeout events:** How many evaluations timed out?

### Recommended Test Matrix

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

### Expected Results

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

### Debugging Tips

If optimization fails:
1. **Verify setup:** Check `error_func(p_true)` is near zero
2. **Check bounds:** Ensure `p_true` is within bounds
3. **Reduce timeout:** May indicate stiff system
4. **Try different distance:** log_L2_norm may help with scaling
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

### Uneven Sampling

Test with non-uniform time sampling (requires modifying `sample_data` call):

```julia
# More points at early times
uneven_times = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

data = sample_data(
    problem, model, outputs, time_interval, p_true, ic, length(uneven_times);
    uneven_sampling = true,
    uneven_sampling_times = uneven_times
)
```

---

## Output Format for globtim

Recommended experiment structure:

```julia
# experiments/test_campaign_2025/config.jl

using Dynamic_objectives

experiments = [
    (
        name = "lv2d_easy",
        model_fn = define_lotka_volterra_2D_model_v3_two_outputs,
        p_true = [1.0, 0.5],
        ic = [1.0, 0.5],
        bounds = [(0.0, 3.0), (0.0, 2.0)],
        time_interval = [0.0, 20.0],
        numpoints = 30,
        distance = L2_norm,
        timeout = nothing
    ),
    # Add more experiments...
]

function create_objective(exp)
    model, params, states, outputs = exp.model_fn()
    return make_error_distance(
        model, outputs, exp.ic, exp.p_true,
        exp.time_interval, exp.numpoints,
        exp.distance, first, nothing;
        return_inf_on_error = true,
        eval_timeout = exp.timeout
    )
end
```

---

## Questions?

For issues or questions:
- Review code documentation in `src/error_metrics.jl`
- Check model definitions in `src/systems/*.jl`
- See `README.md` for basic usage
- See `MIGRATION.md` for integration with globtimcore
