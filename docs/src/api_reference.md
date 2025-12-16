# API Reference

Complete reference for the Dynamic\_objectives package.

## Model Definition Functions

All model definition functions return a tuple:
```julia
(model, parameters, states, measured_quantities)
```

Where:
- `model`: ModelingToolkit ODESystem
- `parameters`: Vector of symbolic parameters
- `states`: Vector of symbolic state variables
- `measured_quantities`: Vector of measurement equations

### Lotka-Volterra Models

```@docs
define_lotka_volterra_2D_model
define_lotka_volterra_2D_model_v2
define_lotka_volterra_2D_model_v3
define_lotka_volterra_2D_model_v3_two_outputs
define_lotka_volterra_3D_model
define_lotka_volterra_3D_model_v2
define_generalized_lotka_volterra_4D
define_constrained_lotka_volterra_4D
```

### DAISY Benchmark Models

```@docs
define_daisy_ex3_model_4D
define_daisy_ex3_model_4D_no_input
```

### Other Models

```@docs
define_fitzhugh_nagumo_3D_model
define_simple_1D_model_locally_identifiable
define_simple_2D_model_locally_identifiable
define_simple_2D_model_locally_identifiable_square
```

---

## Error Metrics

### Main Function

```@docs
make_error_distance
```

### Distance Functions

```@docs
L2_norm
L1_norm
log_L2_norm
```

---

## Data Generation

```@docs
sample_data
```

---

## Display Functions

### Model Display

```@docs
display_model_summary
display_parameters
```

### Data Display

```@docs
display_time_series
display_comparison
```

### Results Display

```@docs
display_error_metrics
display_optimization_result
display_optimization_progress
```

---

## Exported Symbols

The following symbols are exported from `Dynamic_objectives`:

### Model Definitions

- `define_lotka_volterra_2D_model`
- `define_lotka_volterra_2D_model_v2`
- `define_lotka_volterra_2D_model_v3`
- `define_lotka_volterra_2D_model_v3_two_outputs`
- `define_lotka_volterra_3D_model`
- `define_lotka_volterra_3D_model_v2`
- `define_generalized_lotka_volterra_4D`
- `define_constrained_lotka_volterra_4D`
- `define_daisy_ex3_model_4D`
- `define_daisy_ex3_model_4D_no_input`
- `define_fitzhugh_nagumo_3D_model`
- `define_simple_1D_model_locally_identifiable`
- `define_simple_2D_model_locally_identifiable`
- `define_simple_2D_model_locally_identifiable_square`

### Error Functions

- `make_error_distance`
- `L2_norm`
- `L1_norm`
- `log_L2_norm`

### Data Generation

- `sample_data`

### Display Functions

- `display_model_summary`
- `display_parameters`
- `display_time_series`
- `display_comparison`
- `display_error_metrics`
- `display_optimization_result`
- `display_optimization_progress`

---

## Usage Examples

### Creating an Objective Function

```julia
using Dynamic_objectives

# Define model
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()

# Set parameters
p_true = [1.0, 0.5]
ic = [1.0, 0.5]
time_interval = [0.0, 20.0]
numpoints = 30

# Create objective function
error_func = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    numpoints,
    L2_norm,
    first,
    nothing;
    return_inf_on_error = true,
    eval_timeout = 10.0
)

# Test
@assert error_func(p_true) < 1e-6
```

### Using Display Functions

```julia
using Dynamic_objectives

# Define model
model, params, states, outputs = define_daisy_ex3_model_4D()

# Display model summary
display_model_summary(model)

# Display parameters
param_names = [:p1, :p3, :p4, :p6]
p_true = [0.5, 0.3, 0.4, 0.2]
p_test = [0.55, 0.28, 0.42, 0.18]
display_parameters(param_names, [p_true p_test], labels=["True", "Test"])

# Display error metrics
errors = Dict("L2" => 0.045, "L1" => 0.123)
display_error_metrics(errors)
```

### Custom Distance Functions

```julia
# Custom weighted L2 norm
function weighted_L2(y_true, y_test)
    weights = exp.(-0.1 * (0:length(y_true)-1))
    return norm((y_true - y_test) .* weights, 2)
end

# Use with make_error_distance
error_func = make_error_distance(
    model, outputs, ic, p_true, time_interval, numpoints,
    weighted_L2, first, nothing
)
```

---

## See Also

- [Getting Started](getting_started.md) - Installation and basic usage
- [Model Catalog](model_catalog.md) - Complete model list
- [Testing Guide](testing_guide.md) - Recommended configurations
