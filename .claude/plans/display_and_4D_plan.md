# Display Improvements & 4D Testing Plan

## Overview

Standardize display output for critical point quality analysis and prepare for n=4 dimension problems using LV 4D Constrained model.

---

## Phase 1: Display Functions

**Goal**: Add display functions for gradient validation and quality metrics

### Implementation Location
**Package**: `Dynamic_objectives`
**File**: `src/display.jl` (extend existing)

### Functions to Add

```julia
"""
Display gradient norm analysis table
"""
function display_gradient_analysis(
    norms::Vector{Float64};
    tolerance::Float64 = 1e-6,
    title::String = "Gradient Norm Analysis"
)

"""
Display critical point quality summary (combined metrics)
"""
function display_quality_summary(
    refined_results::NamedTuple;
    title::String = "Critical Point Quality"
)

"""
Display degree-by-degree comparison table
"""
function display_degree_comparison(
    degree_results::Vector;
    title::String = "Degree Comparison"
)
```

### Output Format

```
Gradient Norm Analysis
┌──────────────────┬────────────┐
│ Metric           │      Value │
├──────────────────┼────────────┤
│ Valid points     │    78/81   │
│ Invalid points   │     3/81   │
│ Min ||∇f||       │   1.2e-12  │
│ Mean ||∇f||      │   3.4e-09  │
│ Max ||∇f||       │   2.1e-05  │
│ Tolerance        │   1.0e-06  │
└──────────────────┴────────────┘
```

---

## Phase 2: Update Integration Test

**Goal**: Incorporate gradient validation into the 2-stage workflow

### File to Modify
`examples/test_globtim_integration.jl`

### Changes
1. After refinement, compute gradient norms
2. Display gradient validation results
3. Include in final quality summary

---

## Phase 3: 4D Model Setup (LV 4D Constrained)

**Goal**: Prepare and test 4D parameter estimation

### Model Details
- **Model**: `define_constrained_lotka_volterra_4D()`
- **Parameters**: 4 (p_true = [1.0, 0.5, 1.0, 0.5])
- **States**: 4 (ic = [1.0, 0.5, 0.5, 0.5])
- **Bounds**: [(0.0, 3.0), (0.0, 2.0), (0.0, 3.0), (0.0, 2.0)]

### Configuration for 4D
```julia
config = ExperimentParams(
    domain_size = 1.5,
    GN = 10,           # 10^4 = 10,000 grid points
    degree_range = 6:10,  # Lower degrees for 4D
    max_time = 1800.0,    # 30 minutes
    basis = :chebyshev
)

refinement_config = ode_refinement_config(
    max_time_per_point = 30.0,  # More time for 4D
    show_progress = true
)
```

### File to Create
`examples/test_globtim_4D.jl`

---

## Files to Modify/Create

| File | Action |
|------|--------|
| `src/display.jl` | MODIFY (add gradient display) |
| `examples/test_globtim_integration.jl` | MODIFY (add gradient display) |
| `examples/test_globtim_4D.jl` | CREATE |

---

## Success Criteria

1. Display shows valid/invalid count with tolerance
2. 2D test passes with gradient validation
3. 4D test runs and shows progress
4. 4D achieves parameter recovery (may need tuning)
