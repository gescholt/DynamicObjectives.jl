# Per-Model Test Scripts

Individual test scripts for each model, using the full 2-stage optimization pipeline.

## Usage

Run any model test directly:

```bash
cd Dynamic_objectives
julia --project=. examples/models/test_lv_4d_constrained.jl
julia --project=. examples/models/test_daisy_ex3_with_input.jl
julia --project=. examples/models/test_daisy_ex3_no_input.jl
```

## Available Models

| Script | Model | Parameters | Description |
|--------|-------|------------|-------------|
| `test_lv_4d_constrained.jl` | Constrained LV 4D | 4 | Skew-symmetric perturbations |
| `test_daisy_ex3_with_input.jl` | DAISY Ex3 | 4 | Time-dependent input |
| `test_daisy_ex3_no_input.jl` | DAISY Ex3 | 4 | No input variant |

## Output

Each test generates:
- Rich terminal output with progress bars
- `test_results/models/<model_name>/report.md` - Detailed markdown report
- `test_results/models/<model_name>/` - Raw globtim outputs

## Report Contents

Each report includes:
- Model configuration and parameters
- Stage 1 (globtim) results per degree
- Stage 2 (refinement) results
- **Top 5 critical points** with objective values and parameter errors
- Gradient validation statistics
- Parameter recovery metrics
- Timing breakdown

## Customization

Edit the `ModelTestConfig` in each script to adjust:
- `GN` - Grid resolution (higher = more thorough, slower)
- `degree_range` - Polynomial degrees to try
- `max_time` - Max time for Stage 1
- `max_time_per_point` - Refinement time per point

## Shared Utilities

Common functionality is in `model_test_utils.jl`:
- `ModelTestConfig` - Configuration struct
- `run_model_test()` - Main test runner
- `generate_test_report()` - Markdown report generator
