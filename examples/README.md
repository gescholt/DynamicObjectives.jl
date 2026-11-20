# Examples

Self-contained demos. Run each from the package root in an environment that has
DynamicObjectives (and, where noted, Globtim, GlobtimPostProcessing and Optim):

| File | What it shows |
|---|---|
| `quick_start.jl` | One model end to end: objective from an ODE model, Globtim critical points, refinement (needs Globtim, GlobtimPostProcessing). |
| `basic_workflow.jl` | Building and checking an objective by hand: model, parameters, solver, error function. |
| `run_screening_demo.jl` | The `p_true` screening pipeline on DAISY 4D: sweep, probe, rank, save and reload a catalogue (written to `tempdir()`). |
| `exp_candidates/lv4d.toml` | A screening config for `run_screening_from_config`. |
| `configs/lv4d_basic.toml` | An experiment config for `run_experiment_from_config`, on the shipped `paper/catalogue`. |
| `figures/l2_vs_l2squared.jl` | Regenerates the README figure, L2 against L2-squared misfit (needs CairoMakie). |

```bash
julia --project=. examples/run_screening_demo.jl
```
