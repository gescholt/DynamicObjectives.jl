# Examples

Self-contained demos. Run each from the package root in an environment that has
DynamicObjectives and the packages listed for it. The demos print through the
display functions, which are package extensions: tables need PrettyTables, panels
and headers need Term, line plots need UnicodePlots (see the README, "Optional
packages").

| File | What it shows | Also needs |
|---|---|---|
| `quick_start.jl` | One model end to end: objective from an ODE model, Globtim critical points, refinement. | Globtim, GlobtimPostProcessing, PrettyTables, Term |
| `basic_workflow.jl` | Building and checking an objective by hand: model, parameters, solver, error function. | PrettyTables, Term, UnicodePlots |
| `run_screening_demo.jl` | The `p_true` screening pipeline on DAISY 4D: sweep, probe, rank, save and reload a catalogue (written to `tempdir()`). | PrettyTables, Term |
| `exp_candidates/lv4d.toml` | A screening config for `run_screening_from_config`. | |
| `configs/lv4d_basic.toml` | An experiment config for `run_experiment_from_config`, on the shipped `paper/catalogue`. | Globtim, GlobtimPostProcessing, Optim |
| `figures/l2_vs_l2squared.jl` | Regenerates the README figure, L2 against L2-squared misfit. | CairoMakie |

```bash
julia --project=. examples/run_screening_demo.jl
```
