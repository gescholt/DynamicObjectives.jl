"""
Integration utilities for using DynamicObjectives with globtim optimizer.

# Architecture - Standalone Design

DynamicObjectives is **standalone** — Globtim and GlobtimPostProcessing are
**weak dependencies**. The TOML pipeline orchestrator `run_experiment_from_config`
is provided by the `DynamicObjectivesGlobtimExt` package extension, which loads
automatically when both `Globtim` and `GlobtimPostProcessing` are imported.

Integration is achieved through:
- Package extension (ext/DynamicObjectivesGlobtimExt.jl)
- Local dev setup using `Pkg.develop()` (see `setup_dev_packages.jl`)
- File-based data exchange (CSV files)

See ARCHITECTURE.md and REFINEMENT_INTEGRATION_GUIDE.md for complete integration patterns.

# Usage: 1-arg functions work directly with globtim

```julia
# DynamicObjectives creates 1-arg functions
objective = make_error_distance(model, outputs, ic, p_true, [0.0, 20.0], 30,
                                L2_norm; eval_timeout=10.0)

# Works directly with globtim — no wrapper needed
using Globtim
result = run_standard_experiment(
    objective_function = objective, objective_name = "my_model",
    bounds = bounds,
    experiment_config = config, output_dir = "results/my_experiment"
)
```

# TOML-driven pipeline

```julia
using DynamicObjectives
using Globtim                  # ← triggers DynamicObjectivesGlobtimExt
using GlobtimPostProcessing    # ← triggers DynamicObjectivesGlobtimExt

result = run_experiment_from_config("experiments/lv4d_basic.toml")
```

# ForwardDiff Compatibility

DynamicObjectives supports ForwardDiff automatic differentiation through ODE objectives.
Parameter setting uses `SciMLBase.remake` which preserves ForwardDiff.Dual types.

For globtim integration, ForwardDiff-based gradient/Hessian computation is available
but may be slower than FiniteDiff for high-dimensional models due to chunk size scaling.

Created: 2025-11-22
Updated: 2026-02-09 (TOML pipeline via package extension)
"""

# ═══════════════════════════════════════════════════════════════════════════════
# Stub — implemented by DynamicObjectivesGlobtimExt when Globtim is loaded
# ═══════════════════════════════════════════════════════════════════════════════

"""
    run_experiment_from_config(path::String; io::IO=stdout) -> Dict{Symbol, Any}

Run a complete experiment pipeline from a TOML configuration file.

**Requires** `using Globtim`, `using GlobtimPostProcessing`, and `using Optim` to be loaded
(activates the `DynamicObjectivesGlobtimExt` package extension).

See the extension module for full documentation.
"""
function run_experiment_from_config end

"""
    build_experiment_objective(config; io::IO=stdout)
        -> (; objective, bounds, obj_name, entry, p_true)

Resolve a parsed TOML experiment config (`Globtim.ExperimentPipelineConfig`)
into an objective function and domain bounds — the canonical model-resolution
step shared by `run_experiment_from_config` and the per-axis audit /
counterfactual cluster drivers.

**Requires** `using Globtim` (activates the narrow
`DynamicObjectivesGlobtimCoreExt` package extension — GlobtimPostProcessing /
Optim are NOT needed).

See the extension module for full documentation of the resolution semantics.
"""
function build_experiment_objective end

# ═══════════════════════════════════════════════════════════════════════════════
# Solver resolution — shared between core and extension
# ═══════════════════════════════════════════════════════════════════════════════

"""Map solver method name strings to OrdinaryDiffEq solver constructors."""
const _SOLVER_MAP = Dict{String,Function}(
    "Tsit5" => () -> Tsit5(),
    "Vern7" => () -> Vern7(),
    "Vern9" => () -> Vern9(),
    "Rodas5" => () -> Rodas5(),
    "Rosenbrock23" => () -> Rosenbrock23(),
    "TRBDF2" => () -> TRBDF2(),
    "KenCarp4" => () -> KenCarp4(),
    "AutoTsit5" => () -> AutoTsit5(Rosenbrock23()),
)

"""
    _resolve_solver(method::String)

Resolve a solver method name to an OrdinaryDiffEq solver instance.
"""
function _resolve_solver(method::String)
    haskey(_SOLVER_MAP, method) || error(
        "Unknown solver method: \"$method\". " *
        "Known: $(join(sort(collect(keys(_SOLVER_MAP))), ", "))",
    )
    return _SOLVER_MAP[method]()
end
