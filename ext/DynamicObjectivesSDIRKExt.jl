"""
DynamicObjectivesSDIRKExt — the SDIRK stiff solvers TRBDF2 and KenCarp4

Loaded when OrdinaryDiffEqSDIRK is loaded next to DynamicObjectives. It registers
the solver names `"TRBDF2"` and `"KenCarp4"` in `DynamicObjectives._SOLVER_MAP`, so
TOML configs and `_resolve_solver` can name them. Without this extension those two
names raise an error that says which package to load.

Triggered by: using DynamicObjectives; using OrdinaryDiffEqSDIRK
"""
module DynamicObjectivesSDIRKExt

using DynamicObjectives
using OrdinaryDiffEqSDIRK: TRBDF2, KenCarp4

function __init__()
    DynamicObjectives._SOLVER_MAP["TRBDF2"] = () -> TRBDF2()
    DynamicObjectives._SOLVER_MAP["KenCarp4"] = () -> KenCarp4()
    return nothing
end

end # module DynamicObjectivesSDIRKExt
