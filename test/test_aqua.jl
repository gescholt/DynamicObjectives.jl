using Test
using Aqua
using DynamicObjectives

# Reached through DynamicObjectives (which does `using DataStructures`) rather
# than imported directly, so this file does not need DataStructures resolvable in
# whatever environment the tests are driven from.
const BinaryHeap = DynamicObjectives.DataStructures.BinaryHeap

@testset "Aqua.jl Quality Assurance" begin
    # Full Aqua suite. `persistent_tasks` was previously skipped on the theory that
    # ModelingToolkit/SciML globals trip it; that is not true on Julia 1.12.6 /
    # Aqua 0.8.14 — it passes. It is slow here (~2 min: the check spawns a
    # subprocess that loads the whole SciML stack from scratch), but it passes.
    #
    # The one real exclusion is piracy: src/compat_patches.jl defines
    # `Base.empty!(::DataStructures.BinaryHeap)` to work around a DataStructures /
    # OrdinaryDiffEqCore version skew. Re-verified 2026-08-10 — DataStructures
    # still ships no such method (ours is the only one in the method table), and
    # OrdinaryDiffEqCore.reinit! calls it on the integrator-reuse path, so the
    # shim cannot be dropped yet. Scoped to BinaryHeap so piracy stays checked
    # everywhere else; retire this kwarg together with compat_patches.jl.
    #
    # undocumented_names is off by Aqua's own default and stays off for now.
    Aqua.test_all(
        DynamicObjectives;
        piracies = (; treat_as_own = [BinaryHeap]),
    )
end
