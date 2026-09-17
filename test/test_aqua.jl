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
    #
    # persistent_tasks precompiles the package inside a wrapper env that Pkg
    # resolves afresh from the registry. Under Pkg.test that resolve is the test
    # sandbox's and is already precompiled (the whole Aqua suite takes ~47 s on
    # the public mirror's CI). The monorepo's workspace CI job runs this file
    # against profiles/cluster precisely to avoid a fresh resolve (see
    # .github/workflows/ci.yml, test-dynamic_objectives): there the wrapper
    # resolves to versions the workspace never built (OrdinaryDiffEqCore 3.28
    # against the workspace's 3.1, JumpProcesses 9.29 against 9.23 — the pair
    # that job's comment names), precompiles the SciML stack from scratch in a
    # subprocess whose output Aqua discards, and died after 13 min on 2026-09-17
    # (Julia 1.12.7, Aqua 0.8.16). That job sets DYNOBJ_AQUA_PERSISTENT_TASKS=0;
    # everywhere else the check runs.
    persistent_tasks_opt = if get(ENV, "DYNOBJ_AQUA_PERSISTENT_TASKS", "1") == "0"
        @info "Skipping Aqua persistent_tasks — DYNOBJ_AQUA_PERSISTENT_TASKS=0 (workspace CI; the check runs under Pkg.test)"
        false
    else
        # tmax raised from Aqua's 10s default: this package's subprocess loads the
        # whole SciML stack (~2 min), so the default is nowhere near enough and a
        # slow load would be misread as a lingering task. See the 2026-08-13 flake
        # note in the sibling packages' test_aqua.jl.
        (; tmax = 300)
    end
    Aqua.test_all(
        DynamicObjectives;
        piracies = (; treat_as_own = [BinaryHeap]),
        persistent_tasks = persistent_tasks_opt,
    )
end
