"""
    run_aqua_tests(mod)

Run Aqua.jl quality assurance checks on `mod`.

Skips `test_persistent_tasks` to match the globtim convention; ModelingToolkit
+ SciML upstream globals routinely trip that check.
"""
function run_aqua_tests(mod)
    Aqua.test_ambiguities(mod)
    Aqua.test_undefined_exports(mod)
    Aqua.test_unbound_args(mod)
    Aqua.test_stale_deps(mod)
    Aqua.test_deps_compat(mod)
end
