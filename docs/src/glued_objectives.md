# Glued Objectives

`glue`, `glue_catalogue`, and `count_oracle_cps_in_box` build higher-dimensional test
problems with **known** critical-point sets by Cartesian-product gluing: for
`F(z) = f₁(x) + f₂(y)`, the stationary points factorize as `CP(F) = CP(f₁) × CP(f₂)` and
the Hessian is block-diagonal, so the Morse index adds.

The construction, the worked examples, the choice of factor pairs, and the verified
recall results (himmelblau² and camel² — recall = precision = 1.0 at every Morse index)
are documented with the solver they certify:

**→ [Glued Objectives (Globtim docs)](https://gescholt.github.io/Globtim.jl/dev/glued_objectives/)**

For the docstrings of `glue`, `glue_catalogue`, and `count_oracle_cps_in_box`, see the
[API Reference](api.md). Tests: `test/test_glued_objectives.jl` (construction, oracle-CP
gradient-zero check, sub-box counting) and `test/test_glued_recovery.jl` (CI recall
regression).
