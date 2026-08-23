# DAISY Example Models
# Models from DAISY (Differential Algebra for Identifiability of SYstems) benchmark suite

"""
DAISY Example 3 Model (4D with input)

Fixed parameter: p7 = -0.2

System equations:
    dx1/dt = -p1*x1 + x2 + u0
    dx2/dt = p3*x1 - p4*x2 + x3
    dx3/dt = p6*x1 + 0.2*x3
    du0/dt = 1

Measured quantities:
    y1 = x1 + x3
    y2 = x2

Returns:
- model: System
- parameters: [p1, p3, p4, p6]
- states: [x1, x2, x3, u0]
- measured_quantities: [y1 ~ x1 + x3, y2 ~ x2]
"""
function define_daisy_ex3_model_4D()
    @independent_variables t
    @parameters p1 p3 p4 p6
    @variables x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
    D = Differential(t)

    states = [x1, x2, x3, u0]
    parameters = [p1, p3, p4, p6]
    @mtkcompile model = System(
        [
            D(x1) ~ -1 * p1 * x1 + x2 + u0,
            D(x2) ~ p3 * x1 - p4 * x2 + x3,
            D(x3) ~ p6 * x1 + 0.2 * x3,
            D(u0) ~ 1,
        ],
        t,
        states,
        parameters,
    )
    measured_quantities = [y1 ~ x1 + x3, y2 ~ x2]
    model, parameters, states, measured_quantities
end

"""
DAISY Example 3 Model (4D without input)

Fixed parameter: p7 = -0.2

System equations:
    dx1/dt = -p1*x1 + x2
    dx2/dt = p3*x1 - p4*x2 + x3
    dx3/dt = p6*x1 + 0.2*x3

Measured quantities:
    y1 = x1 + x3
    y2 = x2

Returns:
- model: System
- parameters: [p1, p3, p4, p6]
- states: [x1, x2, x3]
- measured_quantities: [y1 ~ x1 + x3, y2 ~ x2]
"""
function define_daisy_ex3_model_4D_no_input()
    @independent_variables t
    @parameters p1 p3 p4 p6
    @variables x1(t) x2(t) x3(t) y1(t) y2(t)
    D = Differential(t)

    states = [x1, x2, x3]
    parameters = [p1, p3, p4, p6]
    @mtkcompile model = System(
        [
            D(x1) ~ -1 * p1 * x1 + x2,
            D(x2) ~ p3 * x1 - p4 * x2 + x3,
            D(x3) ~ p6 * x1 + 0.2 * x3,
        ],
        t,
        states,
        parameters,
    )
    measured_quantities = [y1 ~ x1 + x3, y2 ~ x2]
    model, parameters, states, measured_quantities
end

"""
DAISY Example 3 Model (5D) — the 4D variant with the fifth parameter freed.

The 4D model fixes the original benchmark's fifth parameter at p7 = -1/5,
which appears as the hard-coded `+ 0.2 * x3` term. Here p7 is estimated
(DAISY sign convention: the term is `- p7 * x3`), so the p7 = -0.2 slice of
this model reproduces `define_daisy_ex3_model_4D` exactly — the identity the
DAISY-5D smoke asserts. Defined for the m = 5 dimension-wall measurement
(bead jucj): the paper's Bézout extrapolation at m = 5, d = 8 is 16807 paths.

System equations:
    dx1/dt = -p1*x1 + x2 + u0
    dx2/dt = p3*x1 - p4*x2 + x3
    dx3/dt = p6*x1 - p7*x3
    du0/dt = 1

Measured quantities:  y1 = x1 + x3,  y2 = x2

Returns: model, parameters [p1, p3, p4, p6, p7], states, measured_quantities
"""
function define_daisy_ex3_model_5D()
    @independent_variables t
    @parameters p1 p3 p4 p6 p7
    @variables x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
    D = Differential(t)

    states = [x1, x2, x3, u0]
    parameters = [p1, p3, p4, p6, p7]
    @mtkcompile model = System(
        [
            D(x1) ~ -1 * p1 * x1 + x2 + u0,
            D(x2) ~ p3 * x1 - p4 * x2 + x3,
            D(x3) ~ p6 * x1 - p7 * x3,
            D(u0) ~ 1,
        ],
        t,
        states,
        parameters,
    )
    measured_quantities = [y1 ~ x1 + x3, y2 ~ x2]
    model, parameters, states, measured_quantities
end

"""
DAISY Example 3 Model (6D) — the 5D variant with a sixth parameter freed.

The 4D/5D variants fix the original benchmark's p2 at 1, which appears as the
bare `+ x2` term in dx1. Here p2 is estimated, so the p2 = 1 slice of this
model reproduces `define_daisy_ex3_model_5D` exactly — the identity the
DAISY-6D smoke asserts. Parameter ORDER follows the append convention of the
m-ladder (jucj): the new parameter goes LAST, so p_true/bounds/results extend
the 5D vectors by one coordinate and stay slice-comparable.

System equations:
    dx1/dt = -p1*x1 + p2*x2 + u0
    dx2/dt = p3*x1 - p4*x2 + x3
    dx3/dt = p6*x1 - p7*x3
    du0/dt = 1

Measured quantities:  y1 = x1 + x3,  y2 = x2

Returns: model, parameters [p1, p3, p4, p6, p7, p2], states, measured_quantities
"""
function define_daisy_ex3_model_6D()
    @independent_variables t
    @parameters p1 p3 p4 p6 p7 p2
    @variables x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
    D = Differential(t)

    states = [x1, x2, x3, u0]
    parameters = [p1, p3, p4, p6, p7, p2]
    @mtkcompile model = System(
        [
            D(x1) ~ -1 * p1 * x1 + p2 * x2 + u0,
            D(x2) ~ p3 * x1 - p4 * x2 + x3,
            D(x3) ~ p6 * x1 - p7 * x3,
            D(u0) ~ 1,
        ],
        t,
        states,
        parameters,
    )
    measured_quantities = [y1 ~ x1 + x3, y2 ~ x2]
    model, parameters, states, measured_quantities
end

"""
DAISY Example 3 Model (7D) — the FULL benchmark, all seven parameters freed.

The 6D variant still fixes p5 at 1, which appears as the bare `+ x3` term in
dx2. Here p5 is estimated, so the p5 = 1 slice of this model reproduces
`define_daisy_ex3_model_6D` exactly — the identity the DAISY-7D smoke
asserts. This is the terminus of the m-ladder (4 -> 5 -> 6 -> 7): the
complete original DAISY ex3 with every parameter free. Parameter ORDER
follows the append convention: [p1, p3, p4, p6, p7, p2, p5].

System equations:
    dx1/dt = -p1*x1 + p2*x2 + u0
    dx2/dt = p3*x1 - p4*x2 + p5*x3
    dx3/dt = p6*x1 - p7*x3
    du0/dt = 1

Measured quantities:  y1 = x1 + x3,  y2 = x2

Returns: model, parameters [p1, p3, p4, p6, p7, p2, p5], states, measured_quantities
"""
function define_daisy_ex3_model_7D()
    @independent_variables t
    @parameters p1 p3 p4 p6 p7 p2 p5
    @variables x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
    D = Differential(t)

    states = [x1, x2, x3, u0]
    parameters = [p1, p3, p4, p6, p7, p2, p5]
    @mtkcompile model = System(
        [
            D(x1) ~ -1 * p1 * x1 + p2 * x2 + u0,
            D(x2) ~ p3 * x1 - p4 * x2 + p5 * x3,
            D(x3) ~ p6 * x1 - p7 * x3,
            D(u0) ~ 1,
        ],
        t,
        states,
        parameters,
    )
    measured_quantities = [y1 ~ x1 + x3, y2 ~ x2]
    model, parameters, states, measured_quantities
end
