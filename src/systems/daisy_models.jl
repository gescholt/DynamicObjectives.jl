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
