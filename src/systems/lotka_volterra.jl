# Lotka-Volterra Models
# Various Lotka-Volterra predator-prey and competition models

"""
Generalized 4D Lotka-Volterra model with full interaction matrix.

System equations:
    dx_i/dt = x_i * (a_i + Σ_j(B_ij * x_j))

where:
- a_i are growth rates (can be positive or negative)
- B_ij are interaction coefficients (diagonal = self-regulation, off-diagonal = species interactions)

Parameters (20 total):
- Growth rates: a1, a2, a3, a4
- Interaction matrix: B11, B12, B13, B14, B21, B22, B23, B24, B31, B32, B33, B34, B41, B42, B43, B44

Returns:
- model: ODESystem
- parameters: Vector of 20 parameters [a1, a2, a3, a4, B11, B12, ..., B44]
- states: Vector of 4 state variables [x1, x2, x3, x4]
- measured_quantities: All 4 species are observable [y1 ~ x1, y2 ~ x2, y3 ~ x3, y4 ~ x4]
"""
function define_generalized_lotka_volterra_4D()
    @independent_variables t
    @parameters a1 a2 a3 a4 B11 B12 B13 B14 B21 B22 B23 B24 B31 B32 B33 B34 B41 B42 B43 B44
    @variables x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
    D = Differential(t)

    states = [x1, x2, x3, x4]
    parameters = [a1, a2, a3, a4, B11, B12, B13, B14, B21, B22, B23, B24, B31, B32, B33, B34, B41, B42, B43, B44]

    @mtkbuild model = ODESystem(
        [
            D(x1) ~ x1 * (a1 + B11*x1 + B12*x2 + B13*x3 + B14*x4),
            D(x2) ~ x2 * (a2 + B21*x1 + B22*x2 + B23*x3 + B24*x4),
            D(x3) ~ x3 * (a3 + B31*x1 + B32*x2 + B33*x3 + B34*x4),
            D(x4) ~ x4 * (a4 + B41*x1 + B42*x2 + B43*x3 + B44*x4)
        ],
        t, states, parameters
    )

    measured_quantities = [y1 ~ x1, y2 ~ x2, y3 ~ x3, y4 ~ x4]
    model, parameters, states, measured_quantities
end

"""
Constrained 4D Lotka-Volterra model with skew-symmetric perturbations.

System equations:
    dx_i/dt = x_i * (a_i + Σ_j(B_ij * x_j))

Fixed interaction matrix structure:
    B_fixed = [
        0.2   -0.6  0.0   0.0;
        0.6    0.05 0.0   0.0;
        0.0    0.0  0.2  -0.6;
        0.0    0.0  0.6   0.05
    ]

Perturbed matrix with skew-symmetric constraints:
    B = [
        0.2   -0.6   ε₁    ε₂;
        0.6    0.05  ε₃    ε₄;
       -ε₁    -ε₃   0.2   -0.6;
       -ε₂    -ε₄   0.6    0.05
    ]

Known Parameters (Fixed):
- Growth rates: a = [-0.5, 1.0, -0.5, 1.0]

Unknown Parameters (4 total - to be inferred):
- Epsilon perturbations: eps1, eps2, eps3, eps4

Returns:
- model: ODESystem
- parameters: Vector of 4 parameters [eps1, eps2, eps3, eps4]
- states: Vector of 4 state variables [x1, x2, x3, x4]
- measured_quantities: All 4 species are observable [y1 ~ x1, y2 ~ x2, y3 ~ x3, y4 ~ x4]
"""
function define_constrained_lotka_volterra_4D()
    @independent_variables t
    @parameters eps1 eps2 eps3 eps4
    @variables x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
    D = Differential(t)

    states = [x1, x2, x3, x4]
    parameters = [eps1, eps2, eps3, eps4]

    # Fixed known parameters (growth rates)
    a1, a2, a3, a4 = -0.5, 1.0, -0.5, 1.0

    # Fixed B matrix elements (not perturbed)
    B11, B12, B21, B22 = 0.2, -0.6, 0.6, 0.05
    B33, B34, B43, B44 = 0.2, -0.6, 0.6, 0.05

    @mtkbuild model = ODESystem(
        [
            D(x1) ~ x1 * (a1 + B11*x1 + B12*x2 + eps1*x3 + eps2*x4),
            D(x2) ~ x2 * (a2 + B21*x1 + B22*x2 + eps3*x3 + eps4*x4),
            D(x3) ~ x3 * (a3 + (-eps1)*x1 + (-eps3)*x2 + B33*x3 + B34*x4),
            D(x4) ~ x4 * (a4 + (-eps2)*x1 + (-eps4)*x2 + B43*x3 + B44*x4)
        ],
        t, states, parameters
    )

    measured_quantities = [y1 ~ x1, y2 ~ x2, y3 ~ x3, y4 ~ x4]
    model, parameters, states, measured_quantities
end

"""
3D Lotka-Volterra model (variant 1)

Fixed parameter: c (implicit in equations)

System equations:
    dx1/dt = a*x1 + b*x1*x2
    dx2/dt = b*x1*x2 + c*x2

Measured quantities:
    y1 = x1

Returns:
- model: ODESystem
- parameters: [a, b, c]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lotka_volterra_3D_model()
    @independent_variables t
    @variables x1(t) x2(t) y1(t)
    @parameters a b c
    D = Differential(t)
    params = [a, b, c]
    states = [x1, x2]
    @mtkbuild model = ODESystem(
        [D(x1) ~ a * x1 + b * x1 * x2, D(x2) ~ b * x1 * x2 + c * x2],
        t,
        states,
        params
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
3D Lotka-Volterra model (variant 2)

System equations:
    dx1/dt = a*x1 - b*x1*x2
    dx2/dt = -b*x2 + c*x1*x2

Measured quantities:
    y1 = x1

Returns:
- model: ODESystem
- parameters: [a, b, c]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lotka_volterra_3D_model_v2()
    @independent_variables t
    @variables x1(t) x2(t) y1(t)
    @parameters a b c
    D = Differential(t)
    params = [a, b, c]
    states = [x1, x2]
    @mtkbuild model = ODESystem(
        [D(x1) ~ a * x1 + -b * x1 * x2,
            D(x2) ~ -b * x2 + c * x1 * x2],
        t,
        states,
        params
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
2D Lotka-Volterra model (c = 1)

System equations:
    dx1/dt = a*x1 + b*x1*x2
    dx2/dt = b*x1*x2 + x2

Measured quantities:
    y1 = x1

Returns:
- model: ODESystem
- parameters: [a, b]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lotka_volterra_2D_model()
    @independent_variables t
    @variables x1(t) x2(t) y1(t)
    @parameters a b
    D = Differential(t)
    params = [a, b]
    states = [x1, x2]
    @mtkbuild model = ODESystem(
        [D(x1) ~ a * x1 + b * x1 * x2, D(x2) ~ b * x1 * x2 + x2],
        t,
        states,
        params
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
2D Lotka-Volterra model v3 (c = 0.5)

System equations:
    dx1/dt = a*x1 - b*x1*x2
    dx2/dt = -b*x2 + 0.5*x1*x2

Measured quantities:
    y1 = x1

Returns:
- model: ODESystem
- parameters: [a, b]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lotka_volterra_2D_model_v3()
    @independent_variables t
    @variables x1(t) x2(t) y1(t)
    @parameters a b
    D = Differential(t)
    params = [a, b]
    states = [x1, x2]
    @mtkbuild model = ODESystem(
        [D(x1) ~ a * x1 + -b * x1 * x2,
            D(x2) ~ -b * x2 + 0.5 * x1 * x2],
        t,
        states,
        params
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
2D Lotka-Volterra model v3 with two outputs (c = 0.5)

System equations:
    dx1/dt = a*x1 - b*x1*x2
    dx2/dt = -b*x2 + 0.5*x1*x2

Measured quantities:
    y1 = x1
    y2 = x2

Returns:
- model: ODESystem
- parameters: [a, b]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1, y2 ~ x2]
"""
function define_lotka_volterra_2D_model_v3_two_outputs()
    @independent_variables t
    @variables x1(t) x2(t) y1(t) y2(t)
    @parameters a b
    D = Differential(t)
    params = [a, b]
    states = [x1, x2]
    @mtkbuild model = ODESystem(
        [D(x1) ~ a * x1 + -b * x1 * x2,
            D(x2) ~ -b * x2 + 0.5 * x1 * x2],
        t,
        states,
        params
    )
    outputs = [y1 ~ x1, y2 ~ x2]
    return model, params, states, outputs
end

"""
2D Lotka-Volterra model v2 (c = 0.1)

System equations:
    dx1/dt = a*x1 + b*x1*x2
    dx2/dt = b*x1*x2 + 0.1*x2

Measured quantities:
    y1 = x1

Returns:
- model: ODESystem
- parameters: [a, b]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lotka_volterra_2D_model_v2()
    @independent_variables t
    @variables x1(t) x2(t) y1(t)
    @parameters a b
    D = Differential(t)
    params = [a, b]
    states = [x1, x2]
    @mtkbuild model = ODESystem(
        [
            D(x1) ~ a * x1 + b * x1 * x2,
            D(x2) ~ b * x1 * x2 + 0.1 * x2],
        t,
        states,
        params
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end
