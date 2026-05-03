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
- model: System
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
    parameters = [
        a1,
        a2,
        a3,
        a4,
        B11,
        B12,
        B13,
        B14,
        B21,
        B22,
        B23,
        B24,
        B31,
        B32,
        B33,
        B34,
        B41,
        B42,
        B43,
        B44,
    ]

    @mtkcompile model = System(
        [
            D(x1) ~ x1 * (a1 + B11 * x1 + B12 * x2 + B13 * x3 + B14 * x4),
            D(x2) ~ x2 * (a2 + B21 * x1 + B22 * x2 + B23 * x3 + B24 * x4),
            D(x3) ~ x3 * (a3 + B31 * x1 + B32 * x2 + B33 * x3 + B34 * x4),
            D(x4) ~ x4 * (a4 + B41 * x1 + B42 * x2 + B43 * x3 + B44 * x4),
        ],
        t,
        states,
        parameters,
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
       -0.5   -0.6  0.0   0.0;
        0.6   -0.5  0.0   0.0;
        0.0    0.0 -0.5  -0.6;
        0.0    0.0  0.6  -0.5
    ]

Perturbed matrix with skew-symmetric constraints:
    B = [
       -0.5   -0.6   ε₁    ε₂;
        0.6   -0.5   ε₃    ε₄;
       -ε₁    -ε₃  -0.5  -0.6;
       -ε₂    -ε₄   0.6  -0.5
    ]

Known Parameters (Fixed):
- Growth rates: a = [-0.5, 1.0, -0.5, 1.0]

Unknown Parameters (4 total - to be inferred):
- Epsilon perturbations: eps1, eps2, eps3, eps4

Returns:
- model: System
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
    # Diagonal terms are negative (carrying capacity / self-regulation)
    # to prevent finite-time blow-up. Off-diagonal terms encode predator-prey coupling.
    B11, B12, B21, B22 = -0.5, -0.6, 0.6, -0.5
    B33, B34, B43, B44 = -0.5, -0.6, 0.6, -0.5

    @mtkcompile model = System(
        [
            D(x1) ~ x1 * (a1 + B11 * x1 + B12 * x2 + eps1 * x3 + eps2 * x4),
            D(x2) ~ x2 * (a2 + B21 * x1 + B22 * x2 + eps3 * x3 + eps4 * x4),
            D(x3) ~ x3 * (a3 + (-eps1) * x1 + (-eps3) * x2 + B33 * x3 + B34 * x4),
            D(x4) ~ x4 * (a4 + (-eps2) * x1 + (-eps4) * x2 + B43 * x3 + B44 * x4),
        ],
        t,
        states,
        parameters,
    )

    measured_quantities = [y1 ~ x1, y2 ~ x2, y3 ~ x3, y4 ~ x4]
    model, parameters, states, measured_quantities
end

"""
3D Lotka-Volterra model (variant 1) — globally identifiable

System equations:
    dx1/dt = a*x1 + b*x1*x2
    dx2/dt = b*x1*x2 + c*x2

Measured quantities:
    y1 = x1

This model is structurally globally identifiable: the parameters (a, b, c)
are uniquely determined by the observed output y1 = x1.

Returns:
- model: System
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
    @mtkcompile model = System(
        [D(x1) ~ a * x1 + b * x1 * x2, D(x2) ~ b * x1 * x2 + c * x2],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
3D Lotka-Volterra model — locally identifiable (squared parameterization)

System equations:
    dx1/dt = a²*x1 + b²*x1*x2
    dx2/dt = b²*x1*x2 + c*x2

Measured quantities:
    y1 = x1

This model is structurally **locally** identifiable but not globally identifiable.
The parameters a and b appear only as a² and b², so (a, b) and (-a, -b) produce
identical ODE trajectories. The error landscape therefore has multiple equivalent
global minimizers — useful for testing whether optimization methods can find all
of them.

Returns:
- model: System
- parameters: [a, b, c]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lotka_volterra_3D_model_locally_identifiable()
    @independent_variables t
    @variables x1(t) x2(t) y1(t)
    @parameters a b c
    D = Differential(t)
    params = [a, b, c]
    states = [x1, x2]
    @mtkcompile model = System(
        [D(x1) ~ a^2 * x1 + b^2 * x1 * x2, D(x2) ~ b^2 * x1 * x2 + c * x2],
        t,
        states,
        params,
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
- model: System
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
    @mtkcompile model = System(
        [D(x1) ~ a * x1 + -b * x1 * x2, D(x2) ~ -b * x2 + c * x1 * x2],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
Simple 4D Lotka-Volterra model (2 species, 4 parameters) — NOT identifiable

Paper Example 4.3. Classic Lotka-Volterra predator-prey with all 4 rate
parameters unknown and only one species observed.

System equations:
    dx1/dt = a*x1 + b*x1*x2
    dx2/dt = c*x1*x2 + d*x2

Parameters (4): a, b, c, d
States (2): x1, x2
Output (1): y1 = x1

This model is NOT structurally identifiable (neither globally nor locally).
The paper uses it as a negative result: Globtim should find a continuum/manifold
of equivalent parameters rather than isolated critical points.

Returns:
- model: System
- parameters: [a, b, c, d]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lotka_volterra_4D_simple()
    @independent_variables t
    @variables x1(t) x2(t) y1(t)
    @parameters a b c d
    D = Differential(t)
    params = [a, b, c, d]
    states = [x1, x2]
    @mtkcompile model = System(
        [D(x1) ~ a * x1 + b * x1 * x2, D(x2) ~ c * x1 * x2 + d * x2],
        t,
        states,
        params,
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
- model: System
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
    @mtkcompile model =
        System([D(x1) ~ a * x1 + b * x1 * x2, D(x2) ~ b * x1 * x2 + x2], t, states, params)
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
- model: System
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
    @mtkcompile model = System(
        [D(x1) ~ a * x1 + -b * x1 * x2, D(x2) ~ -b * x2 + 0.5 * x1 * x2],
        t,
        states,
        params,
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
- model: System
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
    @mtkcompile model = System(
        [D(x1) ~ a * x1 + -b * x1 * x2, D(x2) ~ -b * x2 + 0.5 * x1 * x2],
        t,
        states,
        params,
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
- model: System
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
    @mtkcompile model = System(
        [D(x1) ~ a * x1 + b * x1 * x2, D(x2) ~ b * x1 * x2 + 0.1 * x2],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
2D Lotka-Volterra model — SciML benchmark parameterization

System equations:
    dx1/dt =  a1*x1 - 0.6*x1*x2
    dx2/dt = -a2*x2 + 0.8*x1*x2

Measured quantities:
    y1 = x1
    y2 = x2

This parameterization uses the diagonal growth/death rates (a1, a2) as the unknown
parameters, with the off-diagonal interaction coefficients (0.6, 0.8) as known
constants. Both states are observed. This is the parameterization used in the
SciML/DiffEqParamEstim benchmark where local optimizers struggle due to the
landscape geometry on large parameter domains.

Returns:
- model: System
- parameters: [a1, a2]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1, y2 ~ x2]
"""
function define_lotka_volterra_2D_sciml_benchmark()
    @independent_variables t
    @variables x1(t) x2(t) y1(t) y2(t)
    @parameters a1 a2
    D = Differential(t)
    params = [a1, a2]
    states = [x1, x2]
    @mtkcompile model = System(
        [D(x1) ~ a1 * x1 - 0.6 * x1 * x2, D(x2) ~ -a2 * x2 + 0.8 * x1 * x2],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1, y2 ~ x2]
    return model, params, states, outputs
end

"""
LV2D reparameterized 3D model — nonlinear map from (θ₁,θ₂,θ₃) → (a,b)

Takes the LV2D v1 model (which has 12 minima in 2D) and embeds it in 3D via a
nonlinear polynomial map. Each LV2D minimum at (a*,b*) becomes a curved 1D manifold
in θ-space. Different minima produce different curves → multiple non-planar valleys.

Reparameterization (polynomial variant):
    a = θ₁ + θ₃²
    b = θ₂ · θ₃

System equations (same ODE as LV2D v1):
    dx₁/dt = (θ₁ + θ₃²)·x₁ + (θ₂·θ₃)·x₁·x₂
    dx₂/dt = (θ₂·θ₃)·x₁·x₂ + x₂

Parameters (3): θ₁, θ₂, θ₃
States (2): x₁, x₂
Outputs (1): y₁ = x₁

Why it works: θ₃ appears nonlinearly (squared in a, multiplicatively in b),
so level sets cannot be hyperplanes. The known LV2D landscape structure is
preserved but warped into genuinely curved 3D manifolds.

Returns:
- model: System
- parameters: [theta1, theta2, theta3]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lv2d_reparam_3d_model()
    @independent_variables t
    @parameters theta1 theta2 theta3
    @variables x1(t) x2(t) y1(t)
    D = Differential(t)
    params = [theta1, theta2, theta3]
    states = [x1, x2]
    # a = theta1 + theta3^2, b = theta2 * theta3
    @mtkcompile model = System(
        [
            D(x1) ~ (theta1 + theta3^2) * x1 + (theta2 * theta3) * x1 * x2,
            D(x2) ~ (theta2 * theta3) * x1 * x2 + x2,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
Coupled LV2D 3D model — two LV2D subsystems with shared parameter + state coupling.

Constructs a 3D parameter landscape by coupling two LV2D v1 subsystems:
- Subsystem A: (a, b) = (θ₁, θ₃)
- Subsystem B: (a, b) = (α·θ₃, θ₂) with diffusive state coupling ε·(x₁-x₃), ε·(x₂-x₄)

The shared parameter θ₃ and diffusive coupling create genuine 3D non-separability.
Near θ = (0.2, 0.4, 0.4), both subsystems are close to LV2D_paper_1's interesting landscape.

Diffusive coupling keeps subsystems synchronized without adding net energy — unlike
additive feed-forward (ε·x₁), which destabilizes the marginally-stable LV dynamics.

Parameters: θ₁ (a of A), θ₂ (b of B), θ₃ (b of A, determines a of B)
States: x₁, x₂ (subsystem A), x₃, x₄ (subsystem B)
Outputs: y₁ = x₁, y₂ = x₃
"""
function define_coupled_lv2d_3d_model(; epsilon::Float64 = 0.1, alpha::Float64 = 0.5)
    @independent_variables t
    @parameters theta1 theta2 theta3
    @variables x1(t) x2(t) x3(t) x4(t) y1(t) y2(t)
    D = Differential(t)
    params = [theta1, theta2, theta3]
    states = [x1, x2, x3, x4]
    @mtkcompile model = System(
        [
            D(x1) ~ theta1 * x1 + theta3 * x1 * x2,
            D(x2) ~ theta3 * x1 * x2 + x2,
            D(x3) ~ alpha * theta3 * x3 + theta2 * x3 * x4 + epsilon * (x1 - x3),
            D(x4) ~ theta2 * x3 * x4 + x4 + epsilon * (x2 - x4),
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1, y2 ~ x3]
    return model, params, states, outputs
end

"""
LV 3D symmetric — locally identifiable via permutation symmetry

Lotka-Volterra dynamics with parameters through elementary symmetric polynomials.
Parameters a and b appear only via (a+b) and (a·b), so (a,b,c) and (b,a,c) produce
identical trajectories. This creates 2 discrete equivalent global minima.

System equations:
    dx1/dt = (a+b) · x1 + (a·b) · x1·x2
    dx2/dt = (a·b) · x1·x2 + c · x2

The (a·b) product creates hyperbolic level sets (ab=const) which are more curved
than standard LV hyperplanes.

Parameters (3): a, b, c
States (2): x1, x2
Outputs (1): y1 = x1

Non-uniqueness: 2 minima — (a, b, c) ↔ (b, a, c) permutation.

Returns:
- model: System
- parameters: [a, b, c]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_lv_3d_symmetric_model()
    @independent_variables t
    @parameters a b c
    @variables x1(t) x2(t) y1(t)
    D = Differential(t)
    params = [a, b, c]
    states = [x1, x2]
    @mtkcompile model = System(
        [D(x1) ~ (a + b) * x1 + (a * b) * x1 * x2, D(x2) ~ (a * b) * x1 * x2 + c * x2],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
    create_lv2d_localid1d_3d_objective(; alpha=0.3, beta=1.0, ...) -> Function

Create a 3D hybrid objective combining two ODE error surfaces via modulated additive coupling:

    F(a₁, a₂, a₃) = E_LV2D(a₁, a₂) · (1 + α·sin(π·a₃)) + β·E_1D(a₃)

Components:
- E_LV2D: Lotka-Volterra 2D SciML benchmark error (2 params → scalar)
- E_1D: Simple 1D locally-identifiable model error (dx/dt = x + a², 1 param → scalar)

The sin(π·a₃) modulation creates non-separable 3D structure:
- At true 2D params (E_LV2D = 0): a₃ CPs coincide with 1D ODE CPs
- At secondary 2D CPs (E_LV2D > 0): a₃ CPs shift proportionally
- At a₃ = ±1.5 (true 1D params): sin(±3π/2) = ∓1, breaking the ±a symmetry
  → basin at +1.5 is shallower (scale 1-α), basin at -1.5 is deeper (scale 1+α)

Two global minima: (1.0, 1.5, +1.5) and (1.0, 1.5, -1.5), both with F = 0.

# Arguments
- `alpha::Float64 = 0.3`: Modulation strength. Must be < 1 to keep (1 + α·sin) > 0.
- `beta::Float64 = 0.068`: Scale factor for 1D component (calibrated: median E_LV2D / median E_1D).
- `lv2d_ic`, `lv2d_p_true`, `lv2d_time`, `lv2d_numpoints`: LV2D SciML configuration.
- `id1d_ic`, `id1d_p_true`, `id1d_time`, `id1d_numpoints`: Locally-id 1D configuration.

Returns: `f(x::Vector{Float64}) -> Float64` where `x = [a₁, a₂, a₃]`.
"""
function create_lv2d_localid1d_3d_objective(;
    alpha::Float64 = 0.3,
    beta::Float64 = 0.068,
    lv2d_ic::Vector{Float64} = [1.0, 1.0],
    lv2d_p_true::Vector{Float64} = [1.0, 1.5],
    lv2d_time::Vector{Float64} = [0.0, 1.5],
    lv2d_numpoints::Int = 4,
    id1d_ic::Vector{Float64} = [1.0],
    id1d_p_true::Vector{Float64} = [1.5],
    id1d_time::Vector{Float64} = [0.0, 5.0],
    id1d_numpoints::Int = 30,
)
    alpha < 1.0 || error("alpha must be < 1 to ensure (1 + α·sin) > 0 everywhere")

    # Build LV2D SciML objective (2D → scalar)
    lv2d_model, _, _, lv2d_outputs = define_lotka_volterra_2D_sciml_benchmark()
    E_lv2d = make_error_distance(
        lv2d_model,
        lv2d_outputs,
        lv2d_ic,
        lv2d_p_true,
        lv2d_time,
        lv2d_numpoints,
        L2_norm,
        sum;
        return_inf_on_error = true,
    )

    # Build 1D locally-id objective (1D → scalar)
    id1d_model, _, _, id1d_outputs = define_simple_1D_model_locally_identifiable()
    E_1d = make_error_distance(
        id1d_model,
        id1d_outputs,
        id1d_ic,
        id1d_p_true,
        id1d_time,
        id1d_numpoints,
        L2_norm,
        sum;
        return_inf_on_error = true,
    )

    function hybrid(x::Vector{Float64})
        length(x) == 3 || error("Expected 3D input [a₁, a₂, a₃], got $(length(x))D")
        e_2d = E_lv2d(x[1:2])
        e_1d = E_1d([x[3]])
        modulation = 1.0 + alpha * sin(π * x[3])
        return e_2d * modulation + beta * e_1d
    end

    return hybrid
end
