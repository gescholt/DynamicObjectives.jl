# Other Dynamical Systems
# FitzHugh-Nagumo and simple test models for identifiability studies

"""
FitzHugh-Nagumo 3D model

A simplified model of neuronal excitability derived from the Hodgkin-Huxley model.

System equations:
    dV/dt = g*(V - V³/3 + R)
    dR/dt = (1/g)*(V - a + b*R)

Parameters:
- g: coupling parameter
- a: threshold parameter
- b: recovery parameter

Measured quantities:
    y1 = V (voltage)

Returns:
- model: System
- parameters: [g, a, b]
- states: [V, R]
- measured_quantities: [y1 ~ V]
"""
function define_fitzhugh_nagumo_3D_model()
    @independent_variables t
    @parameters g a b
    @variables V(t) R(t) y1(t)
    D = Differential(t)
    states = [V, R]
    params = [g, a, b]
    outputs = [y1 ~ V]
    @mtkcompile model = System(
        [D(V) ~ g * (V - V^3 / 3 + R), D(R) ~ 1 / g * (V - a + b * R)],
        t,
        states,
        params,
    )
    return model, params, states, outputs
end

"""
FitzHugh-Nagumo 3D model — dual output (V + R)

Same ODE dynamics as `define_fitzhugh_nagumo_3D_model`, but both states are observed.
Measuring both voltage V and recovery variable R eliminates the structural
non-identifiability that makes single-output FHN level sets collapse to hyperplanes.

System equations:
    dV/dt = g*(V - V³/3 + R)
    dR/dt = (1/g)*(V - a + b*R)

Parameters:
- g: coupling parameter
- a: threshold parameter
- b: recovery parameter

Measured quantities:
    y1 = V (voltage)
    y2 = R (recovery)

Returns:
- model: System
- parameters: [g, a, b]
- states: [V, R]
- measured_quantities: [y1 ~ V, y2 ~ R]
"""
function define_fitzhugh_nagumo_3D_model_two_outputs()
    @independent_variables t
    @parameters g a b
    @variables V(t) R(t) y1(t) y2(t)
    D = Differential(t)
    states = [V, R]
    params = [g, a, b]
    @mtkcompile model = System(
        [D(V) ~ g * (V - V^3 / 3 + R), D(R) ~ 1 / g * (V - a + b * R)],
        t,
        states,
        params,
    )
    outputs = [y1 ~ V, y2 ~ R]
    return model, params, states, outputs
end

"""
Goodwin oscillator 4D model — Hill function negative feedback loop

Paper Example 4.5. A 3-state biochemical oscillator with Hill function
repression. The model describes mRNA (x1) → protein (x2) → end product (x3)
with x3 repressing x1 production via a Hill function.

System equations:
    dx1/dt = k1 * K^n / (K^n + x3^n) - k2 * x1
    dx2/dt = k3 * x1 - k4 * x2
    dx3/dt = k5 * x2 - k6 * x3

Fixed (known) constants:
    K  = 0.9  (Hill half-maximal constant)
    n  = 10   (Hill exponent)
    k3 = 0.3  (translation rate)
    k6 = 0.5  (end product degradation)

Unknown parameters (4): k1, k2, k4, k5
States (3): x1, x2, x3
Outputs (2): y1 = x1, y2 = x3

Note: The Hill function K^n / (K^n + x3^n) is a rational (non-polynomial)
nonlinearity. This makes the model more challenging for polynomial
approximation methods.

Returns:
- model: System
- parameters: [k1, k2, k4, k5]
- states: [x1, x2, x3]
- measured_quantities: [y1 ~ x1, y2 ~ x3]
"""
function define_goodwin_oscillator_4D()
    @independent_variables t
    @parameters k1 k2 k4 k5
    @variables x1(t) x2(t) x3(t) y1(t) y2(t)
    D = Differential(t)

    params = [k1, k2, k4, k5]
    states = [x1, x2, x3]

    # Fixed constants
    K = 0.9    # Hill half-maximal constant
    n = 10     # Hill exponent
    k3 = 0.3   # translation rate (fixed, known)
    k6 = 0.5   # end product degradation (fixed, known)
    Kn = K^n   # precompute K^10 = 0.9^10

    @mtkcompile model = System(
        [
            D(x1) ~ k1 * Kn / (Kn + x3^n) - k2 * x1,
            D(x2) ~ k3 * x1 - k4 * x2,
            D(x3) ~ k5 * x2 - k6 * x3,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1, y2 ~ x3]
    return model, params, states, outputs
end

"""
Simple 2D model - locally identifiable (product form)

System equations:
    dx1/dt = a*b*x1 + (a + b)

This model is structurally locally identifiable but not globally identifiable
because the parameters a and b appear only in the combinations a*b and a+b.

Measured quantities:
    y1 = x1

Returns:
- model: System
- parameters: [a, b]
- states: [x1]
- measured_quantities: [y1 ~ x1]
"""
function define_simple_2D_model_locally_identifiable()
    @independent_variables t
    @variables x1(t) y1(t)
    @parameters a b
    D = Differential(t)
    params = [a, b]
    states = [x1]
    @mtkcompile model = System([D(x1) ~ a * b * x1 + (a + b)], t, states, params)
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
Simple 2D model - locally identifiable (square form)

System equations:
    dx1/dt = a*x1 + b²

This model is structurally locally identifiable but not globally identifiable
because the parameter b appears only as b².

Measured quantities:
    y1 = x1

Returns:
- model: System
- parameters: [a, b]
- states: [x1]
- measured_quantities: [y1 ~ x1]
"""
function define_simple_2D_model_locally_identifiable_square()
    @independent_variables t
    @variables x1(t) y1(t)
    @parameters a b
    D = Differential(t)
    params = [a, b]
    states = [x1]
    @mtkcompile model = System([D(x1) ~ a * x1 + b^2], t, states, params)
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

"""
Simple 1D model - locally identifiable (square form)

System equations:
    dx1/dt = x1 + a²

This model is structurally locally identifiable but not globally identifiable
because the parameter a appears only as a².

Measured quantities:
    y1 = x1

Returns:
- model: System
- parameters: [a]
- states: [x1]
- measured_quantities: [y1 ~ x1]
"""
function define_simple_1D_model_locally_identifiable()
    @independent_variables t
    @variables x1(t) y1(t)
    @parameters a
    D = Differential(t)
    params = [a]
    states = [x1]
    @mtkcompile model = System([D(x1) ~ x1 + a^2], t, states, params)
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

# ═══════════════════════════════════════════════════════════════════════════════
# Chaotic & nonlinear 3D models for level set landscape exploration
# ═══════════════════════════════════════════════════════════════════════════════

"""
Lorenz 3D model (σ, ρ, β)

The classic Lorenz system of deterministic chaos. Cross-terms x*(ρ-z) and x*y
create strong nonlinear parameter coupling, producing genuinely curved level sets.

System equations:
    dx/dt = σ(y - x)
    dy/dt = x(ρ - z) - y
    dz/dt = xy - βz

Parameters:
- σ (sigma): Prandtl number / coupling
- ρ (rho): Rayleigh number (chaos transition at ρ ≈ 24.74)
- β (beta): geometric factor

Measured quantities:
    y1 = x, y2 = z  (2 of 3 states — avoids over-determination)

Returns:
- model: System
- parameters: [σ, ρ, β]
- states: [x, y, z]
- measured_quantities: [y1 ~ x, y2 ~ z]
"""
function define_lorenz_3D_model()
    @independent_variables t
    @parameters σ ρ β
    @variables x(t) y(t) z(t) y1(t) y2(t)
    D = Differential(t)
    states = [x, y, z]
    params = [σ, ρ, β]
    @mtkcompile model = System(
        [D(x) ~ σ * (y - x), D(y) ~ x * (ρ - z) - y, D(z) ~ x * y - β * z],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x, y2 ~ z]
    return model, params, states, outputs
end

"""
Rössler 3D model (a, b, c)

A prototypical system exhibiting period-doubling cascades and spiral chaos.
The key nonlinearity z*(x - c) creates strong parameter coupling in c.

System equations:
    dx/dt = -y - z
    dy/dt = x + a*y
    dz/dt = b + z*(x - c)

Parameters:
- a: controls spiral stability
- b: offset in z dynamics
- c: controls period-doubling cascade (key nonlinearity via z*(x-c))

Measured quantities:
    y1 = x, y2 = z

Returns:
- model: System
- parameters: [a, b, c]
- states: [x, y, z]
- measured_quantities: [y1 ~ x, y2 ~ z]
"""
function define_rossler_3D_model()
    @independent_variables t
    @parameters a b c
    @variables x(t) y(t) z(t) y1(t) y2(t)
    D = Differential(t)
    states = [x, y, z]
    params = [a, b, c]
    @mtkcompile model =
        System([D(x) ~ -y - z, D(y) ~ x + a * y, D(z) ~ b + z * (x - c)], t, states, params)
    outputs = [y1 ~ x, y2 ~ z]
    return model, params, states, outputs
end

"""
Goodwin oscillator 3D model (k1, k2, k4)

A reduced 3-parameter version of the 4D Goodwin oscillator. Fixes k5=1.5
(known from the 4D model) to create a 3D parameter estimation problem.
The Hill function nonlinearity K^n/(K^n+x3^n) should produce curved level sets.

System equations:
    dx1/dt = k1 * K^n / (K^n + x3^n) - k2 * x1
    dx2/dt = 0.3 * x1 - k4 * x2
    dx3/dt = 1.5 * x2 - 0.5 * x3

Fixed constants:
    K  = 0.9  (Hill half-maximal constant)
    n  = 10   (Hill exponent)
    k3 = 0.3  (translation rate)
    k5 = 1.5  (end product production — fixed from 4D)
    k6 = 0.5  (end product degradation)

Parameters (3): k1, k2, k4
States (3): x1, x2, x3
Outputs (2): y1 = x1, y2 = x3

Returns:
- model: System
- parameters: [k1, k2, k4]
- states: [x1, x2, x3]
- measured_quantities: [y1 ~ x1, y2 ~ x3]
"""
function define_goodwin_oscillator_3D()
    @independent_variables t
    @parameters k1 k2 k4
    @variables x1(t) x2(t) x3(t) y1(t) y2(t)
    D = Differential(t)

    params = [k1, k2, k4]
    states = [x1, x2, x3]

    # Fixed constants (same as 4D model)
    K = 0.9    # Hill half-maximal constant
    n = 10     # Hill exponent
    k3 = 0.3   # translation rate (fixed, known)
    k5 = 1.5   # end product production (fixed from 4D)
    k6 = 0.5   # end product degradation (fixed, known)
    Kn = K^n   # precompute K^10 = 0.9^10

    @mtkcompile model = System(
        [
            D(x1) ~ k1 * Kn / (Kn + x3^n) - k2 * x1,
            D(x2) ~ k3 * x1 - k4 * x2,
            D(x3) ~ k5 * x2 - k6 * x3,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1, y2 ~ x3]
    return model, params, states, outputs
end

"""
Goodwin oscillator 3D model with product observation — y₁ = x₁·x₃

Same ODE as `define_goodwin_oscillator_3D` but observes the product x₁·x₃ instead
of x₁ and x₃ separately. The product observation introduces parameter interactions
through the measurement equation, potentially creating more curved level sets.

System equations (same as Goodwin 3D):
    dx1/dt = k1 * K^n / (K^n + x3^n) - k2 * x1
    dx2/dt = 0.3 * x1 - k4 * x2
    dx3/dt = 1.5 * x2 - 0.5 * x3

Changed output: y₁ = x₁ * x₃  (product observation)

Parameters (3): k1, k2, k4
States (3): x1, x2, x3
Outputs (1): y1 = x1 * x3

Returns:
- model: System
- parameters: [k1, k2, k4]
- states: [x1, x2, x3]
- measured_quantities: [y1 ~ x1 * x3]
"""
function define_goodwin_3d_product_obs_model()
    @independent_variables t
    @parameters k1 k2 k4
    @variables x1(t) x2(t) x3(t) y1(t)
    D = Differential(t)

    params = [k1, k2, k4]
    states = [x1, x2, x3]

    # Fixed constants (same as Goodwin 3D)
    K = 0.9
    n = 10
    k3 = 0.3
    k5 = 1.5
    k6 = 0.5
    Kn = K^n

    @mtkcompile model = System(
        [
            D(x1) ~ k1 * Kn / (Kn + x3^n) - k2 * x1,
            D(x2) ~ k3 * x1 - k4 * x2,
            D(x3) ~ k5 * x2 - k6 * x3,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1 * x3]
    return model, params, states, outputs
end

"""
Rosenzweig-MacArthur 3D model (r, a, K) — predator-prey with Holling type II

A predator-prey model where parameters appear nonlinearly in the Holling type II
functional response a·x/(1 + a·h·x). The denominator makes `a` appear nonlinearly
in the dynamics, creating curved level sets naturally.

System equations:
    dx/dt = r·x·(1 - x/K) - a·x·y/(1 + a·h·x)
    dy/dt = e·a·x·y/(1 + a·h·x) - d·y

Free parameters (3): r (growth rate), a (attack rate), K (carrying capacity)
Fixed constants: h=0.5 (handling time), e=0.6 (conversion efficiency), d=0.3 (death rate)

Measured quantities:
    y1 = x (prey density)

Why it works: `a` appears as both a·x·y in numerator and a·h·x in denominator →
saturating nonlinearity. `K` appears in x/K creating parameter-dependent equilibria.
r·x·(1-x/K) bounds prey dynamics → no blow-up. Classic ecology model with
well-understood bifurcation structure.

Returns:
- model: System
- parameters: [r, a, K]
- states: [x, y]
- measured_quantities: [y1 ~ x]
"""
function define_rosenzweig_macarthur_3d_model()
    @independent_variables t
    @parameters r aa KK
    @variables x(t) yy(t) y1(t)
    D = Differential(t)
    params = [r, aa, KK]
    states = [x, yy]

    # Fixed constants
    h = 0.5     # handling time
    e = 0.6     # conversion efficiency
    d = 0.3     # predator death rate

    @mtkcompile model = System(
        [
            D(x) ~ r * x * (1 - x / KK) - aa * x * yy / (1 + aa * h * x),
            D(yy) ~ e * aa * x * yy / (1 + aa * h * x) - d * yy,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x]
    return model, params, states, outputs
end

# ═══════════════════════════════════════════════════════════════════════════════
# Locally identifiable 3D models — discrete finite symmetry minima
# ═══════════════════════════════════════════════════════════════════════════════

"""
Goodwin oscillator 3D — locally identifiable (squared parameters)

Same Hill function ODE as `define_goodwin_oscillator_3D`, but k1 → k1² and k4 → k4².
The squared parameterization creates sign-flip symmetry: (k1, k2, k4) and (-k1, k2, -k4)
produce identical trajectories. This gives 4 discrete equivalent global minima.

System equations:
    dx1/dt = k1² · K^n / (K^n + x3^n) - k2 · x1
    dx2/dt = 0.3 · x1 - k4² · x2
    dx3/dt = 1.5 · x2 - 0.5 · x3

Fixed constants:
    K  = 0.9  (Hill half-maximal constant)
    n  = 10   (Hill exponent)
    k3 = 0.3  (translation rate)
    k5 = 1.5  (end product production)
    k6 = 0.5  (end product degradation)

Parameters (3): k1, k2, k4
States (3): x1, x2, x3
Outputs (2): y1 = x1, y2 = x3

Non-uniqueness: 4 minima — (±k1, k2, ±k4) sign combinations.

Returns:
- model: System
- parameters: [k1, k2, k4]
- states: [x1, x2, x3]
- measured_quantities: [y1 ~ x1, y2 ~ x3]
"""
function define_goodwin_3d_locally_id_model()
    @independent_variables t
    @parameters k1 k2 k4
    @variables x1(t) x2(t) x3(t) y1(t) y2(t)
    D = Differential(t)

    params = [k1, k2, k4]
    states = [x1, x2, x3]

    # Fixed constants (same as Goodwin 3D)
    K = 0.9
    n = 10
    k3 = 0.3
    k5 = 1.5
    k6 = 0.5
    Kn = K^n

    @mtkcompile model = System(
        [
            D(x1) ~ k1^2 * Kn / (Kn + x3^n) - k2 * x1,
            D(x2) ~ k3 * x1 - k4^2 * x2,
            D(x3) ~ k5 * x2 - k6 * x3,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1, y2 ~ x3]
    return model, params, states, outputs
end

"""
FitzHugh-Nagumo 3D — locally identifiable (squared coupling, dual output)

Same FHN ODE but g → g², with both V and R observed. The squared parameterization
creates sign-flip symmetry: (g, a, b) and (-g, a, b) produce identical trajectories.
Dual output avoids the 60-90% plateau problem of single-output FHN.

System equations:
    dV/dt = g² · (V - V³/3 + R)
    dR/dt = (1/g²) · (V - a + b·R)

Parameters:
- g: coupling parameter (appears as g² and 1/g²)
- a: threshold parameter
- b: recovery parameter

Measured quantities:
    y1 = V (voltage)
    y2 = R (recovery)

Non-uniqueness: 2 minima — (±g, a, b).

Returns:
- model: System
- parameters: [g, a, b]
- states: [V, R]
- measured_quantities: [y1 ~ V, y2 ~ R]
"""
function define_fhn_3d_locally_id_model()
    @independent_variables t
    @parameters g a b
    @variables V(t) R(t) y1(t) y2(t)
    D = Differential(t)
    states = [V, R]
    params = [g, a, b]
    @mtkcompile model = System(
        [D(V) ~ g^2 * (V - V^3 / 3 + R), D(R) ~ (1 / g^2) * (V - a + b * R)],
        t,
        states,
        params,
    )
    outputs = [y1 ~ V, y2 ~ R]
    return model, params, states, outputs
end

"""
Driven FitzHugh-Nagumo 2D family (bead 6pj0) — the control input as a
landscape-shaping lever.

System equations (b = 0.8, I = 0.5 fixed; drive amplitude A and period Td
baked into each registered variant; estimated parameters eps, a):
    dv/dt = v - v³/3 - w + I + A*sin(2π t / Td)
    dw/dt = eps*(v + a - b*w)

Measured quantities: y1 = v (the latent recovery variable w is unobserved).

The autonomous variant (A = 0) has a fringe-multimodal PE landscape over
(eps, a): spike-count aliasing gives ~T-proportional local minima (probe run
20260829T122030_fhn_driven_probe: 53/154/256 raw minima at T = 100/200/400,
fringe spacing halving as the window doubles). A known drive reshapes it:
deep decoys are non-monotone in amplitude (7 -> 1 -> 4 at A = 0/0.15/0.3,
T = 200) — moderate drive nearly convexifies via entrainment, strong drive
re-fragments into tongue lobes with non-smooth boundaries.
"""
function _define_fhn_driven(A::Float64, Td::Float64)
    @independent_variables t
    @parameters epsilon a
    @variables v(t) w(t) y1(t)
    D = Differential(t)
    states = [v, w]
    params = [epsilon, a]
    outputs = [y1 ~ v]
    rhs_v = A == 0.0 ? (0.5) : (0.5 + A * sin(2pi / Td * t))
    @mtkcompile model = System(
        [D(v) ~ v - v^3 / 3 - w + rhs_v,
         D(w) ~ epsilon * (v + a - 0.8 * w)],
        t,
        states,
        params,
    )
    return model, params, states, outputs
end

define_fhn_driven_auto() = _define_fhn_driven(0.0, 25.0)
define_fhn_driven_A015() = _define_fhn_driven(0.15, 25.0)
define_fhn_driven_A030() = _define_fhn_driven(0.3, 25.0)
