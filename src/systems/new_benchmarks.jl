# New Benchmark Models for DynamicObjectives.jl
# SIR, Hindmarsh-Rose, 2-compartment PK, Brusselator, Michaelis-Menten

# ═══════════════════════════════════════════════════════════════════════════════
# Epidemiology — SIR and SEIR compartmental models
# ═══════════════════════════════════════════════════════════════════════════════

"""
SIR model (2 parameters) — Susceptible-Infectious-Recovered

Standard compartmental epidemic model. The total population N = S + I + R is conserved.

System equations:
    dS/dt = -β·S·I / N
    dI/dt =  β·S·I / N - γ·I
    dR/dt =  γ·I

Parameters (2):
- β (beta): transmission rate per contact
- γ (gamma): recovery rate (1/γ = mean infectious period)

Fixed constant:
- N = 1000 (total population, known)

States (3): S, I, R
Output (1): y1 = I (observed incidence / infectious count)

Identifiability: Globally identifiable from I(t) when N is known.
True parameters for verification: β=0.5, γ=0.1 → R₀ = 5.

Returns:
- model: System
- parameters: [β, γ]
- states: [S, I, R]
- measured_quantities: [y1 ~ I]
"""
function define_sir_2d_model()
    @independent_variables t
    @parameters beta gamma
    @variables S(t) I(t) R(t) y1(t)
    D = Differential(t)

    params = [beta, gamma]
    states = [S, I, R]

    N = 1000.0  # total population (known constant)

    @mtkcompile model = System(
        [D(S) ~ -beta * S * I / N, D(I) ~ beta * S * I / N - gamma * I, D(R) ~ gamma * I],
        t,
        states,
        params,
    )
    outputs = [y1 ~ I]
    return model, params, states, outputs
end

"""
SEIR model (3 parameters) — Susceptible-Exposed-Infectious-Recovered

Adds an Exposed (latent) compartment to SIR. The latent period is 1/σ.

System equations:
    dS/dt = -β·S·I / N
    dE/dt =  β·S·I / N - σ·E
    dI/dt =  σ·E - γ·I
    dR/dt =  γ·I

Parameters (3):
- β (beta): transmission rate
- σ (sigma): progression rate from exposed to infectious (1/σ = mean latent period)
- γ (gamma): recovery rate

Fixed constant:
- N = 1000 (total population, known)

States (4): S, E, I, R
Output (1): y1 = I (observed infectious count)

Identifiability: σ and γ can trade off with only I observed; weak identifiability
makes this an interesting test case for globtim's valley-following behavior.
True parameters for verification: β=0.6, σ=0.2, γ=0.1.

Returns:
- model: System
- parameters: [β, σ, γ]
- states: [S, E, I, R]
- measured_quantities: [y1 ~ I]
"""
function define_seir_3d_model()
    @independent_variables t
    @parameters beta sigma gamma
    @variables S(t) E(t) I(t) R(t) y1(t)
    D = Differential(t)

    params = [beta, sigma, gamma]
    states = [S, E, I, R]

    N = 1000.0

    @mtkcompile model = System(
        [
            D(S) ~ -beta * S * I / N,
            D(E) ~ beta * S * I / N - sigma * E,
            D(I) ~ sigma * E - gamma * I,
            D(R) ~ gamma * I,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ I]
    return model, params, states, outputs
end

# ═══════════════════════════════════════════════════════════════════════════════
# Neuroscience — Hindmarsh-Rose (reduced 3-parameter bursting model)
# ═══════════════════════════════════════════════════════════════════════════════

"""
Hindmarsh-Rose 3D model (3 parameters) — neuronal bursting

A simplified 3-parameter variant of the classic Hindmarsh-Rose model for neuronal
bursting. The full model has 10 parameters spanning 4 orders of magnitude; this
reduced version fixes the recovery and adaptation terms to study the core
bursting parameters.

System equations:
    dx/dt = y + a·x² - x³ - z
    dy/dt = 1 - b·x² - y
    dz/dt = c·(x - 1.6) - z

Parameters (3):
- a: fast subsystem gain (controls spiking threshold)
- b: quadratic feedback (controls burst frequency)
- c: slow adaptation rate (controls inter-burst interval)

Fixed constants:
- x_rest = 1.6 (resting potential offset for adaptation current)

States (3): x (membrane potential), y (fast recovery), z (slow adaptation)
Output (1): y1 = x (observed voltage)

True parameters for verification: a=3.0, b=5.0, c=1.0.
At these values the model exhibits regular bursting with ~5 spikes per burst.

Returns:
- model: System
- parameters: [a, b, c]
- states: [x, y, z]
- measured_quantities: [y1 ~ x]
"""
function define_hindmarsh_rose_3d_model()
    @independent_variables t
    @parameters a b c
    @variables x(t) y(t) z(t) y1(t)
    D = Differential(t)

    params = [a, b, c]
    states = [x, y, z]

    x_rest = 1.6  # resting offset

    @mtkcompile model = System(
        [
            D(x) ~ y + a * x^2 - x^3 - z,
            D(y) ~ 1.0 - b * x^2 - y,
            D(z) ~ c * (x - x_rest) - z,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x]
    return model, params, states, outputs
end

# ═══════════════════════════════════════════════════════════════════════════════
# Pharmacokinetics — 2-compartment model
# ═══════════════════════════════════════════════════════════════════════════════

"""
2-compartment pharmacokinetic model (3 parameters) — plasma + tissue

Standard two-compartment PK model describing drug absorption/distribution.
Compartment 1 = central (plasma, observed), Compartment 2 = peripheral (tissue).

System equations:
    dx1/dt = -(k₁₀ + k₁₂)·x1 + k₂₁·x2
    dx2/dt =  k₁₂·x1 - k₂₁·x2

Parameters (3):
- k10: elimination rate from central compartment
- k12: distribution rate (central → peripheral)
- k21: distribution rate (peripheral → central)

States (2): x1 (plasma concentration), x2 (tissue concentration)
Output (1): y1 = x1 (plasma is observed)

Identifiability: Globally identifiable from plasma concentration data alone
for a bolus IV administration (initial condition x1(0) > 0, x2(0) = 0).
True parameters: k10=0.1, k12=0.5, k21=0.3.
The analytical solution is a bi-exponential: x1(t) = A·exp(-αt) + B·exp(-βt),
which makes this an excellent baseline for exact-error verification.

Returns:
- model: System
- parameters: [k10, k12, k21]
- states: [x1, x2]
- measured_quantities: [y1 ~ x1]
"""
function define_pk_2comp_3d_model()
    @independent_variables t
    @parameters k10 k12 k21
    @variables x1(t) x2(t) y1(t)
    D = Differential(t)

    params = [k10, k12, k21]
    states = [x1, x2]

    @mtkcompile model = System(
        [D(x1) ~ -(k10 + k12) * x1 + k21 * x2, D(x2) ~ k12 * x1 - k21 * x2],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end

# ═══════════════════════════════════════════════════════════════════════════════
# Chemistry — Brusselator (autocatalytic oscillator)
# ═══════════════════════════════════════════════════════════════════════════════

"""
Brusselator model (2 parameters) — autocatalytic chemical oscillator

The Brusselator is a classic model of an autocatalytic reaction system that
exhibits sustained oscillations when B > 1 + A².

System equations:
    dx/dt = A + x²·y - (B + 1)·x
    dy/dt = B·x - x²·y

Parameters (2):
- A: feed rate of reactant X
- B: feed rate / autocatalytic gain

States (2): x, y (concentrations of intermediates X and Y)
Output (1): y1 = x

Oscillation condition: B > 1 + A². For verification:
- A = 1.0, B = 3.0 → B = 3 > 1 + 1 = 2, so sustained oscillations occur.
True parameters: A=1.0, B=3.0.

Identifiability: Globally identifiable from x(t) observations.

Returns:
- model: System
- parameters: [A, B]
- states: [x, y]
- measured_quantities: [y1 ~ x]
"""
function define_brusselator_2d_model()
    @independent_variables t
    @parameters A B
    @variables x(t) y(t) y1(t)
    D = Differential(t)

    params = [A, B]
    states = [x, y]

    @mtkcompile model = System(
        [D(x) ~ A + x^2 * y - (B + 1.0) * x, D(y) ~ B * x - x^2 * y],
        t,
        states,
        params,
    )
    outputs = [y1 ~ x]
    return model, params, states, outputs
end

# ═══════════════════════════════════════════════════════════════════════════════
# Biochemistry — Michaelis-Menten enzyme kinetics
# ═══════════════════════════════════════════════════════════════════════════════

"""
Michaelis-Menten model (2 parameters) — enzyme-substrate dynamics

Classic enzyme kinetics model with substrate (S) and enzyme-substrate complex (ES).
The reaction follows E + S ⇌ ES → E + P, with product P not tracked.

System equations:
    dS/dt  = -k1·E_total·S + k2·ES + kcat·ES
    dES/dt =  k1·E_total·S - k2·ES - kcat·ES

Quasi-steady-state reduction: the observable dynamics depend on:
    Vmax = kcat · E_total    (maximum reaction rate)
    Km   = (k2 + kcat) / k1  (Michaelis constant)

For parameter estimation we estimate Vmax and Km directly, which are the
standard PKPD / biochemistry parameters.

Parameters (2):
- Vmax: maximum reaction rate
- Km: substrate concentration at half-maximal rate

Fixed constant:
- E_total = 1.0 (total enzyme concentration, known)

States (2): S (substrate), ES (enzyme-substrate complex)
Output (1): y1 = S (observed substrate depletion)

True parameters: Vmax=2.0, Km=1.0.
The familiar rate law v = Vmax·S/(Km + S) governs the net flux.

Identifiability: Globally identifiable from S(t) data.

Returns:
- model: System
- parameters: [Vmax, Km]
- states: [S, ES]
- measured_quantities: [y1 ~ S]
"""
function define_michaelis_menten_2d_model()
    @independent_variables t
    @parameters Vmax Km
    @variables S(t) ES(t) y1(t)
    D = Differential(t)

    params = [Vmax, Km]
    states = [S, ES]

    E_total = 1.0

    # Reconstruct individual rate constants from Vmax and Km
    # Vmax = kcat * E_total  → kcat = Vmax / E_total
    # Km   = (k2 + kcat) / k1
    # We fix k1 = 1.0, then k2 = Km - kcat
    k1 = 1.0
    kcat = Vmax / E_total
    k2 = Km * k1 - kcat

    @mtkcompile model = System(
        [
            D(S) ~ -k1 * E_total * S + k2 * ES + kcat * ES,
            D(ES) ~ k1 * E_total * S - k2 * ES - kcat * ES,
        ],
        t,
        states,
        params,
    )
    outputs = [y1 ~ S]
    return model, params, states, outputs
end
