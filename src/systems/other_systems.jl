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
- model: ODESystem
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
    @mtkbuild model = ODESystem(
        [D(V) ~ g * (V - V^3 / 3 + R), D(R) ~ 1 / g * (V - a + b * R)],
        t,
        states,
        params
    )
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
- model: ODESystem
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
    @mtkbuild model = ODESystem([D(x1) ~ a * b * x1 + (a + b)], t, states, params)
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
- model: ODESystem
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
    @mtkbuild model = ODESystem([D(x1) ~ a * x1 + b^2], t, states, params)
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
- model: ODESystem
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
    @mtkbuild model = ODESystem([D(x1) ~ x1 + a^2], t, states, params)
    outputs = [y1 ~ x1]
    return model, params, states, outputs
end
