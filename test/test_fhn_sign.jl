# The "FHN" models of this package have the recovery sign reversed relative to
# FitzHugh-Nagumo (Ramsay et al. 2007); define_fitzhugh_nagumo_ramsay_model is the
# correctly signed one. These tests pin both conventions at the catalogue entry
# FHN3D_1, (g, a, b) = (1.109, 1.133, 0.546), (V, R)(0) = (0, 0), t in [0, 20]:
# reversed sign -> monotone runaway, Ramsay sign -> bounded, V settles near 1.4.

using Test
using DynamicObjectives
using DynamicObjectives: Tsit5, ODEProblem
using DynamicObjectives.SciMLBase: solve

function _fhn_trajectory(define_model, p; T = 20.0, n = 201)
    model, params, states, outputs = define_model()
    # bind by the returned state symbols, not by unknowns(model) order
    prob = ODEProblem(model, merge(Dict(states .=> [0.0, 0.0]), Dict(params .=> p)), (0.0, T))
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10, saveat = range(0.0, T; length = n))
    return sol[states[1]], sol[states[2]]
end

@testset "FHN sign conventions" begin
    p = [1.109, 1.133, 0.546]

    V, R = _fhn_trajectory(define_fitzhugh_nagumo_3D_model, p)
    @test all(diff(V) .<= 1e-10)          # monotone decrease, no spiking
    @test isapprox(V[end], -72.32; atol = 0.05)
    @test R[end] < -1.0e5

    V, R = _fhn_trajectory(define_fitzhugh_nagumo_ramsay_model, p)
    @test isapprox(maximum(V), 1.76; atol = 0.01)
    @test isapprox(V[end], 1.399; atol = 0.005)
    @test !all(diff(V) .>= -1e-10)         # rises, overshoots, comes back
end
