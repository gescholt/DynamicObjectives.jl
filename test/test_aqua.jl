using Test
using Aqua
using DynamicObjectives

include("aqua_config.jl")

@testset "Aqua.jl Quality Assurance" begin
    run_aqua_tests(DynamicObjectives)
end
