# Smoke test for DynamicObjectivesGlobtimExt (run_experiment_from_config).
#
# Requires Globtim AND GlobtimPostProcessing on the load path — i.e. a
# workspace profile (profiles/dev, profiles/cluster). The public standalone
# mirror cannot resolve GlobtimPostProcessing (unregistered), so runtests.jl
# only includes this file when the package is resolvable.
#
# Standalone run:
#   julia --project=profiles/dev pkg/DynamicObjectives/test/test_ext_toml_pipeline.jl

using Test
using DynamicObjectives
using Globtim
using GlobtimPostProcessing

@testset "DynamicObjectivesGlobtimExt: TOML pipeline smoke" begin
    ext = Base.get_extension(DynamicObjectives, :DynamicObjectivesGlobtimExt)
    @test ext !== nothing

    mktempdir() do dir
        outdir = joinpath(dir, "out")
        config_path = joinpath(dir, "smoke_2d.toml")
        write(
            config_path,
            """
            [experiment]
            name = "ext_smoke_deuflhard_2d"
            description = "tiny ext smoke — degree 4, GN 6"

            [model]
            analytical_function = "Deuflhard"
            dimension = 2

            [domain]
            bounds = [[-1.2, 1.2], [-1.2, 1.2]]

            [polynomial]
            GN = 6
            degree_range = [4, 2, 4]
            basis = "chebyshev"

            [analysis]
            enabled = true
            gradient_method = "forwarddiff"
            newton_tol = 1e-8
            newton_max_iterations = 50
            hessian_tol = 1e-6
            dedup_fraction = 0.02

            [output]
            dir = "$(outdir)"
            """,
        )

        result = run_experiment_from_config(config_path; io = devnull)

        @test result isa Dict{Symbol,Any}
        for key in (:config, :experiment_result, :degree_results, :output_dir)
            @test haskey(result, key)
        end
        @test !isempty(result[:degree_results])

        # Reproducibility artifacts land in the output dir
        @test isfile(joinpath(result[:output_dir], "experiment_config.toml"))
        @test isfile(joinpath(result[:output_dir], "results_summary.json"))
    end
end
