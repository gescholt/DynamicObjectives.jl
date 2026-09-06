using DynamicObjectives
using DynamicObjectives: DataStructures, ModelingToolkit
import DynamicObjectives: ODEProblem, Tsit5, Vern9, AutoTsit5, Rosenbrock23
using Test

@testset "DynamicObjectives.jl" begin

    # Core functionality tests
    @testset "DAISY LV4D Definition" begin
        model, params, states, outputs = define_daisy_ex3_model_4D()
        @test length(params) == 4
        @test length(states) == 4
        @test length(outputs) == 2
    end

    @testset "Generalized LV4D Definition" begin
        model, params, states, outputs = define_generalized_lotka_volterra_4D()
        @test length(params) == 20  # 4 growth rates + 16 interaction matrix elements
        @test length(states) == 4
        @test length(outputs) == 4
    end

    @testset "Constrained LV4D Definition" begin
        model, params, states, outputs = define_constrained_lotka_volterra_4D()
        @test length(params) == 4  # Only epsilon perturbations
        @test length(states) == 4
        @test length(outputs) == 4
    end

    @testset "Error Function Creation" begin
        model, params, states, outputs = define_daisy_ex3_model_4D()
        p_true = [0.1, 0.2, 0.3, 0.4]
        ic = [1.0, 2.0, 1.0, 1.0]
        time_interval = [0.0, 10.0]

        error_func = make_error_distance(model, outputs, ic, p_true, time_interval, 25)

        # Test at true parameters (should be ~0)
        @test error_func(p_true) < 1e-6

        # Test at perturbed parameters (should be > 0)
        p_test = p_true .+ 0.1
        @test error_func(p_test) > 0
    end

    @testset "Distance Functions" begin
        y_true = [1.0, 2.0, 3.0]
        y_test = [1.1, 2.1, 3.1]

        @test L1_norm(y_true, y_test) > 0
        @test L2_norm(y_true, y_test) > 0
        @test log_L2_norm(y_true, y_test) < 0  # log of small number
    end

    @testset "Display Infrastructure" begin
        # Test DisplayConfig creation
        config = DisplayConfig()
        @test config.use_color == true
        @test config.table_backend == :text
        @test config.precision == 4

        # Test custom config
        custom_config = DisplayConfig(use_color = false, precision = 6)
        @test custom_config.use_color == false
        @test custom_config.precision == 6

        # Test display functions don't error (output is visual, so we just check they run)
        model, params, states, outputs = define_daisy_ex3_model_4D()

        # Suppress output during tests
        original_stdout = stdout
        redirect_stdout(devnull)

        try
            # Test model summary
            @test_nowarn display_model_summary(
                model,
                config = DisplayConfig(use_color = false),
            )

            # Test parameter display
            param_names = [:α, :β, :γ, :δ]
            p_values = [0.1, 0.2, 0.3, 0.4]
            @test_nowarn display_parameters(
                param_names,
                p_values,
                config = DisplayConfig(use_color = false),
            )

            # Test error metrics display
            errors = Dict("L1" => 0.123, "L2" => 0.045)
            @test_nowarn display_error_metrics(
                errors,
                config = DisplayConfig(use_color = false),
            )

            # Test optimization result display
            result = (
                params = p_values,
                error = 0.001,
                iterations = 10,
                converged = true,
                time_elapsed = 1.0,
            )
            @test_nowarn display_optimization_result(
                result,
                config = DisplayConfig(use_color = false),
            )

            # Test time series display
            ic = [1.0, 2.0, 1.0, 1.0]
            time_interval = [0.0, 10.0]
            p_true = [0.1, 0.2, 0.3, 0.4]

            # Bind the initial condition by `unknowns(model)`, which is what
            # `make_error_distance` does; `states` is model_fn's DECLARATION order and the
            # two differ for these models (bead rai5.12).
            problem = ODEProblem(
                model,
                merge(Dict(ModelingToolkit.unknowns(model) .=> ic), Dict(params .=> p_true)),
                time_interval,
            )

            data = sample_data(problem, model, outputs, time_interval, p_true, ic, 10)
            @test_nowarn display_time_series(
                data,
                show_plot = false,
                config = DisplayConfig(use_color = false),
            )
        finally
            redirect_stdout(original_stdout)
        end
    end
end

# Aggregation strategy registry (bead 0iq)
include("test_aggregation_strategies.jl")

# Partial observability helper (bead dds)
include("test_partial_observability.jl")

# TolerantObjective tests
include("test_tolerant_objective.jl")

# Screening tests
include("test_screening.jl")

# Catalogue persistence tests
include("test_catalogue.jl")

# TOML screening config: parse/validate/path-resolution plus one end-to-end
# run_screening_from_config (bead 89rn)
include("test_screening_config.jl")

# Grid-based interestingness scoring tests
include("test_grid_scoring.jl")

# Second difficulty axis — curved multimodality, which interestingness_score
# cannot express because all its weights are positive (bead cbyn.1)
include("test_structure_score.jl")

# End-to-end interestingness pipeline (bead bzf): screen_and_probe →
# score_top_candidates → interestingness_score
include("test_interestingness_end_to_end.jl")

# Candidate-level parallelism for catalogue experiments (bead 1yt)
include("test_catalogue_parallel.jl")

# Threaded Cartesian grid evaluator (bead ychu)
include("test_parallel_eval.jl")

# Cartesian-product gluing for higher-dim test problems (bead zwbs.10.1)
include("test_glued_objectives.jl")

# Recovery regression: globtim enumerates the FULL glued oracle CP set (bead
# zwbs.10.1). Needs the HomotopyContinuation solver extension — present in the
# workspace, absent on the standalone mirror — so guard like test_ext_toml_pipeline.jl.
if Base.identify_package("HomotopyContinuation") !== nothing
    include("test_glued_recovery.jl")
else
    @info "Skipping test_glued_recovery.jl — HomotopyContinuation not resolvable in this environment (public mirror)"
end

# Gradient compatibility tests need ForwardDiff + FiniteDiff, which are NOT
# in DynamicObjectives' [deps] (ForwardDiff support is duck-typed — callers
# bring their own AD backend). The test file exists for manual runs:
#   julia --project=pkg/DynamicObjectives -e 'using Pkg; Pkg.add(["ForwardDiff","FiniteDiff"]); include("test/test_gradient_compatibility.jl")'
# CI should not re-add those deps just to cover one integration test.
# include("test_gradient_compatibility.jl")

# sample_data direct tests
include("test_sample_data.jl")

# create_objective and dimension() tests
include("test_create_objective.jl")

# Experiment utilities tests (build_bounds, CandidateResult, print_candidate_summary_table)
include("test_experiment_utils.jl")

# Validate all registered models
include("test_all_registered_models.jl")

# globtim integration tests (comprehensive test suite)
include("test_globtim_integration.jl")

# DynamicObjectivesGlobtimExt smoke (needs GlobtimPostProcessing, which the
# public standalone mirror cannot resolve — unregistered; workspace runs cover it)
if Base.identify_package("GlobtimPostProcessing") !== nothing
    include("test_ext_toml_pipeline.jl")
else
    @info "Skipping test_ext_toml_pipeline.jl — GlobtimPostProcessing not resolvable in this environment (public mirror)"
end

# Aqua.jl quality assurance (bead eti8)
include("test_aqua.jl")
