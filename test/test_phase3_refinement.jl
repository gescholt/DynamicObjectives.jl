"""
Test suite for Phase 3 globtimpostprocessing refinement

Tests the local refinement of critical points found by globtimcore using
gradient-free optimization methods.

Test structure:
- Phase 3.1: Configuration and setup
- Phase 3.2: Basic refinement functionality
- Phase 3.3: Refinement accuracy and improvement
- Phase 3.4: Error handling and edge cases
"""

using Test
using Dynamic_objectives
using LinearAlgebra
using CSV
using DataFrames

# Test helper functions
function create_mock_critical_points(output_dir, objective, bounds; n_points=5, degree=8)
    """Create mock critical points CSV for testing"""
    mkpath(output_dir)

    # Generate grid points
    dimension = length(bounds)
    grid_ranges = [range(b[1], b[2], length=n_points) for b in bounds]

    # Sample points and evaluate
    points = []
    for indices in Iterators.product([1:n_points for _ in 1:dimension]...)
        point = [grid_ranges[i][indices[i]] for i in 1:dimension]
        value = objective(point)
        if isfinite(value)
            push!(points, (point, value))
        end
    end

    # Sort by objective value
    sort!(points, by = x -> x[2])

    # Take best points
    best_points = points[1:min(n_points, length(points))]

    # Create DataFrame
    df = DataFrame()
    df.point_id = 1:length(best_points)

    for i in 1:dimension
        df[!, "p$i"] = [p[1][i] for p in best_points]
    end

    df.objective_value = [p[2] for p in best_points]

    # Write CSV
    csv_path = joinpath(output_dir, "critical_points_deg_$(degree).csv")
    CSV.write(csv_path, df)

    return csv_path, best_points
end

@testset "Phase 3: Refinement Integration" begin

    # ==============================================================================
    # Phase 3.1: Configuration and Setup
    # ==============================================================================

    @testset "Phase 3.1: RefinementConfig" begin

        @testset "ode_refinement_config default values" begin
            config = ode_refinement_config()

            @test config isa RefinementConfig
            @test config.max_time_per_point == 60.0
            @test config.optimization_method == :nelder_mead
            @test config.max_iterations == 1000
            @test config.x_tol == 1e-6
            @test config.f_tol == 1e-6
            @test config.distance_metric == :L2
            @test config.verbose == true
        end

        @testset "ode_refinement_config custom values" begin
            config = ode_refinement_config(
                max_time_per_point = 120.0,
                max_iterations = 500,
                x_tol = 1e-4,
                verbose = false
            )

            @test config.max_time_per_point == 120.0
            @test config.max_iterations == 500
            @test config.x_tol == 1e-4
            @test config.verbose == false
        end
    end

    # ==============================================================================
    # Phase 3.2: Basic Refinement Functionality
    # ==============================================================================

    @testset "Phase 3.2: Basic Refinement" begin

        @testset "Refine 2D LV model" begin
            # Create test model
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]
            bounds = [(0.0, 3.0), (0.0, 2.0)]

            # Create objective function
            objective = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm, first, nothing;
                return_inf_on_error = true,
                eval_timeout = 10.0
            )

            # Create mock critical points
            output_dir = mktempdir()
            csv_path, mock_points = create_mock_critical_points(
                output_dir, objective, bounds, n_points=3
            )

            @test isfile(csv_path)

            # Configure refinement
            config = ode_refinement_config(
                max_time_per_point = 30.0,
                verbose = false
            )

            # Run refinement
            result = refine_experiment_results(output_dir, objective, config)

            # Verify results structure
            @test haskey(result, :n_raw)
            @test haskey(result, :n_converged)
            @test haskey(result, :mean_improvement)
            @test haskey(result, :best_raw_value)
            @test haskey(result, :best_refined_value)
            @test haskey(result, :refined_points)
            @test haskey(result, :best_refined_idx)
            @test haskey(result, :refinement_details)

            # Basic sanity checks
            @test result[:n_raw] == 3
            @test result[:n_converged] >= 0
            @test result[:n_converged] <= result[:n_raw]
            @test result[:mean_improvement] >= 1.0  # Should improve or stay same
            @test result[:best_refined_value] <= result[:best_raw_value]  # Should improve
            @test length(result[:refined_points]) == result[:n_raw]

            # Clean up
            rm(output_dir, recursive=true)
        end

        @testset "Output files generated" begin
            # Create test model
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]
            bounds = [(0.0, 3.0), (0.0, 2.0)]

            objective = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm, first, nothing;
                return_inf_on_error = true
            )

            # Create mock data
            output_dir = mktempdir()
            create_mock_critical_points(output_dir, objective, bounds, n_points=2)

            config = ode_refinement_config(verbose = false)
            result = refine_experiment_results(output_dir, objective, config)

            # Check output files exist
            @test isfile(joinpath(output_dir, "critical_points_refined_deg_8.csv"))
            @test isfile(joinpath(output_dir, "refinement_comparison_deg_8.csv"))
            @test isfile(joinpath(output_dir, "refinement_summary.json"))

            # Clean up
            rm(output_dir, recursive=true)
        end
    end

    # ==============================================================================
    # Phase 3.3: Refinement Accuracy and Improvement
    # ==============================================================================

    @testset "Phase 3.3: Refinement Accuracy" begin

        @testset "Improvement over raw critical points" begin
            # Create simple 2D model
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]
            bounds = [(0.0, 3.0), (0.0, 2.0)]

            objective = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm, first, nothing;
                return_inf_on_error = true
            )

            # Create mock critical points with some offset from true params
            output_dir = mktempdir()

            # Manually create points near true params
            near_true_points = [
                [1.1, 0.6],  # Close to true
                [1.2, 0.4],  # Close to true
                [0.9, 0.55]  # Close to true
            ]

            values = [objective(p) for p in near_true_points]

            df = DataFrame(
                point_id = 1:3,
                p1 = [p[1] for p in near_true_points],
                p2 = [p[2] for p in near_true_points],
                objective_value = values
            )

            CSV.write(joinpath(output_dir, "critical_points_deg_8.csv"), df)

            # Refine
            config = ode_refinement_config(
                max_time_per_point = 30.0,
                verbose = false
            )

            result = refine_experiment_results(output_dir, objective, config)

            # Verify improvement
            @test result[:best_refined_value] < result[:best_raw_value]
            @test result[:mean_improvement] > 1.0

            # Verify convergence
            @test result[:n_converged] > 0

            # Clean up
            rm(output_dir, recursive=true)
        end

        @testset "Near-optimal points converge quickly" begin
            # Points very close to true params should converge fast
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]

            objective = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm, first, nothing;
                return_inf_on_error = true
            )

            output_dir = mktempdir()

            # Points very close to optimum
            near_optimal = [[1.01, 0.51]]
            values = [objective(p) for p in near_optimal]

            df = DataFrame(
                point_id = [1],
                p1 = [near_optimal[1][1]],
                p2 = [near_optimal[1][2]],
                objective_value = values
            )

            CSV.write(joinpath(output_dir, "critical_points_deg_8.csv"), df)

            config = ode_refinement_config(
                max_time_per_point = 30.0,
                verbose = false
            )

            result = refine_experiment_results(output_dir, objective, config)

            # Should converge
            @test result[:n_converged] == 1

            # Should be very close to true params
            best_params = result[:refined_points][result[:best_refined_idx]]
            recovery_error = norm(best_params .- p_true) / norm(p_true)
            @test recovery_error < 0.05  # Within 5%

            # Clean up
            rm(output_dir, recursive=true)
        end
    end

    # ==============================================================================
    # Phase 3.4: Error Handling and Edge Cases
    # ==============================================================================

    @testset "Phase 3.4: Error Handling" begin

        @testset "Timeout handling" begin
            # Create model with potential timeout issues
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]

            objective = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm, first, nothing;
                return_inf_on_error = true,
                eval_timeout = 5.0
            )

            output_dir = mktempdir()
            create_mock_critical_points(output_dir, objective, [(0.0, 3.0), (0.0, 2.0)], n_points=2)

            # Use very short timeout
            config = ode_refinement_config(
                max_time_per_point = 0.1,  # Very short
                verbose = false
            )

            result = refine_experiment_results(output_dir, objective, config)

            # Should complete without errors (even if timeouts occur)
            @test result[:n_raw] == 2

            # Clean up
            rm(output_dir, recursive=true)
        end

        @testset "No critical points found error" begin
            # Empty directory should error gracefully
            output_dir = mktempdir()

            objective = x -> sum(x.^2)
            config = ode_refinement_config(verbose = false)

            @test_throws ErrorException refine_experiment_results(output_dir, objective, config)

            # Clean up
            rm(output_dir, recursive=true)
        end

        @testset "Invalid objective values" begin
            # Test handling of Inf/NaN objective values
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]

            objective = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm, first, nothing;
                return_inf_on_error = true
            )

            output_dir = mktempdir()

            # Create points that might cause issues
            bad_points = [
                [100.0, 100.0],  # Extreme values
                [-10.0, -10.0]   # Negative (may cause ODE failure)
            ]

            values = [objective(p) for p in bad_points]

            df = DataFrame(
                point_id = 1:2,
                p1 = [p[1] for p in bad_points],
                p2 = [p[2] for p in bad_points],
                objective_value = values
            )

            CSV.write(joinpath(output_dir, "critical_points_deg_8.csv"), df)

            config = ode_refinement_config(
                max_time_per_point = 30.0,
                verbose = false
            )

            # Should handle gracefully without crashing
            result = refine_experiment_results(output_dir, objective, config)

            @test result[:n_raw] == 2
            @test isfinite(result[:n_converged])

            # Clean up
            rm(output_dir, recursive=true)
        end
    end

    # ==============================================================================
    # Success Criteria
    # ==============================================================================

    @testset "Success Criteria Summary" begin
        @testset "Refinement improves accuracy" begin
            # Overall test: refinement should consistently improve results
            model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
            p_true = [1.0, 0.5]
            ic = [1.0, 0.5]
            bounds = [(0.0, 3.0), (0.0, 2.0)]

            objective = make_error_distance(
                model, outputs, ic, p_true,
                [0.0, 20.0], 30,
                L2_norm, first, nothing;
                return_inf_on_error = true
            )

            output_dir = mktempdir()
            create_mock_critical_points(output_dir, objective, bounds, n_points=5)

            config = ode_refinement_config(verbose = false)
            result = refine_experiment_results(output_dir, objective, config)

            # Success criteria from POSTPROCESSING_REQUIREMENTS.md:
            # - Mean improvement > 1.0 (refinement helps)
            # - Best refined value <= best raw value (monotonic improvement)
            # - Convergence rate > 0 (at least some points converge)

            @test result[:mean_improvement] >= 1.0
            @test result[:best_refined_value] <= result[:best_raw_value]
            @test result[:n_converged] > 0

            # Clean up
            rm(output_dir, recursive=true)
        end
    end
end
