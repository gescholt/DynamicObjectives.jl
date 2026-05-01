# Tests for create_objective and dimension() from catalogue.jl

@testset "create_objective and dimension" begin

    entry = CatalogueEntry(
        name = "test_obj",
        description = "test",
        model_fn = define_lotka_volterra_2D_model_v3,
        p_true = [1.0, 0.5],
        ic = [1.0, 0.5],
        bounds = [(0.5, 2.0), (0.1, 1.0)],
        time_interval = [0.0, 20.0],
        numpoints = 15,
    )

    # ── dimension ────────────────────────────────────────────────────────

    @testset "dimension returns length of p_true" begin
        @test Dynamic_objectives.dimension(entry) == 2
    end

    @testset "dimension for 1D entry" begin
        entry_1d = CatalogueEntry(
            name = "test_1d",
            description = "test",
            model_fn = define_simple_1D_model_locally_identifiable,
            p_true = [1.5],
            ic = [1.0],
            bounds = [(-3.0, 3.0)],
            time_interval = [0.0, 5.0],
            numpoints = 15,
        )
        @test Dynamic_objectives.dimension(entry_1d) == 1
    end

    # ── create_objective ─────────────────────────────────────────────────

    @testset "returns a callable" begin
        obj = create_objective(entry)
        @test obj isa Function
    end

    @testset "near-zero at p_true" begin
        obj = create_objective(entry)
        @test obj(collect(entry.p_true)) < 1e-6
    end

    @testset "positive away from p_true" begin
        obj = create_objective(entry)
        @test obj([1.5, 0.3]) > 0.01
    end

    @testset "returns finite or Inf for any params (no exceptions)" begin
        obj = create_objective(entry)
        # create_objective uses return_inf_on_error=true — should never throw
        result = obj([100.0, 100.0])
        @test result isa Number
        @test result == Inf || isfinite(result)
    end

    @testset "deterministic" begin
        obj = create_objective(entry)
        p = [1.3, 0.7]
        @test obj(p) == obj(p)
    end

    @testset "separate instances are independent" begin
        obj1 = create_objective(entry)
        obj2 = create_objective(entry)
        p = [1.2, 0.6]
        @test obj1(p) == obj2(p)
    end
end
