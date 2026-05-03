# Tests for data_generation.jl — sample_data direct tests

@testset "sample_data" begin

    # Use a simple fast model for all tests
    model, params, states, outputs = define_lotka_volterra_2D_model_v3()
    p_true = [1.0, 0.5]
    ic = [1.0, 0.5]
    time_interval = [0.0, 10.0]

    # Build the ODEProblem the same way make_error_distance does
    problem = ODEProblem(
        model,
        merge(Dict(states .=> ic), Dict(params .=> p_true)),
        time_interval,
    )

    # ── Basic functionality ──────────────────────────────────────────────

    @testset "returns OrderedDict with t key" begin
        data = sample_data(problem, model, outputs, time_interval, p_true, ic, 15)
        @test data isa DataStructures.OrderedDict
        @test haskey(data, "t")
        @test length(data["t"]) == 15
    end

    @testset "correct number of timepoints" begin
        for n in [5, 10, 25, 50]
            data = sample_data(problem, model, outputs, time_interval, p_true, ic, n)
            @test length(data["t"]) == n
            for (key, vals) in data
                key == "t" && continue
                @test length(vals) == n
            end
        end
    end

    @testset "time spans correct interval" begin
        data = sample_data(problem, model, outputs, time_interval, p_true, ic, 20)
        t = data["t"]
        @test first(t) ≈ 0.0
        @test last(t) ≈ 10.0
    end

    @testset "values are finite" begin
        data = sample_data(problem, model, outputs, time_interval, p_true, ic, 15)
        for (key, vals) in data
            @test all(isfinite, vals)
        end
    end

    @testset "deterministic — same params give same data" begin
        data1 = sample_data(problem, model, outputs, time_interval, p_true, ic, 15)
        data2 = sample_data(problem, model, outputs, time_interval, p_true, ic, 15)
        for key in keys(data1)
            @test data1[key] == data2[key]
        end
    end

    @testset "different params give different data" begin
        data1 = sample_data(problem, model, outputs, time_interval, p_true, ic, 15)
        data2 = sample_data(problem, model, outputs, time_interval, [1.5, 0.3], ic, 15)
        # At least one variable should differ
        differs = false
        for key in keys(data1)
            key == "t" && continue
            if data1[key] != data2[key]
                differs = true
                break
            end
        end
        @test differs
    end

    # ── Noise injection ──────────────────────────────────────────────────

    @testset "noise injection changes data" begin
        data_clean = sample_data(problem, model, outputs, time_interval, p_true, ic, 15)
        data_noisy = sample_data(
            problem,
            model,
            outputs,
            time_interval,
            p_true,
            ic,
            15;
            inject_noise = true,
            mean_noise = 0.0,
            stddev_noise = 0.1,
        )
        # Noise should make data differ (with very high probability)
        differs = false
        for key in keys(data_clean)
            key == "t" && continue
            if data_clean[key] != data_noisy[key]
                differs = true
                break
            end
        end
        @test differs
    end

    # ── Custom solver and tolerances ─────────────────────────────────────

    @testset "custom solver" begin
        data = sample_data(
            problem,
            model,
            outputs,
            time_interval,
            p_true,
            ic,
            15;
            solver = Tsit5(),
            abstol = 1e-6,
            reltol = 1e-6,
        )
        @test length(data["t"]) == 15
        @test all(isfinite, data["t"])
    end

    # ── Uneven sampling ──────────────────────────────────────────────────

    @testset "uneven sampling times" begin
        custom_times = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
        data = sample_data(
            problem,
            model,
            outputs,
            time_interval,
            p_true,
            ic,
            length(custom_times);
            uneven_sampling = true,
            uneven_sampling_times = custom_times,
        )
        @test collect(data["t"]) ≈ custom_times
    end

    @testset "uneven sampling — length mismatch errors" begin
        @test_throws ErrorException sample_data(
            problem,
            model,
            outputs,
            time_interval,
            p_true,
            ic,
            10;
            uneven_sampling = true,
            uneven_sampling_times = [0.0, 1.0],
        )
    end

    @testset "uneven sampling — empty times errors" begin
        @test_throws ErrorException sample_data(
            problem,
            model,
            outputs,
            time_interval,
            p_true,
            ic,
            5;
            uneven_sampling = true,
            uneven_sampling_times = Float64[],
        )
    end
end
