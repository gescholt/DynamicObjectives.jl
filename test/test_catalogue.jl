# Tests for CatalogueEntry JSONL persistence

@testset "Catalogue Persistence" begin

    # ── Function Registries ──────────────────────────────────────────────
    @testset "Function Registries" begin
        @testset "Registries populated at load" begin
            # Lower bounds only — the exact name list is enforced by the
            # "All model functions registered" testset below. Using `==`
            # here forced every new model/distance/aggregation to also
            # bump this count, which goes stale silently (CI only catches
            # it when Pkg.test's sandbox doesn't blow up first).
            @test length(MODEL_REGISTRY) >= 18
            @test length(DISTANCE_REGISTRY) >= 3
            @test length(AGGREGATION_REGISTRY) >= 4
        end

        @testset "All model functions registered" begin
            expected_models = [
                "define_daisy_ex3_model_4D",
                "define_daisy_ex3_model_4D_no_input",
                "define_generalized_lotka_volterra_4D",
                "define_constrained_lotka_volterra_4D",
                "define_lotka_volterra_4D_simple",
                "define_lotka_volterra_3D_model",
                "define_lotka_volterra_3D_model_locally_identifiable",
                "define_lotka_volterra_3D_model_v2",
                "define_lotka_volterra_2D_model",
                "define_lotka_volterra_2D_model_v2",
                "define_lotka_volterra_2D_model_v3",
                "define_lotka_volterra_2D_model_v3_two_outputs",
                "define_lotka_volterra_2D_sciml_benchmark",
                "define_fitzhugh_nagumo_3D_model",
                "define_goodwin_oscillator_4D",
                "define_simple_2D_model_locally_identifiable",
                "define_simple_2D_model_locally_identifiable_square",
                "define_simple_1D_model_locally_identifiable",
            ]
            for name in expected_models
                @test haskey(MODEL_REGISTRY, name)
            end
        end

        @testset "Distance functions registered" begin
            @test haskey(DISTANCE_REGISTRY, "L1_norm")
            @test haskey(DISTANCE_REGISTRY, "L2_norm")
            @test haskey(DISTANCE_REGISTRY, "log_L2_norm")
            @test DISTANCE_REGISTRY["L2_norm"] === L2_norm
        end

        @testset "Aggregation functions registered" begin
            @test haskey(AGGREGATION_REGISTRY, "sum")
            @test haskey(AGGREGATION_REGISTRY, "maximum")
            @test haskey(AGGREGATION_REGISTRY, "mean")
            @test haskey(AGGREGATION_REGISTRY, "minimum")
            @test AGGREGATION_REGISTRY["sum"] === sum
        end

        @testset "Custom registration" begin
            my_dist(a, b) = abs(a - b)
            register_distance!("my_dist", my_dist)
            @test DISTANCE_REGISTRY["my_dist"] === my_dist
            # Clean up
            delete!(DISTANCE_REGISTRY, "my_dist")
        end

        @testset "Unregistered function error" begin
            unregistered_fn = x -> x^2
            entry = CatalogueEntry(
                name = "bad",
                description = "unregistered model_fn",
                model_fn = unregistered_fn,
                p_true = [1.0],
                ic = [1.0],
                bounds = [(0.0, 1.0)],
            )
            @test_throws ErrorException Dynamic_objectives._entry_to_dict(entry)
        end

        @testset "Unknown name on load error" begin
            bad_dict = Dict{String,Any}(
                "name" => "bad",
                "description" => "test",
                "model_fn" => "nonexistent_model_function",
                "p_true" => [1.0],
                "ic" => [1.0],
                "bounds" => [[0.0, 1.0]],
                "time_interval" => [0.0, 20.0],
                "numpoints" => 30,
                "distance_function" => "L2_norm",
                "aggregate_distances" => "sum",
                "eval_timeout" => nothing,
            )
            @test_throws ErrorException Dynamic_objectives._dict_to_entry(bad_dict)
        end
    end

    # ── Round-trip serialization ─────────────────────────────────────────
    @testset "Round-trip Serialization" begin
        tmpfile = tempname() * ".jsonl"

        entry = CatalogueEntry(
            name = "LV_2D_v3_roundtrip",
            description = "Round-trip test entry",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = [1.0, 0.5],
            ic = [1.0, 0.5],
            bounds = [(0.0, 3.0), (0.0, 2.0)],
            time_interval = [0.0, 20.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = 10.0,
        )

        try
            save_catalogue(tmpfile, [entry])
            loaded = load_catalogue(tmpfile)

            @test length(loaded) == 1
            e = loaded[1]
            @test e.name == "LV_2D_v3_roundtrip"
            @test e.description == "Round-trip test entry"
            @test e.model_fn === define_lotka_volterra_2D_model_v3
            @test e.p_true == [1.0, 0.5]
            @test e.ic == [1.0, 0.5]
            @test e.bounds == [(0.0, 3.0), (0.0, 2.0)]
            @test e.time_interval == [0.0, 20.0]
            @test e.numpoints == 30
            @test e.distance_function === L2_norm
            @test e.aggregate_distances === sum
            @test e.eval_timeout == 10.0
        finally
            rm(tmpfile, force=true)
        end
    end

    @testset "Round-trip with defaults" begin
        tmpfile = tempname() * ".jsonl"

        entry = CatalogueEntry(
            name = "LV_2D_defaults",
            description = "Uses all defaults",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = [1.0, 0.5],
            ic = [1.0, 0.5],
            bounds = [(0.0, 3.0), (0.0, 2.0)],
        )

        try
            save_catalogue(tmpfile, [entry])
            loaded = load_catalogue(tmpfile)

            e = loaded[1]
            @test e.time_interval == [0.0, 20.0]
            @test e.numpoints == 30
            @test e.distance_function === L2_norm
            @test e.aggregate_distances === sum
            @test e.eval_timeout === nothing
        finally
            rm(tmpfile, force=true)
        end
    end

    @testset "Round-trip all model functions" begin
        tmpfile = tempname() * ".jsonl"

        # Create one entry per registered model to verify all names survive round-trip
        entries = CatalogueEntry[]
        for (name, fn) in MODEL_REGISTRY
            push!(entries, CatalogueEntry(
                name = name,
                description = "auto-test for $name",
                model_fn = fn,
                p_true = [0.1],  # dummy
                ic = [1.0],
                bounds = [(0.0, 1.0)],
            ))
        end

        try
            save_catalogue(tmpfile, entries)
            loaded = load_catalogue(tmpfile)
            @test length(loaded) == length(entries)

            loaded_names = Set(e.name for e in loaded)
            for entry in entries
                @test entry.name in loaded_names
            end
        finally
            rm(tmpfile, force=true)
        end
    end

    # ── Append mode ──────────────────────────────────────────────────────
    @testset "Append Mode" begin
        tmpfile = tempname() * ".jsonl"

        entry1 = CatalogueEntry(
            name = "entry_1",
            description = "first",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = [1.0, 0.5],
            ic = [1.0, 0.5],
            bounds = [(0.0, 3.0), (0.0, 2.0)],
        )
        entry2 = CatalogueEntry(
            name = "entry_2",
            description = "second",
            model_fn = define_lotka_volterra_3D_model,
            p_true = [1.0, 0.5, 0.3],
            ic = [1.0, 0.5, 0.8],
            bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 1.0)],
        )
        entry3 = CatalogueEntry(
            name = "entry_3",
            description = "third",
            model_fn = define_constrained_lotka_volterra_4D,
            p_true = [0.01, 0.02, -0.01, 0.03],
            ic = [0.8, 1.2, 0.8, 1.2],
            bounds = [(-0.3, 0.3), (-0.3, 0.3), (-0.3, 0.3), (-0.3, 0.3)],
        )

        try
            # Save first two
            save_catalogue(tmpfile, [entry1, entry2])
            @test length(load_catalogue(tmpfile)) == 2

            # Append third
            append_catalogue(tmpfile, entry3)
            loaded = load_catalogue(tmpfile)
            @test length(loaded) == 3
            @test loaded[3].name == "entry_3"
            @test loaded[3].model_fn === define_constrained_lotka_volterra_4D
        finally
            rm(tmpfile, force=true)
        end
    end

    @testset "Append creates file" begin
        tmpfile = tempname() * ".jsonl"

        entry = CatalogueEntry(
            name = "fresh_append",
            description = "append to nonexistent file",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = [1.0, 0.5],
            ic = [1.0, 0.5],
            bounds = [(0.0, 3.0), (0.0, 2.0)],
        )

        try
            @test !isfile(tmpfile)
            append_catalogue(tmpfile, entry)
            @test isfile(tmpfile)
            loaded = load_catalogue(tmpfile)
            @test length(loaded) == 1
            @test loaded[1].name == "fresh_append"
        finally
            rm(tmpfile, force=true)
        end
    end

    # ── Overwrite mode ───────────────────────────────────────────────────
    @testset "Save Overwrites" begin
        tmpfile = tempname() * ".jsonl"

        entry1 = CatalogueEntry(
            name = "original",
            description = "will be overwritten",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = [1.0, 0.5],
            ic = [1.0, 0.5],
            bounds = [(0.0, 3.0), (0.0, 2.0)],
        )
        entry2 = CatalogueEntry(
            name = "replacement",
            description = "new content",
            model_fn = define_lotka_volterra_3D_model,
            p_true = [1.0, 0.5, 0.3],
            ic = [1.0, 0.5, 0.8],
            bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 1.0)],
        )

        try
            save_catalogue(tmpfile, [entry1, entry1, entry1])
            @test length(load_catalogue(tmpfile)) == 3

            save_catalogue(tmpfile, [entry2])
            loaded = load_catalogue(tmpfile)
            @test length(loaded) == 1
            @test loaded[1].name == "replacement"
        finally
            rm(tmpfile, force=true)
        end
    end

    # ── Filtering ────────────────────────────────────────────────────────
    @testset "Load Filtering" begin
        tmpfile = tempname() * ".jsonl"

        entries = [
            CatalogueEntry(name="LV_2D_a", description="a", model_fn=define_lotka_volterra_2D_model_v3,
                p_true=[1.0, 0.5], ic=[1.0, 0.5], bounds=[(0.0, 3.0), (0.0, 2.0)]),
            CatalogueEntry(name="LV_2D_b", description="b", model_fn=define_lotka_volterra_2D_model_v3,
                p_true=[2.0, 0.3], ic=[1.0, 0.5], bounds=[(0.0, 3.0), (0.0, 2.0)]),
            CatalogueEntry(name="LV_3D_a", description="c", model_fn=define_lotka_volterra_3D_model,
                p_true=[1.0, 0.5, 0.3], ic=[1.0, 0.5, 0.8], bounds=[(0.0, 3.0), (0.0, 2.0), (0.0, 1.0)]),
            CatalogueEntry(name="DAISY_a", description="d", model_fn=define_daisy_ex3_model_4D,
                p_true=[0.1, 0.2, 0.3, 0.4], ic=[1.0, 2.0, 1.0, 1.0], bounds=[(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)]),
        ]

        try
            save_catalogue(tmpfile, entries)

            # Filter by model_name prefix
            lv2d = load_catalogue(tmpfile, model_name="LV_2D")
            @test length(lv2d) == 2
            @test all(startswith(e.name, "LV_2D") for e in lv2d)

            lv3d = load_catalogue(tmpfile, model_name="LV_3D")
            @test length(lv3d) == 1

            daisy = load_catalogue(tmpfile, model_name="DAISY")
            @test length(daisy) == 1

            none = load_catalogue(tmpfile, model_name="NONEXISTENT")
            @test length(none) == 0

            # max_entries
            limited = load_catalogue(tmpfile, max_entries=2)
            @test length(limited) == 2

            # Combined filter
            combo = load_catalogue(tmpfile, model_name="LV_2D", max_entries=1)
            @test length(combo) == 1
            @test startswith(combo[1].name, "LV_2D")
        finally
            rm(tmpfile, force=true)
        end
    end

    # ── Summary ──────────────────────────────────────────────────────────
    @testset "Catalogue Summary" begin
        tmpfile = tempname() * ".jsonl"

        entries = [
            CatalogueEntry(name="LV_2D_a", description="a", model_fn=define_lotka_volterra_2D_model_v3,
                p_true=[1.0, 0.5], ic=[1.0, 0.5], bounds=[(0.0, 3.0), (0.0, 2.0)]),
            CatalogueEntry(name="LV_2D_b", description="b", model_fn=define_lotka_volterra_2D_model_v3,
                p_true=[2.0, 0.3], ic=[1.0, 0.5], bounds=[(0.0, 3.0), (0.0, 2.0)]),
            CatalogueEntry(name="LV_3D_a", description="c", model_fn=define_lotka_volterra_3D_model,
                p_true=[1.0, 0.5, 0.3], ic=[1.0, 0.5, 0.8], bounds=[(0.0, 3.0), (0.0, 2.0), (0.0, 1.0)]),
        ]

        try
            save_catalogue(tmpfile, entries)
            s = catalogue_summary(tmpfile)

            @test s.total == 3
            @test s.per_model["LV_2D_a"] == 1
            @test s.per_model["LV_2D_b"] == 1
            @test s.per_model["LV_3D_a"] == 1
            @test s.dimensions[2] == 2  # two 2D entries
            @test s.dimensions[3] == 1  # one 3D entry
            @test length(s.models) == 3
            @test issorted(s.models)
        finally
            rm(tmpfile, force=true)
        end
    end

    # ── Error handling ───────────────────────────────────────────────────
    @testset "Error Handling" begin
        @testset "Load nonexistent file" begin
            @test_throws ErrorException load_catalogue("/nonexistent/path.jsonl")
        end

        @testset "Summary nonexistent file" begin
            @test_throws ErrorException catalogue_summary("/nonexistent/path.jsonl")
        end

        @testset "Save creates parent dirs" begin
            tmpdir = tempname()
            tmpfile = joinpath(tmpdir, "sub", "dir", "catalogue.jsonl")
            entry = CatalogueEntry(
                name = "mkdir_test",
                description = "test parent dir creation",
                model_fn = define_lotka_volterra_2D_model_v3,
                p_true = [1.0, 0.5],
                ic = [1.0, 0.5],
                bounds = [(0.0, 3.0), (0.0, 2.0)],
            )
            try
                save_catalogue(tmpfile, [entry])
                @test isfile(tmpfile)
                loaded = load_catalogue(tmpfile)
                @test length(loaded) == 1
            finally
                rm(tmpdir, recursive=true, force=true)
            end
        end
    end

    # ── Non-default distance/aggregation ─────────────────────────────────
    @testset "Non-default Functions" begin
        tmpfile = tempname() * ".jsonl"

        entry = CatalogueEntry(
            name = "custom_fns",
            description = "Uses non-default distance and aggregation",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = [1.0, 0.5],
            ic = [1.0, 0.5],
            bounds = [(0.0, 3.0), (0.0, 2.0)],
            distance_function = L1_norm,
            aggregate_distances = maximum,
        )

        try
            save_catalogue(tmpfile, [entry])
            loaded = load_catalogue(tmpfile)
            e = loaded[1]
            @test e.distance_function === L1_norm
            @test e.aggregate_distances === maximum
        finally
            rm(tmpfile, force=true)
        end
    end

    # ── Empty catalogue ──────────────────────────────────────────────────
    @testset "Empty Catalogue" begin
        tmpfile = tempname() * ".jsonl"

        try
            save_catalogue(tmpfile, CatalogueEntry[])
            @test isfile(tmpfile)
            loaded = load_catalogue(tmpfile)
            @test length(loaded) == 0

            s = catalogue_summary(tmpfile)
            @test s.total == 0
            @test isempty(s.per_model)
            @test isempty(s.dimensions)
            @test isempty(s.models)
        finally
            rm(tmpfile, force=true)
        end
    end

    # ── Screening → Catalogue Bridge ──────────────────────────────────────
    @testset "Screening to Catalogue Bridge" begin

        # Build a synthetic ScreeningResult (no ODE solves needed)
        valid_ptrue = [[0.5, 0.3], [1.0, 0.7], [0.2, 0.9], [0.8, 0.1], [0.4, 0.6]]
        # Positional args: elapsed_seconds, n_total, n_valid, n_rejected_solver_failure,
        #                   n_rejected_amplitude, ms_per_candidate, amplitude_min, amplitude_max, amplitude_median
        diag = SweepDiagnostics(0.1, 5, 5, 0, 0, 20.0, 1.0, 5.0, 3.0)
        # Positional args: valid, rejected, pass_rate, diagnostics
        sweep = SweepResult(valid_ptrue, RejectedCandidate[], 1.0, diag)
        probes = [
            ProbeResult(0.0, [1.0, 2.0], [[0.1, 0.1], [0.2, 0.2]], 0.3, 0.5, 2, 0, 10.0),  # idx 1: dynamic_range=0.3
            ProbeResult(0.0, [1.0, 100.0], [[0.1, 0.1], [0.2, 0.2]], 2.0, 50.0, 2, 0, 10.0), # idx 2: dynamic_range=2.0 (best)
            ProbeResult(0.0, [1.0, 10.0], [[0.1, 0.1], [0.2, 0.2]], 1.0, 10.0, 2, 0, 10.0),  # idx 3: dynamic_range=1.0
            ProbeResult(0.0, [1.0, 50.0], [[0.1, 0.1], [0.2, 0.2]], 1.7, 30.0, 2, 0, 10.0),  # idx 4: dynamic_range=1.7
            ProbeResult(0.0, [1.0, 5.0], [[0.1, 0.1], [0.2, 0.2]], 0.7, 3.0, 2, 0, 10.0),    # idx 5: dynamic_range=0.7
        ]
        # Ranked by dynamic_range descending: [2, 4, 3, 5, 1]
        # Positional args: ranked_indices, filtered_out, filter_reasons
        ranking = RankingResult([2, 4, 3, 5, 1], Int[], Symbol[])
        # Positional args: sweep, probes, ranking
        result = ScreeningResult(sweep, probes, ranking)

        ic = [1.0, 0.5]
        bounds = [(0.0, 3.0), (0.0, 2.0)]
        tspan = [0.0, 20.0]

        @testset "Basic conversion" begin
            entries = screening_to_catalogue(
                result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                top_n = 3, name_prefix = "LV2D_test",
            )
            @test length(entries) == 3

            # First entry should be the highest-ranked (index 2 → p_true = [1.0, 0.7])
            @test entries[1].name == "LV2D_test_1"
            @test entries[1].p_true == [1.0, 0.7]
            @test entries[1].model_fn === define_lotka_volterra_2D_model_v3
            @test entries[1].ic == ic
            @test entries[1].bounds == bounds
            @test entries[1].time_interval == tspan

            # Second entry: index 4 → p_true = [0.8, 0.1]
            @test entries[2].name == "LV2D_test_2"
            @test entries[2].p_true == [0.8, 0.1]

            # Third entry: index 3 → p_true = [0.2, 0.9]
            @test entries[3].name == "LV2D_test_3"
            @test entries[3].p_true == [0.2, 0.9]
        end

        @testset "top_n exceeds available" begin
            entries = screening_to_catalogue(
                result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                top_n = 100,
            )
            # Should return all 5 ranked, not error
            @test length(entries) == 5
        end

        @testset "top_n = 1" begin
            entries = screening_to_catalogue(
                result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                top_n = 1,
            )
            @test length(entries) == 1
            @test entries[1].p_true == [1.0, 0.7]  # best candidate
        end

        @testset "Custom kwargs propagate" begin
            entries = screening_to_catalogue(
                result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                top_n = 1,
                description = "custom desc",
                numpoints = 50,
                distance_function = L1_norm,
                aggregate_distances = maximum,
                eval_timeout = 5.0,
            )
            e = entries[1]
            @test e.description == "custom desc"
            @test e.numpoints == 50
            @test e.distance_function === L1_norm
            @test e.aggregate_distances === maximum
            @test e.eval_timeout == 5.0
        end

        @testset "Invalid top_n" begin
            @test_throws ErrorException screening_to_catalogue(
                result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                top_n = 0,
            )
            @test_throws ErrorException screening_to_catalogue(
                result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                top_n = -1,
            )
        end

        @testset "Empty ranking" begin
            empty_ranking = RankingResult(Int[], [1, 2, 3, 4, 5], fill(:low_finite_fraction, 5))
            empty_result = ScreeningResult(sweep, probes, empty_ranking)
            entries = screening_to_catalogue(
                empty_result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                top_n = 5,
            )
            @test length(entries) == 0
        end

        @testset "Round-trip through JSONL" begin
            tmpfile = tempname() * ".jsonl"
            try
                entries = screening_to_catalogue(
                    result, define_lotka_volterra_2D_model_v3, ic, bounds, tspan;
                    top_n = 3, name_prefix = "roundtrip",
                )
                save_catalogue(tmpfile, entries)
                loaded = load_catalogue(tmpfile)

                @test length(loaded) == 3
                @test loaded[1].name == "roundtrip_1"
                @test loaded[1].p_true == [1.0, 0.7]
                @test loaded[1].model_fn === define_lotka_volterra_2D_model_v3
                @test loaded[2].name == "roundtrip_2"
                @test loaded[3].name == "roundtrip_3"
            finally
                rm(tmpfile, force=true)
            end
        end

        @testset "Tuple time_interval converted to Vector" begin
            entries = screening_to_catalogue(
                result, define_lotka_volterra_2D_model_v3, ic, bounds, (0.0, 15.0);
                top_n = 1,
            )
            @test entries[1].time_interval == [0.0, 15.0]
            @test entries[1].time_interval isa Vector{Float64}
        end
    end

end
