"""
Test suite for the TOML screening configuration (src/screening_config.jl).

Covers the parse/validate/resolve layer in full, plus one end-to-end run of the
orchestrator:
- validate_screening_toml — every required section and key, and the bounds shape
- _find_results_root_from — GLOBTIM_RESULTS_ROOT precedence, upward walk, failure
- _resolve_screening_path — absolute, globtim_results/-prefixed, plain relative
- load_screening_config — full parse, defaults for omitted sections, type coercion
- run_screening_from_config — tiny 2D screening end to end into a temp catalogue

Everything writes into mktempdir(); nothing touches the repo's globtim_results/.
"""

using Test
using DynamicObjectives
using DynamicObjectives:
    ScreeningConfig,
    load_screening_config,
    validate_screening_toml,
    run_screening_from_config,
    _resolve_screening_path,
    _find_results_root_from

# Run `f` with GLOBTIM_RESULTS_ROOT set to `val` (or unset when `val === nothing`),
# restoring whatever was there before. The upward-walk tests depend on the variable
# being absent, and the developer running the suite may well have it set.
function _with_results_root(f, val)
    had = haskey(ENV, "GLOBTIM_RESULTS_ROOT")
    old = had ? ENV["GLOBTIM_RESULTS_ROOT"] : ""
    try
        if val === nothing
            delete!(ENV, "GLOBTIM_RESULTS_ROOT")
        else
            ENV["GLOBTIM_RESULTS_ROOT"] = val
        end
        return f()
    finally
        if had
            ENV["GLOBTIM_RESULTS_ROOT"] = old
        else
            delete!(ENV, "GLOBTIM_RESULTS_ROOT")
        end
    end
end

# Smallest dict that passes validation; tests delete from a copy to probe each check.
function _valid_toml_dict()
    return Dict{String,Any}(
        "screening" => Dict{String,Any}("name" => "t"),
        "model" => Dict{String,Any}(
            "model_fn" => "define_lotka_volterra_2D_model_v3_two_outputs",
            "ic" => [0.5, 1.0],
            "bounds" => [[0.5, 2.0], [0.2, 1.0]],
            "time_interval" => [0.0, 10.0],
        ),
        "output" => Dict{String,Any}("catalogue_path" => "/tmp/c.jsonl"),
    )
end

@testset "Screening Config" begin

    # ========================================================================
    # validate_screening_toml
    # ========================================================================

    @testset "validate_screening_toml" begin
        @testset "accepts a minimal valid dict" begin
            @test validate_screening_toml(_valid_toml_dict()) === nothing
        end

        @testset "requires each top-level section" begin
            for section in ("screening", "model", "output")
                d = _valid_toml_dict()
                delete!(d, section)
                @test_throws ErrorException validate_screening_toml(d)
            end
        end

        @testset "requires [screening] name" begin
            d = _valid_toml_dict()
            delete!(d["screening"], "name")
            @test_throws ErrorException validate_screening_toml(d)
        end

        @testset "requires every [model] key" begin
            for key in ("model_fn", "ic", "bounds", "time_interval")
                d = _valid_toml_dict()
                delete!(d["model"], key)
                @test_throws ErrorException validate_screening_toml(d)
            end
        end

        @testset "requires [output] catalogue_path" begin
            d = _valid_toml_dict()
            delete!(d["output"], "catalogue_path")
            @test_throws ErrorException validate_screening_toml(d)
        end

        @testset "rejects malformed bounds" begin
            # Not an array at all
            d = _valid_toml_dict()
            d["model"]["bounds"] = 3.0
            @test_throws ErrorException validate_screening_toml(d)

            # Element is not a pair
            d = _valid_toml_dict()
            d["model"]["bounds"] = [[0.5, 2.0], [0.2]]
            @test_throws ErrorException validate_screening_toml(d)

            # Element is a triple
            d = _valid_toml_dict()
            d["model"]["bounds"] = [[0.5, 2.0, 3.0]]
            @test_throws ErrorException validate_screening_toml(d)

            # Element is a scalar
            d = _valid_toml_dict()
            d["model"]["bounds"] = [1.0, 2.0]
            @test_throws ErrorException validate_screening_toml(d)
        end
    end

    # ========================================================================
    # _find_results_root_from
    # ========================================================================

    @testset "_find_results_root_from" begin
        @testset "GLOBTIM_RESULTS_ROOT wins when it exists" begin
            mktempdir() do dir
                root = joinpath(dir, "explicit_root")
                mkpath(root)
                # A globtim_results/ sitting right next to the start dir must NOT
                # be preferred over the environment variable.
                mkpath(joinpath(dir, "globtim_results"))

                got = _with_results_root(root) do
                    _find_results_root_from(dir)
                end
                @test got == abspath(root)
            end
        end

        @testset "GLOBTIM_RESULTS_ROOT set but missing is an error" begin
            mktempdir() do dir
                missing_root = joinpath(dir, "does_not_exist")
                @test_throws ErrorException _with_results_root(missing_root) do
                    _find_results_root_from(dir)
                end
            end
        end

        @testset "walks up to the nearest globtim_results/" begin
            mktempdir() do dir
                results = joinpath(dir, "globtim_results")
                mkpath(results)
                nested = joinpath(dir, "experiments", "cluster", "deep")
                mkpath(nested)

                got = _with_results_root(nothing) do
                    _find_results_root_from(nested)
                end
                @test got == abspath(results)
            end
        end

        @testset "finds it in the start dir itself" begin
            mktempdir() do dir
                results = joinpath(dir, "globtim_results")
                mkpath(results)
                got = _with_results_root(nothing) do
                    _find_results_root_from(dir)
                end
                @test got == abspath(results)
            end
        end

        @testset "errors when there is nothing to find" begin
            mktempdir() do dir
                # No globtim_results/ here or anywhere up to the filesystem root.
                @test_throws ErrorException _with_results_root(nothing) do
                    _find_results_root_from(dir)
                end
            end
        end
    end

    # ========================================================================
    # _resolve_screening_path
    # ========================================================================

    @testset "_resolve_screening_path" begin
        @testset "absolute paths pass through untouched" begin
            p = joinpath(homedir(), "somewhere", "cat.jsonl")
            @test _resolve_screening_path(p, "/any/config/dir") == p
        end

        @testset "globtim_results/ prefix resolves against the results root" begin
            mktempdir() do dir
                results = joinpath(dir, "globtim_results")
                mkpath(results)
                config_dir = joinpath(dir, "experiments", "cluster")
                mkpath(config_dir)

                got = _with_results_root(nothing) do
                    _resolve_screening_path("globtim_results/lv4d.jsonl", config_dir)
                end
                @test got == joinpath(abspath(results), "lv4d.jsonl")
            end
        end

        @testset "other relative paths resolve against config_dir" begin
            @test _resolve_screening_path("out/cat.jsonl", "/cfg") ==
                  joinpath("/cfg", "out/cat.jsonl")
            @test _resolve_screening_path("cat.jsonl", "/cfg") ==
                  joinpath("/cfg", "cat.jsonl")
        end

        @testset "a path merely containing globtim_results is not special-cased" begin
            # Only the leading prefix triggers root resolution.
            @test _resolve_screening_path("nested/globtim_results/c.jsonl", "/cfg") ==
                  joinpath("/cfg", "nested/globtim_results/c.jsonl")
        end
    end

    # ========================================================================
    # load_screening_config
    # ========================================================================

    @testset "load_screening_config" begin
        @testset "missing file is an error" begin
            @test_throws ErrorException load_screening_config(
                joinpath(mktempdir(), "nope.toml"),
            )
        end

        @testset "parses a fully-specified config" begin
            mktempdir() do dir
                path = joinpath(dir, "full.toml")
                cat_path = joinpath(dir, "out.jsonl")
                write(
                    path,
                    """
                    [screening]
                    name = "full_test"
                    description = "every key set"

                    [model]
                    model_fn = "define_lotka_volterra_2D_model_v3_two_outputs"
                    ic = [0.5, 1.0]
                    bounds = [[0.5, 2.0], [0.2, 1.0]]
                    time_interval = [0.0, 10.0]

                    [sweep]
                    n_candidates = 42
                    margin = 0.25
                    sampling = "sobol"
                    threshold = 1e5
                    numpoints_screen = 11

                    [probe]
                    n_probes = 7
                    numpoints_probe = 13
                    distance_function = "L1_norm"

                    [solver]
                    method = "Vern9"
                    abstol = 1e-8
                    reltol = 1e-9

                    [ranking]
                    min_finite_fraction = 0.75
                    max_noise_ratio = 3.5

                    [output]
                    catalogue_path = "$cat_path"
                    name_prefix = "PFX"
                    top_n = 4
                    description = "cat description"
                    numpoints = 17
                    aggregate_distances = "mean"
                    eval_timeout = 12.5
                    """,
                )

                c = load_screening_config(path)

                @test c isa ScreeningConfig
                @test c.name == "full_test"
                @test c.description == "every key set"

                @test c.model_fn == "define_lotka_volterra_2D_model_v3_two_outputs"
                @test c.ic == [0.5, 1.0]
                @test c.bounds == [(0.5, 2.0), (0.2, 1.0)]
                @test c.bounds isa Vector{Tuple{Float64,Float64}}
                @test c.time_interval == [0.0, 10.0]

                @test c.n_candidates == 42
                @test c.margin == 0.25
                @test c.sampling === :sobol
                @test c.threshold == 1e5
                @test c.numpoints_screen == 11

                @test c.n_probes == 7
                @test c.numpoints_probe == 13
                @test c.distance_function == "L1_norm"

                @test c.solver_method == "Vern9"
                @test c.solver_abstol == 1e-8
                @test c.solver_reltol == 1e-9

                @test c.min_finite_fraction == 0.75
                @test c.max_noise_ratio == 3.5

                @test c.catalogue_path == cat_path
                @test c.name_prefix == "PFX"
                @test c.top_n == 4
                @test c.catalogue_description == "cat description"
                @test c.catalogue_numpoints == 17
                @test c.aggregate_distances == "mean"
                @test c.eval_timeout == 12.5
            end
        end

        @testset "omitted sections fall back to documented defaults" begin
            mktempdir() do dir
                path = joinpath(dir, "minimal.toml")
                cat_path = joinpath(dir, "out.jsonl")
                # Only the three required sections; no [sweep], [probe], [solver],
                # [ranking], and only catalogue_path under [output].
                write(
                    path,
                    """
                    [screening]
                    name = "minimal"

                    [model]
                    model_fn = "define_lotka_volterra_2D_model_v3_two_outputs"
                    ic = [0.5, 1.0]
                    bounds = [[0.5, 2.0], [0.2, 1.0]]
                    time_interval = [0.0, 10.0]

                    [output]
                    catalogue_path = "$cat_path"
                    """,
                )

                c = load_screening_config(path)

                @test c.description == ""

                @test c.n_candidates == 1000
                @test c.margin == 0.1
                @test c.sampling === :random
                @test c.threshold == 1e6
                @test c.numpoints_screen == 20

                @test c.n_probes == 10
                @test c.numpoints_probe == 30
                @test c.distance_function == "L2_norm"

                @test c.solver_method == "Tsit5"
                @test c.solver_abstol == 1e-4
                @test c.solver_reltol == 1e-4

                @test c.min_finite_fraction == 0.5
                @test c.max_noise_ratio == Inf

                @test c.name_prefix == "screened"
                @test c.top_n == 10
                @test c.catalogue_numpoints == 30
                @test c.aggregate_distances == "sum"
                # Absent eval_timeout is nothing, not a sentinel number.
                @test c.eval_timeout === nothing
            end
        end

        @testset "integer-valued TOML floats are coerced" begin
            mktempdir() do dir
                path = joinpath(dir, "ints.toml")
                cat_path = joinpath(dir, "out.jsonl")
                # ic / bounds / time_interval written as TOML integers, and
                # margin / abstol likewise — all must land as Float64.
                write(
                    path,
                    """
                    [screening]
                    name = "ints"

                    [model]
                    model_fn = "define_lotka_volterra_2D_model_v3_two_outputs"
                    ic = [1, 2]
                    bounds = [[0, 2], [0, 1]]
                    time_interval = [0, 10]

                    [sweep]
                    margin = 1

                    [solver]
                    abstol = 1

                    [output]
                    catalogue_path = "$cat_path"
                    eval_timeout = 30
                    """,
                )

                c = load_screening_config(path)

                @test c.ic == [1.0, 2.0]
                @test c.ic isa Vector{Float64}
                @test c.bounds == [(0.0, 2.0), (0.0, 1.0)]
                @test c.time_interval == [0.0, 10.0]
                @test c.margin === 1.0
                @test c.solver_abstol === 1.0
                @test c.eval_timeout === 30.0
            end
        end

        @testset "catalogue_path is resolved at parse time" begin
            mktempdir() do dir
                results = joinpath(dir, "globtim_results")
                mkpath(results)
                config_dir = joinpath(dir, "configs")
                mkpath(config_dir)
                path = joinpath(config_dir, "rel.toml")
                write(
                    path,
                    """
                    [screening]
                    name = "rel"

                    [model]
                    model_fn = "define_lotka_volterra_2D_model_v3_two_outputs"
                    ic = [0.5, 1.0]
                    bounds = [[0.5, 2.0], [0.2, 1.0]]
                    time_interval = [0.0, 10.0]

                    [output]
                    catalogue_path = "globtim_results/rel_cat.jsonl"
                    """,
                )

                c = _with_results_root(nothing) do
                    load_screening_config(path)
                end
                @test c.catalogue_path == joinpath(abspath(results), "rel_cat.jsonl")
                @test isabspath(c.catalogue_path)
            end
        end

        @testset "invalid TOML content is rejected by validation" begin
            mktempdir() do dir
                path = joinpath(dir, "bad.toml")
                # [output] present but catalogue_path missing.
                write(
                    path,
                    """
                    [screening]
                    name = "bad"

                    [model]
                    model_fn = "define_lotka_volterra_2D_model_v3_two_outputs"
                    ic = [0.5, 1.0]
                    bounds = [[0.5, 2.0], [0.2, 1.0]]
                    time_interval = [0.0, 10.0]

                    [output]
                    name_prefix = "X"
                    """,
                )
                @test_throws ErrorException load_screening_config(path)
            end
        end
    end

    # ========================================================================
    # run_screening_from_config — end to end
    # ========================================================================

    @testset "run_screening_from_config" begin
        @testset "tiny 2D screening writes a catalogue" begin
            mktempdir() do dir
                path = joinpath(dir, "tiny.toml")
                cat_path = joinpath(dir, "tiny_catalogue.jsonl")
                # Deliberately small: 8 candidates, 2 probes, coarse sampling.
                # This is a wiring test for the orchestrator, not a screening run.
                write(
                    path,
                    """
                    [screening]
                    name = "tiny"
                    description = "orchestrator wiring check"

                    [model]
                    model_fn = "define_lotka_volterra_2D_model_v3_two_outputs"
                    ic = [0.5, 1.0]
                    bounds = [[0.5, 2.0], [0.2, 1.0]]
                    time_interval = [0.0, 10.0]

                    [sweep]
                    n_candidates = 8
                    numpoints_screen = 10

                    [probe]
                    n_probes = 2
                    numpoints_probe = 10

                    [output]
                    catalogue_path = "$cat_path"
                    name_prefix = "TINY"
                    top_n = 2
                    numpoints = 10
                    """,
                )

                result = run_screening_from_config(path; io = devnull, verbose = false)

                @test result isa Dict{Symbol,Any}
                @test issetequal(
                    keys(result),
                    (:config, :screening_result, :entries, :catalogue_path),
                )
                @test result[:config] isa ScreeningConfig
                @test result[:config].name == "tiny"
                @test result[:catalogue_path] == cat_path

                entries = result[:entries]
                @test entries isa Vector{CatalogueEntry}
                # LV 2D over these bounds is well behaved, so the sweep should
                # yield candidates and the catalogue should exist on disk.
                @test !isempty(entries)
                @test length(entries) <= 2          # top_n
                @test isfile(cat_path)

                for e in entries
                    @test startswith(e.name, "TINY")
                    @test length(e.p_true) == 2
                    @test all(isfinite, e.p_true)
                end

                # Round-trips through the catalogue reader.
                reloaded = load_catalogue(cat_path)
                @test length(reloaded) == length(entries)
                @test [e.name for e in reloaded] == [e.name for e in entries]
            end
        end

        @testset "writes nothing but still returns when no candidate ranks" begin
            mktempdir() do dir
                path = joinpath(dir, "empty.toml")
                cat_path = joinpath(dir, "empty_catalogue.jsonl")
                # A threshold this small rejects every trajectory (LV amplitudes are
                # order 1), so the sweep yields no valid candidates, `ranked` comes
                # back empty and the early-return branch is taken. It has to be
                # positive — is_trajectory_bounded asserts that — so 0.0 will not do.
                write(
                    path,
                    """
                    [screening]
                    name = "empty"

                    [model]
                    model_fn = "define_lotka_volterra_2D_model_v3_two_outputs"
                    ic = [0.5, 1.0]
                    bounds = [[0.5, 2.0], [0.2, 1.0]]
                    time_interval = [0.0, 10.0]

                    [sweep]
                    n_candidates = 4
                    numpoints_screen = 10
                    threshold = 1e-12

                    [probe]
                    n_probes = 1
                    numpoints_probe = 10

                    [output]
                    catalogue_path = "$cat_path"
                    """,
                )

                result = run_screening_from_config(path; io = devnull, verbose = false)

                @test isempty(result[:entries])
                @test result[:catalogue_path] == cat_path
                # The early return must not create a catalogue file.
                @test !isfile(cat_path)
            end
        end
    end
end
