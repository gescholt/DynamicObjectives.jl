# Optional dependencies (package extensions)
#
# Terminal display needs PrettyTables / Term / UnicodePlots, and the SDIRK solvers
# TRBDF2 / KenCarp4 need OrdinaryDiffEqSDIRK. runtests.jl loads all of them, so
# the "package alone" behaviour can only be observed in a fresh process: the
# subprocess below runs `using DynamicObjectives` and nothing else, checks that
# the extensions stay unloaded, that each moved function raises an error naming
# its package, and that the objective/screening path works; it then loads each
# trigger and checks that the feature appears. It prints one `name<TAB>value`
# line per check, and the parent turns each line into a @test.

using Test

const _OPTIONAL_DEPS_CODE = raw"""
using DynamicObjectives

out(k, v) = println("CHECK\t", k, "\t", v)
errmsg(f) = try
    f()
    "NO ERROR"
catch e
    replace(sprint(showerror, e), '\n' => ' ')
end
loaded(ext) = Base.get_extension(DynamicObjectives, ext) !== nothing

for ext in (:DynamicObjectivesPrettyTablesExt, :DynamicObjectivesTermExt,
            :DynamicObjectivesUnicodePlotsExt, :DynamicObjectivesSDIRKExt)
    out("alone_loaded_$(ext)", loaded(ext))
end
for pkg in ("PrettyTables", "Term", "UnicodePlots", "OrdinaryDiffEqSDIRK")
    out("alone_module_$(pkg)", any(p.name == pkg for p in keys(Base.loaded_modules)))
end

ts = DynamicObjectives.OrderedDict("t" => [0.0, 1.0], "x" => [1.0, 2.0])
out("stub_display_results", errmsg(() -> display_results(["a" => 1])))
out("stub_display_parameters", errmsg(() -> display_parameters([:a], [1.0])))
out("stub_display_time_series_table", errmsg(() -> display_time_series(ts; show_plot = false)))
out("stub_display_section", errmsg(() -> display_section("x")))
out("stub_display_model_summary",
    errmsg(() -> display_model_summary(first(define_lotka_volterra_2D_model_v3()))))
out("stub_progress_bar", errmsg(() -> progress_bar(3)))
out("stub_stop_spinner", errmsg(() -> stop_spinner!(nothing)))
out("stub_display_comparison", errmsg(() -> display_comparison(ts, ts)))
out("solver_TRBDF2", errmsg(() -> DynamicObjectives._resolve_solver("TRBDF2")))
out("solver_KenCarp4", errmsg(() -> DynamicObjectives._resolve_solver("KenCarp4")))
out("export_TRBDF2", isdefined(Main, :TRBDF2))
out("format_duration", format_duration(75))

# Core objective and screening path, with no display package and no SDIRK
entry = CatalogueEntry(
    name = "optional_deps",
    description = "test",
    model_fn = define_lotka_volterra_2D_model_v3,
    p_true = [1.0, 0.5],
    ic = [1.0, 0.5],
    bounds = [(0.5, 2.0), (0.1, 1.0)],
    time_interval = [0.0, 20.0],
    numpoints = 15,
)
obj = create_objective(entry)
out("objective_at_p_true_small", obj(collect(entry.p_true)) < 1e-6)
model = first(define_lotka_volterra_2D_model_v3())
bounded = is_trajectory_bounded(model, [1.0, 0.5], [1.0, 0.5], [0.0, 10.0];
                                solver = DynamicObjectives._resolve_solver("Tsit5"))
out("screening_bounded", bounded.bounded)

# Loading each trigger brings its feature in
using PrettyTables
out("pt_loaded", loaded(:DynamicObjectivesPrettyTablesExt))
out("pt_display_results", errmsg(() -> redirect_stdout(() -> display_results(["a" => 1]), devnull)))
# Table drawn, plots still need UnicodePlots
out("pt_time_series_plot", errmsg(() -> redirect_stdout(() -> display_time_series(ts), devnull)))
using Term
out("term_loaded", loaded(:DynamicObjectivesTermExt))
out("term_display_section", errmsg(() -> redirect_stdout(() -> display_section("x"), devnull)))
out("term_progress", errmsg(() -> redirect_stdout(() -> finish_progress!(progress_bar(3)), devnull)))
using UnicodePlots
out("up_loaded", loaded(:DynamicObjectivesUnicodePlotsExt))
out("up_display_comparison", errmsg(() -> redirect_stdout(() -> display_comparison(ts, ts), devnull)))
out("up_time_series", errmsg(() -> redirect_stdout(() -> display_time_series(ts), devnull)))
using OrdinaryDiffEqSDIRK
out("sdirk_loaded", loaded(:DynamicObjectivesSDIRKExt))
out("sdirk_TRBDF2", typeof(DynamicObjectives._resolve_solver("TRBDF2")) <: OrdinaryDiffEqSDIRK.TRBDF2)
out("sdirk_KenCarp4", typeof(DynamicObjectives._resolve_solver("KenCarp4")) <: OrdinaryDiffEqSDIRK.KenCarp4)
"""

@testset "Optional dependencies (package extensions)" begin
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $_OPTIONAL_DEPS_CODE`
    output = IOBuffer()
    proc = run(pipeline(ignorestatus(cmd); stdout = output, stderr = output))
    text = String(take!(output))
    @test success(proc)
    success(proc) || println(text)

    checks = Dict{String,String}()
    for line in split(text, '\n')
        parts = split(line, '\t')
        length(parts) == 3 && parts[1] == "CHECK" && (checks[parts[2]] = parts[3])
    end
    get_check(k) = get(checks, k, "MISSING CHECK $k")

    @testset "package alone: extensions and their packages not loaded" begin
        for ext in ("PrettyTables", "Term", "UnicodePlots", "SDIRK")
            @test get_check("alone_loaded_DynamicObjectives$(ext)Ext") == "false"
        end
        for pkg in ("PrettyTables", "Term", "UnicodePlots", "OrdinaryDiffEqSDIRK")
            @test get_check("alone_module_$(pkg)") == "false"
        end
    end

    @testset "package alone: stubs name the package to load" begin
        for k in ("stub_display_results", "stub_display_parameters",
                  "stub_display_time_series_table")
            @test occursin("requires PrettyTables.jl", get_check(k))
            @test occursin("using PrettyTables", get_check(k))
        end
        for k in ("stub_display_section", "stub_display_model_summary",
                  "stub_progress_bar", "stub_stop_spinner")
            @test occursin("requires Term.jl", get_check(k))
        end
        @test occursin("requires UnicodePlots.jl", get_check("stub_display_comparison"))
        for s in ("TRBDF2", "KenCarp4")
            msg = get_check("solver_$s")
            @test occursin("solver \"$s\" requires OrdinaryDiffEqSDIRK.jl", msg)
            @test occursin("using OrdinaryDiffEqSDIRK", msg)
        end
        @test get_check("export_TRBDF2") == "false"
        @test get_check("format_duration") == "01:15"
    end

    @testset "package alone: objective and screening path work" begin
        @test get_check("objective_at_p_true_small") == "true"
        @test get_check("screening_bounded") == "true"
    end

    @testset "after loading the triggers" begin
        @test get_check("pt_loaded") == "true"
        @test get_check("pt_display_results") == "NO ERROR"
        @test occursin("requires UnicodePlots.jl", get_check("pt_time_series_plot"))
        @test get_check("term_loaded") == "true"
        @test get_check("term_display_section") == "NO ERROR"
        @test get_check("term_progress") == "NO ERROR"
        @test get_check("up_loaded") == "true"
        @test get_check("up_display_comparison") == "NO ERROR"
        @test get_check("up_time_series") == "NO ERROR"
        @test get_check("sdirk_loaded") == "true"
        @test get_check("sdirk_TRBDF2") == "true"
        @test get_check("sdirk_KenCarp4") == "true"
    end
end
