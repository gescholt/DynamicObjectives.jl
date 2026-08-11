"""
Test suite for evaluate_grid_threaded (src/parallel_eval.jl).

Uses pure closures rather than ODE objectives — the point here is the grid
evaluator's own contract, not the dynamics:
- values match a single-threaded reference, and land at the right indices
- output shape follows grid_axes for 1D/2D/3D
- the factory is called once per outer-axis slice, and instances stay independent
- failures and non-finite results both become Inf
- InterruptException escapes instead of being swallowed as Inf
"""

using Test
using DynamicObjectives
using DynamicObjectives: evaluate_grid_threaded

# Single-threaded reference: same Cartesian sweep, no tasks.
function _reference_grid(f, grid_axes)
    sizes = Tuple(length(ax) for ax in grid_axes)
    out = fill(Inf, sizes)
    for idx in CartesianIndices(sizes)
        p = [grid_axes[d][idx[d]] for d in 1:length(grid_axes)]
        val = f(p)
        out[idx] = isfinite(val) ? val : Inf
    end
    return out
end

# Threads.@threads wraps a throwing body in TaskFailedException / CompositeException,
# and the nesting depends on scheduling, so unwrap rather than assert a shape.
function _unwrap_exception(e)
    while true
        if e isa TaskFailedException
            e = e.task.exception
        elseif e isa CompositeException
            isempty(e.exceptions) && return e
            e = e.exceptions[1]
        else
            return e
        end
    end
end

@testset "evaluate_grid_threaded" begin

    # ========================================================================
    # Basic correctness
    # ========================================================================

    @testset "matches single-threaded reference — 2D" begin
        f = p -> p[1]^2 + 2 * p[2]
        axes = [collect(range(-1.0, 1.0, length = 7)), collect(range(0.0, 3.0, length = 5))]

        got = evaluate_grid_threaded(() -> f, axes; show_progress = false)
        want = _reference_grid(f, axes)

        @test size(got) == (7, 5)
        @test got == want
    end

    @testset "values land at the right indices" begin
        # Encode the index pair in the value so a transposed or shifted write is caught.
        axes = [[10.0, 20.0, 30.0], [1.0, 2.0]]
        got = evaluate_grid_threaded(() -> (p -> p[1] + p[2]), axes; show_progress = false)

        @test got[1, 1] == 11.0
        @test got[1, 2] == 12.0
        @test got[2, 1] == 21.0
        @test got[3, 2] == 32.0
    end

    @testset "shape follows grid_axes — 1D, 3D" begin
        ax1 = [collect(range(0.0, 1.0, length = 4))]
        @test size(evaluate_grid_threaded(() -> sum, ax1; show_progress = false)) == (4,)

        ax3 = [[0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0, 2.0, 3.0]]
        got3 = evaluate_grid_threaded(() -> sum, ax3; show_progress = false)
        @test size(got3) == (3, 2, 4)
        @test got3 == _reference_grid(sum, ax3)
    end

    @testset "single-point grid" begin
        got = evaluate_grid_threaded(() -> (p -> 42.0), [[0.0]]; show_progress = false)
        @test got == fill(42.0, 1)
    end

    @testset "rejects empty grid_axes" begin
        @test_throws ArgumentError evaluate_grid_threaded(
            () -> sum,
            Vector{Vector{Float64}}();
            show_progress = false,
        )
    end

    # ========================================================================
    # Factory pattern
    # ========================================================================

    @testset "factory is called once per outer-axis slice" begin
        n_calls = Threads.Atomic{Int}(0)
        make_obj = function ()
            Threads.atomic_add!(n_calls, 1)
            return sum
        end

        axes = [collect(1.0:6.0), collect(1.0:3.0)]
        evaluate_grid_threaded(make_obj, axes; show_progress = false)

        # One instance per index along the first (parallelized) axis.
        @test n_calls[] == 6
    end

    @testset "instances are independent" begin
        # Each objective counts only its own calls. If the factory handed back one
        # shared instance, a single counter would see all 12 evaluations rather
        # than each of the four seeing 3.
        #
        # No locking needed: the factory loop is sequential, and each slot is then
        # written only by its own outer-axis iteration, which runs on one thread.
        counters = Int[]
        make_obj = function ()
            push!(counters, 0)
            my_slot = length(counters)
            return function (p)
                counters[my_slot] += 1
                return sum(p)
            end
        end

        axes = [collect(1.0:4.0), collect(1.0:3.0)]
        evaluate_grid_threaded(make_obj, axes; show_progress = false)

        @test length(counters) == 4
        @test all(==(3), counters)
    end

    # ========================================================================
    # Error handling
    # ========================================================================

    @testset "thrown errors become Inf" begin
        # Fails on the second outer slice only; the rest must still be computed.
        f = p -> p[1] == 2.0 ? error("boom") : p[1] + p[2]
        axes = [[1.0, 2.0, 3.0], [10.0, 20.0]]

        got = evaluate_grid_threaded(() -> f, axes; show_progress = false)

        @test got[2, 1] == Inf
        @test got[2, 2] == Inf
        @test got[1, 1] == 11.0
        @test got[3, 2] == 23.0
    end

    @testset "non-finite results become Inf" begin
        vals = Dict(1.0 => NaN, 2.0 => -Inf, 3.0 => Inf, 4.0 => 7.0)
        got = evaluate_grid_threaded(
            () -> (p -> vals[p[1]]),
            [[1.0, 2.0, 3.0, 4.0]];
            show_progress = false,
        )

        # NaN and -Inf are normalized to +Inf, not passed through.
        @test got[1] == Inf
        @test got[2] == Inf
        @test got[3] == Inf
        @test got[4] == 7.0
        @test !any(isnan, got)
    end

    @testset "InterruptException propagates" begin
        f = p -> throw(InterruptException())
        axes = [[1.0, 2.0], [1.0]]

        caught = try
            evaluate_grid_threaded(() -> f, axes; show_progress = false)
            nothing
        catch e
            _unwrap_exception(e)
        end

        @test caught isa InterruptException
    end

    # ========================================================================
    # Progress reporting
    # ========================================================================

    @testset "show_progress does not change results" begin
        f = p -> p[1] * p[2]
        axes = [[1.0, 2.0, 3.0], [4.0, 5.0]]

        quiet = evaluate_grid_threaded(() -> f, axes; show_progress = false)
        loud = redirect_stdout(devnull) do
            evaluate_grid_threaded(() -> f, axes; show_progress = true)
        end

        @test loud == quiet
    end

    @testset "periodic progress reporter runs" begin
        # The reporter is an @async task that only prints while evaluations are
        # still outstanding. On a fast grid the main loop finishes before the task
        # is ever scheduled, so the periodic branch needs an objective slow enough
        # to yield — hence the sleep. "ETA" appears only in the periodic line, not
        # in the final summary, so it is what distinguishes the two.
        f = function (p)
            sleep(0.05)
            return p[1]
        end
        axes = [[1.0, 2.0], [1.0, 2.0]]

        # redirect_stdout needs a real stream, not an IOBuffer — go via a temp file.
        got = Ref{Any}(nothing)
        printed = mktemp() do path, io
            redirect_stdout(io) do
                got[] = evaluate_grid_threaded(() -> f, axes; show_progress = true)
            end
            flush(io)
            return read(path, String)
        end
        got = got[]

        @test got == [1.0 1.0; 2.0 2.0]
        @test occursin("ETA", printed)
        @test occursin("finite values", printed)   # final summary line
    end
end
