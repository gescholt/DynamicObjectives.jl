#!/usr/bin/env julia
#
# LV-2D Landscape Explorer
#
# Interactive 2D level-set visualization of ODE parameter estimation landscapes.
# Cmd+click (macOS) or Ctrl+click (Linux) to place candidate critical points,
# then "Refine" to run trust-region Newton and classify via Hessian eigenvalues.
# Newton iteration traces are shown as dashed lines.
# Press 'g' to toggle gradient method between FiniteDiff and ForwardDiff —
# the current method is shown in the plot title and in refinement output.
#
# "Find CPs" uses Globtim: fits a Chebyshev polynomial to the landscape,
# then finds ALL critical points via homotopy continuation. The raw
# polynomial CPs are displayed on the contour — Cmd+click to select
# individual CPs for Newton refinement on the true objective.
#
# Usage:
#   julia --project=profiles/dev --threads=auto pkg/DynamicObjectives/examples/lv2d_landscape_explorer.jl

# ═══════════════════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════════════════

using DynamicObjectives
using GlobtimPostProcessing: refine_point, refinement_method, NewtonCP, NewtonMinimize
using GlobtimPostProcessing: CriticalPointRefinementResult
using Globtim: TestInput, Constructor, relative_l2_error, solve_and_transform, evaluate
using Printf

const CATALOGUE_PATH =
    joinpath(@__DIR__, "..", "paper", "catalogue", "lv2d_all_catalogue.jsonl")

# ═══════════════════════════════════════════════════════════════════════════════
# Model selection
# ═══════════════════════════════════════════════════════════════════════════════

println("Loading catalogue...")
entries = load_catalogue(CATALOGUE_PATH)

println("\n" * "="^60)
println("  LV-2D Landscape Explorer")
println("="^60)
println("\nAvailable 2D models:\n")
for (i, e) in enumerate(entries)
    p_str = join([@sprintf("%.2f", p) for p in e.p_true], ", ")
    b_str = join([@sprintf("[%.1f,%.1f]", b[1], b[2]) for b in e.bounds], " x ")
    println("  $i) $(rpad(e.name, 16))  p_true=[$p_str]  bounds=$b_str")
end

print("\nSelect model (1-$(length(entries))): ")
idx = parse(Int, strip(readline()))
1 <= idx <= length(entries) || error("Invalid selection: $idx")
entry = entries[idx]

println("\n  Selected: $(entry.name)")
println("  p_true = $(entry.p_true)")
println("  bounds = $(entry.bounds)")

# ═══════════════════════════════════════════════════════════════════════════════
# Objective + threaded grid evaluator
# ═══════════════════════════════════════════════════════════════════════════════

println("\n  Creating objectives...")
objective = create_objective(entry)
refine_objective = create_refinement_objective(entry)
println("  f(p_true)       = $(objective(collect(entry.p_true)))")
if refine_objective !== objective
    println("  f_refine(p_true) = $(refine_objective(collect(entry.p_true)))")
    println("  ⚠  Refinement uses a different objective (non-log distance)")
end

println()
println("  ╔══════════════════════════════════════════════╗")
println("  ║  Gradient method: FiniteDiff (default)       ║")
println("  ║  Press 'g' in plot window to toggle          ║")
println("  ║  FiniteDiff ↔ ForwardDiff                    ║")
println("  ╚══════════════════════════════════════════════╝")
println()

# Threaded grid evaluator: replaces the default single-threaded comprehension
# in interactive_levelset_explorer with parallel ODE evaluation.
function threaded_grid(f, xs, ys)
    make_obj = () -> create_objective(entry)
    evaluate_grid_threaded(make_obj, [collect(Float64, xs), collect(Float64, ys)])
end

# ═══════════════════════════════════════════════════════════════════════════════
# Refine callback (trust-region Newton, toggleable gradient method)
# ═══════════════════════════════════════════════════════════════════════════════

bounds_tuples = [Tuple(b) for b in entry.bounds]

# Toggleable gradient method — press 'g' in the plot window to switch
const _grad_method = Ref{Symbol}(:finitediff)

function _grad_method_label()
    _grad_method[] == :finitediff && return "FiniteDiff"
    _grad_method[] == :forwarddiff && return "ForwardDiff"
    return String(_grad_method[])
end

function refine_callback(
    point::Vector{Float64};
    tol::Float64 = 1e-8,
    accept_tol::Float64 = Inf,
    f_accept_tol::Union{Nothing,Float64} = nothing,
    max_iterations::Int = 100,
    hessian_tol::Float64 = 1e-6,
    hessian_relative_tol::Float64 = 0.0,
    trust_radius_fraction::Float64 = 0.1,
    patience::Int = 10,
    min_improvement_ratio::Float64 = 0.99,
    mode::Symbol = :critical_point,
)
    m = if mode == :minimize
        NewtonMinimize(;
            gradient_method = _grad_method[],
            tol,
            accept_tol,
            f_accept_tol,
            max_iterations,
            hessian_tol,
            hessian_relative_tol,
            trust_radius_fraction,
            patience,
            min_improvement_ratio,
        )
    else
        NewtonCP(;
            gradient_method = _grad_method[],
            tol,
            accept_tol,
            f_accept_tol,
            max_iterations,
            hessian_tol,
            hessian_relative_tol,
            trust_radius_fraction,
            patience,
            min_improvement_ratio,
        )
    end
    refine_point(m, refine_objective, point; bounds = bounds_tuples, trace = true)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Globtim CP finding callback (polynomial + homotopy continuation)
# ═══════════════════════════════════════════════════════════════════════════════

# Cache the TestInput across degree changes — the expensive grid evaluation
# (10K+ ODE solves) only happens once; subsequent degree changes just refit
# the polynomial and re-run HC (fast).
const _cached_TR = Ref{Any}(nothing)

# Cache the latest polynomial for the error window
const _cached_pol = Ref{Any}(nothing)

function find_cps_callback(degree::Int)
    if _cached_TR[] === nothing
        center = [(b[1] + b[2]) / 2 for b in bounds_tuples]
        sample_range = [(b[2] - b[1]) / 2 for b in bounds_tuples]
        println("  Building Chebyshev grid (GN=100, $(101^2) evaluations)...")
        _cached_TR[] = TestInput(
            objective;
            dim = 2,
            center,
            GN = 100,
            sample_range,
            tolerance = nothing,
        )
        println("  Grid built.")
    end

    pol = Constructor(_cached_TR[], degree; basis = :chebyshev)
    _cached_pol[] = pol
    rel_l2 = relative_l2_error(pol)
    cps, solve_time = solve_and_transform(pol, bounds_tuples)

    # Filter to domain
    cps_in = filter(cps) do cp
        all(bounds_tuples[i][1] <= cp[i] <= bounds_tuples[i][2] for i in 1:2)
    end

    # Update the error window if it exists
    if _error_obs[] !== nothing
        update_error_window!(pol, degree, rel_l2)
    end

    return cps_in, rel_l2, solve_time
end

# ═══════════════════════════════════════════════════════════════════════════════
# Polynomial approximation error window
# ═══════════════════════════════════════════════════════════════════════════════

# Observable holding the error grid data — initialized after GLMakie loads.
# Set to `nothing` until the first "Find CPs" call produces a polynomial.
const _error_obs = Ref{Any}(nothing)
const _error_title = Ref{Any}(nothing)
const _error_fig = Ref{Any}(nothing)

"""
Compute log10(|f(x) - poly(x)|) on a grid covering the domain bounds.
Uses threaded ODE evaluation for f(x).
"""
function compute_error_grid(pol, xs, ys)
    nx, ny = length(xs), length(ys)

    # Evaluate f(x) on the grid (threaded — these are ODE solves)
    f_grid = threaded_grid(objective, xs, ys)

    # Evaluate poly(x) on the same grid
    poly_grid = Matrix{Float64}(undef, nx, ny)
    for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
        poly_grid[i, j] = evaluate(pol, [x, y])
    end

    # Compute log10(|f - poly|), clamped to avoid log10(0)
    error_grid = Matrix{Float64}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        err = abs(f_grid[i, j] - poly_grid[i, j])
        error_grid[i, j] = log10(max(err, 1e-20))
    end

    return error_grid
end

"""
Update the error window with a new polynomial (called from find_cps_callback).
"""
function update_error_window!(pol, degree::Int, rel_l2::Float64)
    xs = range(bounds_tuples[1][1], bounds_tuples[1][2]; length = 100)
    ys = range(bounds_tuples[2][1], bounds_tuples[2][2]; length = 100)

    println("  Computing approximation error grid...")
    error_grid = compute_error_grid(pol, xs, ys)

    # Update observables (triggers re-render)
    _error_obs[][] = error_grid
    notify(_error_obs[])
    _error_title[][] = @sprintf("Approximation Error  deg=%d  rel_L2=%.2e", degree, rel_l2)
    notify(_error_title[])
    println("  Error window updated.")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Launch GLMakie explorer
# ═══════════════════════════════════════════════════════════════════════════════

println("\n  Loading GLMakie + GlobtimPlots...")
using GLMakie
GLMakie.activate!(; inline = false)
using GlobtimPlots

bounds_vec = [collect(Float64, b) for b in entry.bounds]

println("  Launching explorer...")
base_title = "$(entry.name)  p_true=$(entry.p_true)"
figs = interactive_levelset_explorer(
    objective,
    bounds_vec;
    refine = refine_callback,
    find_cps = find_cps_callback,
    p_true = entry.p_true,
    compute_grid = threaded_grid,
    enable_click_placement = true,
    resolution = 100,
    n_levels = 30,
    title = "$base_title  [$(_grad_method_label())]",
)

# ═══════════════════════════════════════════════════════════════════════════════
# Create the error visualization window (3rd window)
# ═══════════════════════════════════════════════════════════════════════════════

# Initial empty error grid (populated on first "Find CPs")
xs_err = range(bounds_tuples[1][1], bounds_tuples[1][2]; length = 100)
ys_err = range(bounds_tuples[2][1], bounds_tuples[2][2]; length = 100)
initial_error = fill(NaN, length(xs_err), length(ys_err))

error_data = Observable(initial_error)
error_title_obs = Observable("Approximation Error  (click 'Find CPs' to compute)")

# Store observables in the Refs so find_cps_callback can update them
_error_obs[] = error_data
_error_title[] = error_title_obs

error_fig = Figure(size = (800, 700))
error_ax = Axis(
    error_fig[1, 1];
    xlabel = "p₁",
    ylabel = "p₂",
    title = error_title_obs,
    aspect = DataAspect(),
)

hm = heatmap!(
    error_ax,
    collect(xs_err),
    collect(ys_err),
    error_data;
    colormap = Reverse(:RdYlBu),
    nan_color = :gray90,
)
Colorbar(error_fig[1, 2], hm; label = "log₁₀|f(x) - poly(x)|")

# Mark p_true on the error plot
if entry.p_true !== nothing
    scatter!(
        error_ax,
        [entry.p_true[1]],
        [entry.p_true[2]];
        marker = :star5,
        markersize = 20,
        color = :white,
        strokewidth = 1.5,
        strokecolor = :black,
    )
end

_error_fig[] = error_fig
push!(figs, error_fig)

# ═══════════════════════════════════════════════════════════════════════════════
# Gradient method toggle (press 'g' in plot window)
# ═══════════════════════════════════════════════════════════════════════════════

plot_fig = figs[1]
plot_ax = plot_fig.content[1]  # the Axis in the plot figure

on(events(plot_fig.scene).keyboardbutton) do event
    event.action == Makie.Keyboard.press || return Consume(false)
    if event.key == Makie.Keyboard.g
        # Toggle between :finitediff and :forwarddiff
        _grad_method[] = _grad_method[] == :finitediff ? :forwarddiff : :finitediff
        label = _grad_method_label()
        plot_ax.title[] = "$base_title  [$label]"
        println()
        println("  ╔══════════════════════════════════════════════╗")
        println("  ║  Gradient method → $label$(repeat(" ", 26 - length(label)))║")
        println("  ╚══════════════════════════════════════════════╝")
        println()
        return Consume(true)
    end
    return Consume(false)
end

screens = [display(GLMakie.Screen(), fig) for fig in figs]

println("\n  Explorer ready!")
println("  - Cmd+click (macOS) or Ctrl+click (Linux) to place CP candidates")
println("  - Press 'r' or click 'Refine' to run Newton refinement")
println("  - Press 'g' to toggle gradient method (FiniteDiff ↔ ForwardDiff)")
println("  - Click 'Find CPs' to run Globtim polynomial + HC")
println("  - Press 'v' to refine view at current zoom")
println("  - Press 'l' to toggle log scale")
println("  - Press 'c' to clear CPs")
println("  - The 'Approximation Error' window updates when 'Find CPs' runs")
println("\n  Press Enter to close...")
readline()
for s in screens
    try
        close(s)
    catch
    end
end
println("Done.")
