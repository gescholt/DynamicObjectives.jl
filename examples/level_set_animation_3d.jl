#!/usr/bin/env julia
#
# Level Set Animation for 3-Parameter Dynamical Systems — Interactive Explorer
#
# Loads a catalogue of ODE parameter estimation benchmarks and presents an
# interactive terminal menu for visualizing level sets with different entries
# and adjustable settings.
#
# For each selected entry:
#   Phase 1 — Coarse grid over the tight domain (p_true ± fraction)
#   Phase 2 — Adaptive refinement: denser sub-grids around cells where the
#             objective falls within the target level range [0, level_max]
#   Phase 3 — Low-value augmentation: extra-dense sampling in cells near the
#             global minimum (bottom augment_fraction of the level range)
#   Phase 4 — Interactive GLMakie visualization with slider
#
# Usage:
#   # TOML mode (non-interactive — reads [visualization] + [model] from config):
#   julia --project=profiles/dev --threads=auto pkg/DynamicObjectives/examples/level_set_animation_3d.jl experiments/paper/lv3d_wide.toml
#
#   # Interactive mode (catalogue picker + menu):
#   julia --project=profiles/dev --threads=auto pkg/DynamicObjectives/examples/level_set_animation_3d.jl
#   julia --project=profiles/dev --threads=auto pkg/DynamicObjectives/examples/level_set_animation_3d.jl fhn3d_catalogue.jsonl

# ═══════════════════════════════════════════════════════════════════════════════
# Default settings (mutable — adjustable from the interactive menu)
# ═══════════════════════════════════════════════════════════════════════════════

const CATALOGUE_DIR = joinpath(@__DIR__, "..", "paper", "catalogue")

mutable struct Settings
    n_coarse::Int
    n_refine::Int
    level_max::Float64
    domain_mode::Symbol
    tight_frac::Float64
    level_tol::Float64
    record_animation::Bool
    animation_fps::Int
    animation_duration::Int
    augment_enabled::Bool
    augment_fraction::Float64
    augment_n::Int
    near_zero_enabled::Bool
    near_zero_threshold::Float64
    near_zero_n::Int
end

const SETTINGS = Settings(
    10,       # n_coarse (10³ = 1000 evals; use 20 for high-res)
    2,        # n_refine (light refinement; use 4 for detailed)
    1.0,      # level_max
    :catalogue, # domain_mode (:catalogue for full bounds, :tight for zoomed)
    0.1,      # tight_frac (only used when domain_mode=:tight)
    0.005,    # level_tol
    false,    # record_animation
    30,       # animation_fps
    15,       # animation_duration
    true,     # augment_enabled (extra sampling near low-value regions)
    0.25,     # augment_fraction (fraction of level_max below which cells get extra sampling)
    4,        # augment_n (points per axis in augmented cells)
    true,     # near_zero_enabled (extra-dense sampling near global minimum)
    0.5,      # near_zero_threshold (absolute value threshold for near-zero cells)
    8,         # near_zero_n (points per axis in near-zero cells)
)

# ═══════════════════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════════════════

using DynamicObjectives
using StaticArrays
using Printf: @sprintf

const _viz_loaded = Ref(false)
function _ensure_viz_loaded()
    _viz_loaded[] && return
    println("Loading GLMakie + GlobtimPlots...")
    @eval using GLMakie
    @eval GLMakie.activate!(; inline = false)
    @eval using GlobtimPlots
    _viz_loaded[] = true
end

# ═══════════════════════════════════════════════════════════════════════════════
# Catalogue selection
# ═══════════════════════════════════════════════════════════════════════════════

"""
    list_3d_catalogues() -> Vector{@NamedTuple{file::String, n_entries::Int}}

Scan CATALOGUE_DIR for JSONL files where all entries have exactly 3 parameters.
Returns filenames and entry counts, sorted alphabetically.
"""
function list_3d_catalogues()
    result = @NamedTuple{file::String, n_entries::Int}[]
    for f in sort(readdir(CATALOGUE_DIR))
        endswith(f, ".jsonl") || continue
        path = joinpath(CATALOGUE_DIR, f)
        entries = load_catalogue(path)
        isempty(entries) && continue
        all(e -> length(e.p_true) == 3, entries) || continue
        push!(result, (; file = f, n_entries = length(entries)))
    end
    return result
end

"""
    resolve_catalogue_path(arg::AbstractString) -> String

Resolve a CLI argument to a full catalogue path. Accepts either a full path
or a filename relative to CATALOGUE_DIR.
"""
function resolve_catalogue_path(arg::AbstractString)
    if isfile(arg)
        return arg
    end
    full = joinpath(CATALOGUE_DIR, arg)
    isfile(full) || error("Catalogue not found: $arg (searched $full)")
    return full
end

"""
    select_catalogue_interactive() -> Tuple{Vector, String}

Show a numbered list of 3D-compatible catalogues and prompt user to pick one.
Returns `(entries, catalogue_name)`.
"""
function select_catalogue_interactive()
    cats = list_3d_catalogues()
    isempty(cats) && error("No 3D-compatible catalogues found in $CATALOGUE_DIR")

    println("\n  Available 3D catalogues:")
    for (i, c) in enumerate(cats)
        println("    $i) $(rpad(c.file, 30)) ($(c.n_entries) entries)")
    end
    print("\n  Select catalogue [1]: ")
    input = strip(readline())
    idx = isempty(input) ? 1 : parse(Int, input)
    1 <= idx <= length(cats) || error("Invalid selection: $idx")

    chosen = cats[idx]
    path = joinpath(CATALOGUE_DIR, chosen.file)
    entries = load_catalogue(path)
    name = replace(chosen.file, "_catalogue.jsonl" => "", ".jsonl" => "")
    return entries, name
end

# ═══════════════════════════════════════════════════════════════════════════════
# Evaluation engine
# ═══════════════════════════════════════════════════════════════════════════════

"""
    evaluate_entry(entry, s::Settings)

Run coarse + adaptive refinement evaluation for a single catalogue entry.
Returns `(all_points, all_values)` — flat vectors ready for visualization.
"""
function evaluate_entry(entry, s::Settings)
    name = entry.name
    println("\n  ── $name ──")
    println("     p_true = $(entry.p_true)")

    make_obj = () -> create_objective(entry)

    # Auto-calibrate domain: adjust tight_frac so ~50-95% of coarse grid is in range.
    # Uses a small probe grid (10^3) to be fast, then does the real evaluation at n_coarse.
    frac = s.tight_frac
    if s.domain_mode == :tight
        probe_n = 10
        for attempt in 1:12
            probe_bounds = [(p - frac * abs(p), p + frac * abs(p)) for p in entry.p_true]
            probe_axes =
                [collect(range(lo, hi, length = probe_n)) for (lo, hi) in probe_bounds]
            probe_vals = evaluate_grid_threaded(make_obj, probe_axes; show_progress = false)
            fv = probe_vals[isfinite.(probe_vals)]
            pct = isempty(fv) ? 0.0 : count(v -> 0.0 <= v <= s.level_max, fv) / length(fv)
            if pct < 0.50
                frac *= 0.5   # too few in range → shrink domain
            elseif pct > 0.95
                frac *= 2.0   # too many in range → widen domain
            else
                break
            end
        end
        println("     Auto-calibrated tight_frac: $(@sprintf("%.4f", frac))")
    end

    bounds = if s.domain_mode == :catalogue
        entry.bounds
    elseif s.domain_mode == :tight
        [(p - frac * abs(p), p + frac * abs(p)) for p in entry.p_true]
    else
        error("Unknown domain_mode: $(s.domain_mode)")
    end

    grid_axes = [collect(range(lo, hi, length = s.n_coarse)) for (lo, hi) in bounds]
    for (d, (lo, hi)) in enumerate(bounds)
        println("     Axis $d: [$(@sprintf("%.4f", lo)), $(@sprintf("%.4f", hi))]")
    end

    # Auto-calibrate level_max using a cheap probe grid
    probe_n = 6
    probe_axes = [collect(range(lo, hi, length = probe_n)) for (lo, hi) in bounds]
    println("     Probe: $(probe_n)^3 = $(probe_n^3) evals for level_max calibration...")
    probe_vals = evaluate_grid_threaded(make_obj, probe_axes; show_progress = false)
    probe_finite = sort(filter(v -> isfinite(v) && v >= 0.0, vec(probe_vals)))
    if !isempty(probe_finite)
        n_in_range = count(v -> v <= s.level_max, probe_finite)
        pct = 100.0 * n_in_range / length(probe_finite)
        if pct < 5.0
            # Set level_max to p5 of finite values — focuses on basin near minimum
            # Heavy-tailed objectives need aggressive targeting; p25 is still too wide
            target_idx =
                clamp(round(Int, 0.05 * length(probe_finite)), 1, length(probe_finite))
            new_max = probe_finite[target_idx]
            println(
                "     ⚠ Auto-calibrating level_max: $(s.level_max) → $(@sprintf("%.4g", new_max)) (p5 of probe, min=$(@sprintf("%.4g", probe_finite[1])))",
            )
            s.level_max = new_max
        end
    end

    # Scale near_zero_threshold relative to calibrated level_max
    if s.near_zero_enabled && s.near_zero_threshold > s.level_max
        s.near_zero_threshold = s.level_max * 0.05
        println(
            "     Scaled near_zero_threshold → $(@sprintf("%.4g", s.near_zero_threshold))",
        )
    end

    # Phase 1: coarse grid
    println("     Phase 1: $(s.n_coarse)^3 = $(s.n_coarse^3) coarse evals...")
    values_coarse = evaluate_grid_threaded(make_obj, grid_axes; show_progress = true)

    finite_mask = isfinite.(values_coarse)
    n_finite = count(finite_mask)
    finite_vals = values_coarse[finite_mask]
    n_in_range = count(v -> 0.0 <= v <= s.level_max, finite_vals)
    pct_in_range = 100.0 * n_in_range / max(n_finite, 1)
    println(
        "     Finite: $n_finite/$(length(values_coarse)), in [0,$(s.level_max)]: $n_in_range ($(@sprintf("%.1f", pct_in_range))%)",
    )

    # Phase 2: adaptive refinement
    cell_dims = ntuple(_ -> s.n_coarse - 1, 3)
    interesting_cells = Tuple{Int,Int,Int}[]
    for i in 1:cell_dims[1], j in 1:cell_dims[2], k in 1:cell_dims[3]
        is_interesting = false
        for di in 0:1, dj in 0:1, dk in 0:1
            v = values_coarse[i+di, j+dj, k+dk]
            if isfinite(v) && 0.0 <= v <= s.level_max
                is_interesting = true
                break
            end
        end
        is_interesting && push!(interesting_cells, (i, j, k))
    end

    n_cells = length(interesting_cells)
    println(
        "     Phase 2: $n_cells interesting cells → $(n_cells * s.n_refine^3) refinement evals",
    )

    refine_points = Vector{SVector{3,Float64}}()
    sizehint!(refine_points, n_cells * s.n_refine^3)

    for (ci, cj, ck) in interesting_cells
        x_lo, x_hi = grid_axes[1][ci], grid_axes[1][ci+1]
        y_lo, y_hi = grid_axes[2][cj], grid_axes[2][cj+1]
        z_lo, z_hi = grid_axes[3][ck], grid_axes[3][ck+1]

        rx = range(x_lo, x_hi, length = s.n_refine + 2)[2:(end-1)]
        ry = range(y_lo, y_hi, length = s.n_refine + 2)[2:(end-1)]
        rz = range(z_lo, z_hi, length = s.n_refine + 2)[2:(end-1)]

        for x in rx, y in ry, z in rz
            push!(refine_points, SVector{3,Float64}(x, y, z))
        end
    end

    refine_values = fill(Inf, length(refine_points))

    if !isempty(refine_points)
        n_refine_total = length(refine_points)
        n_per_chunk = max(1, cld(n_refine_total, Threads.nthreads() * 4))
        chunks = [
            i:min(i + n_per_chunk - 1, n_refine_total) for i in 1:n_per_chunk:n_refine_total
        ]
        chunk_objs = [make_obj() for _ in 1:length(chunks)]
        completed = Threads.Atomic{Int}(0)
        t_start = time()

        progress_task = @async begin
            while true
                n_done = completed[]
                n_done >= n_refine_total && break
                elapsed = time() - t_start
                rate = n_done / max(elapsed, 1e-6)
                remaining = (n_refine_total - n_done) / max(rate, 1e-6)
                pct = 100.0 * n_done / n_refine_total
                print(
                    "\r  Refine: $n_done/$n_refine_total ($(@sprintf("%.1f", pct))%) — " *
                    "$(@sprintf("%.0f", rate)) evals/s — ETA $(@sprintf("%.0f", remaining))s   ",
                )
                sleep(2.0)
            end
        end

        Threads.@threads for ci in 1:length(chunks)
            obj = chunk_objs[ci]
            for idx in chunks[ci]
                p = collect(refine_points[idx])
                val = try
                    result = obj(p)
                    isfinite(result) ? result : Inf
                catch e
                    e isa InterruptException && rethrow(e)
                    Inf
                end
                refine_values[idx] = val
                Threads.atomic_add!(completed, 1)
            end
        end

        completed[] = n_refine_total
        sleep(0.1)
        elapsed = time() - t_start
        n_refine_in_range =
            count(v -> isfinite(v) && 0.0 <= v <= s.level_max, refine_values)
        println(
            "\r  Refine: done in $(@sprintf("%.1f", elapsed))s — $n_refine_in_range in range" *
            "                                        ",
        )
    end

    # Phase 3: local augmentation of low-value cells
    augment_points = Vector{SVector{3,Float64}}()
    augment_values = Vector{Float64}()

    if s.augment_enabled
        threshold = s.level_max * s.augment_fraction
        low_cells = filter(interesting_cells) do (ci, cj, ck)
            min_v = Inf
            for di in 0:1, dj in 0:1, dk in 0:1
                v = values_coarse[ci+di, cj+dj, ck+dk]
                if isfinite(v) && v < min_v
                    min_v = v
                end
            end
            min_v <= threshold
        end

        n_low = length(low_cells)
        println(
            "     Phase 3: $n_low low-value cells (threshold=$(@sprintf("%.3f", threshold))) → $(n_low * s.augment_n^3) augmentation evals",
        )

        if n_low > 0
            sizehint!(augment_points, n_low * s.augment_n^3)

            for (ci, cj, ck) in low_cells
                x_lo, x_hi = grid_axes[1][ci], grid_axes[1][ci+1]
                y_lo, y_hi = grid_axes[2][cj], grid_axes[2][cj+1]
                z_lo, z_hi = grid_axes[3][ck], grid_axes[3][ck+1]

                rx = range(x_lo, x_hi, length = s.augment_n + 2)[2:(end-1)]
                ry = range(y_lo, y_hi, length = s.augment_n + 2)[2:(end-1)]
                rz = range(z_lo, z_hi, length = s.augment_n + 2)[2:(end-1)]

                for x in rx, y in ry, z in rz
                    push!(augment_points, SVector{3,Float64}(x, y, z))
                end
            end

            resize!(augment_values, length(augment_points))
            fill!(augment_values, Inf)

            n_aug_total = length(augment_points)
            n_per_chunk = max(1, cld(n_aug_total, Threads.nthreads() * 4))
            aug_chunks =
                [i:min(i + n_per_chunk - 1, n_aug_total) for i in 1:n_per_chunk:n_aug_total]
            aug_objs = [make_obj() for _ in 1:length(aug_chunks)]
            aug_completed = Threads.Atomic{Int}(0)
            t_aug_start = time()

            aug_progress = @async begin
                while true
                    n_done = aug_completed[]
                    n_done >= n_aug_total && break
                    elapsed = time() - t_aug_start
                    rate = n_done / max(elapsed, 1e-6)
                    remaining = (n_aug_total - n_done) / max(rate, 1e-6)
                    pct = 100.0 * n_done / n_aug_total
                    print(
                        "\r  Augment: $n_done/$n_aug_total ($(@sprintf("%.1f", pct))%) — " *
                        "$(@sprintf("%.0f", rate)) evals/s — ETA $(@sprintf("%.0f", remaining))s   ",
                    )
                    sleep(2.0)
                end
            end

            Threads.@threads for ci in 1:length(aug_chunks)
                obj = aug_objs[ci]
                for idx in aug_chunks[ci]
                    p = collect(augment_points[idx])
                    val = try
                        result = obj(p)
                        isfinite(result) ? result : Inf
                    catch e
                        e isa InterruptException && rethrow(e)
                        Inf
                    end
                    augment_values[idx] = val
                    Threads.atomic_add!(aug_completed, 1)
                end
            end

            aug_completed[] = n_aug_total
            sleep(0.1)
            elapsed = time() - t_aug_start
            n_aug_in_range =
                count(v -> isfinite(v) && 0.0 <= v <= s.level_max, augment_values)
            println(
                "\r  Augment: done in $(@sprintf("%.1f", elapsed))s — $n_aug_in_range in range" *
                "                                        ",
            )
        end
    end

    # Phase 3b: near-zero densification (absolute threshold, not relative)
    nz_points = Vector{SVector{3,Float64}}()
    nz_values = Vector{Float64}()

    if s.near_zero_enabled
        # Identify cells where any corner value is below the absolute threshold
        nz_cells = filter(interesting_cells) do (ci, cj, ck)
            min_v = Inf
            for di in 0:1, dj in 0:1, dk in 0:1
                v = values_coarse[ci+di, cj+dj, ck+dk]
                if isfinite(v) && v < min_v
                    min_v = v
                end
            end
            min_v <= s.near_zero_threshold
        end

        # Also check augmented/refined points: find cells containing any sub-threshold value
        # Build a set of cell indices that already qualify
        nz_cell_set = Set(nz_cells)
        for (pt, v) in Iterators.flatten((
            zip(refine_points, refine_values),
            zip(augment_points, augment_values),
        ))
            isfinite(v) && v <= s.near_zero_threshold || continue
            # Map point back to coarse cell
            ci = clamp(searchsortedlast(grid_axes[1], pt[1]), 1, s.n_coarse - 1)
            cj = clamp(searchsortedlast(grid_axes[2], pt[2]), 1, s.n_coarse - 1)
            ck = clamp(searchsortedlast(grid_axes[3], pt[3]), 1, s.n_coarse - 1)
            cell = (ci, cj, ck)
            if cell ∉ nz_cell_set
                push!(nz_cell_set, cell)
                push!(nz_cells, cell)
            end
        end

        n_nz = length(nz_cells)
        println(
            "     Phase 3b: $n_nz near-zero cells (threshold=$(@sprintf("%.3f", s.near_zero_threshold))) → $(n_nz * s.near_zero_n^3) near-zero evals",
        )

        if n_nz > 0
            sizehint!(nz_points, n_nz * s.near_zero_n^3)

            for (ci, cj, ck) in nz_cells
                x_lo, x_hi = grid_axes[1][ci], grid_axes[1][ci+1]
                y_lo, y_hi = grid_axes[2][cj], grid_axes[2][cj+1]
                z_lo, z_hi = grid_axes[3][ck], grid_axes[3][ck+1]

                rx = range(x_lo, x_hi, length = s.near_zero_n + 2)[2:(end-1)]
                ry = range(y_lo, y_hi, length = s.near_zero_n + 2)[2:(end-1)]
                rz = range(z_lo, z_hi, length = s.near_zero_n + 2)[2:(end-1)]

                for x in rx, y in ry, z in rz
                    push!(nz_points, SVector{3,Float64}(x, y, z))
                end
            end

            resize!(nz_values, length(nz_points))
            fill!(nz_values, Inf)

            n_nz_total = length(nz_points)
            n_per_chunk = max(1, cld(n_nz_total, Threads.nthreads() * 4))
            nz_chunks =
                [i:min(i + n_per_chunk - 1, n_nz_total) for i in 1:n_per_chunk:n_nz_total]
            nz_objs = [make_obj() for _ in 1:length(nz_chunks)]
            nz_completed = Threads.Atomic{Int}(0)
            t_nz_start = time()

            nz_progress = @async begin
                while true
                    n_done = nz_completed[]
                    n_done >= n_nz_total && break
                    elapsed = time() - t_nz_start
                    rate = n_done / max(elapsed, 1e-6)
                    remaining = (n_nz_total - n_done) / max(rate, 1e-6)
                    pct = 100.0 * n_done / n_nz_total
                    print(
                        "\r  Near-zero: $n_done/$n_nz_total ($(@sprintf("%.1f", pct))%) — " *
                        "$(@sprintf("%.0f", rate)) evals/s — ETA $(@sprintf("%.0f", remaining))s   ",
                    )
                    sleep(2.0)
                end
            end

            Threads.@threads for ci in 1:length(nz_chunks)
                obj = nz_objs[ci]
                for idx in nz_chunks[ci]
                    p = collect(nz_points[idx])
                    val = try
                        result = obj(p)
                        isfinite(result) ? result : Inf
                    catch e
                        e isa InterruptException && rethrow(e)
                        Inf
                    end
                    nz_values[idx] = val
                    Threads.atomic_add!(nz_completed, 1)
                end
            end

            nz_completed[] = n_nz_total
            sleep(0.1)
            elapsed = time() - t_nz_start
            n_nz_in_range = count(v -> isfinite(v) && 0.0 <= v <= s.level_max, nz_values)
            println(
                "\r  Near-zero: done in $(@sprintf("%.1f", elapsed))s — $n_nz_in_range in range" *
                "                                        ",
            )
        end
    end

    # Merge coarse + refined + augmented + near-zero
    n_total_est =
        s.n_coarse^3 + length(refine_points) + length(augment_points) + length(nz_points)
    all_points = Vector{SVector{3,Float64}}()
    all_values = Vector{Float64}()
    sizehint!(all_points, n_total_est)
    sizehint!(all_values, n_total_est)

    for k in 1:s.n_coarse, j in 1:s.n_coarse, i in 1:s.n_coarse
        push!(
            all_points,
            SVector{3,Float64}(grid_axes[1][i], grid_axes[2][j], grid_axes[3][k]),
        )
        v = values_coarse[i, j, k]
        push!(all_values, isfinite(v) ? v : NaN)
    end
    for (pt, v) in zip(refine_points, refine_values)
        push!(all_points, pt)
        push!(all_values, isfinite(v) ? v : NaN)
    end
    for (pt, v) in zip(augment_points, augment_values)
        push!(all_points, pt)
        push!(all_values, isfinite(v) ? v : NaN)
    end
    for (pt, v) in zip(nz_points, nz_values)
        push!(all_points, pt)
        push!(all_values, isfinite(v) ? v : NaN)
    end

    n_total_in_range = count(v -> !isnan(v) && 0.0 <= v <= s.level_max, all_values)
    println(
        "     Total: $(length(all_points)) points, $n_total_in_range in [0,$(s.level_max)]",
    )

    return all_points, all_values
end

# ═══════════════════════════════════════════════════════════════════════════════
# Visualization
# ═══════════════════════════════════════════════════════════════════════════════

"""
    show_visualization(entry, all_points, all_values, s::Settings) -> screen

Build and display the level set figure. Returns the GLMakie screen handle.
"""
function show_visualization(entry, all_points, all_values, s::Settings)
    _ensure_viz_loaded()
    return Base.invokelatest(_show_visualization_impl, entry, all_points, all_values, s)
end

function _show_visualization_impl(entry, all_points, all_values, s::Settings)
    n_total = length(all_points)
    grid_sv = reshape(collect(all_points), (n_total, 1, 1))
    values_viz = reshape(collect(all_values), (n_total, 1, 1))
    z_range = (0.0, s.level_max)

    viz_params = GlobtimPlots.VisualizationParameters(
        point_tolerance = s.level_tol,
        point_window = 2 * s.level_tol,
        fig_size = (1200, 900),
    )

    fig = create_level_set_visualization(values_viz, grid_sv, nothing, z_range, viz_params)
    Label(
        fig[0, 1],
        "$(entry.name): p_true = $(round.(entry.p_true, digits=3))",
        fontsize = 20,
        tellwidth = false,
    )

    screen = display(GLMakie.Screen(), fig)

    if s.record_animation
        output_file = "$(lowercase(entry.name))_level_sets.mp4"
        println("  Recording animation to $output_file ...")
        create_level_set_animation(
            values_viz,
            grid_sv,
            nothing,
            z_range,
            viz_params;
            fps = s.animation_fps,
            duration = s.animation_duration,
            output_file = output_file,
        )
        println("  Animation saved!")
    end

    return screen, fig
end

# ═══════════════════════════════════════════════════════════════════════════════
# Interactive menu helpers
# ═══════════════════════════════════════════════════════════════════════════════

function print_entry_table(entries)
    println()
    for (i, e) in enumerate(entries)
        p_str = join([@sprintf("%.2f", p) for p in e.p_true], ", ")
        println("  $i) $(rpad(e.name, 8))  p_true = [$p_str]")
    end
    println()
end

function print_settings(s::Settings)
    println("  n_coarse         = $(s.n_coarse)      (coarse grid points per axis)")
    println(
        "  n_refine         = $(s.n_refine)       (refinement points per axis per cell)",
    )
    println("  level_max        = $(s.level_max)    (upper bound of level range)")
    println("  domain_mode      = $(s.domain_mode)   (:tight or :catalogue)")
    println("  tight_frac       = $(s.tight_frac)     (p_true ± frac*|p_true|)")
    println("  level_tol        = $(s.level_tol)     (tolerance for level set display)")
    println(
        "  augment_enabled  = $(s.augment_enabled)  (extra sampling near low-value regions)",
    )
    println(
        "  augment_fraction = $(s.augment_fraction)  (fraction of level_max for threshold)",
    )
    println(
        "  augment_n        = $(s.augment_n)       (points per axis in augmented cells)",
    )
    println(
        "  near_zero_enabled   = $(s.near_zero_enabled)  (extra-dense sampling near global minimum)",
    )
    println("  near_zero_threshold = $(s.near_zero_threshold)  (absolute value threshold)")
    println(
        "  near_zero_n         = $(s.near_zero_n)       (points per axis in near-zero cells)",
    )
end

function print_help()
    println("""
  Commands:
    <number>            Plot entry by number (e.g. "3")
    list                Show all catalogue entries
    catalogue           Switch to a different catalogue
    settings            Show current settings
    set                 Change a setting interactively
    quit / q            Exit
    help / ?            Show this help""")
end

"""
    run_set_submenu!(s::Settings)

Interactive submenu for adjusting settings. Walks through each parameter
and asks whether to change it.
"""
function run_set_submenu!(s::Settings)
    println("\n  ── Adjust Settings ──")
    println("  Press Enter to keep current value, or type a new value.\n")

    fields = [
        (:n_coarse, "Coarse grid per axis", Int),
        (:n_refine, "Refine grid per axis/cell", Int),
        (:level_max, "Level range upper bound", Float64),
        (:tight_frac, "Tight domain fraction", Float64),
        (:level_tol, "Level set tolerance", Float64),
        (:domain_mode, "Domain mode (tight/catalogue)", Symbol),
        (:augment_enabled, "Augment low-value cells", Bool),
        (:augment_fraction, "Augment threshold (frac of level_max)", Float64),
        (:augment_n, "Augment points per axis/cell", Int),
        (:near_zero_enabled, "Near-zero densification", Bool),
        (:near_zero_threshold, "Near-zero threshold (absolute)", Float64),
        (:near_zero_n, "Near-zero points per axis/cell", Int),
    ]

    for (field, label, T) in fields
        current = getfield(s, field)
        print("  $label [$(current)]: ")
        input = strip(readline())
        if !isempty(input)
            try
                new_val = if T == Symbol
                    Symbol(input)
                elseif T == Int
                    parse(Int, input)
                elseif T == Float64
                    parse(Float64, input)
                elseif T == Bool
                    lowercase(input) in ("true", "1", "yes") ? true :
                    lowercase(input) in ("false", "0", "no") ? false :
                    error("expected true/false")
                end
                setfield!(s, field, new_val)
                println("    → $(field) = $(new_val)")
            catch e
                println("    Invalid input, keeping $(current)")
            end
        end
    end

    println("\n  Updated settings:")
    print_settings(s)
end

# ═══════════════════════════════════════════════════════════════════════════════
# TOML config support
# ═══════════════════════════════════════════════════════════════════════════════

"""
    load_from_toml(toml_path::String) -> (entry, Settings)

Load visualization settings and catalogue entry from a TOML experiment config.
Requires `[model]` (catalogue_path + entry_name) and optionally `[visualization]`.
"""
function load_from_toml(toml_path::String)
    config = Globtim.load_experiment_config(toml_path)

    config.catalogue_path !== nothing ||
        error("TOML must use catalogue mode ([model] catalogue_path + entry_name)")
    config.entry_name !== nothing || error("TOML must specify [model] entry_name")

    entries = load_catalogue(config.catalogue_path)
    idx = findfirst(e -> e.name == config.entry_name, entries)
    idx !== nothing ||
        error("Entry '$(config.entry_name)' not found in $(config.catalogue_path)")
    entry = entries[idx]

    s = Settings(
        config.viz_n_coarse,
        config.viz_n_refine,
        config.viz_level_max,
        config.viz_domain_mode,
        config.viz_tight_frac,
        config.viz_level_tol,
        config.viz_record_animation,
        config.viz_animation_fps,
        config.viz_animation_duration,
        config.viz_augment_enabled,
        config.viz_augment_fraction,
        config.viz_augment_n,
        config.viz_near_zero_enabled,
        config.viz_near_zero_threshold,
        config.viz_near_zero_n,
    )

    return entry, s
end

# ═══════════════════════════════════════════════════════════════════════════════
# Main — TOML mode (non-interactive) or interactive mode
# ═══════════════════════════════════════════════════════════════════════════════

if !isempty(ARGS) && endswith(ARGS[1], ".toml")
    # ── TOML mode: load config, evaluate, visualize, done ──
    using Globtim
    entry, settings = load_from_toml(ARGS[1])
    println()
    println("="^70)
    println("  Level Set Explorer — TOML: $(ARGS[1])")
    println("  Entry: $(entry.name)  p_true = $(entry.p_true)")
    println("="^70)
    print_settings(settings)
    println()

    all_points, all_values = evaluate_entry(entry, settings)
    screen, _ = show_visualization(entry, all_points, all_values, settings)

    println("\n  Visualization active. Close the window or press Ctrl-C to exit.")
    # Keep alive until window is closed
    try
        wait(screen)
    catch e
        e isa InterruptException || rethrow(e)
    end

else
    # ── Interactive mode: catalogue picker + menu loop ──
    catalogue_name = ""
    if !isempty(ARGS)
        cat_path = resolve_catalogue_path(ARGS[1])
        entries = load_catalogue(cat_path)
        catalogue_name =
            replace(basename(cat_path), "_catalogue.jsonl" => "", ".jsonl" => "")
    else
        entries, catalogue_name = select_catalogue_interactive()
    end

    println()
    println("="^70)
    println("  Level Set Explorer — $(catalogue_name) — $(length(entries)) entries loaded")
    println("="^70)
    print_entry_table(entries)
    println("  Current settings:")
    print_settings(SETTINGS)
    println()
    print_help()

    active_screen = nothing  # track open GLMakie window

    while true
        print("\nlevelset> ")
        input = strip(readline())
        isempty(input) && continue

        cmd = lowercase(input)

        if cmd in ("quit", "q", "exit")
            if active_screen !== nothing
                try
                    close(active_screen)
                catch
                end
            end
            println("Goodbye.")
            break

        elseif cmd in ("help", "?")
            print_help()

        elseif cmd == "list"
            print_entry_table(entries)

        elseif cmd == "settings"
            print_settings(SETTINGS)

        elseif cmd == "set"
            run_set_submenu!(SETTINGS)

        elseif cmd == "catalogue"
            new_entries, new_name = select_catalogue_interactive()
            global entries = new_entries
            global catalogue_name = new_name
            println("\n  Switched to $(catalogue_name) — $(length(entries)) entries")
            print_entry_table(entries)

        else
            # Try to parse as entry number
            idx = tryparse(Int, cmd)
            if idx !== nothing && 1 <= idx <= length(entries)
                local entry = entries[idx]

                # Close previous window
                if active_screen !== nothing
                    try
                        close(active_screen)
                    catch
                    end
                    global active_screen = nothing
                end

                # Evaluate and display
                local all_points, all_values = evaluate_entry(entry, SETTINGS)
                global active_screen
                active_screen, _ =
                    show_visualization(entry, all_points, all_values, SETTINGS)

                println("\n  Showing: $(entry.name)  p_true = $(entry.p_true)")
                println("  (type another number, 'set' to adjust, or 'quit')")
            else
                println("  Unknown command: '$input'. Type 'help' for options.")
            end
        end
    end
end
