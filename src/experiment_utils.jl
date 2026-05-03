"""
Shared experiment utilities for multi-candidate Globtim validation scripts.

Extracted from duplicated code in run_lv4d_globtim.jl and run_daisy4d_multi_ptrue.jl.
"""

"""
    CandidateResult

Per-candidate, per-degree result from a Globtim validation run.
"""
struct CandidateResult
    name::String
    p_true::Vector{Float64}
    degree::Int
    n_cps::Int
    raw_obj::Float64
    raw_recovery::Float64
    l2_error::Float64
    rel_l2::Float64
    capture_1pct::Float64
    capture_5pct::Float64
end

"""
    build_bounds(p_true, radius) -> Vector{Tuple{Float64,Float64}}

Build symmetric parameter bounds centered on `p_true` with half-width `radius`.
"""
function build_bounds(p_true, radius::Real)
    return [(Float64(p - radius), Float64(p + radius)) for p in p_true]
end

"""
    build_bounds(p_true, radii::AbstractVector{<:Real}) -> Vector{Tuple{Float64,Float64}}

Build anisotropic parameter bounds centered on `p_true` with per-dimension half-widths.

# Example
```julia
build_bounds([0.1, 0.2, 0.3], [0.01, 0.05, 0.02])
# => [(0.09, 0.11), (0.15, 0.25), (0.28, 0.32)]
```
"""
function build_bounds(p_true, radii::AbstractVector{<:Real})
    length(radii) == length(p_true) || error(
        "radii length ($(length(radii))) must match p_true length ($(length(p_true)))",
    )
    return [(Float64(p - r), Float64(p + r)) for (p, r) in zip(p_true, radii)]
end

"""
    print_candidate_summary_table(all_results, entries; title, radius, degree_range, name_width=20, p_true_precision=4)

Print a formatted summary table of Globtim candidate results.

Prints three sections:
1. Full per-candidate, per-degree table
2. Best per candidate (lowest raw recovery error)
3. Aggregate statistics

Returns the `best_per_candidate` vector for further use.

# Arguments
- `all_results::Vector{CandidateResult}`: All results across candidates and degrees
- `entries`: Vector of candidate entries (must have `.name` field)
- `title::String`: Title for the summary header
- `radius::Real`: Domain radius used
- `degree_range`: Degree range used (for header display)
- `name_width::Int=20`: Column width for candidate name
- `p_true_precision::Int=4`: Decimal places for p_true display
"""
function print_candidate_summary_table(
    all_results::Vector{CandidateResult},
    entries;
    title::String,
    radius::Real,
    degree_range,
    name_width::Int = 20,
    p_true_precision::Int = 4,
)
    _pad(s, w) = rpad(s, w)

    println()
    println("="^120)
    println("  $title")
    @printf(
        "  Radius: %.4f, %d candidates, degrees %s\n",
        radius,
        length(entries),
        string(collect(degree_range))
    )
    println("="^120)
    println()

    # ── Full per-candidate, per-degree table ──
    println(
        rpad("#", 4),
        "  ",
        _pad("Name", name_width),
        "  ",
        lpad("Deg", 4),
        "  ",
        lpad("#CPs", 5),
        "  ",
        lpad("Raw f(x)", 11),
        "  ",
        lpad("Raw Recov", 12),
        "  ",
        lpad("L2 Error", 10),
        "  ",
        lpad("Rel L2", 10),
        "  ",
        lpad("Cap@1%", 8),
        "  ",
        lpad("Cap@5%", 8),
    )
    println("-"^120)

    for r in all_results
        raw_str = isnan(r.raw_obj) ? "N/A" : @sprintf("%.3e", r.raw_obj)
        rec_str = isnan(r.raw_recovery) ? "N/A" : @sprintf("%.3e", r.raw_recovery)
        rel_str = isnan(r.rel_l2) ? "N/A" : @sprintf("%.3e", r.rel_l2)
        c1_str = isnan(r.capture_1pct) ? "---" : @sprintf("%.0f%%", 100 * r.capture_1pct)
        c5_str = isnan(r.capture_5pct) ? "---" : @sprintf("%.0f%%", 100 * r.capture_5pct)

        ci = findfirst(e -> e.name == r.name, entries)
        println(
            rpad(string(ci), 4),
            "  ",
            _pad(r.name, name_width),
            "  ",
            lpad(string(r.degree), 4),
            "  ",
            lpad(string(r.n_cps), 5),
            "  ",
            lpad(raw_str, 11),
            "  ",
            lpad(rec_str, 12),
            "  ",
            lpad(@sprintf("%.3e", r.l2_error), 10),
            "  ",
            lpad(rel_str, 10),
            "  ",
            lpad(c1_str, 8),
            "  ",
            lpad(c5_str, 8),
        )
    end

    # ── Best per candidate ──
    println()
    println("─"^120)
    println("  Best per candidate (lowest raw recovery error across degrees):")
    println("─"^120)
    println(
        rpad("#", 4),
        "  ",
        _pad("Name", name_width),
        "  ",
        lpad("Deg", 4),
        "  ",
        lpad("Raw Recov", 12),
        "  ",
        lpad("Raw f(x)", 12),
        "  ",
        lpad("Cap@1%", 8),
        "  ",
        lpad("Cap@5%", 8),
        "  ",
        "p_true",
    )
    println("-"^120)

    best_per_candidate = CandidateResult[]
    p_fmt = "%.$(p_true_precision)f"
    for (ci, entry) in enumerate(entries)
        rows = filter(r -> r.name == entry.name && !isnan(r.raw_recovery), all_results)
        if isempty(rows)
            println(
                rpad(string(ci), 4),
                "  ",
                _pad(entry.name, name_width),
                "  ",
                lpad("---", 4),
                "  ",
                lpad("---", 12),
                "  ",
                lpad("---", 12),
                "  ",
                lpad("---", 8),
                "  ",
                lpad("---", 8),
            )
            continue
        end
        best = rows[argmin([r.raw_recovery for r in rows])]
        push!(best_per_candidate, best)
        raw_str = isnan(best.raw_obj) ? "N/A" : @sprintf("%.3e", best.raw_obj)
        rec_str = @sprintf("%.3e", best.raw_recovery)
        c1_str =
            isnan(best.capture_1pct) ? "---" : @sprintf("%.0f%%", 100 * best.capture_1pct)
        c5_str =
            isnan(best.capture_5pct) ? "---" : @sprintf("%.0f%%", 100 * best.capture_5pct)
        p_str = "[" * join([@sprintf("%.4f", x) for x in best.p_true], ", ") * "]"
        println(
            rpad(string(ci), 4),
            "  ",
            _pad(best.name, name_width),
            "  ",
            lpad(string(best.degree), 4),
            "  ",
            lpad(rec_str, 12),
            "  ",
            lpad(raw_str, 12),
            "  ",
            lpad(c1_str, 8),
            "  ",
            lpad(c5_str, 8),
            "  ",
            p_str,
        )
    end

    # ── Aggregate statistics ──
    if !isempty(best_per_candidate)
        recoveries = [r.raw_recovery for r in best_per_candidate if !isnan(r.raw_recovery)]
        if !isempty(recoveries)
            println()
            println("─"^120)
            println("  Aggregate Statistics (best degree per candidate):")
            println("─"^120)
            @printf("  Candidates tested:   %d\n", length(entries))
            @printf("  With valid recovery: %d\n", length(recoveries))
            @printf("  Mean raw recovery:   %.3e\n", sum(recoveries) / length(recoveries))
            @printf(
                "  Median raw recovery: %.3e\n",
                sort(recoveries)[div(length(recoveries) + 1, 2)]
            )
            @printf("  Min raw recovery:    %.3e\n", minimum(recoveries))
            @printf("  Max raw recovery:    %.3e\n", maximum(recoveries))

            cap5_rates =
                [r.capture_5pct for r in best_per_candidate if !isnan(r.capture_5pct)]
            if !isempty(cap5_rates)
                n_full_cap = count(r -> r == 1.0, cap5_rates)
                @printf(
                    "  Cap@5%% = 100%%:       %d/%d candidates\n",
                    n_full_cap,
                    length(cap5_rates)
                )
            end

            cap1_rates =
                [r.capture_1pct for r in best_per_candidate if !isnan(r.capture_1pct)]
            if !isempty(cap1_rates)
                n_full_cap1 = count(r -> r == 1.0, cap1_rates)
                @printf(
                    "  Cap@1%% = 100%%:       %d/%d candidates\n",
                    n_full_cap1,
                    length(cap1_rates)
                )
            end
        end
    end

    println()
    return best_per_candidate
end
