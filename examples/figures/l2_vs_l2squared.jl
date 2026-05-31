# l2_vs_l2squared.jl
# Regenerate the L2 vs L2² illustrative figure for DynamicObjectives docs.
#
# Run with:
#   julia --project=profiles/viz pkg/DynamicObjectives/examples/figures/l2_vs_l2squared.jl
#
# Output: pkg/DynamicObjectives/docs/src/assets/l2_vs_l2squared.png
#
# Setup mirrors the LV2D_paper_1 catalogue entry
# (globtim_results/lv2d_catalogue.jsonl).

using DynamicObjectives
using DynamicObjectives: L2_norm, L2_squared, make_error_distance
using DynamicObjectives: define_lotka_volterra_2D_sciml_benchmark
using CairoMakie

const ASSET_DIR = abspath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets"))
isdir(ASSET_DIR) || mkpath(ASSET_DIR)

# We use the SciML-benchmark LV2D variant because it has TWO observable outputs
# (y1 = x1, y2 = x2). With both species observed the parameters (a1, a2) are
# globally identifiable, so the objective near p* is an isolated minimum
# rather than a valley — the contour panels then show the cone (L2) vs the
# smooth bowl (L2²) directly, not just on a 1-D slice.

println("Building L2 vs L2² figure (LV2D SciML benchmark, 2 outputs) …")

model, _params, _states, outputs = define_lotka_volterra_2D_sciml_benchmark()
ic            = [1.0, 1.0]
p_true        = [1.0, 1.5]
time_interval = [0.0, 5.0]
numpoints     = 25

obj_L2   = make_error_distance(model, outputs, ic, p_true, time_interval, numpoints, L2_norm)
obj_L2sq = make_error_distance(model, outputs, ic, p_true, time_interval, numpoints, L2_squared)

# Contour grid: tight window centred on p_true
δ      = 0.18
ngrid  = 80
a_axis = range(p_true[1] - δ, p_true[1] + δ; length = ngrid)
b_axis = range(p_true[2] - δ, p_true[2] + δ; length = ngrid)

println("  evaluating $(ngrid*ngrid) grid points × 2 metrics …")
Z_L2   = [obj_L2([a, b])   for a in a_axis, b in b_axis]
Z_L2sq = [obj_L2sq([a, b]) for a in a_axis, b in b_axis]
println("  L2  range: $(round(minimum(Z_L2);   sigdigits=3)) … $(round(maximum(Z_L2);   sigdigits=3))")
println("  L2² range: $(round(minimum(Z_L2sq); sigdigits=3)) … $(round(maximum(Z_L2sq); sigdigits=3))")

# 1D slice along a-axis (b fixed at b*) — widen to 1.5× δ so the V/U shape reads
δ_slice  = 1.5 * δ
nslice   = 200
a_slice  = collect(range(p_true[1] - δ_slice, p_true[1] + δ_slice; length = nslice))
y_L2     = [obj_L2([a, p_true[2]])   for a in a_slice]
y_L2sq   = [obj_L2sq([a, p_true[2]]) for a in a_slice]

# ─── Figure ──────────────────────────────────────────────────────────────────
fig = Figure(size = (1200, 760), fontsize = 14)

Label(fig[0, 1:4],
    "LV2D parameter-estimation objective near p* = $(p_true)   ·   SciML benchmark, 2 outputs (fully identifiable)";
    fontsize = 17, font = :bold, halign = :center)

# Top row: paired contour panels
ax_L2 = Axis(fig[1, 1];
    title  = "L₂  distance  ‖y_obs − y(p)‖₂   (cone at p*)",
    xlabel = "a₁", ylabel = "a₂",
    aspect = DataAspect(),
)
cf_L2 = contourf!(ax_L2, a_axis, b_axis, Z_L2;
    levels = 20, colormap = :viridis)
contour!(ax_L2, a_axis, b_axis, Z_L2;
    levels = 20, color = (:white, 0.25), linewidth = 0.5)
scatter!(ax_L2, [p_true[1]], [p_true[2]];
    marker = :star5, color = :red, markersize = 22,
    strokecolor = :black, strokewidth = 1, label = "p*")
Colorbar(fig[1, 2], cf_L2; label = "‖·‖₂", width = 14)

ax_L2sq = Axis(fig[1, 3];
    title  = "L₂²  distance  ‖y_obs − y(p)‖₂²   (smooth at p*)",
    xlabel = "a₁", ylabel = "a₂",
    aspect = DataAspect(),
)
cf_L2sq = contourf!(ax_L2sq, a_axis, b_axis, Z_L2sq;
    levels = 20, colormap = :viridis)
contour!(ax_L2sq, a_axis, b_axis, Z_L2sq;
    levels = 20, color = (:white, 0.25), linewidth = 0.5)
scatter!(ax_L2sq, [p_true[1]], [p_true[2]];
    marker = :star5, color = :red, markersize = 22,
    strokecolor = :black, strokewidth = 1, label = "p*")
Colorbar(fig[1, 4], cf_L2sq; label = "‖·‖₂²", width = 14)

# Bottom row: 1D slice spanning the full width
ax_slice = Axis(fig[2, 1:4];
    title  = "1-D slice through p*  (a₂ ≡ $(p_true[2]))",
    xlabel = "a₁", ylabel = "distance to data",
)
lines!(ax_slice, a_slice, y_L2;
    color = :crimson, linewidth = 2.8,
    label = "L₂   (Lipschitz cone — gradient discontinuous at p*)")
lines!(ax_slice, a_slice, y_L2sq;
    color = :steelblue, linewidth = 2.8,
    label = "L₂²  (C² quadratic — admits polynomial approximation)")
vlines!(ax_slice, [p_true[1]];
    color = :black, linestyle = :dash, linewidth = 1)
text!(ax_slice, p_true[1], 0.0;
    text = "  p*", align = (:left, :bottom),
    color = :black, fontsize = 12)
axislegend(ax_slice; position = :ct, framevisible = true,
    backgroundcolor = (:white, 0.85), labelsize = 12)

Label(fig[3, 1:4],
    "L₂ is Lipschitz at p* — gradient discontinuity violates the c ≥ max(3, βn+1) " *
    "smoothness assumption of Thm 1.\n" *
    "L₂² is C² at p* and shares the same minimiser — opt-in via " *
    "[model].distance_function_override = \"L2_squared\" in the TOML config.";
    fontsize = 11, color = :gray30, halign = :center)

rowgap!(fig.layout, 1, 4)
rowgap!(fig.layout, 2, 8)

out_path = joinpath(ASSET_DIR, "l2_vs_l2squared.png")
save(out_path, fig; px_per_unit = 2)
println("→ wrote $(out_path)")
