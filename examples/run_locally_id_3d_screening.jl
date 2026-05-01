# Run with: julia --project=profiles/dev --threads=auto pkg/Dynamic_objectives/examples/run_locally_id_3d_screening.jl
#
# Screen 3 locally-identifiable 3D ODE models designed for discrete finite minima.
#
# These models use squared parameters / symmetric polynomials to create sign-flip
# or permutation symmetry, guaranteeing multiple discrete equivalent global minima
# (unlike globally-identifiable models which have a single basin).
#
#   1. Goodwin3D-LI (k1², k2, k4²) — 4 discrete minima from ±k1, ±k4 sign flips
#   2. FHN3D-LI (g², a, b) — 2 discrete minima from ±g sign flip, dual output
#   3. LV3D-Sym (a+b, a·b, c) — 2 discrete minima from (a,b) ↔ (b,a) permutation
#
# Pipeline per model:
#   1. screen_and_probe(200 candidates, 10 probes) — discovery
#   2. score_landscape_grid(points_per_dim=15) — grid-based metrics on best candidate
#   3. Print comparison table with interestingness metrics
#   4. Save top entries to globtim_results/ catalogues

using Printf
using Dynamic_objectives

# ═══════════════════════════════════════════════════════════════════════════════
# Model specifications
# ═══════════════════════════════════════════════════════════════════════════════

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))

struct ModelSpec
    name::String
    define_fn::Function
    ic::Vector{Float64}
    bounds::Vector{Tuple{Float64,Float64}}
    tspan::Vector{Float64}
    solver::Any
    param_names::Vector{String}
end

const MODELS = [
    ModelSpec(
        "Goodwin3D_LI",
        define_goodwin_3d_locally_id_model,
        [0.5, 0.5, 0.5],
        [(-3.0, 3.0), (0.1, 2.0), (-2.0, 2.0)],  # k1 symmetric, k2 positive, k4 symmetric
        [0.0, 30.0],
        AutoTsit5(Rosenbrock23()),
        ["k1", "k2", "k4"],
    ),
    ModelSpec(
        "FHN3D_LI",
        define_fhn_3d_locally_id_model,
        [0.0, 0.0],
        [(-2.0, 2.0), (0.0, 2.0), (0.0, 2.0)],    # g symmetric, a/b positive
        [0.0, 20.0],
        AutoTsit5(Rosenbrock23()),
        ["g", "a", "b"],
    ),
    ModelSpec(
        "LV3D_Sym",
        define_lv_3d_symmetric_model,
        [0.3, 0.6],
        [(-1.0, 1.5), (-1.0, 1.5), (-0.5, 1.5)],   # a, b, c
        [0.0, 1.0],
        Tsit5(),
        ["a", "b", "c"],
    ),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Screening parameters
# ═══════════════════════════════════════════════════════════════════════════════

const N_CANDIDATES = 200
const N_PROBES = 10
const GRID_PPD = 15
const TOP_N_CATALOGUE = 5

# ═══════════════════════════════════════════════════════════════════════════════
# Run screening + grid scoring for each model
# ═══════════════════════════════════════════════════════════════════════════════

println("=" ^ 100)
println("  Locally-Identifiable 3D Models — Discrete Symmetry Minima")
println("  Screening: $N_CANDIDATES candidates, $N_PROBES probes per model")
println("  Grid scoring: $(GRID_PPD)^3 = $(GRID_PPD^3) points per model")
println("=" ^ 100)
println()

# Store results for comparison table
struct ModelResult
    spec::ModelSpec
    screen_result::ScreeningResult
    n_valid::Int
    pass_rate::Float64
    screen_time::Float64
    grid_score::Union{GridScoreResult,Nothing}
    grid_p_true::Union{Vector{Float64},Nothing}
end

results = ModelResult[]

for spec in MODELS
    println("─" ^ 100)
    println("  Model: $(spec.name)")
    @printf("  IC: %s   Bounds: %s\n", string(spec.ic), string(spec.bounds))
    @printf("  Tspan: %s   Params: %s\n", string(spec.tspan), join(spec.param_names, ", "))
    println("─" ^ 100)
    println()

    # Step 1: Screen
    model, params, states, outputs = spec.define_fn()

    println("  [1/2] Screening ($N_CANDIDATES candidates, $N_PROBES probes)...")
    screen_start = time()
    screen_result = screen_and_probe(
        model,
        outputs,
        spec.ic,
        spec.bounds,
        spec.tspan;
        n_candidates = N_CANDIDATES,
        n_probes = N_PROBES,
        solver = spec.solver,
        verbose = false,
    )
    screen_time = time() - screen_start

    n_valid = screen_result.sweep.diagnostics.n_valid
    pass_rate = screen_result.sweep.pass_rate
    ranked = screen_result.ranking.ranked_indices

    @printf(
        "  Pass rate: %d/%d (%.1f%%)  Ranked: %d  Time: %.1fs\n",
        n_valid,
        N_CANDIDATES,
        100.0 * pass_rate,
        length(ranked),
        screen_time
    )

    if n_valid == 0
        println("  ⚠ No valid candidates — skipping grid scoring")
        push!(
            results,
            ModelResult(spec, screen_result, 0, 0.0, screen_time, nothing, nothing),
        )
        println()
        continue
    end

    # Step 2: Grid score on best candidate
    best_idx = ranked[1]
    best_p = screen_result.sweep.valid[best_idx]
    @printf("  Best p_true: [%s]\n", join([@sprintf("%.4f", x) for x in best_p], ", "))

    println("  [2/2] Grid scoring ($(GRID_PPD)^3 = $(GRID_PPD^3) points)...")
    grid_start = time()
    gs = score_landscape_grid(
        model,
        outputs,
        spec.ic,
        best_p,
        spec.bounds,
        spec.tspan;
        points_per_dim = GRID_PPD,
        solver = spec.solver,
    )
    grid_time = time() - grid_start

    @printf(
        "  Grid time: %.1fs  Finite: %d/%d (%.1f%%)\n",
        grid_time,
        gs.n_finite,
        gs.n_total,
        100.0 * gs.n_finite / gs.n_total
    )
    @printf(
        "  Local minima: %d  Dynamic range: %.2f  Basin: %.1f%%  Deceptive: %d\n",
        gs.n_local_minima,
        gs.dynamic_range,
        100.0 * gs.basin_fraction,
        gs.n_deceptive
    )
    @printf(
        "  Dir. variation: [%s]  Curvature: %.3f\n",
        join([@sprintf("%.3f", v) for v in gs.directional_variation], ", "),
        gs.curvature_score
    )
    @printf(
        "  Basin depth: %.3f  Conditioning: %.1f  Ruggedness: %.3f  Plateau: %.2f\n",
        gs.basin_depth_ratio,
        gs.conditioning,
        gs.ruggedness,
        gs.plateau_fraction
    )
    @printf("  Interestingness: %.3f\n", interestingness_score(gs))
    println()

    push!(
        results,
        ModelResult(spec, screen_result, n_valid, pass_rate, screen_time, gs, best_p),
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Comparison table
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("=" ^ 100)
println("  COMPARISON TABLE — Locally Identifiable 3D Models")
println("=" ^ 100)
println()

@printf(
    "  %-14s  %9s  %8s  %10s  %8s  %8s  %10s  %10s  %7s  %6s  %6s  %6s  %12s\n",
    "Model",
    "pass_rate",
    "%finite",
    "dyn_range",
    "n_minima",
    "basin_%",
    "n_deceptive",
    "curvature",
    "depth",
    "cond",
    "rugg",
    "plat",
    "interesting."
)
println("  " * "-" ^ 138)

for r in results
    if r.grid_score === nothing
        @printf(
            "  %-14s  %8.1f%%  %8s  %10s  %8s  %8s  %10s  %10s  %7s  %6s  %6s  %6s  %12s\n",
            r.spec.name,
            100.0 * r.pass_rate,
            "—",
            "—",
            "—",
            "—",
            "—",
            "—",
            "—",
            "—",
            "—",
            "—",
            "—"
        )
    else
        gs = r.grid_score
        @printf(
            "  %-14s  %8.1f%%  %7.1f%%  %10.2f  %8d  %7.1f%%  %10d  %10.3f  %7.3f  %6.1f  %6.3f  %6.2f  %12.3f\n",
            r.spec.name,
            100.0 * r.pass_rate,
            100.0 * gs.n_finite / gs.n_total,
            gs.dynamic_range,
            gs.n_local_minima,
            100.0 * gs.basin_fraction,
            gs.n_deceptive,
            gs.curvature_score,
            gs.basin_depth_ratio,
            gs.conditioning,
            gs.ruggedness,
            gs.plateau_fraction,
            interestingness_score(gs)
        )
    end
end
println()

# Rank by interestingness
scored = [
    (i, interestingness_score(r.grid_score)) for
    (i, r) in enumerate(results) if r.grid_score !== nothing
]
sort!(scored; by = x -> -x[2])

if !isempty(scored)
    winner_idx = scored[1][1]
    winner = results[winner_idx]
    @printf("  Winner: %s (interestingness = %.3f)\n", winner.spec.name, scored[1][2])
    println()
end

# ═══════════════════════════════════════════════════════════════════════════════
# Save catalogue entries for all models with valid results
# ═══════════════════════════════════════════════════════════════════════════════

println("─" ^ 100)
println("  Saving catalogues")
println("─" ^ 100)

const CATALOGUE_NAMES = Dict(
    "Goodwin3D_LI" => "goodwin3d_locally_id",
    "FHN3D_LI" => "fhn3d_locally_id",
    "LV3D_Sym" => "lv3d_symmetric",
)

for r in results
    r.n_valid == 0 && continue

    spec = r.spec
    cat_name = CATALOGUE_NAMES[spec.name]
    cat_path = joinpath(REPO_ROOT, "globtim_results", "$(cat_name)_catalogue.jsonl")

    entries = screening_to_catalogue(
        r.screen_result,
        spec.define_fn,
        spec.ic,
        spec.bounds,
        spec.tspan;
        top_n = TOP_N_CATALOGUE,
        name_prefix = spec.name,
        description = "Locally-identifiable 3D model: $(spec.name) ($(join(spec.param_names, ", "))). " *
                      "Squared/symmetric parameterization for discrete minima.",
    )

    save_catalogue(cat_path, entries)
    @printf("  %s: %d entries → %s\n", spec.name, length(entries), cat_path)
end

println()
println("=" ^ 100)
println("  DONE — Run grid scoring to verify discrete minima structure:")
println(
    "  julia --project=profiles/dev --threads=auto experiments/sandbox/run_grid_scoring.jl",
)
println("=" ^ 100)
println()
