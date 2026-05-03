"""
TOML Screening Configuration — Candidate Discovery Pipeline

Parses TOML config files into `ScreeningConfig` and orchestrates the
screening pipeline: sweep → probe → rank → catalogue.

This is **Pipeline 1** (candidate discovery). Pipeline 2 (CP recovery via Globtim)
is handled by `run_experiment_from_config()` in the DynamicObjectivesGlobtimExt
package extension.

# Two-Pipeline Architecture
```
Pipeline 1: Candidate Discovery (this module)
  screening TOML → run_screening_from_config() → catalogue JSONL

Pipeline 2: CP Recovery (DynamicObjectivesGlobtimExt)
  experiment TOML → run_experiment_from_config() → Globtim results
```

The catalogue JSONL file is the handoff point between the two pipelines.
"""

# ═══════════════════════════════════════════════════════════════════════════════
# ScreeningConfig
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ScreeningConfig

Flat configuration struct parsed from a TOML screening file.

# Sections (matching TOML layout)
- `[screening]`: name, description
- `[model]`: model_fn (registry key), ic, bounds, time_interval
- `[sweep]`: n_candidates, margin, sampling, threshold, numpoints_screen
- `[probe]`: n_probes, numpoints_probe, distance_function
- `[solver]`: method, abstol, reltol
- `[ranking]`: min_finite_fraction, max_noise_ratio
- `[output]`: catalogue_path, name_prefix, top_n, description, numpoints,
              aggregate_distances, eval_timeout
"""
struct ScreeningConfig
    # [screening]
    name::String
    description::String

    # [model]
    model_fn::String                           # MODEL_REGISTRY key
    ic::Vector{Float64}
    bounds::Vector{Tuple{Float64,Float64}}
    time_interval::Vector{Float64}

    # [sweep]
    n_candidates::Int
    margin::Float64
    sampling::Symbol
    threshold::Float64
    numpoints_screen::Int

    # [probe]
    n_probes::Int
    numpoints_probe::Int
    distance_function::String                  # DISTANCE_REGISTRY key

    # [solver]
    solver_method::String
    solver_abstol::Float64
    solver_reltol::Float64

    # [ranking]
    min_finite_fraction::Float64
    max_noise_ratio::Float64

    # [output]
    catalogue_path::String                     # absolute (resolved at parse time)
    name_prefix::String
    top_n::Int
    catalogue_description::String
    catalogue_numpoints::Int
    aggregate_distances::String                # AGGREGATION_REGISTRY key
    eval_timeout::Union{Nothing,Float64}
end

# ═══════════════════════════════════════════════════════════════════════════════
# Path Resolution (same pattern as globtim/src/config_loader.jl)
# ═══════════════════════════════════════════════════════════════════════════════

const _SCREENING_RESULTS_PREFIX = "globtim_results/"  # always forward slash — TOML files use "/" on all platforms

"""
    _resolve_screening_path(p, config_dir) -> String

Resolve a relative path from a screening TOML config to an absolute path.

- Absolute paths pass through unchanged
- Paths starting with `"globtim_results/"` are resolved by walking up from
  `config_dir` to find the repo root (directory containing `globtim_results/`)
- Other relative paths resolve relative to `config_dir`
"""
function _resolve_screening_path(p::AbstractString, config_dir::AbstractString)
    isabspath(p) && return p

    if startswith(p, _SCREENING_RESULTS_PREFIX)
        # Walk up from config_dir to find repo root containing globtim_results/
        results_root = _find_results_root_from(config_dir)
        return joinpath(results_root, p[length(_SCREENING_RESULTS_PREFIX)+1:end])
    end

    return joinpath(config_dir, p)
end

"""
    _find_results_root_from(start_dir) -> String

Walk up from `start_dir` to find a directory containing `globtim_results/`.
Returns the absolute path to that `globtim_results/` directory.

Respects the `GLOBTIM_RESULTS_ROOT` environment variable (same as
PathManager.get_results_root() in Globtim) for consistency between
the screening and experiment pipelines.
"""
function _find_results_root_from(start_dir::AbstractString)::String
    # Respect GLOBTIM_RESULTS_ROOT if set (same precedence as PathManager)
    if haskey(ENV, "GLOBTIM_RESULTS_ROOT")
        root = abspath(ENV["GLOBTIM_RESULTS_ROOT"])
        isdir(root) ||
            error("GLOBTIM_RESULTS_ROOT is set but directory does not exist: $root")
        return root
    end

    # Walk up from start_dir to find repo root containing globtim_results/
    current = abspath(start_dir)
    max_iterations = 20
    for _ in 1:max_iterations
        candidate = joinpath(current, "globtim_results")
        if isdir(candidate)
            return candidate
        end
        parent = dirname(current)
        parent == current && break
        current = parent
    end
    error("""
        Cannot find globtim_results/ directory.
        Searched upward from: $start_dir

        Either:
        1. Run from within the repository, or
        2. Set the GLOBTIM_RESULTS_ROOT environment variable, or
        3. Use an absolute path for catalogue_path in your TOML config
        """)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    validate_screening_toml(d::Dict)

Validate that a parsed TOML dict has all required sections and keys
for a screening config.
"""
function validate_screening_toml(d::Dict)
    # Required sections
    for section in ("screening", "model", "output")
        haskey(d, section) || error("Missing required TOML section: [$section]")
    end

    # Required keys in [screening]
    scr = d["screening"]
    haskey(scr, "name") || error("[screening] missing required key: name")

    # Required keys in [model]
    mod = d["model"]
    for key in ("model_fn", "ic", "bounds", "time_interval")
        haskey(mod, key) || error("[model] missing required key: $key")
    end

    # Validate bounds is array of pairs
    bounds = mod["bounds"]
    bounds isa AbstractVector || error("[model] bounds must be an array of [lo, hi] pairs")
    for (i, pair) in enumerate(bounds)
        (pair isa AbstractVector && length(pair) == 2) ||
            error("[model] bounds[$i] must be a [lo, hi] pair, got: $pair")
    end

    # Required keys in [output]
    out = d["output"]
    haskey(out, "catalogue_path") || error("[output] missing required key: catalogue_path")

    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Parser
# ═══════════════════════════════════════════════════════════════════════════════

"""
    load_screening_config(path::String) -> ScreeningConfig

Parse a TOML screening configuration file and return a validated
`ScreeningConfig`.

# Example
```julia
config = load_screening_config("DynamicObjectives/examples/exp_candidates/lv4d.toml")
config.model_fn      # => "define_constrained_lotka_volterra_4D"
config.n_candidates   # => 200
config.catalogue_path # => "/abs/path/to/globtim_results/lv4d_catalogue.jsonl"
```
"""
function load_screening_config(path::String)
    isfile(path) || error("Screening config file not found: $path")

    d = TOML.parsefile(path)
    validate_screening_toml(d)

    scr = d["screening"]
    mod = d["model"]
    swp = get(d, "sweep", Dict())
    prb = get(d, "probe", Dict())
    sol = get(d, "solver", Dict())
    rnk = get(d, "ranking", Dict())
    out = d["output"]

    # Parse bounds: [[lo, hi], ...] → Vector{Tuple{Float64,Float64}}
    bounds = [
        Tuple{Float64,Float64}((Float64(pair[1]), Float64(pair[2]))) for
        pair in mod["bounds"]
    ]

    # Resolve catalogue_path
    config_dir = dirname(abspath(path))
    catalogue_path = _resolve_screening_path(String(out["catalogue_path"]), config_dir)

    return ScreeningConfig(
        # [screening]
        String(scr["name"]),
        String(get(scr, "description", "")),
        # [model]
        String(mod["model_fn"]),
        Float64.(mod["ic"]),
        bounds,
        Float64.(mod["time_interval"]),
        # [sweep]
        Int(get(swp, "n_candidates", 1000)),
        Float64(get(swp, "margin", 0.1)),
        Symbol(get(swp, "sampling", "random")),
        Float64(get(swp, "threshold", 1e6)),
        Int(get(swp, "numpoints_screen", 20)),
        # [probe]
        Int(get(prb, "n_probes", 10)),
        Int(get(prb, "numpoints_probe", 30)),
        String(get(prb, "distance_function", "L2_norm")),
        # [solver]
        String(get(sol, "method", "Tsit5")),
        Float64(get(sol, "abstol", 1e-4)),
        Float64(get(sol, "reltol", 1e-4)),
        # [ranking]
        Float64(get(rnk, "min_finite_fraction", 0.5)),
        Float64(get(rnk, "max_noise_ratio", Inf)),
        # [output]
        catalogue_path,
        String(get(out, "name_prefix", "screened")),
        Int(get(out, "top_n", 10)),
        String(get(out, "description", "Auto-discovered via run_screening_from_config")),
        Int(get(out, "numpoints", 30)),
        String(get(out, "aggregate_distances", "sum")),
        haskey(out, "eval_timeout") ? Float64(out["eval_timeout"]) : nothing,
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

"""
    run_screening_from_config(path::String; io::IO=stdout, verbose::Bool=true)
        -> Dict{Symbol, Any}

Run the full screening pipeline from a TOML configuration file.

Executes: model resolution → screen_and_probe → screening_to_catalogue → save_catalogue

No Globtim dependency required — this is purely DynamicObjectives.

# Returns
A Dict with keys:
- `:config` — the parsed `ScreeningConfig`
- `:screening_result` — the `ScreeningResult` from `screen_and_probe`
- `:entries` — `Vector{CatalogueEntry}` written to the catalogue
- `:catalogue_path` — absolute path to the saved catalogue JSONL

# Example
```julia
using DynamicObjectives

result = run_screening_from_config("DynamicObjectives/examples/exp_candidates/lv4d.toml")
result[:entries]         # Vector{CatalogueEntry}
result[:catalogue_path]  # "/abs/path/to/globtim_results/lv4d_catalogue.jsonl"
```
"""
function run_screening_from_config(path::String; io::IO = stdout, verbose::Bool = true)
    config = load_screening_config(path)

    println(io, "="^100)
    println(io, "  Screening: $(config.name)")
    config.description != "" && println(io, "  $(config.description)")
    println(io, "="^100)
    println(io)

    # ── 1. Resolve model from MODEL_REGISTRY ──
    model_fn = _resolve_function(config.model_fn, MODEL_REGISTRY, "model")
    model, _params, _states, outputs = model_fn()
    println(io, "  Model: $(config.model_fn)")
    println(io, "  IC: $(config.ic)")
    println(io, "  Bounds: $(config.bounds)")
    println(io, "  Time interval: $(config.time_interval)")
    println(io)

    # ── 2. Resolve solver ──
    solver = _resolve_solver(config.solver_method)
    println(
        io,
        "  Solver: $(config.solver_method) (abstol=$(config.solver_abstol), reltol=$(config.solver_reltol))",
    )

    # ── 3. Resolve distance function ──
    distance_fn = _resolve_function(config.distance_function, DISTANCE_REGISTRY, "distance")

    # ── 4. Run screen_and_probe ──
    println(io)
    println(io, "─"^100)
    println(
        io,
        "  screen_and_probe: $(config.n_candidates) candidates, $(config.n_probes) probes",
    )
    println(io, "─"^100)
    println(io)

    screen_start = time()
    screening_result = screen_and_probe(
        model,
        outputs,
        config.ic,
        config.bounds,
        config.time_interval;
        n_candidates = config.n_candidates,
        n_probes = config.n_probes,
        margin = config.margin,
        sampling = config.sampling,
        threshold = config.threshold,
        numpoints_screen = config.numpoints_screen,
        numpoints_probe = config.numpoints_probe,
        distance_function = distance_fn,
        solver = solver,
        abstol = config.solver_abstol,
        reltol = config.solver_reltol,
        min_finite_fraction = config.min_finite_fraction,
        max_noise_ratio = config.max_noise_ratio,
        verbose = verbose,
    )
    screen_time = time() - screen_start

    sweep = screening_result.sweep
    ranked = screening_result.ranking.ranked_indices

    println(io)
    println(io, "─"^100)
    println(
        io,
        "  Screening complete: $(length(ranked)) ranked / $(length(sweep.valid)) valid / $(sweep.diagnostics.n_total) total",
    )
    @printf(io, "  Pass rate: %.1f%%  |  Time: %.1fs\n", sweep.pass_rate * 100, screen_time)
    println(io, "─"^100)

    if isempty(ranked)
        println(io)
        println(io, "  WARNING: No viable candidates found. Try:")
        println(io, "    - Widening [model] bounds")
        println(io, "    - Increasing [sweep] n_candidates")
        println(io, "    - Lowering [ranking] min_finite_fraction")
        println(io)
        return Dict{Symbol,Any}(
            :config => config,
            :screening_result => screening_result,
            :entries => CatalogueEntry[],
            :catalogue_path => config.catalogue_path,
        )
    end

    # ── 5. Print top candidates ──
    println(io)
    n_show = min(config.top_n, length(ranked))
    println(io, "  Top $n_show candidates:")
    for i in 1:n_show
        idx = ranked[i]
        p = sweep.valid[idx]
        probe = screening_result.probes[idx]
        @printf(
            io,
            "    %d. p_true = [%s]  dynamic_range=%.2f  obj@true=%.2e\n",
            i,
            join([@sprintf("%.6f", x) for x in p], ", "),
            probe.dynamic_range,
            probe.objective_at_true
        )
    end
    println(io)

    # ── 6. Convert to catalogue entries ──
    aggregate_fn =
        _resolve_function(config.aggregate_distances, AGGREGATION_REGISTRY, "aggregation")

    entries = screening_to_catalogue(
        screening_result,
        model_fn,
        config.ic,
        config.bounds,
        config.time_interval;
        top_n = config.top_n,
        name_prefix = config.name_prefix,
        description = config.catalogue_description,
        numpoints = config.catalogue_numpoints,
        distance_function = distance_fn,
        aggregate_distances = aggregate_fn,
        eval_timeout = config.eval_timeout,
    )

    # ── 7. Save catalogue ──
    save_catalogue(config.catalogue_path, entries)
    println(io, "  Saved $(length(entries)) entries to: $(config.catalogue_path)")
    println(io)

    return Dict{Symbol,Any}(
        :config => config,
        :screening_result => screening_result,
        :entries => entries,
        :catalogue_path => config.catalogue_path,
    )
end
