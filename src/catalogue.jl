"""
    Catalogue — Typed entries for ODE model benchmark configurations

Provides `CatalogueEntry` struct with canonical field names for defining
ODE parameter estimation benchmarks, `create_objective` for converting
a catalogue entry into an optimization-ready error function, JSONL
persistence for saving/loading discovered configurations, and a bridge
from `ScreeningResult` to `CatalogueEntry` for the discovery pipeline.

Created: 2026-02-06 (Phase 7 of API naming standardization)
Extended: 2026-02-06 (JSONL persistence)
Extended: 2026-02-06 (screening → catalogue bridge)
"""

# ═══════════════════════════════════════════════════════════════════════════════
# Function Registries — bidirectional name ↔ function mapping for serialization
# ═══════════════════════════════════════════════════════════════════════════════

"""Global registry mapping string names to model-definition functions."""
const MODEL_REGISTRY = Dict{String,Function}()

"""Global registry mapping string names to distance functions."""
const DISTANCE_REGISTRY = Dict{String,Function}()

"""Global registry mapping string names to aggregation functions."""
const AGGREGATION_REGISTRY = Dict{String,Function}()

"""
    register_model!(name::String, fn::Function)

Register a model-definition function under `name` for catalogue serialization.
The function must accept no arguments and return `(model, params, states, outputs)`.

# Example
```julia
register_model!("define_lotka_volterra_2D_model_v3", define_lotka_volterra_2D_model_v3)
```
"""
function register_model!(name::String, fn::Function)
    MODEL_REGISTRY[name] = fn
    return nothing
end

"""
    register_distance!(name::String, fn::Function)

Register a distance function under `name` for catalogue serialization.

# Example
```julia
register_distance!("L2_norm", L2_norm)
```
"""
function register_distance!(name::String, fn::Function)
    DISTANCE_REGISTRY[name] = fn
    return nothing
end

"""
    register_aggregation!(name::String, fn::Function)

Register an aggregation function under `name` for catalogue serialization.

# Example
```julia
register_aggregation!("sum", sum)
```
"""
function register_aggregation!(name::String, fn::Function)
    AGGREGATION_REGISTRY[name] = fn
    return nothing
end

"""
    _resolve_function(name::String, registry::Dict{String,Function}, kind::String) -> Function

Look up a function by name in the given registry. Raises an error if not found,
listing all registered names for debugging.
"""
function _resolve_function(name::String, registry::Dict{String,Function}, kind::String)
    haskey(registry, name) && return registry[name]
    registered = sort(collect(keys(registry)))
    error("Unknown $kind function: \"$name\". Registered: $(join(registered, ", "))")
end

"""
    _function_name(fn::Function, registry::Dict{String,Function}, kind::String) -> String

Reverse-lookup: find the registered name for a function. Raises an error if
the function is not registered.
"""
function _function_name(fn::Function, registry::Dict{String,Function}, kind::String)
    for (name, registered_fn) in registry
        registered_fn === fn && return name
    end
    error(
        "$kind function is not registered for serialization. " *
        "Register it with register_$(lowercase(kind))!(name, fn) before saving.",
    )
end

"""
    _init_registries!()

Populate the default registries with all known model, distance, and aggregation
functions. Called once during module initialization.
"""
function _init_registries!()
    # Distance functions
    register_distance!("L1_norm", L1_norm)
    register_distance!("L2_norm", L2_norm)
    register_distance!("L2_squared", L2_squared)
    register_distance!("log_L2_norm", log_L2_norm)

    # Aggregation functions
    register_aggregation!("sum", sum)
    register_aggregation!("maximum", maximum)
    register_aggregation!("mean", mean)
    register_aggregation!("minimum", minimum)

    # Model functions — all define_* from systems/
    register_model!("define_daisy_ex3_model_4D", define_daisy_ex3_model_4D)
    register_model!("define_daisy_ex3_model_5D", define_daisy_ex3_model_5D)
    register_model!(
        "define_daisy_ex3_model_4D_no_input",
        define_daisy_ex3_model_4D_no_input,
    )
    register_model!(
        "define_generalized_lotka_volterra_4D",
        define_generalized_lotka_volterra_4D,
    )
    register_model!(
        "define_constrained_lotka_volterra_4D",
        define_constrained_lotka_volterra_4D,
    )
    register_model!("define_lotka_volterra_4D_simple", define_lotka_volterra_4D_simple)
    register_model!("define_lotka_volterra_3D_model", define_lotka_volterra_3D_model)
    register_model!(
        "define_lotka_volterra_3D_model_locally_identifiable",
        define_lotka_volterra_3D_model_locally_identifiable,
    )
    register_model!("define_lotka_volterra_3D_model_v2", define_lotka_volterra_3D_model_v2)
    register_model!("define_lotka_volterra_2D_model", define_lotka_volterra_2D_model)
    register_model!("define_lotka_volterra_2D_model_v2", define_lotka_volterra_2D_model_v2)
    register_model!("define_lotka_volterra_2D_model_v3", define_lotka_volterra_2D_model_v3)
    register_model!(
        "define_lotka_volterra_2D_model_v3_two_outputs",
        define_lotka_volterra_2D_model_v3_two_outputs,
    )
    register_model!(
        "define_lotka_volterra_2D_sciml_benchmark",
        define_lotka_volterra_2D_sciml_benchmark,
    )
    register_model!("define_fitzhugh_nagumo_3D_model", define_fitzhugh_nagumo_3D_model)
    register_model!(
        "define_fitzhugh_nagumo_3D_model_two_outputs",
        define_fitzhugh_nagumo_3D_model_two_outputs,
    )
    register_model!("define_goodwin_oscillator_4D", define_goodwin_oscillator_4D)
    register_model!(
        "define_simple_2D_model_locally_identifiable",
        define_simple_2D_model_locally_identifiable,
    )
    register_model!(
        "define_simple_2D_model_locally_identifiable_square",
        define_simple_2D_model_locally_identifiable_square,
    )
    register_model!(
        "define_simple_1D_model_locally_identifiable",
        define_simple_1D_model_locally_identifiable,
    )
    register_model!("define_lorenz_3D_model", define_lorenz_3D_model)
    register_model!("define_rossler_3D_model", define_rossler_3D_model)
    register_model!("define_goodwin_oscillator_3D", define_goodwin_oscillator_3D)
    register_model!("define_coupled_lv2d_3d_model", define_coupled_lv2d_3d_model)
    register_model!("define_lv2d_reparam_3d_model", define_lv2d_reparam_3d_model)
    register_model!(
        "define_rosenzweig_macarthur_3d_model",
        define_rosenzweig_macarthur_3d_model,
    )
    register_model!(
        "define_goodwin_3d_product_obs_model",
        define_goodwin_3d_product_obs_model,
    )
    register_model!(
        "define_goodwin_3d_locally_id_model",
        define_goodwin_3d_locally_id_model,
    )
    register_model!("define_fhn_3d_locally_id_model", define_fhn_3d_locally_id_model)
    register_model!("define_lv_3d_symmetric_model", define_lv_3d_symmetric_model)

    # New benchmark models (epidemiology, neuroscience, PK, chemistry, biochemistry)
    register_model!("define_sir_2d_model", define_sir_2d_model)
    register_model!("define_seir_3d_model", define_seir_3d_model)
    register_model!("define_hindmarsh_rose_3d_model", define_hindmarsh_rose_3d_model)
    register_model!("define_pk_2comp_3d_model", define_pk_2comp_3d_model)
    register_model!("define_brusselator_2d_model", define_brusselator_2d_model)
    register_model!("define_michaelis_menten_2d_model", define_michaelis_menten_2d_model)

    # 4-D trophic model (cbyn.4 — 3-species RMA food chain)
    register_model!(
        "define_rosenzweig_macarthur_4d_model",
        define_rosenzweig_macarthur_4d_model,
    )

    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# CatalogueEntry struct
# ═══════════════════════════════════════════════════════════════════════════════

"""
    CatalogueEntry

A typed, immutable description of an ODE parameter estimation benchmark.

# Fields
- `name::String`: Short identifier (e.g., "LV_2D_v3")
- `description::String`: Human-readable description
- `model_fn::Function`: Returns `(model, params, states, outputs)` when called with no arguments
- `p_true::Vector{Float64}`: True parameter values (generates reference data)
- `ic::Vector{Float64}`: Initial conditions for ODE integration
- `bounds::Vector{Tuple{Float64,Float64}}`: Parameter search bounds as (lo, hi) tuples
- `time_interval::Vector{Float64}`: Integration time span `[t_start, t_end]`
- `numpoints::Int`: Number of time points to sample
- `distance_function::Function`: Per-output distance metric (e.g., `L2_norm`)
- `aggregate_distances::Function`: How to combine multi-output errors (e.g., `sum`)
- `eval_timeout::Union{Nothing,Float64}`: Timeout in seconds per evaluation (nothing = no timeout)

# Example
```julia
entry = CatalogueEntry(
    name = "LV_2D_v3",
    description = "2D Lotka-Volterra with c=0.5",
    model_fn = define_lotka_volterra_2D_model_v3,
    p_true = [1.0, 0.5],
    ic = [1.0, 0.5],
    bounds = [(0.0, 3.0), (0.0, 2.0)],
    time_interval = [0.0, 20.0],
    numpoints = 30,
    distance_function = L2_norm,
)

# Create objective function
objective = create_objective(entry)
objective([1.1, 0.6])  # returns Float64
```
"""
Base.@kwdef struct CatalogueEntry
    name::String
    description::String
    model_fn::Function
    p_true::Vector{Float64}
    ic::Vector{Float64}
    bounds::Vector{Tuple{Float64,Float64}}
    time_interval::Vector{Float64} = [0.0, 20.0]
    numpoints::Int = 30
    distance_function::Function = L2_norm
    aggregate_distances::Function = sum
    eval_timeout::Union{Nothing,Float64} = nothing
    refinement_distance_function::Union{Nothing,Function} = nothing
end

"""
    create_objective(entry::CatalogueEntry) -> Function

Create an optimization-ready error function from a catalogue entry.

Returns a single-argument function `f(p::Vector{Float64}) -> Float64`
suitable for use with globtim's `run_standard_experiment`.

# Example
```julia
entry = CatalogueEntry(name="LV_2D_v3", ...)
objective = create_objective(entry)
objective([1.0, 0.5])  # → 0.0 (at true parameters)
```
"""
function create_objective(entry::CatalogueEntry)
    model, _, _, outputs = entry.model_fn()

    return make_error_distance(
        model,
        outputs,
        entry.ic,
        entry.p_true,
        entry.time_interval,
        entry.numpoints,
        entry.distance_function,
        entry.aggregate_distances;
        return_inf_on_error = true,
        eval_timeout = entry.eval_timeout,
    )
end

"""
    create_tolerant_objective(entry::CatalogueEntry; abstol=1e-10, reltol=1e-10, solver=Vern9()) -> TolerantObjective

Create a `TolerantObjective` from a catalogue entry. Unlike `create_objective` (which
returns a plain closure), the returned object supports `set_tolerance!` for switching
solver tolerances between phases (e.g., coarse grid evaluation vs fine refinement).

# Example
```julia
entry = load_catalogue("data/catalogue.jsonl")[1]
obj = create_tolerant_objective(entry; abstol=1e-4, reltol=1e-4)
obj([1.0, 0.5])          # evaluate
set_tolerance!(obj, 1e-8) # switch to tight tolerances
```
"""
function create_tolerant_objective(
    entry::CatalogueEntry;
    abstol::Real = 1e-10,
    reltol::Real = 1e-10,
    solver = Vern9(),
)
    model, _, _, outputs = entry.model_fn()

    return TolerantObjective(
        model,
        outputs,
        entry.ic,
        entry.p_true,
        entry.time_interval,
        entry.numpoints,
        entry.distance_function,
        entry.aggregate_distances;
        return_inf_on_error = true,
        eval_timeout = entry.eval_timeout,
        solver = solver,
        abstol = abstol,
        reltol = reltol,
    )
end

"""
    create_refinement_objective(entry::CatalogueEntry) -> Function

Create a refinement-specific objective from a catalogue entry. If the entry has a
`refinement_distance_function` set (e.g. raw `L2_norm` instead of `log_L2_norm`),
uses that for the objective. Otherwise falls back to `create_objective`.

This separation allows display objectives (e.g. `log_L2_norm` for visualization)
to differ from refinement objectives (e.g. `L2_norm` for Newton convergence).
"""
function create_refinement_objective(entry::CatalogueEntry)
    dist_fn = entry.refinement_distance_function
    if dist_fn === nothing
        return create_objective(entry)
    end
    model, _, _, outputs = entry.model_fn()
    return make_error_distance(
        model,
        outputs,
        entry.ic,
        entry.p_true,
        entry.time_interval,
        entry.numpoints,
        dist_fn,
        entry.aggregate_distances;
        return_inf_on_error = true,
        eval_timeout = entry.eval_timeout,
    )
end

"""
    _non_log_distance(fn::Function) -> Union{Nothing, Function}

If `fn` is a registered `log_*` distance function, return the corresponding
non-log base function. Returns `nothing` if `fn` is not a log variant or
the base function is not registered.
"""
function _non_log_distance(fn::Function)::Union{Nothing,Function}
    for (name, registered_fn) in DISTANCE_REGISTRY
        registered_fn === fn || continue
        startswith(name, "log_") || return nothing
        base_name = name[5:end]
        haskey(DISTANCE_REGISTRY, base_name) || return nothing
        return DISTANCE_REGISTRY[base_name]
    end
    return nothing
end

"""
    dimension(entry::CatalogueEntry) -> Int

Number of parameters to estimate (length of `p_true`).
"""
dimension(entry::CatalogueEntry) = length(entry.p_true)

# ═══════════════════════════════════════════════════════════════════════════════
# JSONL Serialization — convert CatalogueEntry ↔ Dict for JSON persistence
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _entry_to_dict(entry::CatalogueEntry) -> Dict{String,Any}

Convert a `CatalogueEntry` to a JSON-serializable `Dict`. Function fields
are stored as their registered string names.
"""
function _entry_to_dict(entry::CatalogueEntry)
    d = Dict{String,Any}(
        "name" => entry.name,
        "description" => entry.description,
        "model_fn" => _function_name(entry.model_fn, MODEL_REGISTRY, "model"),
        "p_true" => entry.p_true,
        "ic" => entry.ic,
        "bounds" => [[lo, hi] for (lo, hi) in entry.bounds],
        "time_interval" => entry.time_interval,
        "numpoints" => entry.numpoints,
        "distance_function" =>
            _function_name(entry.distance_function, DISTANCE_REGISTRY, "distance"),
        "aggregate_distances" => _function_name(
            entry.aggregate_distances,
            AGGREGATION_REGISTRY,
            "aggregation",
        ),
        "eval_timeout" => entry.eval_timeout,
    )
    if entry.refinement_distance_function !== nothing
        d["refinement_distance_function"] = _function_name(
            entry.refinement_distance_function,
            DISTANCE_REGISTRY,
            "distance",
        )
    end
    return d
end

"""
    _dict_to_entry(d::Dict) -> CatalogueEntry

Reconstruct a `CatalogueEntry` from a parsed JSON `Dict`. Function fields
are resolved from their registered string names.

Throws an error if any function name is not found in the registries.
"""
function _dict_to_entry(d::Dict)
    dist_fn =
        _resolve_function(String(d["distance_function"]), DISTANCE_REGISTRY, "distance")

    # Resolve refinement_distance_function: explicit in dict → auto-detect from log_ prefix → nothing
    ref_dist =
        if haskey(d, "refinement_distance_function") &&
           d["refinement_distance_function"] !== nothing
            _resolve_function(
                String(d["refinement_distance_function"]),
                DISTANCE_REGISTRY,
                "distance",
            )
        else
            _non_log_distance(dist_fn)
        end

    return CatalogueEntry(
        name = String(d["name"]),
        description = String(d["description"]),
        model_fn = _resolve_function(String(d["model_fn"]), MODEL_REGISTRY, "model"),
        p_true = Float64.(d["p_true"]),
        ic = Float64.(d["ic"]),
        bounds = [
            Tuple{Float64,Float64}((Float64(pair[1]), Float64(pair[2]))) for
            pair in d["bounds"]
        ],
        time_interval = Float64.(d["time_interval"]),
        numpoints = Int(d["numpoints"]),
        distance_function = dist_fn,
        aggregate_distances = _resolve_function(
            String(d["aggregate_distances"]),
            AGGREGATION_REGISTRY,
            "aggregation",
        ),
        eval_timeout = d["eval_timeout"] === nothing ? nothing : Float64(d["eval_timeout"]),
        refinement_distance_function = ref_dist,
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# JSONL Persistence — save / append / load / summarize catalogue files
# ═══════════════════════════════════════════════════════════════════════════════

"""
    save_catalogue(path::AbstractString, entries::Vector{CatalogueEntry})

Write all entries to a JSONL file, **overwriting** any existing content.
Each line is one JSON object representing a `CatalogueEntry`.

Creates parent directories if they don't exist.

# Example
```julia
save_catalogue("data/catalogue.jsonl", [entry1, entry2, entry3])
```
"""
function save_catalogue(path::AbstractString, entries::Vector{CatalogueEntry})
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        for entry in entries
            d = _entry_to_dict(entry)
            JSON3.write(io, d)
            write(io, '\n')
        end
    end
    return nothing
end

"""
    append_catalogue(path::AbstractString, entry::CatalogueEntry)

Append a single entry to a JSONL catalogue file. Creates the file (and
parent directories) if it doesn't exist.

# Example
```julia
append_catalogue("data/catalogue.jsonl", new_entry)
```
"""
function append_catalogue(path::AbstractString, entry::CatalogueEntry)
    mkpath(dirname(abspath(path)))
    open(path, "a") do io
        d = _entry_to_dict(entry)
        JSON3.write(io, d)
        write(io, '\n')
    end
    return nothing
end

"""
    load_catalogue(path::AbstractString; model_name=nothing, max_entries=nothing) -> Vector{CatalogueEntry}

Load entries from a JSONL catalogue file. Optionally filter by model name
and/or limit the number of entries returned.

Throws an error if the file does not exist.

# Keyword Arguments
- `model_name::Union{Nothing,String}`: If set, only return entries whose `name`
  starts with this prefix (case-sensitive).
- `max_entries::Union{Nothing,Int}`: If set, return at most this many entries.

# Example
```julia
# Load all entries
entries = load_catalogue("data/catalogue.jsonl")

# Load only LV_2D entries
lv2d = load_catalogue("data/catalogue.jsonl", model_name="LV_2D")

# Load first 10
top10 = load_catalogue("data/catalogue.jsonl", max_entries=10)
```
"""
function load_catalogue(
    path::AbstractString;
    model_name::Union{Nothing,String} = nothing,
    max_entries::Union{Nothing,Int} = nothing,
)
    isfile(path) || error("Catalogue file not found: $path")

    entries = CatalogueEntry[]
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        d = JSON3.read(stripped, Dict{String,Any})
        entry = _dict_to_entry(d)

        # Apply model_name filter
        if model_name !== nothing && !startswith(entry.name, model_name)
            continue
        end

        push!(entries, entry)

        # Apply max_entries limit
        if max_entries !== nothing && length(entries) >= max_entries
            break
        end
    end

    return entries
end

"""
    catalogue_summary(path::AbstractString) -> NamedTuple

Compute summary statistics for a JSONL catalogue file.

# Returns
A NamedTuple with fields:
- `total::Int`: Total number of entries
- `per_model::Dict{String,Int}`: Count of entries per model name
- `dimensions::Dict{Int,Int}`: Count of entries per parameter dimension
- `models::Vector{String}`: Sorted list of unique model names

# Example
```julia
s = catalogue_summary("data/catalogue.jsonl")
println("Total entries: \$(s.total)")
println("Models: \$(s.models)")
println("Per-model: \$(s.per_model)")
```
"""
function catalogue_summary(path::AbstractString)
    isfile(path) || error("Catalogue file not found: $path")

    total = 0
    per_model = Dict{String,Int}()
    dimensions = Dict{Int,Int}()

    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        d = JSON3.read(stripped, Dict{String,Any})

        total += 1
        name = String(d["name"])
        per_model[name] = get(per_model, name, 0) + 1

        dim = length(d["p_true"])
        dimensions[dim] = get(dimensions, dim, 0) + 1
    end

    return (
        total = total,
        per_model = per_model,
        dimensions = dimensions,
        models = sort(collect(keys(per_model))),
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Screening → Catalogue Bridge
# ═══════════════════════════════════════════════════════════════════════════════

"""
    screening_to_catalogue(result, model_fn, ic, bounds, time_interval; kwargs...) -> Vector{CatalogueEntry}

Convert the top-ranked candidates from a `ScreeningResult` into `CatalogueEntry`
objects ready for persistence. The screening pipeline discovers good `p_true`
values; this function packages them with the context needed to reconstruct
objectives later.

# Arguments
- `result::ScreeningResult`: Output from `screen_and_probe`
- `model_fn::Function`: Zero-arg model factory (e.g., `define_daisy_ex3_model_4D`).
  Must be registered in `MODEL_REGISTRY` for JSONL serialization.
- `ic::Vector{Float64}`: Initial conditions used during screening
- `bounds::Vector{Tuple{Float64,Float64}}`: Parameter search bounds
- `time_interval`: Integration time span `[t_start, t_end]`

# Keyword Arguments
- `top_n::Int = 10`: Number of top-ranked candidates to convert
- `name_prefix::String = "screened"`: Prefix for auto-generated entry names
  (entries named `"{prefix}_1"`, `"{prefix}_2"`, etc.)
- `description::String = "Auto-discovered via screen_and_probe"`: Description
  for all entries
- `numpoints::Int = 30`: Number of time points for objective evaluation
- `distance_function::Function = L2_norm`: Per-output distance metric
- `aggregate_distances::Function = sum`: Multi-output aggregation
- `eval_timeout::Union{Nothing,Float64} = nothing`: Per-evaluation timeout

# Returns
`Vector{CatalogueEntry}` of length `min(top_n, length(result.ranking.ranked_indices))`.

# Example
```julia
result = screen_and_probe(model, outputs, ic, bounds, tspan; n_candidates=200)
entries = screening_to_catalogue(result, define_daisy_ex3_model_4D, ic, bounds, tspan;
    top_n=5, name_prefix="DAISY_4D")
save_catalogue("discovered.jsonl", entries)
```
"""
function screening_to_catalogue(
    result::ScreeningResult,
    model_fn::Function,
    ic::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}},
    time_interval;
    top_n::Int = 10,
    name_prefix::String = "screened",
    description::String = "Auto-discovered via screen_and_probe",
    numpoints::Int = 30,
    distance_function::Function = L2_norm,
    aggregate_distances::Function = sum,
    eval_timeout::Union{Nothing,Float64} = nothing,
)
    top_n > 0 || error("top_n must be positive, got $top_n")

    ranked = result.ranking.ranked_indices
    n = min(top_n, length(ranked))

    entries = CatalogueEntry[]
    sizehint!(entries, n)

    for i = 1:n
        idx = ranked[i]
        p_true = result.sweep.valid[idx]

        push!(
            entries,
            CatalogueEntry(
                name = "$(name_prefix)_$i",
                description = description,
                model_fn = model_fn,
                p_true = p_true,
                ic = ic,
                bounds = bounds,
                time_interval = collect(Float64, time_interval),
                numpoints = numpoints,
                distance_function = distance_function,
                aggregate_distances = aggregate_distances,
                eval_timeout = eval_timeout,
            ),
        )
    end

    return entries
end
