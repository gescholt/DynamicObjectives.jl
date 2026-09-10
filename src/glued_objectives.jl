# glued_objectives.jl
#
# Generate higher-dim test problems by Cartesian-product gluing of lower-dim
# objectives. The key property:
#
#     F(x, y) := f1(x) + f2(y)
#     ⇒ ∇F = 0  iff  ∇f1 = 0 AND ∇f2 = 0
#     ⇒ CP(F) = CP(f1) × CP(f2)        (exact Cartesian product)
#     ⇒ H_F = blockdiag(H_f1, H_f2)    (block-diagonal Hessian)
#     ⇒ index(p, q) = index_f1(p) + index_f2(q)
#
# This gives us 4D (or generally n+m-D) test problems where:
# - The full critical-point set is known analytically
# - Each minimum corresponds to a (min, min) factor pair
# - Saddle topology is decomposable

# ── Construction ─────────────────────────────────────────────────────────────

"""
    glue(f1, bounds1, f2, bounds2; cps1=nothing, cps2=nothing) -> NamedTuple

Glue two objectives by independent-coordinate sum:

    F([x; y]) := f1(x) + f2(y)

# Arguments
- `f1, f2`: Each `f_i` must accept a `Vector{Float64}` of length `dim(f_i)`.
- `bounds1, bounds2`: `Vector{Tuple{Float64,Float64}}` for each factor.

# Keyword arguments
- `cps1, cps2`: Optional vectors of `Vector{Float64}` critical points for each
  factor. When both are supplied, the returned NamedTuple includes
  `oracle_cps` = the analytic Cartesian product (the full CP set of `F`).

# Returns
NamedTuple with:
- `F`: the glued objective `F(z) = f1(z[1:n1]) + f2(z[n1+1:end])`
- `bounds`: concatenated `[bounds1..., bounds2...]`
- `dim`: `length(bounds1) + length(bounds2)`
- `n1, n2`: factor dimensions
- `oracle_cps::Union{Vector{Vector{Float64}}, Nothing}`: present when both
  factor CP sets supplied; `nothing` otherwise.

# Example
```julia
f1(x) = (x[1]^2 - 1)^2 + x[2]^2          # 2D quartic, 2 minima at (±1, 0)
f2(y) = (y[1] - 0.3)^2 + (y[2] + 0.5)^2  # 2D bowl, 1 minimum at (0.3, -0.5)
g = glue(f1, [(-2.0, 2.0), (-1.0, 1.0)], f2, [(-1.0, 1.0), (-1.0, 1.0)];
         cps1=[[1.0, 0.0], [-1.0, 0.0], [0.0, 0.0]],   # 2 min + 1 saddle
         cps2=[[0.3, -0.5]])                            # 1 min
g.dim == 4
length(g.oracle_cps) == 3  # 3 × 1 = 3
```
"""
function glue(
    f1,
    bounds1::Vector{Tuple{Float64,Float64}},
    f2,
    bounds2::Vector{Tuple{Float64,Float64}};
    cps1::Union{Nothing,Vector{Vector{Float64}}} = nothing,
    cps2::Union{Nothing,Vector{Vector{Float64}}} = nothing,
)
    n1 = length(bounds1)
    n2 = length(bounds2)
    @assert n1 >= 1 && n2 >= 1 "glue requires nonempty bounds for both factors"

    # Closure-captured indices for fast slicing without allocations.
    F = let n1 = n1, f1 = f1, f2 = f2
        function _F(z::AbstractVector)
            x = view(z, 1:n1)
            y = view(z, (n1+1):length(z))
            return f1(collect(x)) + f2(collect(y))
        end
    end

    bounds = vcat(bounds1, bounds2)

    oracle_cps = if cps1 !== nothing && cps2 !== nothing
        result = Vector{Float64}[]
        for p in cps1, q in cps2
            push!(result, vcat(p, q))
        end
        result
    else
        nothing
    end

    return (; F, bounds, dim = n1 + n2, n1, n2, oracle_cps)
end

# ── Convenience: catalogue-based gluing ──────────────────────────────────────

"""
    glue_catalogue(name1::String, cat_path1::String,
                   name2::String, cat_path2::String;
                   cps1=nothing, cps2=nothing) -> NamedTuple

Load two catalogue entries and glue them into a single objective. Each
returned `objective` is a `TolerantObjective` evaluated independently on its
own coordinate slice — same construction as `glue`, but the factors come
from the catalogue's recorded models instead of arbitrary callables.

The wide-domain rule still applies (see project memory): only glue entries
whose individual catalogue bounds are tight enough that their own ODE solves
don't blow up. For `lv2d × lv2d`, this is fine; for two wide 4D entries you
should expect Inf-region pruning to dominate.

# Returns
NamedTuple with `F`, `bounds`, `dim`, `n1`, `n2`, `oracle_cps`,
`p_true_glued = vcat(p_true_1, p_true_2)`, and `factor_names`.
"""
function glue_catalogue(
    name1::String,
    cat_path1::String,
    name2::String,
    cat_path2::String;
    cps1::Union{Nothing,Vector{Vector{Float64}}} = nothing,
    cps2::Union{Nothing,Vector{Vector{Float64}}} = nothing,
    solver_method::String = "Tsit5",
    abstol::Float64 = 1e-4,
    reltol::Float64 = 1e-4,
)
    e1 = first(filter(e -> e.name == name1, load_catalogue(cat_path1; model_name = name1)))
    e2 = first(filter(e -> e.name == name2, load_catalogue(cat_path2; model_name = name2)))

    obj1 = _build_objective_for_glue(e1, solver_method, abstol, reltol)
    obj2 = _build_objective_for_glue(e2, solver_method, abstol, reltol)

    g = glue(obj1, e1.bounds, obj2, e2.bounds; cps1 = cps1, cps2 = cps2)

    return (;
        g.F,
        g.bounds,
        g.dim,
        g.n1,
        g.n2,
        g.oracle_cps,
        p_true_glued = vcat(e1.p_true, e2.p_true),
        factor_names = (name1, name2),
    )
end

function _build_objective_for_glue(entry, solver_method, abstol, reltol)
    model_mtk, _, _, outputs = entry.model_fn()
    return TolerantObjective(
        model_mtk,
        outputs,
        entry.ic,
        entry.p_true,
        entry.time_interval,
        entry.numpoints,
        entry.distance_function,
        entry.aggregate_distances,
        nothing;
        solver = _resolve_solver(solver_method),
        abstol = abstol,
        reltol = reltol,
    )
end

# ── Brute-force CP-count verification helper ─────────────────────────────────

"""
    count_oracle_cps_in_box(oracle_cps, bounds) -> Int

Number of analytic oracle CPs that fall inside a given hyperrectangle.
Useful for sanity-checking subdivision recovery against the full oracle.
"""
function count_oracle_cps_in_box(
    oracle_cps::Vector{Vector{Float64}},
    bounds::Vector{Tuple{Float64,Float64}},
)
    n = 0
    for p in oracle_cps
        if all(bounds[d][1] <= p[d] <= bounds[d][2] for d in 1:length(p))
            n += 1
        end
    end
    return n
end
