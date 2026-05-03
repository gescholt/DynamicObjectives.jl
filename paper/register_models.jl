#!/usr/bin/env julia
# paper/register_models.jl — Register paper ODE models in catalogue
#
# Creates JSONL catalogue entries for models referenced in paper experiments.
# Idempotent: skips entries that already exist.
#
# Pre-built catalogue files are already shipped in paper/catalogue/.
# Only run this script if you need to regenerate them.
#
# Usage:
#   julia --project=. paper/register_models.jl

using DynamicObjectives

# ─────────────────────────────────────────────────────────────────────────────
# Configuration: define all entries to register
# ─────────────────────────────────────────────────────────────────────────────

const ENTRIES_TO_REGISTER = [
    # ── LV-2D (c = 1) ──────────────────────────────────────────────────────
    # ODE:  dx1/dt = a*x1 + b*x1*x2,  dx2/dt = b*x1*x2 + x2
    # Observed: y1 = x1
    (
        path = "paper/catalogue/lv2d_catalogue.jsonl",
        entry = CatalogueEntry(
            name = "LV2D_1",
            description = "2D Lotka-Volterra with c=1, params (a,b), 1 output",
            model_fn = define_lotka_volterra_2D_model,
            p_true = [0.5, -0.3],
            ic = [1.0, 0.5],
            bounds = [(-1.0, 2.0), (-1.0, 0.0)],
            time_interval = [0.0, 10.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = nothing,
        ),
    ),

    # ── LV-2D-v3 (c = 0.5) — time sensitivity study (Sec 5.2) ─────────────
    # ODE:  dx1/dt = a*x1 - b*x1*x2,  dx2/dt = -b*x2 + 0.5*x1*x2
    # Observed: y1 = x1
    (
        path = "paper/catalogue/lv2d_catalogue.jsonl",
        entry = CatalogueEntry(
            name = "LV2D_v3_1",
            description = "2D Lotka-Volterra with c=0.5, params (a,b), 1 output — time sensitivity study",
            model_fn = define_lotka_volterra_2D_model_v3,
            p_true = [1.0, 0.5],
            ic = [1.0, 0.5],
            bounds = [(0.0, 3.0), (0.0, 2.0)],
            time_interval = [0.0, 10.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = nothing,
        ),
    ),

    # ── LV-3D locally identifiable (a², b², c) ─────────────────────────────
    # ODE:  dx1/dt = a²*x1 + b²*x1*x2,  dx2/dt = b²*x1*x2 + c*x2
    # Observed: y1 = x1
    # Squared parameterization: (a,b) and (-a,-b) give identical dynamics
    (
        path = "paper/catalogue/lv3d_catalogue.jsonl",
        entry = CatalogueEntry(
            name = "LV3D_1",
            description = "3D LV locally identifiable (a²,b²,c), 1 output, multiple equivalent minimizers",
            model_fn = define_lotka_volterra_3D_model_locally_identifiable,
            p_true = [0.5, -0.3, 0.2],
            ic = [1.0, 0.5],
            bounds = [(-1.0, 2.0), (-1.0, 0.0), (-1.0, 1.0)],
            time_interval = [0.0, 10.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = nothing,
        ),
    ),

    # ── LV-4D Simple (non-identifiable negative result) ───────────────────
    # ODE:  dx1/dt = a*x1 + b*x1*x2,  dx2/dt = c*x1*x2 + d*x2
    # Observed: y1 = x1
    # NOT identifiable — paper Example 4.3 negative result
    (
        path = "paper/catalogue/lv4d_simple_catalogue.jsonl",
        entry = CatalogueEntry(
            name = "LV4D_simple_1",
            description = "Simple 2-species LV with 4 params (a,b,c,d), 1 output — NOT identifiable (paper negative result)",
            model_fn = define_lotka_volterra_4D_simple,
            p_true = [0.2, 0.4, 0.7, 0.9],
            ic = [0.3, 0.6],
            bounds = [(-0.5, 1.5), (-0.5, 1.5), (-0.5, 1.5), (-0.5, 1.5)],
            time_interval = [0.0, 1.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = nothing,
        ),
    ),

    # ── Goodwin oscillator 4D ──────────────────────────────────────────────
    # ODE:  dx1/dt = k1*K^n/(K^n+x3^n) - k2*x1
    #       dx2/dt = 0.3*x1 - k4*x2
    #       dx3/dt = k5*x2 - 0.5*x3
    # Observed: y1=x1, y2=x3
    (
        path = "paper/catalogue/goodwin4d_catalogue.jsonl",
        entry = CatalogueEntry(
            name = "Goodwin4D_1",
            description = "Goodwin oscillator 3-state, 4 unknown params (k1,k2,k4,k5), 2 outputs — Hill function repression",
            model_fn = define_goodwin_oscillator_4D,
            # Auto-discovered via screening; original defaults were [1.0, 0.5, 0.5, 0.3]
            p_true = [0.8088, 1.2467, 1.2985, 1.5309],
            ic = [0.5, 0.5, 0.5],
            bounds = [(0.1, 3.0), (0.1, 2.0), (0.1, 2.0), (0.1, 2.0)],
            time_interval = [0.0, 20.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = nothing,
        ),
    ),

    # ── Locally-identifiable 1D (a^2) ───────────────────────────────────────
    # ODE:  dx1/dt = x1 + a^2
    # Observed: y1 = x1
    # Non-globally-identifiable: +/-a give same output
    (
        path = "paper/catalogue/locally_id_catalogue.jsonl",
        entry = CatalogueEntry(
            name = "LocalID_1D_1",
            description = "1D locally-identifiable model, dx/dt = x + a^2, param (a) — non-globally-identifiable",
            model_fn = define_simple_1D_model_locally_identifiable,
            p_true = [1.5],
            ic = [1.0],
            bounds = [(-3.0, 3.0)],
            time_interval = [0.0, 5.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = nothing,
        ),
    ),

    # ── Locally-identifiable 2D square (b^2) ────────────────────────────────
    # ODE:  dx1/dt = a*x1 + b^2
    # Observed: y1 = x1
    # Non-globally-identifiable: +/-b give same output
    (
        path = "paper/catalogue/locally_id_catalogue.jsonl",
        entry = CatalogueEntry(
            name = "LocalID_2D_sq_1",
            description = "2D locally-identifiable model, dx/dt = a*x + b^2, params (a,b) — non-globally-identifiable",
            model_fn = define_simple_2D_model_locally_identifiable_square,
            p_true = [0.5, 2.0],
            ic = [1.0],
            bounds = [(-1.0, 2.0), (-3.0, 3.0)],
            time_interval = [0.0, 5.0],
            numpoints = 30,
            distance_function = L2_norm,
            aggregate_distances = sum,
            eval_timeout = nothing,
        ),
    ),
]

# ─────────────────────────────────────────────────────────────────────────────
# Registration logic (idempotent)
# ─────────────────────────────────────────────────────────────────────────────

function main()
    println("="^72)
    println("Register paper models in catalogue")
    println("="^72)

    registered = 0
    skipped = 0

    for spec in ENTRIES_TO_REGISTER
        path = spec.path
        entry = spec.entry

        # Check if entry already exists
        if isfile(path)
            existing = load_catalogue(path)
            existing_names = Set(e.name for e in existing)
            if entry.name in existing_names
                println("  SKIP  $(entry.name) — already in $(path)")
                skipped += 1
                continue
            end
        end

        # Append (creates file + dirs if needed)
        append_catalogue(path, entry)
        println("  ADD   $(entry.name) -> $(path)")
        println("        p_true = $(entry.p_true)")
        println("        ic     = $(entry.ic)")
        println("        bounds = $(entry.bounds)")
        registered += 1
    end

    println()
    println("Done: $registered registered, $skipped skipped")
    println("="^72)
end

main()
