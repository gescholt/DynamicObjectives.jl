# Test Campaign Configurations
# Import this file to get pre-configured test cases

"""
Test configurations for the original 14 benchmark models.
See docs/src/model_catalog.md for the full catalog of all 29 models.
"""

using DynamicObjectives

# ==============================================================================
# PHASE 1: EASY (2 parameters)
# ==============================================================================

easy_configs = CatalogueEntry[
    CatalogueEntry(
        name = "LV_2D_v1",
        model_fn = define_lotka_volterra_2D_model,
        p_true = [0.5, -0.3],
        ic = [1.0, 0.5],
        bounds = [(-1.0, 2.0), (-1.0, 0.0)],
        time_interval = [0.0, 20.0],
        numpoints = 30,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "Classic 2D Lotka-Volterra, 1 output",
    ),
    CatalogueEntry(
        name = "LV_2D_v2",
        model_fn = define_lotka_volterra_2D_model_v2,
        p_true = [0.5, -0.3],
        ic = [1.0, 0.5],
        bounds = [(-1.0, 2.0), (-1.0, 0.0)],
        time_interval = [0.0, 20.0],
        numpoints = 30,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "2D LV with c=0.1, 1 output",
    ),
    CatalogueEntry(
        name = "LV_2D_v3",
        model_fn = define_lotka_volterra_2D_model_v3,
        p_true = [1.0, 0.5],
        ic = [1.0, 0.5],
        bounds = [(0.0, 3.0), (0.0, 2.0)],
        time_interval = [0.0, 20.0],
        numpoints = 30,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "2D LV with c=0.5, 1 output",
    ),
    CatalogueEntry(
        name = "LV_2D_v3_2outputs",
        model_fn = define_lotka_volterra_2D_model_v3_two_outputs,
        p_true = [1.0, 0.5],
        ic = [1.0, 0.5],
        bounds = [(0.0, 3.0), (0.0, 2.0)],
        time_interval = [0.0, 20.0],
        numpoints = 30,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "2D LV with full observability (2 outputs)",
    ),
]

# ==============================================================================
# PHASE 2: MEDIUM (3-4 parameters)
# ==============================================================================

medium_configs = CatalogueEntry[
    CatalogueEntry(
        name = "LV_3D_locally_id",
        model_fn = define_lotka_volterra_3D_model_locally_identifiable,
        p_true = [0.5, -0.3, 0.2],
        ic = [1.0, 0.5],
        bounds = [(-1.0, 2.0), (-1.0, 0.0), (-1.0, 1.0)],
        time_interval = [0.0, 20.0],
        numpoints = 40,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "3D LV locally identifiable (a²,b²,c), multiple equivalent minimizers",
    ),
    CatalogueEntry(
        name = "LV_3D_v2",
        model_fn = define_lotka_volterra_3D_model_v2,
        p_true = [1.0, 0.5, 0.3],
        ic = [1.0, 0.5],
        bounds = [(0.0, 3.0), (0.0, 2.0), (0.0, 1.0)],
        time_interval = [0.0, 20.0],
        numpoints = 40,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "3D LV variant 2, standard formulation",
    ),
    CatalogueEntry(
        name = "DAISY_Ex3_with_input",
        model_fn = define_daisy_ex3_model_4D,
        p_true = [0.5, 0.3, 0.4, 0.2],
        ic = [1.0, 0.5, 0.3, 0.0],
        bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)],
        time_interval = [0.0, 20.0],
        numpoints = 50,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "DAISY benchmark with time-dependent input",
    ),
    CatalogueEntry(
        name = "DAISY_Ex3_no_input",
        model_fn = define_daisy_ex3_model_4D_no_input,
        p_true = [0.5, 0.3, 0.4, 0.2],
        ic = [1.0, 0.5, 0.3],
        bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)],
        time_interval = [0.0, 20.0],
        numpoints = 50,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "DAISY benchmark without input",
    ),
    CatalogueEntry(
        name = "LV_4D_Constrained",
        model_fn = define_constrained_lotka_volterra_4D,
        p_true = [0.01, 0.02, -0.01, 0.03],
        ic = [0.8, 1.2, 0.8, 1.2],
        bounds = [(-0.1, 0.1), (-0.1, 0.1), (-0.1, 0.1), (-0.1, 0.1)],
        time_interval = [0.0, 30.0],
        numpoints = 60,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "4D LV with skew-symmetric perturbations",
    ),
]

# ==============================================================================
# PHASE 3: HARD (identifiability issues or oscillatory)
# ==============================================================================

hard_configs = CatalogueEntry[
    CatalogueEntry(
        name = "FitzHugh_Nagumo",
        model_fn = define_fitzhugh_nagumo_3D_model,
        p_true = [0.8, 0.7, 0.8],
        ic = [0.0, 0.0],
        bounds = [(0.1, 2.0), (0.0, 2.0), (0.0, 2.0)],
        time_interval = [0.0, 50.0],
        numpoints = 100,
        distance_function = L2_norm,
        eval_timeout = 10.0,
        description = "Neuronal excitability model, oscillatory",
    ),
    CatalogueEntry(
        name = "Simple_1D_square",
        model_fn = define_simple_1D_model_locally_identifiable,
        p_true = [1.5],
        ic = [1.0],
        bounds = [(-3.0, 3.0)],
        time_interval = [0.0, 5.0],
        numpoints = 20,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "Non-identifiable: ±a give same output",
    ),
    CatalogueEntry(
        name = "Simple_2D_product",
        model_fn = define_simple_2D_model_locally_identifiable,
        p_true = [2.0, 3.0],
        ic = [1.0],
        bounds = [(0.1, 5.0), (0.1, 5.0)],
        time_interval = [0.0, 5.0],
        numpoints = 20,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "Non-identifiable: only a*b and a+b observable",
    ),
    CatalogueEntry(
        name = "Simple_2D_square",
        model_fn = define_simple_2D_model_locally_identifiable_square,
        p_true = [0.5, 2.0],
        ic = [1.0],
        bounds = [(-1.0, 2.0), (-3.0, 3.0)],
        time_interval = [0.0, 5.0],
        numpoints = 20,
        distance_function = L2_norm,
        eval_timeout = nothing,
        description = "Non-identifiable: ±b give same output",
    ),
]

# ==============================================================================
# PHASE 4: VERY HARD (high-dimensional)
# ==============================================================================

very_hard_configs = CatalogueEntry[CatalogueEntry(
    name = "LV_4D_Generalized_20D",
    model_fn = define_generalized_lotka_volterra_4D,
    p_true = [
        1.0,
        1.0,
        1.0,
        1.0,  # Growth rates
        -0.5,
        -0.1,
        -0.1,
        -0.1,
        -0.1,
        -0.5,
        -0.1,
        -0.1,
        -0.1,
        -0.1,
        -0.5,
        -0.1,
        -0.1,
        -0.1,
        -0.1,
        -0.5,
    ],
    ic = [0.5, 0.5, 0.5, 0.5],
    bounds = [
        # Growth rates
        (0.0, 3.0),
        (0.0, 3.0),
        (0.0, 3.0),
        (0.0, 3.0),
        # Interaction matrix (diagonal negative, off-diagonal can be positive/negative)
        (-2.0, 0.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-2.0, 0.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-2.0, 0.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
        (-2.0, 0.0),
    ],
    time_interval = [0.0, 15.0],
    numpoints = 40,
    distance_function = log_L2_norm,  # Note: using log scale for wide range
    eval_timeout = 15.0,
    description = "20D full interaction matrix - ultimate challenge",
),]

# ==============================================================================
# CONVENIENCE COLLECTIONS
# ==============================================================================

# All configs
all_configs = vcat(easy_configs, medium_configs, hard_configs, very_hard_configs)

# By phase
configs_by_phase = Dict(
    "EASY" => easy_configs,
    "MEDIUM" => medium_configs,
    "HARD" => hard_configs,
    "VERY_HARD" => very_hard_configs,
)

# Quick access by name
configs_by_name = Dict(c.name => c for c in all_configs)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

"""
    get_config(name::String)

Get a specific configuration by name
"""
function get_config(name::String)
    return configs_by_name[name]
end

"""
    get_configs_by_phase(phase::String)

Get all configurations for a specific phase
"""
function get_configs_by_phase(phase::String)
    return configs_by_phase[uppercase(phase)]
end

"""
    list_configs()

List all available configurations
"""
function list_configs()
    println("Available test configurations:")
    println()

    for phase in ["EASY", "MEDIUM", "HARD", "VERY_HARD"]
        configs = configs_by_phase[phase]
        println("$phase ($(length(configs)) models):")
        for c in configs
            nparams = dimension(c)
            println("  $(c.name) ($nparams params) - $(c.description)")
        end
        println()
    end
end

# Export for convenience
export easy_configs, medium_configs, hard_configs, very_hard_configs
export all_configs, configs_by_phase, configs_by_name
export get_config, get_configs_by_phase, list_configs
