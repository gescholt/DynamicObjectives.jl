module Dynamic_objectives

using ModelingToolkit
using StaticArrays
using OrdinaryDiffEq
using DataStructures
using LinearAlgebra
using Logging
using PrettyTables
using Term: Panel, apply_style
using UnicodePlots

# Include submodules
include("systems/daisy_models.jl")
include("systems/lotka_volterra.jl")
include("systems/other_systems.jl")
include("data_generation.jl")
include("error_metrics.jl")
include("support.jl")
include("display.jl")
include("globtim_integration.jl")
include("globtim_postprocessing.jl")

# System definitions - DAISY models
export define_daisy_ex3_model_4D,
    define_daisy_ex3_model_4D_no_input

# System definitions - Lotka-Volterra variants
export define_generalized_lotka_volterra_4D,
    define_constrained_lotka_volterra_4D,
    define_lotka_volterra_3D_model,
    define_lotka_volterra_3D_model_v2,
    define_lotka_volterra_2D_model,
    define_lotka_volterra_2D_model_v2,
    define_lotka_volterra_2D_model_v3,
    define_lotka_volterra_2D_model_v3_two_outputs

# System definitions - Other models
export define_fitzhugh_nagumo_3D_model,
    define_simple_2D_model_locally_identifiable,
    define_simple_2D_model_locally_identifiable_square,
    define_simple_1D_model_locally_identifiable

# Data generation
export sample_data

# Error metrics
export make_error_distance,
    L1_norm,
    L2_norm,
    log_L2_norm

# Display infrastructure
export DisplayConfig,
    set_display_config!,
    display_model_summary,
    display_parameters,
    display_time_series,
    display_comparison,
    display_error_metrics,
    display_optimization_result,
    display_optimization_progress

# Globtim integration
export create_globtim_objective,
    run_globtim_optimization

# Globtim postprocessing (Phase 3)
export RefinementConfig,
    ode_refinement_config,
    refine_experiment_results

end # module Dynamic_objectives
