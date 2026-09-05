module DynamicObjectives

using ModelingToolkit
using StaticArrays
using OrdinaryDiffEqTsit5: Tsit5, AutoTsit5
using OrdinaryDiffEqVerner: Vern7, Vern9
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas5
using OrdinaryDiffEqSDIRK: TRBDF2, KenCarp4
using SciMLBase
using DataStructures
using LinearAlgebra
using JSON3
using Logging
using Printf
using PrettyTables
using Statistics: mean
using Term: Panel, apply_style
using TOML
using UnicodePlots

# Include submodules
include("compat_patches.jl")
include("systems/daisy_models.jl")
include("systems/lotka_volterra.jl")
include("systems/other_systems.jl")
include("systems/new_benchmarks.jl")
include("data_generation.jl")
include("error_metrics.jl")
include("screening.jl")
include("catalogue.jl")
include("display.jl")
include("experiment_utils.jl")
include("globtim_integration.jl")
include("screening_config.jl")
include("grid_scoring.jl")
include("parallel_eval.jl")
include("experiment_harness.jl")
include("glued_objectives.jl")

# System definitions - DAISY models
export define_daisy_ex3_model_4D, define_daisy_ex3_model_4D_no_input
export define_daisy_ex3_model_5D
export define_daisy_ex3_model_6D
export define_daisy_ex3_model_7D
export define_daisy_ex3_model_8D

# System definitions - Lotka-Volterra variants
export define_generalized_lotka_volterra_4D,
    define_constrained_lotka_volterra_4D,
    define_lotka_volterra_4D_simple,
    define_lotka_volterra_3D_model,
    define_lotka_volterra_3D_model_locally_identifiable,
    define_lotka_volterra_3D_model_v2,
    define_lotka_volterra_2D_model,
    define_lotka_volterra_2D_model_v2,
    define_lotka_volterra_2D_model_v3,
    define_lotka_volterra_2D_model_v3_two_outputs,
    define_lotka_volterra_2D_sciml_benchmark,
    define_coupled_lv2d_3d_model,
    define_lv2d_reparam_3d_model,
    define_lv_3d_symmetric_model,
    create_lv2d_localid1d_3d_objective

# System definitions - Other models
export define_fhn_driven_auto, define_fhn_driven_A015, define_fhn_driven_A030
export define_fhn_driven3_auto
export define_fhn_coupled_A030_k000, define_fhn_coupled_A030_k002
export define_fhn_coupled_A030_k005, define_fhn_coupled_A030_k010
export define_fitzhugh_nagumo_3D_model,
    define_fitzhugh_nagumo_3D_model_two_outputs,
    define_goodwin_oscillator_4D,
    define_goodwin_oscillator_5D,
    define_goodwin_oscillator_6D,
    define_goodwin_oscillator_4D_hill2,
    define_goodwin_oscillator_4D_hill4,
    define_goodwin_oscillator_4D_hill6,
    define_goodwin_oscillator_3D,
    define_goodwin_3d_product_obs_model,
    define_rosenzweig_macarthur_3d_model,
    define_lorenz_3D_model,
    define_rossler_3D_model,
    define_goodwin_3d_locally_id_model,
    define_fhn_3d_locally_id_model,
    define_simple_2D_model_locally_identifiable,
    define_simple_2D_model_locally_identifiable_square,
    define_simple_1D_model_locally_identifiable,
    # New benchmark models (epidemiology, neuroscience, PK, chemistry, biochemistry)
    define_sir_2d_model,
    define_seir_3d_model,
    define_hindmarsh_rose_3d_model,
    define_pk_2comp_3d_model,
    define_brusselator_2d_model,
    define_michaelis_menten_2d_model,
    define_mm_chain_3d_model,
    define_rosenzweig_macarthur_4d_model

# Data generation
export sample_data

# Error metrics
export make_error_distance,
    L1_norm,
    L2_norm,
    L2_squared,
    log_L2_norm,
    # Aggregation strategy registry (bead 0iq)
    AGGREGATION_STRATEGIES,
    resolve_aggregation,
    # Partial observability helper (bead dds)
    make_partial_observability_distance,
    # Tolerant objective wrapper (mutable solver/tolerance settings)
    TolerantObjective,
    set_tolerance!,
    set_solver!,
    get_tolerance,
    get_solver,
    # Re-exported ODE solvers (needed for TolerantObjective / make_error_distance,
    # ext module, tests, and examples)
    Tsit5,
    Vern7,
    Vern9,
    AutoTsit5,
    Rosenbrock23,
    Rodas5,
    TRBDF2,
    KenCarp4,
    # Re-exported SciMLBase symbols (needed by tests and examples)
    ODEProblem

# Screening — functions
export is_trajectory_bounded, sweep_p_true, probe_landscape, rank_probes, screen_and_probe

# Screening — result types
export BoundednessResult,
    RejectedCandidate,
    SweepDiagnostics,
    SweepResult,
    ProbeResult,
    RankingResult,
    ScreeningResult

# Display infrastructure
export DisplayConfig,
    set_display_config!,
    display_model_summary,
    display_parameters,
    display_time_series,
    display_comparison,
    display_error_metrics,
    display_optimization_result,
    display_optimization_progress,
    display_section,
    display_subsection,
    display_results,
    display_gradient_analysis,
    display_quality_summary,
    display_degree_comparison,
    # Progress bars and spinners
    ProgressBar,
    progress_bar,
    update_progress!,
    finish_progress!,
    with_progress,
    format_duration,
    Spinner,
    spinner,
    spin!,
    stop_spinner!,
    with_spinner,
    display_step,
    # Critical points display
    display_top_critical_points,
    format_critical_points_markdown

# Catalogue — struct + objective creation
export CatalogueEntry,
    create_objective, create_tolerant_objective, create_refinement_objective, dimension

# Catalogue — function registries
export register_model!,
    register_distance!,
    register_aggregation!,
    MODEL_REGISTRY,
    DISTANCE_REGISTRY,
    AGGREGATION_REGISTRY

# Catalogue — JSONL persistence
export save_catalogue, append_catalogue, load_catalogue, catalogue_summary

# Catalogue — screening bridge
export screening_to_catalogue

# Experiment utilities (shared across validation scripts)
export CandidateResult, build_bounds, print_candidate_summary_table

# TOML-driven experiment pipeline
export run_experiment_from_config
export build_experiment_objective

# TOML-driven screening pipeline
export ScreeningConfig, load_screening_config, run_screening_from_config

# Grid-based interestingness scoring
export GridScoreResult,
    score_landscape_grid,
    interestingness_score,
    score_top_candidates,
    print_grid_score,
    structure_score,
    difficulty_profile

# Parallel grid evaluation (thread-safe factory pattern)
export evaluate_grid_threaded

# Candidate-level parallelism for catalogue experiments (bead 1yt)
export run_catalogue_experiments

# Cartesian-product gluing for higher-dim test problems with known CPs (bead zwbs.10.1)
export glue, glue_catalogue, count_oracle_cps_in_box

# Initialize function registries with all known model/distance/aggregation functions
_init_registries!()

end # module DynamicObjectives
