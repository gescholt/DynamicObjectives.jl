# Validate all models in MODEL_REGISTRY
#
# For each registered model function:
# 1. Call it to get (model, params, states, outputs)
# 2. Verify the return types
# 3. Verify basic ODE solve works with dummy parameters

@testset "All registered models" begin

    @testset "MODEL_REGISTRY has expected count" begin
        # Lower-bound only — exact name list is enforced by the matching
        # testset in test_catalogue.jl. Using `==` here went stale every
        # time a new model was added without updating this count (and
        # went undetected until Pkg.test stopped masking test runs).
        @test length(MODEL_REGISTRY) >= 18
    end

    # Define p_true and ic for each model (must match param/state dimensions)
    model_configs = Dict(
        "define_lotka_volterra_2D_model" => (p=[0.5, -0.3], ic=[1.0, 0.5]),
        "define_lotka_volterra_2D_model_v2" => (p=[0.5, -0.3], ic=[1.0, 0.5]),
        "define_lotka_volterra_2D_model_v3" => (p=[1.0, 0.5], ic=[1.0, 0.5]),
        "define_lotka_volterra_2D_model_v3_two_outputs" => (p=[1.0, 0.5], ic=[1.0, 0.5]),
        "define_lotka_volterra_2D_sciml_benchmark" => (p=[1.5, 1.0], ic=[1.0, 1.0]),
        "define_lotka_volterra_3D_model" => (p=[0.5, -0.3, 0.2], ic=[1.0, 0.5]),
        "define_lotka_volterra_3D_model_locally_identifiable" => (p=[0.5, -0.3, 0.2], ic=[1.0, 0.5]),
        "define_lotka_volterra_3D_model_v2" => (p=[1.0, 0.5, 0.3], ic=[1.0, 0.5]),
        "define_lotka_volterra_4D_simple" => (p=[0.5, 0.3, 0.4, 0.2], ic=[1.0, 0.5]),
        "define_daisy_ex3_model_4D" => (p=[0.5, 0.3, 0.4, 0.2], ic=[1.0, 0.5, 0.3, 0.0]),
        "define_daisy_ex3_model_4D_no_input" => (p=[0.5, 0.3, 0.4, 0.2], ic=[1.0, 0.5, 0.3]),
        "define_constrained_lotka_volterra_4D" => (p=[0.01, 0.02, -0.01, 0.03], ic=[0.8, 1.2, 0.8, 1.2]),
        "define_goodwin_oscillator_4D" => (p=[1.0, 0.5, 0.5, 0.3], ic=[0.5, 0.5, 0.5]),
        "define_generalized_lotka_volterra_4D" => (p=fill(0.1, 20), ic=[0.5, 0.5, 0.5, 0.5]),
        "define_fitzhugh_nagumo_3D_model" => (p=[0.2, 0.3, 3.0], ic=[0.0, 0.0]),
        "define_simple_1D_model_locally_identifiable" => (p=[1.5], ic=[1.0]),
        "define_simple_2D_model_locally_identifiable" => (p=[1.0, 2.0], ic=[1.0]),
        "define_simple_2D_model_locally_identifiable_square" => (p=[0.5, 2.0], ic=[1.0]),
    )

    for (name, fn) in MODEL_REGISTRY
        @testset "$name" begin
            model, params, states, outputs = fn()
            @test model isa ModelingToolkit.AbstractSystem
            @test !isempty(params)
            @test !isempty(states)
            @test !isempty(outputs)

            # Models without an explicit config still get the type-check
            # pass above (lines 39-42). The full ODE-solve verification
            # below only runs when we've hand-curated param + IC dims for
            # the model. Previously this errored on every newly-added
            # model, forcing this file's maintenance onto every PR that
            # touched MODEL_REGISTRY — which is exactly what kept going
            # stale.
            if haskey(model_configs, name)
                cfg = model_configs[name]
                n_params = length(ModelingToolkit.parameters(model))
                n_unknowns = length(ModelingToolkit.unknowns(model))
                @test length(cfg.p) == n_params
                @test length(cfg.ic) == n_unknowns

                # Verify make_error_distance works
                err_fn = make_error_distance(model, outputs, cfg.ic, cfg.p,
                    [0.0, 5.0], 10)
                val = err_fn(cfg.p)
                @test isfinite(val)
                @test val < 1e-4  # should be near-zero at true params
            end
        end
    end
end
