# Basic Workflow Example
# A simple, practical example showing typical usage of DynamicObjectives

using DynamicObjectives

println("\n" * "="^60)
println("Basic Parameter Estimation Workflow")
println("="^60 * "\n")

# 1. Define the ODE system
println("Step 1: Define ODE model")
model, params, states, outputs = define_daisy_ex3_model_4D()
display_model_summary(model)

# 2. Set up the problem
println("\nStep 2: Configure problem parameters")
param_names = [:α, :β, :γ, :δ]
p_true = [0.1, 0.2, 0.3, 0.4]  # True parameter values
p_initial = [0.15, 0.25, 0.35, 0.45]  # Initial guess
ic = [1.0, 2.0, 1.0, 1.0]  # Initial conditions
time_interval = [0.0, 10.0]  # Time span
num_points = 25  # Number of sampling points

display_parameters(param_names, [p_true p_initial], labels = ["True", "Initial Guess"])

# 3. Generate reference data (in real use, this would be experimental data)
println("\nStep 3: Generate reference data")
using DynamicObjectives: ODEProblem
problem =
    # Bind the initial condition by `unknowns(model)`, which is what
    # `make_error_distance` does; `states` is model_fn's DECLARATION order and the
    # two differ for these models.
    ODEProblem(model, merge(Dict(ModelingToolkit.unknowns(model) .=> ic), Dict(params .=> p_true)), time_interval)

data_reference = sample_data(problem, model, outputs, time_interval, p_true, ic, num_points)

println("Generated $(num_points) time points")
display_time_series(data_reference, show_plot = true)

# 4. Create error/objective function
println("\nStep 4: Create objective function")
error_func = make_error_distance(
    model,
    outputs,
    ic,
    p_true,
    time_interval,
    num_points,
    L2_norm,  # Distance metric
    first,    # Aggregate function (use first output)
    nothing;  # No noise
    return_inf_on_error = true,
    eval_timeout = 5.0,
)

# 5. Evaluate at different parameter values
println("\nStep 5: Evaluate objective at different parameters")
error_at_true = error_func(p_true)
error_at_initial = error_func(p_initial)

error_comparison = Dict(
    "Error at True Params" => error_at_true,
    "Error at Initial Guess" => error_at_initial,
    "Ratio" => error_at_initial / (error_at_true + 1e-10),
)

display_error_metrics(error_comparison)

# 6. Visualize the difference
println("\nStep 6: Compare predictions")
data_initial =
    sample_data(problem, model, outputs, time_interval, p_initial, ic, num_points)

display_comparison(
    data_reference,
    data_initial,
    labels = ["True Parameters", "Initial Guess"],
)

# 7. Summary
println("\n" * "="^60)
println("Workflow Complete!")
println("="^60)
println("""
Next steps for parameter estimation:
  1. Use an optimization algorithm (e.g., Optimization.jl)
  2. Pass 'error_func' as the objective
  3. Start from 'p_initial' and search for parameters that minimize error
  4. Use display functions to monitor progress and visualize results

Example optimization setup:
  using Optimization, OptimizationOptimJL

  # Define optimization problem
  optf = OptimizationFunction(error_func, Optimization.AutoForwardDiff())
  prob = OptimizationProblem(optf, p_initial)

  # Solve with callback for progress display
  callback = (p, l, iter) -> begin
      display_optimization_progress(iter, p, l)
      return false
  end

  sol = solve(prob, BFGS(), callback=callback)

  # Display final results
  result = (
      params=sol.u,
      error=sol.objective,
      iterations=sol.iterations,
      converged=sol.retcode == :success,
      time_elapsed=sol.solve_time
  )
  display_optimization_result(result)
""")
println("="^60 * "\n")

# ==============================================================================
# Step 8: Solver Comparison for DAISY Ex3
# ==============================================================================

println("\n" * "="^60)
println("Step 8: ODE Solver Comparison for DAISY Ex3")
println("="^60 * "\n")

println("DAISY Ex3 may exhibit mild stiffness depending on parameters.")
println("Comparing explicit vs implicit-explicit solvers:\n")

using DynamicObjectives: Vern9, Vern7, Tsit5, AutoTsit5, Rosenbrock23
using Printf

solver_configs = [
    (Tsit5(), "Tsit5 (explicit RK5/4)"),
    (Vern7(), "Vern7 (explicit RK7/6)"),
    (Vern9(), "Vern9 (explicit RK9/8)"),
    (AutoTsit5(Rosenbrock23()), "AutoTsit5(Rosenbrock23) (auto-stiff)"),
]

for tol in [1e-6, 1e-10]
    println("  --- Tolerance = $(tol) ---")
    for (slvr, name) in solver_configs
        obj_s = make_error_distance(
            model,
            outputs,
            ic,
            p_true,
            time_interval,
            num_points,
            L2_norm;
            return_inf_on_error = true,
            solver = slvr,
            abstol = tol,
            reltol = tol,
        )

        # Benchmark: evaluate at a few points
        n_bench = 10
        t_start = time()
        for _ in 1:n_bench
            obj_s(p_initial)
        end
        avg_ms = 1000 * (time() - t_start) / n_bench

        val_at_true = obj_s(p_true)
        val_at_test = obj_s(p_initial)

        @printf(
            "    %-35s: val@true=%.2e  val@test=%.6e  avg=%.1fms\n",
            name,
            val_at_true,
            val_at_test,
            avg_ms
        )
    end
    println()
end

println("Note: If AutoTsit5(Rosenbrock23) is significantly slower but no more accurate,")
println("the system is not stiff at these parameters. Use Tsit5 or Vern7 for speed.")
println("\n" * "="^60 * "\n")
