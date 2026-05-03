# Screening
# Fast pre-screening functions for p_true candidate validation.
# These use coarse ODE solves to reject parameter vectors that produce
# divergent, explosive, or pathological trajectories before expensive
# grid evaluations.

# ============================================================================
# Result types
# ============================================================================

"""
    BoundednessResult

Result of trajectory boundedness pre-screen (`is_trajectory_bounded`).

# Fields
- `bounded::Bool`: true if trajectory stays within threshold and solver succeeds
- `max_amplitude::Float64`: maximum |x_i(t)| across all states and time points (Inf on failure)
- `solver_success::Bool`: true if ODE solver returned successfully
"""
struct BoundednessResult
    bounded::Bool
    max_amplitude::Float64
    solver_success::Bool
end

"""
    RejectedCandidate

A candidate p_true that failed the boundedness pre-screen.

# Fields
- `p_true::Vector{Float64}`: The rejected parameter vector
- `reason::Symbol`: `:solver_failure` or `:amplitude_exceeded`
- `max_amplitude::Float64`: Maximum amplitude observed (Inf on solver failure)
"""
struct RejectedCandidate
    p_true::Vector{Float64}
    reason::Symbol
    max_amplitude::Float64
end

"""
    SweepDiagnostics

Timing and amplitude statistics from a `sweep_p_true` run.

# Fields
- `elapsed_seconds::Float64`: Total wall time for the sweep
- `n_total::Int`: Total number of candidates generated
- `n_valid::Int`: Number of candidates that passed the pre-screen
- `n_rejected_solver_failure::Int`: Number rejected due to ODE solver failure
- `n_rejected_amplitude::Int`: Number rejected due to amplitude exceeding threshold
- `ms_per_candidate::Float64`: Average wall time per candidate in milliseconds
- `amplitude_min::Float64`: Minimum max-amplitude among valid candidates (NaN if none)
- `amplitude_max::Float64`: Maximum max-amplitude among valid candidates (NaN if none)
- `amplitude_median::Float64`: Median max-amplitude among valid candidates (NaN if none)
"""
struct SweepDiagnostics
    elapsed_seconds::Float64
    n_total::Int
    n_valid::Int
    n_rejected_solver_failure::Int
    n_rejected_amplitude::Int
    ms_per_candidate::Float64
    amplitude_min::Float64
    amplitude_max::Float64
    amplitude_median::Float64
end

"""
    SweepResult

Result of batch candidate screening (`sweep_p_true`).

# Fields
- `valid::Vector{Vector{Float64}}`: p_true vectors that passed the pre-screen
- `rejected::Vector{RejectedCandidate}`: Rejected candidates with reasons
- `pass_rate::Float64`: Fraction of candidates that passed (0.0 to 1.0)
- `diagnostics::SweepDiagnostics`: Timing and amplitude statistics
"""
struct SweepResult
    valid::Vector{Vector{Float64}}
    rejected::Vector{RejectedCandidate}
    pass_rate::Float64
    diagnostics::SweepDiagnostics
end

"""
    ProbeResult

Result of landscape probing for a single p_true (`probe_landscape`).

# Fields
- `objective_at_true::Float64`: Objective value at p_true (sanity check, should be ~0)
- `probe_values::Vector{Float64}`: Objective values at random probe points
- `probe_points::Vector{Vector{Float64}}`: The random points used for probing
- `dynamic_range::Float64`: log10(max/min) of finite nonzero probe values (NaN if <2 finite)
- `variance::Float64`: Variance of finite probe values (NaN if <2 finite)
- `n_finite::Int`: Number of probes that returned finite values
- `n_inf::Int`: Number of probes that returned Inf (ODE failure regions)
- `construction_time_ms::Float64`: Time to build the objective function in milliseconds
"""
struct ProbeResult
    objective_at_true::Float64
    probe_values::Vector{Float64}
    probe_points::Vector{Vector{Float64}}
    dynamic_range::Float64
    variance::Float64
    n_finite::Int
    n_inf::Int
    construction_time_ms::Float64
end

"""
    RankingResult

Result of ranking probe results (`rank_probes`).

# Fields
- `ranked_indices::Vector{Int}`: Indices into the probes vector, sorted best-first
- `filtered_out::Vector{Int}`: Indices that failed quality filters
- `filter_reasons::Vector{Symbol}`: Reason for each filtered-out index
  (`:low_finite_fraction` or `:high_noise_floor`)
"""
struct RankingResult
    ranked_indices::Vector{Int}
    filtered_out::Vector{Int}
    filter_reasons::Vector{Symbol}
end

"""
    ScreeningResult

Result of the full screening pipeline (`screen_and_probe`).

# Fields
- `sweep::SweepResult`: Raw output from `sweep_p_true`
- `probes::Vector{ProbeResult}`: Probe results, parallel to `sweep.valid`
- `ranking::RankingResult`: Ranked indices and filter results
"""
struct ScreeningResult
    sweep::SweepResult
    probes::Vector{ProbeResult}
    ranking::RankingResult
end

# ============================================================================
# Trajectory boundedness pre-screen
# ============================================================================

"""
    is_trajectory_bounded(model, ic, p_true, time_interval; kwargs...) -> BoundednessResult

Fast pre-screen: solve the ODE at coarse tolerance and check that all state
variables stay bounded. Returns a `BoundednessResult` struct.

Uses cheap Tsit5() at 1e-4 tolerance by default — this is a fast reject filter,
not a precision tool.

# Arguments
- `model::ModelingToolkit.AbstractSystem`: ODE system (already compiled via @mtkcompile)
- `ic::Vector{Float64}`: Initial conditions for all state variables
- `p_true::Vector{Float64}`: Candidate parameter vector to screen
- `time_interval`: `[t_start, t_end]` for simulation

# Keyword Arguments
- `threshold::Float64 = 1e6`: Maximum allowed absolute value for any state variable
- `numpoints::Int = 20`: Number of time points to sample (coarse is fine for screening)
- `solver = Tsit5()`: ODE solver (cheap default for fast screening)
- `abstol::Float64 = 1e-4`: Absolute tolerance (coarse for speed)
- `reltol::Float64 = 1e-4`: Relative tolerance (coarse for speed)

# Returns
`BoundednessResult` with fields:
- `bounded::Bool`: true if trajectory stays within threshold and solver succeeds
- `max_amplitude::Float64`: maximum |x_i(t)| across all states and time points (Inf on failure)
- `solver_success::Bool`: true if ODE solver returned successfully

# Examples
```julia
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
ic = [0.5, 1.0]
p_true = [1.0, 0.5]
result = is_trajectory_bounded(model, ic, p_true, [0.0, 10.0])
result.bounded      # true
result.max_amplitude  # e.g. 2.3
```
"""
function is_trajectory_bounded(
    model::ModelingToolkit.AbstractSystem,
    ic::Vector{Float64},
    p_true::Vector{Float64},
    time_interval;
    threshold::Float64 = 1e6,
    numpoints::Int = 20,
    solver = Tsit5(),
    abstol::Float64 = 1e-4,
    reltol::Float64 = 1e-4,
)
    @assert length(p_true) == length(ModelingToolkit.parameters(model)) (
        "Parameter vector length mismatch: got $(length(p_true)), " *
        "expected $(length(ModelingToolkit.parameters(model)))"
    )
    @assert length(ic) == length(ModelingToolkit.unknowns(model)) (
        "Initial conditions length mismatch: got $(length(ic)), " *
        "expected $(length(ModelingToolkit.unknowns(model)))"
    )
    @assert length(time_interval) == 2 "Time interval must be [start_time, end_time]"
    @assert time_interval[2] > time_interval[1] "End time must be greater than start time"
    @assert threshold > 0 "Threshold must be positive"
    @assert numpoints > 0 "Number of points must be positive"

    # Construct ODEProblem (same pattern as make_error_distance)
    problem = ODEProblem(
        model,
        merge(
            Dict(ModelingToolkit.unknowns(model) .=> ic),
            Dict(ModelingToolkit.parameters(model) .=> p_true),
        ),
        time_interval,
    )

    sampling_times = range(time_interval[1], time_interval[2], length = numpoints)

    # Solve with coarse tolerances, suppress all warnings
    local solution
    try
        solution = Logging.with_logger(Logging.NullLogger()) do
            ModelingToolkit.solve(
                problem,
                solver,
                saveat = sampling_times;
                abstol = abstol,
                reltol = reltol,
                verbose = false,
                maxiters = 1_000_000,
            )
        end
    catch e
        # ODE solver threw an exception — trajectory is not solvable
        if e isa InterruptException
            rethrow(e)
        end
        return BoundednessResult(false, Inf, false)
    end

    # Check solver return code (SciMLBase.ReturnCode.T enum)
    if solution.retcode != SciMLBase.ReturnCode.Success
        return BoundednessResult(false, Inf, false)
    end

    # Extract all state variable values and compute max absolute amplitude.
    # solution.u is a Vector of state vectors at each time point.
    max_amp = 0.0
    for u in solution.u
        for val in u
            abs_val = abs(val)
            if isnan(abs_val) || isinf(abs_val)
                return BoundednessResult(false, Inf, true)
            end
            if abs_val > max_amp
                max_amp = abs_val
            end
        end
    end

    bounded = max_amp < threshold
    return BoundednessResult(bounded, max_amp, true)
end

# ============================================================================
# Candidate generation
# ============================================================================

"""
    _generate_candidates(bounds, n, margin, sampling) -> Vector{Vector{Float64}}

Generate `n` candidate parameter vectors within `bounds`, respecting `margin`.

Internal function — not exported.

# Sampling strategies
- `:random` — uniform random within margined bounds
- `:grid` — uniform grid with `ceil(n^(1/dim))` points per dimension, truncated to `n`
"""
function _generate_candidates(
    bounds::Vector{Tuple{Float64,Float64}},
    n::Int,
    margin::Float64,
    sampling::Symbol,
)::Vector{Vector{Float64}}
    dim = length(bounds)

    # Compute margined ranges: shrink each dimension by margin on each side
    lo = [b[1] + margin * (b[2] - b[1]) for b in bounds]
    hi = [b[2] - margin * (b[2] - b[1]) for b in bounds]

    if sampling == :random
        return [[lo[d] + rand() * (hi[d] - lo[d]) for d in 1:dim] for _ in 1:n]
    elseif sampling == :grid
        n_per_dim = max(2, ceil(Int, n^(1.0 / dim)))
        grids = [collect(range(lo[d], hi[d], length = n_per_dim)) for d in 1:dim]
        # Cartesian product
        candidates = Vector{Vector{Float64}}()
        for combo in Iterators.product(grids...)
            push!(candidates, collect(combo))
            if length(candidates) >= n
                break
            end
        end
        return candidates[1:min(n, length(candidates))]
    else
        error("Unknown sampling strategy: :$sampling. Use :random or :grid.")
    end
end

# ============================================================================
# Sweep and filter
# ============================================================================

"""
    sweep_p_true(model, ic, bounds, time_interval; kwargs...) -> SweepResult

Generate N candidate p_true vectors, filter each through `is_trajectory_bounded`,
and return the valid set with pass-rate diagnostics.

This is the batch pipeline for discovering parameter vectors that produce
well-behaved ODE trajectories — a prerequisite for building interesting
objective functions.

# Arguments
- `model::ModelingToolkit.AbstractSystem`: ODE system (compiled via @mtkcompile)
- `ic::Vector{Float64}`: Initial conditions for all state variables
- `bounds::Vector{Tuple{Float64,Float64}}`: Per-parameter (lo, hi) bounds
- `time_interval`: `[t_start, t_end]` for trajectory simulation

# Keyword Arguments
- `n_candidates::Int = 1000`: Number of candidate p_true vectors to generate
- `margin::Float64 = 0.1`: Safety margin from bounds (0.1 = stay within inner 80%)
- `sampling::Symbol = :random`: Sampling strategy (:random or :grid)
- `threshold::Float64 = 1e6`: Max amplitude for trajectory boundedness
- `numpoints::Int = 20`: Time points for coarse ODE solve
- `solver = Tsit5()`: ODE solver for screening
- `abstol::Float64 = 1e-4`: Absolute tolerance for screening
- `reltol::Float64 = 1e-4`: Relative tolerance for screening
- `verbose::Bool = false`: Print progress and summary

# Returns
`SweepResult` with fields:
- `valid::Vector{Vector{Float64}}`: p_true vectors that passed the pre-screen
- `rejected::Vector{RejectedCandidate}`: Rejected candidates with reasons
- `pass_rate::Float64`: fraction of candidates that passed (0.0 to 1.0)
- `diagnostics::SweepDiagnostics`: timing, amplitude statistics, candidate counts

# Examples
```julia
model, params, states, outputs = define_lotka_volterra_2D_model_v3_two_outputs()
ic = [0.5, 1.0]
bounds = [(0.5, 2.0), (0.2, 1.0)]

result = sweep_p_true(model, ic, bounds, [0.0, 10.0]; n_candidates=100)
println("Pass rate: \$(result.pass_rate)")          # e.g. 0.95
println("Found \$(length(result.valid)) valid p_true vectors")
println("Rejected \$(length(result.rejected)) candidates")

# Use valid p_true for building objectives
for p in result.valid[1:3]
    obj = make_error_distance(model, outputs, ic, p, [0.0, 10.0], 30, L2_norm)
    println("  p=\$p  f(p)=\$(obj(p))")
end
```
"""
function sweep_p_true(
    model::ModelingToolkit.AbstractSystem,
    ic::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}},
    time_interval;
    n_candidates::Int = 1000,
    margin::Float64 = 0.1,
    sampling::Symbol = :random,
    threshold::Float64 = 1e6,
    numpoints::Int = 20,
    solver = Tsit5(),
    abstol::Float64 = 1e-4,
    reltol::Float64 = 1e-4,
    verbose::Bool = false,
)
    n_params = length(ModelingToolkit.parameters(model))
    @assert length(bounds) == n_params (
        "Bounds length mismatch: got $(length(bounds)), expected $n_params parameters"
    )
    @assert length(ic) == length(ModelingToolkit.unknowns(model)) (
        "Initial conditions length mismatch: got $(length(ic)), " *
        "expected $(length(ModelingToolkit.unknowns(model)))"
    )
    @assert n_candidates > 0 "n_candidates must be positive"
    @assert 0.0 <= margin < 0.5 "margin must be in [0.0, 0.5)"
    @assert all(b -> b[2] > b[1], bounds) "All bounds must have lo < hi"

    # Generate candidates
    candidates = _generate_candidates(bounds, n_candidates, margin, sampling)

    # Filter through is_trajectory_bounded
    valid = Vector{Vector{Float64}}()
    rejected = Vector{RejectedCandidate}()
    amplitudes = Float64[]

    t_start = time()

    for (i, p) in enumerate(candidates)
        result = is_trajectory_bounded(
            model,
            ic,
            p,
            time_interval;
            threshold = threshold,
            numpoints = numpoints,
            solver = solver,
            abstol = abstol,
            reltol = reltol,
        )

        if result.bounded
            push!(valid, p)
            push!(amplitudes, result.max_amplitude)
        else
            reason = result.solver_success ? :amplitude_exceeded : :solver_failure
            push!(rejected, RejectedCandidate(p, reason, result.max_amplitude))
        end

        if verbose && (i % max(1, n_candidates ÷ 10) == 0 || i == length(candidates))
            n_done = i
            n_valid = length(valid)
            rate = n_valid / n_done
            elapsed = time() - t_start
            ms_per = elapsed / n_done * 1000
            println(
                "  [$n_done/$(length(candidates))] " *
                "$n_valid valid ($(round(100*rate, digits=1))%) " *
                "$(round(ms_per, digits=1)) ms/candidate",
            )
        end
    end

    elapsed = time() - t_start
    n_total = length(candidates)
    n_valid = length(valid)
    n_solver_fail = count(r -> r.reason == :solver_failure, rejected)
    n_amplitude = count(r -> r.reason == :amplitude_exceeded, rejected)

    diagnostics = SweepDiagnostics(
        elapsed,
        n_total,
        n_valid,
        n_solver_fail,
        n_amplitude,
        elapsed / n_total * 1000,
        isempty(amplitudes) ? NaN : minimum(amplitudes),
        isempty(amplitudes) ? NaN : maximum(amplitudes),
        isempty(amplitudes) ? NaN : sort(amplitudes)[max(1, length(amplitudes) ÷ 2)],
    )

    if verbose
        println()
        println(
            "  Sweep complete: $n_valid/$n_total passed ($(round(100 * n_valid/n_total, digits=1))%)",
        )
        println(
            "  Rejected: $n_solver_fail solver failures, $n_amplitude amplitude exceeded",
        )
        println(
            "  Time: $(round(elapsed, digits=2))s ($(round(diagnostics.ms_per_candidate, digits=1)) ms/candidate)",
        )
        if !isempty(amplitudes)
            println(
                "  Amplitude range: $(round(diagnostics.amplitude_min, digits=2)) — $(round(diagnostics.amplitude_max, digits=2))",
            )
        end
    end

    pass_rate = n_total > 0 ? n_valid / n_total : 0.0

    return SweepResult(valid, rejected, pass_rate, diagnostics)
end

# ============================================================================
# Landscape probing
# ============================================================================

"""
    probe_landscape(model, outputs, ic, p_true, bounds, time_interval; kwargs...) -> ProbeResult

Lightweight landscape probe: construct an objective function at coarse tolerances
and evaluate at random points to quickly estimate landscape structure.

This is a fast filter (~100ms per candidate at coarse tolerance) for identifying
potentially interesting p_true values before expensive full-precision evaluation.

# Arguments
- `model::ModelingToolkit.AbstractSystem`: ODE system (compiled via @mtkcompile)
- `outputs::Vector{ModelingToolkit.Equation}`: Measured quantities
- `ic::Vector{Float64}`: Initial conditions for all state variables
- `p_true::Vector{Float64}`: Candidate parameter vector (must pass is_trajectory_bounded)
- `bounds::Vector{Tuple{Float64,Float64}}`: Per-parameter (lo, hi) bounds for probe points
- `time_interval`: `[t_start, t_end]` for simulation

# Keyword Arguments
- `n_probes::Int = 10`: Number of random probe points to evaluate
- `numpoints::Int = 30`: Time points for ODE solve (used in make_error_distance)
- `distance_function = L2_norm`: Distance metric for objective
- `solver = Tsit5()`: ODE solver (cheap default for fast probing)
- `abstol::Float64 = 1e-4`: Absolute tolerance (coarse for speed)
- `reltol::Float64 = 1e-4`: Relative tolerance (coarse for speed)

# Returns
`ProbeResult` with fields:
- `objective_at_true::Float64`: Objective at p_true (should be ~0, sanity check)
- `probe_values::Vector{Float64}`: Objective values at random probe points
- `probe_points::Vector{Vector{Float64}}`: The random points used
- `dynamic_range::Float64`: log10(max/min) of finite nonzero probe values (NaN if <2 finite)
- `variance::Float64`: Variance of finite probe values (NaN if <2 finite)
- `n_finite::Int`: Number of probes that returned finite values
- `n_inf::Int`: Number of probes that returned Inf (ODE failure regions)
- `construction_time_ms::Float64`: Time to build the objective function

# Examples
```julia
model, params, states, outputs = define_daisy_ex3_model_4D()
ic = [1.0, 0.5, 0.3, 0.0]
p_true = [0.5, 0.3, 0.4, 0.2]
bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)]

result = probe_landscape(model, outputs, ic, p_true, bounds, [0.0, 20.0])
println("Objective at true: \$(result.objective_at_true)")
println("Dynamic range: \$(result.dynamic_range)")
println("Finite probes: \$(result.n_finite)/\$(length(result.probe_values))")
```
"""
function probe_landscape(
    model::ModelingToolkit.AbstractSystem,
    outputs::Vector{ModelingToolkit.Equation},
    ic::Vector{Float64},
    p_true::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}},
    time_interval;
    n_probes::Int = 10,
    numpoints::Int = 30,
    distance_function = L2_norm,
    solver = Tsit5(),
    abstol::Float64 = 1e-4,
    reltol::Float64 = 1e-4,
)
    n_params = length(ModelingToolkit.parameters(model))
    @assert length(bounds) == n_params (
        "Bounds length mismatch: got $(length(bounds)), expected $n_params parameters"
    )
    @assert length(p_true) == n_params (
        "p_true length mismatch: got $(length(p_true)), expected $n_params parameters"
    )
    @assert n_probes > 0 "n_probes must be positive"

    # Construct objective at coarse tolerances
    t_construct_start = time()
    local error_func
    try
        error_func = make_error_distance(
            model,
            outputs,
            ic,
            p_true,
            time_interval,
            numpoints,
            distance_function;
            return_inf_on_error = true,
            solver = solver,
            abstol = abstol,
            reltol = reltol,
        )
    catch e
        if e isa InterruptException
            rethrow(e)
        end
        # Construction failed — model/p_true combination cannot produce reference data
        error("probe_landscape: failed to construct objective for p_true=$p_true: $e")
    end
    construction_time_ms = (time() - t_construct_start) * 1000

    # Evaluate at p_true (sanity check: should be ~0)
    objective_at_true = error_func(p_true)

    # Generate random probe points within bounds
    dim = length(bounds)
    probe_points = [
        [bounds[d][1] + rand() * (bounds[d][2] - bounds[d][1]) for d in 1:dim] for
        _ in 1:n_probes
    ]

    # Evaluate at each probe point
    probe_values = Float64[]
    for p in probe_points
        val = error_func(p)
        push!(probe_values, val)
    end

    # Compute summary statistics from finite probe values
    finite_vals = filter(v -> isfinite(v) && v > 0, probe_values)
    n_finite = count(isfinite, probe_values)
    n_inf = count(v -> !isfinite(v), probe_values)

    dynamic_range = if length(finite_vals) >= 2
        log10(maximum(finite_vals) / minimum(finite_vals))
    else
        NaN
    end

    variance = if length(finite_vals) >= 2
        mean_val = sum(finite_vals) / length(finite_vals)
        sum((v - mean_val)^2 for v in finite_vals) / (length(finite_vals) - 1)
    else
        NaN
    end

    return ProbeResult(
        objective_at_true,
        probe_values,
        probe_points,
        dynamic_range,
        variance,
        n_finite,
        n_inf,
        construction_time_ms,
    )
end

# ============================================================================
# Ranking
# ============================================================================

"""
    rank_probes(probes; min_finite_fraction, max_noise_ratio) -> RankingResult

Rank a vector of probe results by landscape interestingness.

Applies quality filters, then sorts by dynamic_range (descending) with
variance as tiebreaker.

# Arguments
- `probes::Vector{<:ProbeResult}`: Results from `probe_landscape`, one per candidate

# Keyword Arguments
- `min_finite_fraction::Float64 = 0.5`: Reject candidates where fewer than this
  fraction of probes returned finite values
- `max_noise_ratio::Float64 = Inf`: Reject candidates where `objective_at_true`
  exceeds `max_noise_ratio * median(finite_probes)`. Set to Inf to disable.

# Returns
`RankingResult` with fields:
- `ranked_indices::Vector{Int}`: Indices into `probes`, sorted best-first
- `filtered_out::Vector{Int}`: Indices that failed quality filters
- `filter_reasons::Vector{Symbol}`: Reason for each filtered-out index
  (:low_finite_fraction or :high_noise_floor)

# Ranking criteria
1. **Filter**: n_finite >= min_finite_fraction * n_probes AND noise check
2. **Sort**: dynamic_range descending (higher = more landscape structure)
3. **Tiebreak**: variance descending (more uniformly varied > one spike)

Candidates with NaN dynamic_range (< 2 finite probes) are sorted to the end
of the ranked list rather than filtered out — they pass the finite fraction
check but lack enough data for a reliable dynamic_range.

# Examples
```julia
probes = [probe_landscape(model, outputs, ic, p, bounds, tspan) for p in valid_pts]
ranking = rank_probes(probes)
best_idx = ranking.ranked_indices[1]
println("Best candidate: p=\$(valid_pts[best_idx])")
println("Dynamic range: \$(probes[best_idx].dynamic_range)")
```
"""
function rank_probes(
    probes::AbstractVector{<:ProbeResult};
    min_finite_fraction::Float64 = 0.5,
    max_noise_ratio::Float64 = Inf,
)
    @assert 0.0 <= min_finite_fraction <= 1.0 "min_finite_fraction must be in [0, 1]"
    @assert max_noise_ratio > 0.0 "max_noise_ratio must be positive"

    ranked_indices = Int[]
    filtered_out = Int[]
    filter_reasons = Symbol[]

    for (i, probe) in enumerate(probes)
        n_total = probe.n_finite + probe.n_inf

        # Filter 1: sufficient finite probes
        if n_total > 0 && probe.n_finite / n_total < min_finite_fraction
            push!(filtered_out, i)
            push!(filter_reasons, :low_finite_fraction)
            continue
        end

        # Filter 2: noise floor check
        if isfinite(max_noise_ratio) && probe.n_finite >= 1
            finite_vals = sort(filter(v -> isfinite(v) && v > 0, probe.probe_values))
            if !isempty(finite_vals)
                median_val = finite_vals[max(1, length(finite_vals) ÷ 2)]
                if median_val > 0 && probe.objective_at_true > max_noise_ratio * median_val
                    push!(filtered_out, i)
                    push!(filter_reasons, :high_noise_floor)
                    continue
                end
            end
        end

        push!(ranked_indices, i)
    end

    # Sort: dynamic_range descending, variance as tiebreaker
    sort!(
        ranked_indices,
        by = i -> begin
            dr = probes[i].dynamic_range
            v = probes[i].variance
            # NaN sorts to end: use -Inf as sentinel
            (isnan(dr) ? -Inf : dr, isnan(v) ? -Inf : v)
        end,
        rev = true,
    )

    return RankingResult(ranked_indices, filtered_out, filter_reasons)
end

# ============================================================================
# Batch orchestrator
# ============================================================================

"""
    screen_and_probe(model, outputs, ic, bounds, time_interval; kwargs...) -> ScreeningResult

Batch pipeline: generate candidates -> screen for bounded trajectories ->
probe landscape structure -> rank by interestingness.

Composes `sweep_p_true`, `probe_landscape`, and `rank_probes` into a single
call. For more control, call these functions individually.

# Arguments
- `model::ModelingToolkit.AbstractSystem`: ODE system (compiled via @mtkcompile)
- `outputs::Vector{ModelingToolkit.Equation}`: Measured quantities
- `ic::Vector{Float64}`: Initial conditions for all state variables
- `bounds::Vector{Tuple{Float64,Float64}}`: Per-parameter (lo, hi) bounds
- `time_interval`: `[t_start, t_end]` for simulation

# Keyword Arguments

## Sweep parameters (passed to `sweep_p_true`)
- `n_candidates::Int = 1000`: Number of candidate p_true vectors to generate
- `margin::Float64 = 0.1`: Safety margin from bounds for candidate generation
- `sampling::Symbol = :random`: Sampling strategy (:random or :grid)
- `threshold::Float64 = 1e6`: Max amplitude for trajectory boundedness

## Probe parameters (passed to `probe_landscape`)
- `n_probes::Int = 10`: Number of random probe points per candidate
- `numpoints_screen::Int = 20`: Time points for screening ODE solve
- `numpoints_probe::Int = 30`: Time points for probing ODE solve
- `distance_function = L2_norm`: Distance metric for objective
- `solver = Tsit5()`: ODE solver for both screening and probing
- `abstol::Float64 = 1e-4`: Absolute tolerance
- `reltol::Float64 = 1e-4`: Relative tolerance

## Ranking parameters (passed to `rank_probes`)
- `min_finite_fraction::Float64 = 0.5`: Min fraction of finite probes to pass filter
- `max_noise_ratio::Float64 = Inf`: Max objective_at_true / median(probes) ratio

## Output control
- `verbose::Bool = false`: Print progress and summary

# Returns
`ScreeningResult` with fields:
- `sweep::SweepResult`: Raw output from `sweep_p_true`
- `probes::Vector{ProbeResult}`: Probe results, parallel to `sweep.valid`
- `ranking::RankingResult`: Output from `rank_probes` (ranked_indices, filtered_out, filter_reasons)

# Examples
```julia
model, params, states, outputs = define_daisy_ex3_model_4D()
ic = [1.0, 0.5, 0.3, 0.0]
bounds = [(0.0, 2.0), (0.0, 1.0), (0.0, 1.0), (-0.5, 0.5)]

result = screen_and_probe(model, outputs, ic, bounds, [0.0, 20.0];
    n_candidates=100, n_probes=10, verbose=true)

# Top 5 most interesting candidates
for i in result.ranking.ranked_indices[1:min(5, end)]
    p = result.sweep.valid[i]
    dr = result.probes[i].dynamic_range
    println("p=\$(round.(p, digits=3))  dynamic_range=\$(round(dr, digits=2))")
end
```
"""
function screen_and_probe(
    model::ModelingToolkit.AbstractSystem,
    outputs::Vector{ModelingToolkit.Equation},
    ic::Vector{Float64},
    bounds::Vector{Tuple{Float64,Float64}},
    time_interval;
    # Sweep parameters
    n_candidates::Int = 1000,
    margin::Float64 = 0.1,
    sampling::Symbol = :random,
    threshold::Float64 = 1e6,
    # Probe parameters
    n_probes::Int = 10,
    numpoints_screen::Int = 20,
    numpoints_probe::Int = 30,
    distance_function = L2_norm,
    solver = Tsit5(),
    abstol::Float64 = 1e-4,
    reltol::Float64 = 1e-4,
    # Ranking parameters
    min_finite_fraction::Float64 = 0.5,
    max_noise_ratio::Float64 = Inf,
    # Output control
    verbose::Bool = false,
)
    # Phase 1: Sweep
    if verbose
        println("Phase 1: Screening $n_candidates candidates...")
    end
    sweep = sweep_p_true(
        model,
        ic,
        bounds,
        time_interval;
        n_candidates = n_candidates,
        margin = margin,
        sampling = sampling,
        threshold = threshold,
        numpoints = numpoints_screen,
        solver = solver,
        abstol = abstol,
        reltol = reltol,
        verbose = verbose,
    )

    n_valid = length(sweep.valid)
    if verbose
        println()
        println("Phase 2: Probing $n_valid valid candidates ($n_probes probes each)...")
    end

    # Phase 2: Probe each valid candidate
    probes = Vector{ProbeResult}(undef, n_valid)
    t_probe_start = time()
    for (i, p) in enumerate(sweep.valid)
        probes[i] = probe_landscape(
            model,
            outputs,
            ic,
            p,
            bounds,
            time_interval;
            n_probes = n_probes,
            numpoints = numpoints_probe,
            distance_function = distance_function,
            solver = solver,
            abstol = abstol,
            reltol = reltol,
        )

        if verbose && (i % max(1, n_valid ÷ 10) == 0 || i == n_valid)
            elapsed = time() - t_probe_start
            ms_per = elapsed / i * 1000
            println("  [$i/$n_valid] $(round(ms_per, digits=1)) ms/candidate")
        end
    end
    t_probe_total = time() - t_probe_start

    # Phase 3: Rank
    ranking = rank_probes(
        probes;
        min_finite_fraction = min_finite_fraction,
        max_noise_ratio = max_noise_ratio,
    )

    if verbose
        n_ranked = length(ranking.ranked_indices)
        n_filtered = length(ranking.filtered_out)
        println()
        println("Phase 3: Ranking")
        println("  $n_ranked candidates ranked, $n_filtered filtered out")
        if n_ranked > 0
            best = ranking.ranked_indices[1]
            println("  Best dynamic_range: $(round(probes[best].dynamic_range, digits=2))")
            println(
                "  Probe time: $(round(t_probe_total, digits=2))s " *
                "($(round(t_probe_total/n_valid*1000, digits=1)) ms/candidate)",
            )
        end
        println(
            "  Total time: $(round(sweep.diagnostics.elapsed_seconds + t_probe_total, digits=2))s",
        )
    end

    return ScreeningResult(sweep, probes, ranking)
end
