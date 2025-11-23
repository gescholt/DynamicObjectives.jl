# globtim Integration Objectives

## Overview

This document outlines the objectives for integrating the globtim optimizer with the Dynamic_objectives benchmark suite for systematic testing of time parameter estimation problems.

## Primary Objectives

### 1. Validate globtim Performance on Parameter Estimation

**Goal:** Test globtim's ability to recover true parameters from time series data across varying problem difficulties.

**Success Metrics:**
- Parameter recovery accuracy: `||p_best - p_true|| < 0.01` for identifiable problems
- Convergence success rate: >90% on EASY problems, >70% on MEDIUM, >50% on HARD
- Objective value at solution: `error_func(p_best) < 1e-6` for perfect recovery

### 2. Characterize Performance vs. Problem Characteristics

**Goal:** Understand how globtim performs across different problem types.

**Dimensions to Test:**
- **Parameter count**: 1D → 2D → 3D → 4D → 20D
- **Observability**: Partial (1 output) vs. Full (all states measured)
- **System type**: Linear vs. Nonlinear vs. Oscillatory
- **Identifiability**: Globally identifiable vs. locally identifiable

**Expected Insights:**
- Scaling behavior with dimensionality
- Impact of partial observability on convergence
- Handling of oscillatory dynamics (FitzHugh-Nagumo)
- Behavior on non-identifiable problems (should find equivalent parameters)

### 3. Benchmark Function Evaluation Efficiency

**Goal:** Measure computational cost of parameter estimation problems.

**Metrics to Track:**
- Function evaluations to convergence
- Wall-clock time per evaluation
- Total optimization time
- Timeout events (for hard problems)

**Target Efficiency:**
- EASY problems: <500 evaluations
- MEDIUM problems: <2000 evaluations
- HARD problems: <5000 evaluations
- Minimize timeout events

### 4. Test Robustness to Problem Configuration

**Goal:** Verify globtim handles various objective function configurations.

**Configurations to Test:**
- **Distance metrics**: L2 (default), L1 (robust), log-L2 (wide-range)
- **Timeout settings**: None, 5s, 10s, 15s
- **Time series length**: 20, 50, 100 points
- **Time horizons**: Short (5s), Medium (20s), Long (50s)

**Expected Outcomes:**
- Consistent convergence across distance metrics
- Graceful handling of timeouts (return Inf)
- Performance-accuracy tradeoffs with time series length

### 5. Identify Problem Difficulty Ladder

**Goal:** Create a validated difficulty ordering of benchmark problems.

**Phases (as hypothesized):**
1. **EASY**: LV 2D models (2 params, well-conditioned)
2. **MEDIUM**: LV 3D, DAISY 4D (3-4 params, moderate)
3. **HARD**: FitzHugh-Nagumo, identifiability tests (oscillatory, non-identifiable)
4. **VERY HARD**: 20D Generalized LV (high-dimensional)

**Validation:**
- Confirm ordering via success rates and evaluation counts
- Identify any surprises (e.g., supposedly easy problems that are actually hard)

## Secondary Objectives

### 6. Establish Best Practices for globtim Usage

**Goal:** Document recommended settings for parameter estimation problems.

**Questions to Answer:**
- What are optimal globtim settings for 2D vs. 4D vs. 20D problems?
- When should timeout be enabled?
- Which distance metric works best for different problem types?
- How many time points are sufficient for accurate estimation?

### 7. Create Performance Baseline for Future Comparisons

**Goal:** Generate reference results for comparing optimization algorithms.

**Deliverables:**
- Performance report for all 13 models
- Success rates by difficulty phase
- Function evaluation statistics
- Timing benchmarks

**Use Cases:**
- Compare globtim versions (e.g., v1.0 vs v2.0)
- Compare against other optimizers
- Track improvements over time

### 8. Identify globtim Strengths and Weaknesses

**Goal:** Understand where globtim excels and where it struggles.

**Analysis:**
- Best performance: Which problem classes?
- Struggles: Which problems fail to converge?
- Unique capabilities: What can globtim do that others can't?
- Limitations: Where does globtim need improvement?

## Technical Requirements

### Integration Points

**1. Objective Function Interface:**
```julia
# globtim receives:
error_func::Function           # objective(p) → scalar error
bounds::Vector{Tuple{Float64, Float64}}  # [(lb₁, ub₁), (lb₂, ub₂), ...]

# globtim returns:
result = (
    minimizer = p_best,        # best parameters found
    minimum = error_best,      # best error value
    converged = true/false,    # convergence flag
    iterations = n,            # iterations used
    f_calls = m               # function evaluations
)
```

**2. Expected Behavior:**
- Handle `Inf` returns gracefully (from timeout or ODE failures)
- Respect parameter bounds strictly
- Work with 1D to 20D problems
- Converge to global optimum (or equivalent for non-identifiable problems)

### Performance Expectations

**Baseline Requirements:**
- Reliably find true parameters on EASY problems (>95% success)
- Handle 20D problem without crashing (may not converge perfectly)
- Complete EASY problems in <1 minute each
- Complete MEDIUM problems in <5 minutes each

**Stretch Goals:**
- >50% success on HARD identifiability problems
- Find equivalent parameters on non-identifiable problems
- Competitive with state-of-the-art global optimizers

## Testing Workflow

### Phase 1: Single Model Validation (Week 1)
1. Integrate globtim into `quick_start.jl`
2. Test on 1-2 EASY models
3. Verify parameter recovery
4. Debug any integration issues

### Phase 2: Systematic Easy Testing (Week 2)
1. Test all 4 EASY models (2D LV variants)
2. Measure success rates and efficiency
3. Document best settings
4. Generate preliminary report

### Phase 3: Progressive Difficulty (Week 3)
1. Test MEDIUM models (3-4D)
2. Test HARD models (oscillatory, non-identifiable)
3. Attempt 20D stress test
4. Analyze scaling behavior

### Phase 4: Comprehensive Evaluation (Week 4)
1. Run full suite (all 13 models)
2. Test configuration variations
3. Generate final performance report
4. Document findings and recommendations

## Success Criteria

### Minimum Viable Success
- ✅ globtim successfully integrates with Dynamic_objectives
- ✅ Recovers parameters on ≥3 EASY models
- ✅ Completes execution without crashes
- ✅ Produces interpretable results

### Full Success
- ✅ >90% success rate on EASY models
- ✅ >70% success rate on MEDIUM models
- ✅ Handles all 13 models without crashes
- ✅ Complete performance characterization

### Exceptional Success
- ✅ >50% success on HARD non-identifiable models
- ✅ Finds equivalent parameters on non-identifiable problems
- ✅ Competitive evaluation counts vs. other optimizers
- ✅ Successfully handles 20D problem

## Deliverables

### Code
- ✅ Working integration in all example scripts
- ✅ Documented globtim wrapper function
- ✅ Configuration templates for different problem types

### Documentation
- ✅ Integration guide (how to use globtim with Dynamic_objectives)
- ✅ Performance report (results on all 13 models)
- ✅ Best practices document (recommended settings)

### Data
- ✅ Results files with detailed metrics
- ✅ Success rates by phase
- ✅ Timing benchmarks
- ✅ Convergence curves (if available)

## Open Questions

1. **What is globtim's optimization algorithm?**
   - Global search strategy?
   - Handles box constraints natively?
   - Gradient-free or gradient-based?

2. **What are globtim's configuration options?**
   - Max evaluations
   - Convergence tolerances
   - Population size (if evolutionary)
   - Other hyperparameters?

3. **What output does globtim provide?**
   - Just final solution?
   - Convergence history?
   - Multiple local optima?
   - Uncertainty estimates?

4. **Special features to leverage?**
   - Parallel evaluation?
   - Warm-start from previous runs?
   - Adaptive parameter tuning?
   - Built-in timeout handling?

## Risk Factors

### Technical Risks
- **Integration complexity**: globtim API may not match expected interface
- **Performance**: globtim may be too slow for ODE-based objectives
- **Scaling**: May not handle 20D problem effectively
- **Convergence**: May struggle with non-convex, multi-modal landscapes

### Mitigation Strategies
- Start with simplest models to validate integration
- Use timeout to prevent hanging on hard problems
- Lower expectations for 20D problem (stress test only)
- Document failures as learning opportunities

## Timeline Estimate

- **Week 1**: Integration and validation (1-2 EASY models)
- **Week 2**: Easy model testing (4 models)
- **Week 3**: Medium and hard testing (9 models)
- **Week 4**: Analysis and documentation

**Total**: ~4 weeks for comprehensive evaluation

## Next Actions

1. **Understand globtim API**: Review documentation/examples
2. **Create wrapper function**: Adapt globtim to expected interface
3. **Test integration**: Run `quick_start.jl` with globtim
4. **Validate on easy model**: Verify parameter recovery
5. **Scale up**: Progress through difficulty phases

---

*Last Updated: 2025-11-20*
*Package: Dynamic_objectives v0.1.0*
*Models: 13 time parameter estimation benchmarks*
