# Repository Cleanup Plan

## Current Status: CLEAN ✅

The Dynamic_objectives repository is already well-organized with no significant orphaned files.

## Repository Structure (Current)

```
Dynamic_objectives/
├── .git/                           # Git metadata
├── .gitignore                      # Git ignore rules
├── Project.toml                    # ✅ Package manifest
├── README.md                       # ✅ Main documentation
├── MIGRATION.md                    # ✅ Migration guide from DynamicalSystems.jl
├── TESTING_GUIDE.md                # ✅ Comprehensive testing guide
├── MODEL_CATALOG.md                # ✅ Quick reference for all models
├── src/
│   ├── Dynamic_objectives.jl       # ✅ Main module
│   ├── data_generation.jl          # ✅ sample_data function
│   ├── error_metrics.jl            # ✅ make_error_distance + distance functions
│   ├── support.jl                  # ⚠️  EMPTY PLACEHOLDER - candidate for removal
│   └── systems/
│       ├── daisy_models.jl         # ✅ DAISY benchmarks (2 models)
│       ├── lotka_volterra.jl       # ✅ LV variants (8 models)
│       └── other_systems.jl        # ✅ FitzHugh-Nagumo + identifiability tests (3 models)
└── test/
    ├── runtests.jl                 # ✅ Basic test suite
    └── validate_all_models.jl      # ✅ NEW - Comprehensive validation script
```

## Identified Issues

### Minor Issues

1. **Empty Placeholder File**
   - `src/support.jl` - Only contains comments, no actual code
   - Included in `src/Dynamic_objectives.jl` line 16
   - **Action:** Remove if not needed, or add a TODO comment explaining future purpose

## Recommended Actions

### Option 1: Remove Empty Placeholder (Recommended)

If you don't need `support.jl` immediately:

1. Remove the include statement from `src/Dynamic_objectives.jl`
2. Delete `src/support.jl`
3. Document in comments where to add support structures if needed later

### Option 2: Keep Placeholder with Better Documentation

If you plan to add content later:

1. Add a clear TODO comment explaining what will go there
2. Add a note in README.md about planned features
3. Keep the file but improve the placeholder comment

## Proposed Changes

### 1. Handle support.jl

**Remove it:**
```diff
--- a/src/Dynamic_objectives.jl
+++ b/src/Dynamic_objectives.jl
@@ -13,7 +13,6 @@ include("systems/other_systems.jl")
 include("data_generation.jl")
 include("error_metrics.jl")
-include("support.jl")

 # System definitions - DAISY models
```

Then delete the file.

**OR keep it with better docs:**
```julia
# src/support.jl
# Support Structures
# Additional support types and utilities

# TODO: Future additions planned:
# - EllipseSupport: For constrained parameter spaces (needed for ellipsoid constraints in globtim)
# - BoundedDomain: For domain-specific parameter bounds
# - ParameterTransforms: For log-space or other transformations
#
# Current status: Not yet implemented
# Leave empty for now - no exports from this file
```

### 2. Add CHANGELOG.md (Optional)

Track changes over time:
```bash
# Create CHANGELOG.md
echo "# Changelog

## [0.1.0] - 2025-11-19

### Added
- Initial release with 13 ODE benchmark models
- Comprehensive testing documentation (TESTING_GUIDE.md)
- Model catalog (MODEL_CATALOG.md)
- Migration guide from DynamicalSystems.jl
- Timeout mechanism for robust optimization
- Three distance metrics: L1, L2, log-L2

### Documentation
- Complete testing strategy for globtim
- Recommended configurations for all models
- Four-phase testing campaign
" > CHANGELOG.md
```

### 3. Add LICENSE File (Recommended)

Currently says "MIT License" in README but no LICENSE file:
```bash
# Add MIT License file
# (Copy standard MIT license with your details)
```

### 4. Add .github/ Directory (Optional - For Future)

If you plan to use GitHub Actions or issue templates:
```
.github/
├── workflows/
│   └── test.yml           # CI/CD for running tests
└── ISSUE_TEMPLATE/
    └── bug_report.md      # Template for bug reports
```

## Files That Are Fine As-Is

These files are all necessary and well-organized:

- ✅ All source files in `src/` (except support.jl)
- ✅ All documentation files (README, MIGRATION, TESTING_GUIDE, MODEL_CATALOG)
- ✅ Test files
- ✅ Project.toml
- ✅ .gitignore

## Migration/Archive Strategy

Since there are no legacy files to archive, this section focuses on future maintenance:

### For Future Deprecated Code

Create an `archive/` directory structure:
```
archive/
├── README.md                    # Explains what's archived and why
└── YYYY-MM-DD_description/      # Date-stamped directories
    └── [archived files]
```

Example:
```
archive/
├── README.md
└── 2025-11-19_old_daisy_models/
    └── deprecated_daisy_v1.jl
```

### Archive Guidelines

1. **Never delete** - move to `archive/` instead
2. **Date-stamp** archive directories: `YYYY-MM-DD_description/`
3. **Document** in archive/README.md:
   - What was archived
   - Why it was archived
   - When it was archived
   - Where new code lives (if replaced)
4. **Add to .gitignore** if you don't want archives in version control

## Next Steps

### Immediate (Do Now)

1. ✅ Run validation script: `julia --project=. test/validate_all_models.jl`
2. ✅ Run basic tests: `julia --project=. -e 'using Pkg; Pkg.test()'`
3. ⚠️  Decide on `support.jl`: Remove or improve documentation

### Short Term (This Week)

1. Add LICENSE file
2. Add CHANGELOG.md
3. Consider adding examples/ directory with simple usage examples
4. Test with globtim on at least 3-5 models

### Long Term (As Needed)

1. Add CI/CD if moving to GitHub/GitLab
2. Set up documentation generation with Documenter.jl
3. Add benchmarking suite for performance tracking
4. Consider publishing to Julia General registry

## Conclusion

**Current State:** Clean repository, minimal cleanup needed

**Priority Actions:**
1. Remove or document `support.jl` properly
2. Run validation scripts to ensure everything works
3. Start testing campaign with globtim

The repository is in excellent shape for your testing campaign!
