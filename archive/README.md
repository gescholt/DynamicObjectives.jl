# Archived Documentation

This directory contains documentation that has been superseded or completed. Files are preserved here for historical reference.

## Directory Structure

### completed_planning/
Documents from planning and migration phases that are now complete.

- **CLEANUP_PLAN.md**: Initial repository cleanup plan (Nov 2025)
  - **Status**: COMPLETED ✅
  - **Archived**: 2025-11-23
  - **Reason**: Cleanup objectives achieved; repository is clean

- **MIGRATION.md**: DynamicalSystems.jl migration guide
  - **Status**: COMPLETED ✅
  - **Archived**: 2025-11-23
  - **Reason**: All code migrated from include() pattern to proper package

- **INTEGRATION_OBJECTIVES.md**: Pre-Phase 3 integration objectives
  - **Status**: ACHIEVED ✅
  - **Archived**: 2025-11-23
  - **Reason**: Phase 3 integration complete; objectives met

- **ORGANIZATIONAL_REVIEW.md**: Documentation reorganization plan
  - **Status**: COMPLETED ✅
  - **Archived**: 2025-11-27
  - **Reason**: All recommended reorganization completed

### superseded_docs/
Documentation superseded by updated versions.

- **GLOBTIM_INTEGRATION.md**: Pre-Phase 2 integration guide
  - **Superseded by**: INTEGRATION.md
  - **Archived**: 2025-11-23
  - **Reason**: References deprecated wrapper functions; Phase 2 architecture uses direct 1-arg functions
  - **Note**: Some content may be useful for historical reference

## What to Use Instead

| Archived File | Current Alternative |
|--------------|---------------------|
| CLEANUP_PLAN.md | (No replacement needed - work complete) |
| MIGRATION.md | (No replacement needed - migration complete) |
| INTEGRATION_OBJECTIVES.md | ARCHITECTURE.md, TESTING_GUIDE.md |
| docs/GLOBTIM_INTEGRATION.md | INTEGRATION.md |

## Restoration

If you need to restore any archived file:
```bash
# From the archive/ directory
cp completed_planning/FILENAME.md ../FILENAME.md
```

## Maintenance Policy

**Never delete - always archive**
- Preserves institutional knowledge
- Maintains historical context
- Allows easy restoration if needed

**When to archive:**
- Planning documents after completion
- Documentation after major refactoring
- Guides that reference deprecated APIs
- Examples using old workflow patterns

---

**Last Updated**: 2025-11-27
**Archived during**: Phase 3 post-integration cleanup
