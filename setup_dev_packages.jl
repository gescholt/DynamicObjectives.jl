#!/usr/bin/env julia
# Setup script for DynamicObjectives experiment pipeline
#
# Adds Globtim and GlobtimPostProcessing as development dependencies.
# These packages must be available either as:
#   - Sibling directories (../globtim, ../globtimpostprocessing) in a monorepo layout
#   - Already installed in your Julia environment
#
# Usage:
#   julia setup_dev_packages.jl

using Pkg

println("="^72)
println("DynamicObjectives — Package Setup")
println("="^72)
println()

Pkg.activate(@__DIR__)

# ── Locate and add Globtim ────────────────────────────────────────────────

globtim_path = joinpath(@__DIR__, "..", "globtim")
if isdir(globtim_path)
    println("Found local Globtim at: $globtim_path")
    Pkg.develop(path = globtim_path)
else
    error(
        "Globtim not found at $globtim_path\n" *
        "Either:\n" *
        "  1. Clone Globtim as a sibling directory: git clone <globtim-url> ../globtim\n" *
        "  2. Or add it manually: julia -e 'using Pkg; Pkg.add(url=\"<globtim-url>\")'",
    )
end

# ── Locate and add GlobtimPostProcessing ──────────────────────────────────

gpp_path = joinpath(@__DIR__, "..", "globtimpostprocessing")
if isdir(gpp_path)
    println("Found local GlobtimPostProcessing at: $gpp_path")
    Pkg.develop(path = gpp_path)
else
    error(
        "GlobtimPostProcessing not found at $gpp_path\n" *
        "Either:\n" *
        "  1. Clone as a sibling directory: git clone <gpp-url> ../globtimpostprocessing\n" *
        "  2. Or add it manually: julia -e 'using Pkg; Pkg.add(url=\"<gpp-url>\")'",
    )
end

# ── Resolve and verify ────────────────────────────────────────────────────

println()
println("Resolving dependencies...")
Pkg.resolve()
println("Precompiling...")
Pkg.precompile()

println()
println("Verifying packages load...")
@eval using Globtim
println("  Globtim: OK")
@eval using GlobtimPostProcessing
println("  GlobtimPostProcessing: OK")

println()
println("="^72)
println("Setup complete. Run experiments with:")
println("  julia --project=. paper/run_experiment.jl paper/configs/deuflhard_2d.toml")
println("="^72)
