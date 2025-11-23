#!/usr/bin/env julia
"""
Setup script for Dynamic_objectives integration testing

This script sets up local development versions of globtimcore and
globtimpostprocessing so that integration tests can find them.

Usage:
    julia setup_dev_packages.jl

This only needs to be run ONCE. After that, the packages will always
be available in the Dynamic_objectives environment.
"""

using Pkg

println("=" ^ 80)
println("Dynamic_objectives - Integration Package Setup")
println("=" ^ 80)
println()

# Activate Dynamic_objectives environment
println("Activating Dynamic_objectives environment...")
Pkg.activate(@__DIR__)
println("✅ Environment activated")
println()

# Dev globtimcore
println("Setting up globtimcore (local dev version)...")
globtimcore_path = joinpath(@__DIR__, "..", "globtimcore")
if isdir(globtimcore_path)
    Pkg.develop(path=globtimcore_path)
    println("✅ globtimcore dev version added")
else
    @error "globtimcore not found at: $globtimcore_path"
    exit(1)
end
println()

# Dev globtimpostprocessing
println("Setting up globtimpostprocessing (local dev version)...")
globtimpost_path = joinpath(@__DIR__, "..", "globtimpostprocessing")
if isdir(globtimpost_path)
    Pkg.develop(path=globtimpost_path)
    println("✅ globtimpostprocessing dev version added")
else
    @error "globtimpostprocessing not found at: $globtimpost_path"
    exit(1)
end
println()

# Verify packages can be loaded
println("Verifying packages can be loaded...")
try
    using Globtim
    println("  ✅ Globtim loaded successfully")
catch e
    @error "Failed to load Globtim" exception=e
    exit(1)
end

try
    using GlobtimPostProcessing
    println("  ✅ GlobtimPostProcessing loaded successfully")
catch e
    @error "Failed to load GlobtimPostProcessing" exception=e
    exit(1)
end
println()

println("=" ^ 80)
println("✅ Setup Complete!")
println("=" ^ 80)
println()
println("You can now run integration tests:")
println("  ./run_tests.sh")
println("  julia --project=. -e 'using Pkg; Pkg.test()'")
println()
println("Or use refinement pipeline in your scripts:")
println("  using Globtim")
println("  using GlobtimPostProcessing")
println()
