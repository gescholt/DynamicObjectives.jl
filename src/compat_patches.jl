# Compat patches for DataStructures / OrdinaryDiffEqCore version skew (vemc, 2026-05).
#
# In the current Manifest, OrdinaryDiffEqCore.reinit! invokes
# `empty!(::DataStructures.BinaryHeap)` but DataStructures.jl no longer ships a
# matching method. BinaryHeap stores its payload as `.valtree::Vector`, so the
# fix is a one-line forward: empty the underlying vector.
#
# Validated by experiments/sandbox/spike_integrator_reuse.jl. Remove this file
# once either DataStructures.jl restores the method or OrdinaryDiffEqCore stops
# calling it.

@assert isdefined(DataStructures, :BinaryHeap) "DataStructures.BinaryHeap not defined; remove this patch"

if !hasmethod(Base.empty!, Tuple{DataStructures.BinaryHeap})
    Base.empty!(h::DataStructures.BinaryHeap) = (empty!(h.valtree); h)
end
