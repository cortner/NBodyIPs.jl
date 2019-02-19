using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot

# ----------- Abstract Supertype for pure NBodyFunctions --------------

"""
`NBodyFunction` : abstract supertype of all "pure" N-body functions.
concrete subtypes must implement

* `bodyorder`
* `evaluate`
* `evaluate_d`
"""
abstract type NBodyFunction{N, DT} <: AbstractCalculator end

"""
`NBodyDescriptor`: abstract supertype for different descriptors
of N-body configurations.
"""
abstract type NBodyDescriptor end

"""
`NullDesc` : a descriptor that contains no information => used for
subtyping when an NBodyFunction subtype does not have a descriptor
(traits would be nice right now)
"""
struct NullDesc <: NBodyDescriptor end

"""
`NBSiteDescriptor`: abstract supertype for descriptors that start from
a site-based formulation. I.e., the underlying symmetry is with the
site fixed.
"""
abstract type NBSiteDescriptor <: NBodyDescriptor end

struct NullSiteDesc <: NBSiteDescriptor end

"""
`NBClusterDescriptor`: abstract supertype for descriptors that start from
a cluster-based formulation. I.e. the underlying symmetry allows
permutation of arbitrary atoms within a cluster.

WARNING: For now, the `NBClusterDescriptor` is treated as a special
case of an `NBSiteDescriptor` with some suitable hacks to skip
repeated simplices. For the future we need a "bespoke" neighbourlist
loop for this case.
"""
abstract type NBClusterDescriptor <: NBSiteDescriptor end

struct SpaceTransform{FT, FDT, VT}
   id::String
   f::FT
   f_d::FDT
   valid::VT
end

SpaceTransform(id, f, f_d) = SpaceTransform(id, f, f_d, Val(Symbol(id)))


struct Cutoff{FT, DFT}
   sym::Symbol
   params::Vector{Float64}
   f::FT
   f_d::DFT
   rcut::Float64
end


"""
`NBodyIP` : wraps `NBodyFunction`s or similar into a JuLIP calculator, defining
* `site_energies`
* `energy`
* `forces`
* `virial`
* `cutoff`

Use `load_ip` to load from a file (normally `jld2` or `json`)
"""
mutable struct NBodyIP <: AbstractCalculator
   components::Vector{AbstractCalculator}
end

export BondAngleDesc

struct BondAngleDesc{TT <: SpaceTransform, TC <: Cutoff} <: NBSiteDescriptor
   transform::TT
   cutoff::TC
end


export BondLengthDesc

struct BondLengthDesc{TT <: SpaceTransform, TC <: Cutoff} <: NBSiteDescriptor
   transform::TT
   cutoff::TC
end


export ClusterBLDesc

struct ClusterBLDesc{TT <: SpaceTransform, TC <: Cutoff} <: NBClusterDescriptor
   transform::TT
   cutoff::TC
end


"""
`mutable struct OneBody{T}  <: NBodyFunction{1}`

this should not normally be constructed by a user, but instead E0 should be
passed to the relevant lsq functions, which will construct it.
"""
mutable struct OneBody{T} <: NBodyFunction{1, NullDesc}
   E0::T
end

@pot OneBody
