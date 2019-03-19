using JuLIP: AbstractCalculator
using JuLIP.Potentials: @pot
import Base: hash

# ----------- Abstract Supertype for pure NBodyFunctions --------------

"""
an artificial type to allow dispatch but our different packages
don't need to know about it => this avoids a superpackage of abstract
type definitions.

This `::BASIS` is e.g. used to create hashs of basis functions that
don't include the parameters but only the basis function class
"""
const BASIS = Val{:basis}


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


abstract type SpaceTransform end

function transform end
function transform_d end

include("transforms.jl")


abstract type NBCutoff end

# prototypes for cutoffs
function fcut end
function fcut_d end

"""
`hash(::BASIS, b) -> UInt`: provides a description of the basis function
`b` such that, two basis functions `b1, b2` can be combined into one
if and only if their "basis-hash" match.
"""
hash(::BASIS, C::NBCutoff) = hash(C)

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
export NBodyIP

export BondAngleDesc

struct BondAngleDesc{TT <: SpaceTransform, TC <: NBCutoff} <: NBSiteDescriptor
   transform::TT
   cutoff::TC
end


export BondLengthDesc

struct BondLengthDesc{TT <: SpaceTransform, TC <: NBCutoff} <: NBSiteDescriptor
   transform::TT
   cutoff::TC
end


export ClusterBLDesc

struct ClusterBLDesc{TT <: SpaceTransform, TC <: NBCutoff} <: NBClusterDescriptor
   transform::TT
   cutoff::TC
end

transform(D::NBodyDescriptor) = D.transform
transform(V::NBodyFunction) = transform(descriptor(V))
