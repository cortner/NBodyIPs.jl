
using StaticArrays

using JuLIP: AbstractCalculator,
             Atoms


import Base: Dict,
             ==

import JuLIP: cutoff,
              energy,
              forces,
              site_energies,
              virial

export NBodyIP,
       fast,
       load_ip,
       save_ip,
       energy,
       forces,
       site_energies,
       virial

# ----------- some auxiliaries we will use throughout -------------

"""
`push_str!(ex::Vector{Expr}, s::String)` : parse the string and
append it to the expressions.
"""
push_str!(ex::Vector{Expr}, s::String) = push!(ex, parse(s))


"""
`append_str!(ex::Vector{Expr}, s::Vector{String})` : parse the strings and
append them to the expressions.
"""
append_str!(ex::Vector{Expr}, s::Vector{String}) = append!(ex, parse.(s))

"""
`bo2edges(N)` : bodyorder-to-edges
"""
bo2edges(N::Integer) = (N * (N-1)) ÷ 2
bo2edges(::Val{N}) where {N} = (N * (N-1)) ÷ 2

"""
`edges2bo(M)`: "edges-to-bodyorder", an internal function that translates
the number of edges in a simplex into the body-order
"""
edges2bo(M::Integer) = (M <= 0) ? 1 : round(Int, 0.5 + sqrt(0.25 + 2 * M))


# ----------- some generic functions that we like to have available globally
#             to overload as needed

"""
some measure of degree - could be but need not be polymomial degree; this
is basis-specific
"""
function degree end


"""
`basisid`: returns an identifier (string) specifying which basis functions
may be combined into a single function.
"""
function basisid end

"""
`combinebasis`:  if `basis::Vector{AbstractCalculator}` with identical
`basisid`, then `combinebasis(basis, coeffs)` should return a new calculator
that combines all basis function into a single IP.
"""
function combinebasis end

# prototypes for the invariants
function invariants end
function invariants_d end
function invariants_ed end

# prototypes for obtaining the descriptor
function descriptor end

function evaluate_I end
function evaluate_I_d end
function evaluate_I_ed end


# ----------- Abstract Supertype for pure NBodyFunctions --------------

"""
`NBodyFunction` : abstract supertype of all "pure" N-body functions.
concrete subtypes must implement

* `bodyorder`
* `evaluate`
* `evaluate_d`
"""
abstract type NBodyFunction{N} <: AbstractCalculator end

bodyorder(V::NBodyFunction{N}) where {N} = N

"""
`AbstractDescriptor`: abstract supertype for different descriptors
of N-body configurations.
"""
abstract type AbstractDescriptor end

"""
`NBSiteDescriptor`: abstract supertype for descriptors that start from
a site-based formulation.
"""
abstract type NBSiteDescriptor <: AbstractDescriptor end


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

==(V1::NBodyIP, V2::NBodyIP) = V1.components == V2.components
cutoff(V::NBodyIP) = maximum( cutoff.(V.components) )
energy(V::NBodyIP, at::Atoms) = sum( energy(Vn, at)  for Vn in V.components )
forces(V::NBodyIP, at::Atoms) = sum( forces(Vn, at)  for Vn in V.components )
virial(V::NBodyIP, at::Atoms) = sum( virial(Vn, at)  for Vn in V.components )

"""
turn a potentially slow representation of an IP into a fast one,
by switching to a different representation.
"""
fast(IP::NBodyIP) = NBodyIP( fast.(IP.components) )

# construct an NBodyIP from a basis
function NBodyIP(basis, coeffs)
   components = AbstractCalculator[]
   tps = typeof.(basis)
   for tp in unique(tps)
      # find all basis functions that have the same type, which in particular
      # incorporated the body-order
      Itp = find(tps .== tp)
      # construct a new basis function by combining all of them in one
      # (this assumes that we have NBody types)
      V_N = combinebasis([basis[Itp]...], coeffs[Itp])
      push!(components, V_N)
   end
   return NBodyIP(components)
end

# functionality for pure NBodyFunctions for computing
#  * site_energies
#  * energy
#  * forces
#  * virial
include("eval_nbody.jl")

# IO of NBodyIPs
include("io.jl")

# space transforms and cutoffs
include("aux.jl")




# ======= experimental ============

function evaluate(V::NBodyFunction{N}, r::SVector{M}) where {N, M}
   D = descriptor(V)
   return evaluate_I(V, invariants(D, r)) * fcut(D, r)
end

function evaluate_d(V::NBodyFunction{N}, r::SVector{M}) where {N, M}
   D = descriptor(V)
   fc, fc_d = fcut_d(D, r)
   Vn, Vn_d = evaluate_I_ed(V, invariants_ed(D, r))
   return fc * Vn_d + fc_d * Vn
end
