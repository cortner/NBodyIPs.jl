
using StaticArrays

using JuLIP: AbstractCalculator,
             Atoms,


import Base:Dict

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


# ----------- Abstract Supertype for pure NBodyFunctions --------------

"""
`NBodyFunction` : abstract supertype of all "pure" N-body functions.
concrete subtypes must implement

* `bodyorder`
* `evaluate`
* `evaluate_d`
"""
abstract type NBodyFunction{N} <: AbstractCalculator end

bodyorder(V::NBodyFunction{N}) = N



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

cutoff(V::NBodyIP) = maximum( cutoff.(V.orders) )
energy(V::NBodyIP, at::Atoms) = sum( energy(Vn, at)  for Vn in V.orders )
forces(V::NBodyIP, at::Atoms) = sum( forces(Vn, at)  for Vn in V.orders )
virial(V::NBodyIP, at::Atoms) = sum( virial(Vn, at)  for Vn in V.orders )

"""
turn a potentially slow representation of an IP into a fast one,
by switching to a different representation.
"""
fast(IP::NBodyIP) = NBodyIP( fast.(IP.orders) )

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
      V_N = combine_basis([basis[Itp]...], coeffs[Itp])
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
