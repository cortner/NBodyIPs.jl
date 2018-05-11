
using StaticArrays

using JuLIP: AbstractCalculator, Atoms, neighbourlist
using NeighbourLists: nbodies, maptosites!, maptosites_d!, virial!
using JuLIP.Potentials: evaluate, evaluate_d

import JuLIP: cutoff, energy, forces, site_energies, virial

export NBodyIP,
       bodyorder,
       fast

"""
`NBodyFunction` : abstract supertype of all "pure" N-body functions.
concrete subtypes must implement

* `bodyorder`
* `evaluate`
* `evaluate_d`
"""
abstract type NBodyFunction{N} <: AbstractCalculator end

# prototypes of function defined on `NBodyFunction`
function bodyorder end

function site_energies(V::NBodyFunction, at::Atoms{T}) where {T}
   nlist = neighbourlist(at, cutoff(V))
   return maptosites!(r -> evaluate(V, r),
                      zeros(T, length(at)),
                      nbodies(bodyorder(V), nlist))
end

energy(V::NBodyFunction, at::Atoms) =
      sum_kbn(site_energies(V, at))

function forces(V::NBodyFunction, at::Atoms{T}) where {T}
   nlist = neighbourlist(at, cutoff(V))
   return scale!(maptosites_d!(r -> evaluate_d(V, r),
                 zeros(SVector{3, T}, length(at)),
                 nbodies(bodyorder(V), nlist)), -1)
end

function virial(V::NBodyFunction, at::Atoms{T}) where {T}
   nlist = neighbourlist(at, cutoff(V))
   temp = @MMatrix zeros(3,3)
   virial!(r -> evaluate_d(V, r), temp, nbodies(bodyorder(V), nlist))
   return SMatrix(temp)
end

# ------ special treatment of 1-body functions

site_energies(V::NBodyFunction{1}, at::Atoms) =
      fill(V(), length(at))

forces(V::NBodyFunction{1}, at::Atoms{T}) where {T} =
      zeros(SVector{3, T}, length(at))

virial(V::NBodyFunction{1}, at::Atoms{T}) where {T} =
      zero(SMatrix{3, 3, T})

"""
`NBodyIP` : wraps `NBodyFunction`s into a JuLIP calculator, defining
`energy`, `forces` and `cutoff`.

TODO: `stress`, `site_energies`, etc.
"""
struct NBodyIP <: AbstractCalculator
   orders::Vector{NBodyFunction}
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
