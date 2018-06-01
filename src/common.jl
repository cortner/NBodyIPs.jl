
using StaticArrays, ForwardDiff

using JuLIP: AbstractCalculator, Atoms, neighbourlist, @D, JVec
using NeighbourLists: nbodies, maptosites!, maptosites_d!, virial!, max_neigs


import JuLIP.Potentials: evaluate, evaluate_d, evaluate_dd
import JuLIP: cutoff, energy, forces, site_energies, virial, stress

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


include("eval_nbody.jl")



function site_energies(V::NBodyFunction{N}, at::Atoms{T}) where {N, T}
   Es = zeros(T, length(at))
   for (i, j, r, R) in sites(at, cutoff(V))
      Es[i] = eval_site_nbody!(Val(N), R, cutoff(V),
                               ((out, s, _,_1,_2) -> out + evaluate(V, s)),
                               zero(T), nothing)
   end
   return Es/N
end



# this is probably already in JuLIP??? if not, it should be moved to JuLIP??
energy(V::NBodyFunction, at::Atoms) =
      sum_kbn(site_energies(V, at))


# this appears to be a nice generic implementation of forces with a
# temporary array => move this to JuLIP!
function forces(V::NBodyFunction{N}, at::Atoms{T}) where {N, T}
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   for (i, j, r, R) in sites(nlist)
      dVsite .*= 0.0
      eval_site_nbody!(
            Val(N), R, cutoff(V),
            (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(V, s)/N, J, S),
            dVsite, nothing )
      # write site energy gradient into forces
      for n = 1:length(j)
         F[j[n]] -= dVsite[n]
         F[i] += dVsite[n]
      end
   end
   return F
end



function virial(V::NBodyFunction{N}, at::Atoms{T}) where {N, T}
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   S = @SMatrix zeros(3,3)
   dVsite = zeros(JVec{T}, maxneigs)
   for (i, j, r, R) in sites(nlist)
      dVsite .*= 0.0
      # eval_site_d!(dVsite, V, R)
      eval_site_nbody!(
            Val(N), R, cutoff(V),
            (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(V, s)/N, J, S),
            dVsite, nothing )
      S += JuLIP.Potentials.site_virial(dVsite, R)
   end
   return S
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
`energy`, `forces`, `virial` and `cutoff`.

TODO: `site_energies`, etc.
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


# some simplified access functions

evaluate(V::NBodyFunction{2}, r::Number) = evaluate(V, SVector(r))

evaluate_d(V::NBodyFunction{2}, r::Number) = evaluate_d(V, SVector(r))[1]

evaluate_dd(V::NBodyFunction{2}, r::Number) =
      ((@D V(r+1e-5)) - (@D V(r-1e-5))) / 1e-5

evaluate(V::NBodyFunction{3}, r1::Number, r2::Number, r3::Number) =
      evaluate(V, SVector(r1, r2, r3))


# TODO: implement direct n-body access for NBodyIP


# ========================= assembly support for LSQ system ====================

# For assembling the LSQ system efficiently we need a way to evaluate all basis
# functions of the same body-order at the same time. Otherwise we would be
# re-computing the invariants many many times, which is very expensive.
# To achieve this we just wrap all basis functions of a body-order into
# a new type `NBBasis` which evaluates to a long vector
#
# at the moment, it seems we need to hard-code this to the Polys
# sub-module, but it would be good if this can be fixed, so we keep this
# "interface" here.

function evaluate_many! end
function evaluate_many_d! end


energy(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ energy(b, at) for b in B ]

function energy(B::AbstractVector{TB}, at::Atoms{T}
                ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   temp = zeros(T, length(B))
   E = zeros(T, length(B))
   for (i, j, r, R) in sites(nlist)
      # evaluate all the site energies at the same time
      # for each simples, write the nB energies into temp
      # then add them to E, which is just passed through all the
      # various loops, so no need to update it here again
      eval_site_nbody!(Val(N), R, rcut,
                       (out, s, S, _1, temp) -> (out .+= evaluate_many!(temp, B, s)),
                       E, temp)
   end
   # rescale to account for permutations
   E ./= N
   return E
end



function _acc_manyfrcs(B, dVsite, s, S, J, temp)
   dV = evaluate_many_d!(temp, B, s)
   for ib = 1:length(dVsite)
      _grad_len2pos!(dVsite[ib], dV[ib], J, S)
   end
   return dVsite
end

function forces(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nedges = (N*(N-1))÷2
   nB = length(B)
   # forces
   F =      [ zeros(JVec{T}, length(at)) for n = 1:nB ]
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:nB ]
   # n-body gradients
   dV =     [ zeros(T, nedges)      for n = 1:nB ]
   # temporary arrays to compute the site gradients
   temp = ( zeros(T, nB),
            zeros(T, nedges, nB),
            zeros(T, nedges, nB),
            dV )
   accum_fun = let B=B
      (out, s, S, J, temp) -> _acc_manyfrcs(B, out, s, S, J, temp)
   end

   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:nB
         dVsite[n] .*= 0.0
      end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut, accum_fun, dVsite, temp)
      # write it into the force vectors
      for ib = 1:nB, n = 1:length(j)
         F[ib][j[n]] -= dVsite[ib][n]/N
         F[ib][i] += dVsite[ib][n]/N
      end
   end
   return F
end

forces(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ forces(b, at) for b in B ]


function virial(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nedges = (N*(N-1))÷2
   # virials (main output)
   S = fill((@SMatrix zeros(3,3)), length(B))
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:length(B) ]
   # n-body gradients
   dV =     [ zeros(T, nedges)      for n = 1:length(B) ]
   # temporary arrays to compute the site gradients
   temp = ( zeros(T, length(B)),
            zeros(T, nedges, length(B)),
            zeros(T, nedges, length(B)),
            dV )
   accum_fun = let B=B
      (out, s, S, J, temp) -> _acc_manyfrcs(B, out, s, S, J, temp)
   end

   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:length(dVsite)
         dVsite[n] .*= 0.0
      end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut, accum_fun, dVsite, temp)
      # update the virials
      for iB = 1:length(B)
         S[iB] += JuLIP.Potentials.site_virial(dVsite[iB], R) / N
      end
   end
   return S
end


virial(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ virial(b, at) for b in B ]
