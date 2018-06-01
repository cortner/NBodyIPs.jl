
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
                               ((out, s, S, _1, _2) -> out + evaluate(V, s)/N),
                               zero(T), nothing)
   end
   return Es
end


# REVISIT using maptosites!
# site_energies(V::NBodyFunction, at::Atoms{T}) where {T} =
#    maptosites!( (i,j,r,R) -> eval_site(V, R),
#                 zeros(T, length(at)),
#                 sites(at, cutoff(V)) )


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
      # eval_site_d!(dVsite, V, R)
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


# function forces(V::NBodyFunction, at::Atoms{T}) where {T}
#    nlist = neighbourlist(at, cutoff(V))
#    return scale!(maptosites_d!(r -> evaluate_d(V, r),
#                  zeros(SVector{3, T}, length(at)),
#                  nbodies(bodyorder(V), nlist)), -1)
# end



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


# function virial(V::NBodyFunction, at::Atoms{T}) where {T}
#    nlist = neighbourlist(at, cutoff(V))
#    temp = @MMatrix zeros(3,3)
#    virial!(r -> evaluate_d(V, r), temp, nbodies(bodyorder(V), nlist))
#    return SMatrix(temp)
# end

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

_alloc_svec(T::Type, ::Val{N}) where {N} = zero(SVector{N, T})
_alloc_svec(T::Type, N::Integer) = _alloc_svec(T, Val(N))

_alloc_smat(T::Type, ::Val{N}, ::Val{M}) where {N, M} = zero(SMatrix{N, M, T})
_alloc_smat(T::Type, N, M) = _alloc_smat(T, Val(N), Val(M))

_alloc_mvec(T::Type, ::Val{N}) where {N} = zero(MVector{N, T})
_alloc_mvec(T::Type, N::Integer) = _alloc_mvec(T, Val(N))

_alloc_mmat(T::Type, ::Val{N}, ::Val{M}) where {N, M} = zero(MMatrix{N, M, T})
_alloc_mmat(T::Type, N, M) = _alloc_mmat(T, Val(N), Val(M))


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
                       (out, s, S, _1, temp) -> (out .+= evaluate_many!(temp, B, s)/N),
                       E, temp)
   end
   return E
end


function _many_grad_len2pos!(dVsite, dV, J, S)
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
   # forces
   F =      [ zeros(JVec{T}, length(at)) for n = 1:length(B) ]
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:length(B) ]
   # n-body gradients
   dV =     [ zeros(T, nedges)      for n = 1:length(B) ]
   # temporary arrays to compute the site gradients
   temp = ( zeros(T, length(B)),
            zeros(T, nedges, length(B)),
            zeros(T, nedges, length(B)),
            dV )
   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:length(dVsite)
         dVsite[n] .*= 0.0
      end
      # fill dVsite
      eval_site_nbody!(
            Val(N), R, rcut,
            (out, s, S, J, temp) -> _many_grad_len2pos!(out,
                                       evaluate_many_d!(temp, B, s)/N, J, S),
            dVsite, temp)
      # write it into the force vectors
      for ib = 1:length(B), n = 1:length(j)
         F[ib][j[n]] -= dVsite[ib][n]
         F[ib][i] += dVsite[ib][n]
      end
   end
   return F
end

forces(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ forces(b, at) for b in B ]


function virial(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   nlist = neighbourlist(at, cutoff(B[1]))
   z2 = _alloc_svec(T, length(B))
   temp = ( _alloc_mvec(T, length(B)),
            _alloc_mmat(T, (N*(N-1))÷2, length(B)),
            _alloc_mmat(T, (N*(N-1))÷2, length(B)),
            _alloc_mvec(typeof(z2), (N*(N-1))÷2)
          )
   out = fill((@SMatrix zeros(3,3)), length(B))
   virial!( r -> evaluate_many_d!(temp, B, r),
            out,
            nbodies(N, nlist) )
   # stress = - virial(c, a) / det(defm(a))
   # scale!(out, -1/det(defm(at)))
   return out
end

virial(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ virial(b, at) for b in B ]


#
# teach NeighbourLists.jl how to assemble collections of forces
#
import NeighbourLists._m2s_mul_

@generated function _m2s_mul_(X::SVector{M,T}, S::SVector{N,T}) where {M,N,T}
   exprs = Expr[]
   for i = 1:M
      push_str!(exprs, "xS_$i = X[$i] * S")
   end
   coll = "p = @SVector ["
   for i = 1:M
      coll *= "xS_$i, "
   end
   coll *= "]"
   push_str!(exprs, coll)

   quote
      $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, exprs...))
      return p
   end
end


#
# teach NeighbourLists.jl how to assemble collections of stresses
#
import NeighbourLists._inc_stress_!

function _inc_stress_!(out::Vector{T1}, s::Float64, df::SVector{N,Float64}, S::SVector{3,Float64}
         ) where T1 <: SMatrix{3,3,Float64} where N
   # @assert length(out) == length(df)
   # error("stop here")
   for n = 1:length(df)
      out[n] -= (s*df[n]) * (S * S')
   end
end
