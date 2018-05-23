
using StaticArrays, ForwardDiff

using JuLIP: AbstractCalculator, Atoms, neighbourlist, @D
using NeighbourLists: nbodies, maptosites!, maptosites_d!, virial!


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


# generics

evaluate(V::NBodyFunction{2}, r::Number) = evaluate(V, SVector(r))

evaluate_d(V::NBodyFunction{2}, r::Number) = evaluate_d(V, SVector(r))[1]

evaluate_dd(V::NBodyFunction{2}, r::Number) =
      ((@D V(r+1e-5)) - (@D V(r-1e-5))) / 1e-5

evaluate(V::NBodyFunction{3}, r1::Number, r2::Number, r3::Number) =
      evaluate(V, SVector(r1, r2, r3))


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

function energy(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   # @assert isleaftype{TB}
   nlist = neighbourlist(at, cutoff(B[1]))
   z = _alloc_svec(T, length(B))
   temp = _alloc_mvec(T, length(B))
   return maptosites!(r -> evaluate_many!(temp, B, r),
                      [ copy(z) for _ = 1:length(at) ],
                      nbodies(N, nlist)) |> sum
end

function forces(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   # @assert isleaftype{TB}
   nlist = neighbourlist(at, cutoff(B[1]))
   z = _alloc_svec(JVec{T}, length(B))
   z2 = _alloc_svec(T, length(B))
   temp = ( _alloc_mvec(T, length(B)),
            _alloc_mmat(T, (N*(N-1))÷2, length(B)),
            _alloc_mmat(T, (N*(N-1))÷2, length(B)),
            _alloc_mvec(typeof(z2), (N*(N-1))÷2)
          )
   Fpre = maptosites_d!(r -> evaluate_many_d!(temp, B, r),
                      [ copy(z) for _ = 1:length(at) ],
                      nbodies(N, nlist))
   F = [ [ -Fpre[i][j] for i = 1:length(Fpre) ]  for j = 1:length(B) ]
   return F
end

function stress(B::AbstractVector{TB}, at::Atoms{T}
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
   scale!(out, -1/det(defm(at)))
   return out
end


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
