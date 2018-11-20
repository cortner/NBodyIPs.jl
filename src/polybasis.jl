
module PolyBasis

using StaticArrays
import NBodyIPs

using NBodyIPs:        BondLengthDesc,
                       BondAngleDesc,
                       edges2bo,
                       bo2edges

using NBodyIPs.Polys:  Tup,
                       VecTup,
                       NBPoly

export blpolys,
       bapolys,
       nbpolys

# hack to get the correct tdegrees function
tdegrees(desc::BondLengthDesc, vN) = NBodyIPs.BLInvariants.tdegrees(vN)
tdegrees(desc::BondAngleDesc, vN) = NBodyIPs.BAInvariants.tdegrees(vN)

rtdegrees(desc::BondAngleDesc, vN) = NBodyIPs.BAInvariants.rtdegrees(vN)

"""
compute the total degree of the polynomial represented by α.
Note that `M = K-1` where `K` is the tuple length while
`M` is the number of edges.
"""
function tdegree(desc, α)
   K = length(α)
   degs1, degs2 = tdegrees(desc, Val(edges2bo(K-1)))
   # primary invariants
   d = sum(α[j] * degs1[j] for j = 1:K-1)
   # secondary invariants
   d += degs2[1+α[end]]
   return d
end

"""
compute the total rt degree of the polynomial represented by α
(only for bond-angles).
Note that `M = K-1` where `K` is the tuple length while
`M` is the number of edges.
"""
function rtdegree(desc::BondAngleDesc, α)
   K = length(α)
   degs1, degs2 = rtdegrees(desc, Val(edges2bo(K-1)))
   # primary invariants
   dr = sum(α[j] * degs1[j][1] for j = 1:K-1)
   dt = sum(α[j] * degs1[j][2] for j = 1:K-1)
   # secondary invariants
   dr += degs2[1+α[end]][1]
   dt += degs2[1+α[end]][2]
   return dr,dt
end


"""
`gen_tuples(N, deg; tuplebound = ...)` : generates a list of tuples, where
each tuple defines a basis function. Use `gen_basis` to convert the tuples
into a basis, or use `gen_basis` directly.

* `N` : body order
* `deg` : maximal degree
* `tuplebound` : a function that takes a tuple as an argument and returns
`true` if that tuple should be in the basis and `false` if not. The default is
`α -> (0 < tdegree(α) <= deg)` i.e. the standard monomial degree. (note this is
the degree w.r.t. lengths, not w.r.t. invariants!) The `tuplebound` function
must be **monotone**, that is, `α,β` are tuples with `all(α .≦ β)` then
`tuplebound(β) == true` must imply that also `tuplebound(α) == true`.
"""
gen_tuples(desc, N, deg; tuplebound = (α -> (0 < tdegree(desc, α) <= deg))) =
   gen_tuples(desc, Val(N), Val(bo2edges(Val(N))+1), deg, tuplebound)

function gen_tuples(desc, vN::Val{N}, vK::Val{K}, deg, tuplebound) where {N, K}
   A = Tup{K}[]
   degs1, degs2 = tdegrees(desc, vN)

   α = @MVector zeros(Int, K)
   α[1] = 1
   lastinc = 1

   while true
      admit_tuple = false
      if α[end] <= length(degs2)-1
         if tuplebound(α)
            admit_tuple = true
         end
      end
      if admit_tuple
         push!(A, SVector(α).data)
         α[1] += 1
         lastinc = 1
      else
         if lastinc == K
            return A
         end
         α[1:lastinc] = 0
         α[lastinc+1] += 1
         lastinc += 1
      end
   end
   error("I shouldn't be here!")
end




nbpolys(N::Integer, desc, tdeg; kwargs...) =
   nbpolys(gen_tuples(desc, N, tdeg; kwargs...), desc)

nbpolys(ts::VecTup, desc) = [NBPoly(t, 1.0, desc) for t in ts]

"""
`blpolys(N, trans, cutoff, tdeg; tuplebound = ...) -> Vector{<:NBPoly}`

generates a basis set of Bond-length many-body functions, with
* body-order `N`
* transformation `trans`
* cutoff mechanism described by `cutoff`
* maximal total degree `tdeg`

For more flexibility the kwargs `tuplebound` can be used, use with care.
"""
blpolys(N::Integer, trans::String, cutoff::String, deg; kwargs...) =
   nbpolys(N, BondLengthDesc(trans, cutoff), deg; kwargs...)

"""
`bapolys(N, trans, cutoff, tdeg; tuplebound = ...) -> Vector{<:NBPoly}`

generates a basis set of Bond-angle many-body functions, with
* body-order `N`
* transformation `trans`
* cutoff mechanism described by `cutoff`
* maximal total degree `tdeg`

For more flexibility the kwargs `tuplebound` can be used, use with care.
"""
bapolys(N::Integer, trans::String, cutoff::String, deg; kwargs...) =
   nbpolys(N, BondAngleDesc(trans, cutoff), deg; kwargs...)

# bapolys(N::Integer, trans::String, cutoff::String, (degr,degt); kwargs...) =
#    nbpolys(gen_tuples(desc, N, tdeg; tuplebound = rt_tuplebound), desc)


# --------------------------
# Using RT degrees to generate bond-angles polynomials with given R and theta degrees
# --------------------------


function rt_tuplebound(α,degr,degt,desc::BondAngleDesc)
   dr,dt = rtdegree(desc, α)
   return (dr <= degr)&&(dt <= degt)
end


function gen_tuples_rt(desc::BondAngleDesc, vN::Val{N}, vK::Val{K}, degr, degt;  tuplebound = (α -> rt_tuplebound(α,degr,degt,desc) )) where {N, K}
   A = Tup{K}[]
   degs1, degs2 = rtdegrees(desc, vN)

   α = @MVector zeros(Int, K)
   α[1] = 1
   lastinc = 1

   while true
      admit_tuple = false
      if α[end] <= length(degs2)-1
         if tuplebound(α)
            admit_tuple = true
         end
      end
      if admit_tuple
         push!(A, SVector(α).data)
         α[1] += 1
         lastinc = 1
      else
         if lastinc == K
            return A
         end
         α[1:lastinc] = 0
         α[lastinc+1] += 1
         lastinc += 1
      end
   end
   error("I shouldn't be here!")
end

gen_tuples_rt(desc::BondAngleDesc, N, degr, degt) =
   gen_tuples_rt(desc,
                 Val(N),  Val(bo2edges(Val(N))+1),
                 degr, degt)


nbpolys(N::Integer, desc, degr, degt; kwargs...) =
        nbpolys(gen_tuples_rt(desc, N, degr, degt; kwargs...), desc)



end
