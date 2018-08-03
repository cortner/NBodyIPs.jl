
"""
`module BLPolys`

The exported symbols are
* `BLDictionary`: collects all the information about a basis
* `NBody`: an N-body function wrapped into a JuLIP calculator
* `poly_basis` : generate a basis of N-body functions

## Usage

## Notation

* `N` : N-body
* `M` : number of edges, i.e., M = N (N-1)/2
* `K` : length of the tuples defining polynomials, K = M+1
"""
module BLPolys

import StaticPolynomials

using JuLIP, NeighbourLists, StaticArrays

using NBodyIPs: NBodyFunction,
                SpaceTransform,
                Cutoff

using NBodyIPs.FastPolys: fpoly,
                          fpoly_d

import Base: length,
             Dict,
             ==

import JuLIP: cutoff

import JuLIP.Potentials: evaluate,
                         evaluate_d

import NBodyIPs: bodyorder,
                 fast,
                 evaluate_many!,
                 evaluate_many_d!,
                 degree,
                 combine_basis

const Tup{M} = NTuple{M, Int}
const VecTup{M} = Vector{NTuple{M, Int}}

export BLNBody,
       BLDictionary,
       bl_basis


# import the raw invariant polynomials
include("blinvariants.jl")

using BLInvariants: invariants,
                    invariants_d,
                    invariants_ed,
                    corners,
                    tdegrees,
                    bo2edges,
                    edges2bo

# --------------- BLDictionary -------------

"""
`struct BLDictionary` : specifies all details about the basis functions

## Constructor

`BLDictionary(ftrans, fcut)`, where `ftrans`, `fcut` specify in one of several
ways how the dictionary is defined. If
the radius is not provided then `BLDictionary` will try to infer it from the
`fcut` argument. E.g.,
```julia
D = BLDictionary("r -> 1/r", (:cos, rcut1, rcut2))
```

Known symbols for the cutoff are
```
[:cos, :sw, :spline, :square, :cos2s]
```

## Developer Doc: Methods associated with a `D::BLDictionary`:

* `invariants`, `invariants_d`, `invariants_ed`: compute invariants and jacobian
   in transformed coordinates defined by `D.transform`
* `evaluate`, `evaluate_d`: evaluate the (univariate) basis
   function associated with this dictionary; at the moment only
   standard polynomials are admitted
* `fcut, fcut_d`: evulate the cut-off function and its derivative / gradient
   when interpreted as a product of cut-off functions
"""
struct BLDictionary{TT, TC}
   transform::TT
   cutoff::TC
end

function ==(D1::BLDictionary, D2::BLDictionary)
   return (D1.transform == D2.transform) && (D1.cutoff == D2.cutoff)
end

Dict(D::BLDictionary) = Dict("__id__" => "BLDictionary",
                             "transform" => Dict(D.transform),
                             "cutoff" => Dict(D.cutoff))

BLDictionary(D::Dict) = BLDictionary( SpaceTransform(D["transform"]),
                                      Cutoff(D["cutoff"]) )

Base.convert(::Val{:BLDictionary}, D::Dict) = BLDictionary(D)


BLDictionary(transform::String, cutoff::Union{String, Tuple}) =
         BLDictionary(SpaceTransform(transform), Cutoff(cutoff))


@inline transform(D::BLDictionary, r::Number) = transform.f(r)
@inline transform_d(D::BLDictionary, r::Number) = transform.f_d(r)
@inline fcut(D::BLDictionary, r::Number) = D.cutoff.f(r)
@inline fcut_d(D::BLDictionary, r::Number) = D.cutoff.f_d(r)
@inline cutoff(D::BLDictionary) = D.cutoff.rcut



# ==================================================================
#           Polynomials of Invariants
# ==================================================================


@pot struct BLNBody{N, M, T, TD} <: NBodyFunction{N}
   t::VecTup{M}               # tuples M = #edges + 1
   c::Vector{T}               # coefficients
   D::TD                      # BLDictionary (or nothing)
   valN::Val{N}               # encodes that this is an N-body function
end

==(V1::BLNBody, V2::BLNBody) = ( (V1.t == V2.t) && (V1.c == V2.c) && (V1.D == V2.D) )


"""
`struct BLNBody`  (N-Body Basis Function)

A struct storing the information for a pure N-body potential, i.e., containing
*only* terms of a specific body-order. Several `BLNBody`s can be
combined into an interatomic potential via `NBodyIP`.

### Fields

* `t::Vector{NTuple{M,TI}}` : list of M-tuples containing basis function information
e.g., if M = 3, α = t[1] is a 3-vector then this corresponds to the basis function
`f[α[1]](Q[1]) * f[α[2]](Q[2]) * f[α[3]](Q[3])` where `Q` are the 3-body invariants.

* `c::Vector{T}`: vector of coefficients for the basis functions

* `D`: a `BLDictionary`
"""
BLNBody

# standard constructor (N can be inferred)
BLNBody(t::VecTup{K}, c, D) where {K} = BLNBody(t, c, D, Val(edges2bo(K-1)))

# BLNBody made from a single basis function rather than a collection
BLNBody(t::Tup, c, D) = BLNBody([t], [c], D)

# collect multiple basis functions represented as BLNBody's into a single BLNBody
# (for performance reasons)
BLNBody(B::Vector{TB}, c, D) where {TB <: BLNBody} =
      BLNBody([b.t[1] for b in B], c .* [b.c[1] for b in B], D)

# 1-body term (on-site energy)
BLNBody(c::Float64) =
      BLNBody([Tup{0}()], [c], nothing, Val(1))

# number of basis functions which this term is made from
length(V::BLNBody) = length(V.t)

cutoff(V::BLNBody) = cutoff(V.D)

function match_dictionary(V::BLNBody, V1::BLNBody)
   if V.D != V1.D
      if V.D.s != V1.D.s
         warn("matching two non-matching dictionaries!")
      end
   end
   return BLNBody(V.t, V.c, V1.D, V.valN)
end

combine_basis(basis::AbstractVector{TV}, coeffs) where {TV <: BLNBody} =
      BLNBody(basis, coeffs, basis[1].D)

function degree(V::BLNBody)
   if length(V) == 1
      return tdegree(V.t[1])
   end
   error("`degree` is only defined for `BLNBody` basis functions, length == 1")
end

degree(V::BLNBody{1}) = 0

function Base.info(B::Vector{T}; indent = 2) where T <: BLNBody
   ind = repeat(" ", indent)
   println(ind * "body-order = $(bodyorder(B[1]))")
   println(ind * "    length = $(length(B))")
   if bodyorder(B[1]) > 1
      println(ind * " transform = $(B[1].D.s[1])")
      println(ind * "    cutoff = $(B[1].D.s[2])")
   end
end

# -------------- Infrastructure to read/write BLNBody  --------


Dict(V::BLNBody{N}) where {N} = Dict( "__id__" => "BLNBody",
                                      "t" => V.t,
                                      "c" => V.c,
                                      "D" => Dict(V.D),
                                      "N" => N )

function BLNBody(D::Dict)
   N = D["N"]
   t = [ tuple(ti...) for ti in D["t"] ]
   c = Vector{Float64}(D["c"])
   if N == 1
      return BLNBody(t, c, nothing, Val(1))
   end
   return BLNBody(t, c, BLDictionary(D["D"]), Val(N))
end

Base.convert(::Val{:BLNBody}, D::Dict) = BLNBody(D)



# ==================================================================
#    StBLNBody
# ==================================================================

@pot struct StBLNBody{N, TD, TP} <: NBodyFunction{N}
   D::TD       # BLDictionary (or nothing)
   P::TP       # a static polynomial
   valN::Val{N}
end

"""
`struct StBLNBody`  (N-Body Bond-length Polynomial)

fast evaluation of the outer polynomial using `StaticPolynomials`
"""
StBLNBody


function StBLNBody(V::BLNBody{N}) where {N}
   M = bo2edges(N)  # number of edges for body-order N
   I1, I2 = BLInvariants.invariants(V.D, SVector(ones(M)...))  # how many invariants
   nI1 = length(I1)
   ninvariants = length(I1) + length(I2)
   nmonomials = length(V.c)
   # generate the exponents for the StaticPolynomial
   exps = zeros(Int, ninvariants, nmonomials)
   for (i, α) in enumerate(V.t)  # i = 1:nmonomials
      for (j, a) in enumerate(α[1:end-1])   #  ∏ I1[j]^α[j]
         exps[j, i] = a    # I1[j]^α[j]
      end
      exps[nI1+1+α[end], i] = 1   # I2[α[end]] * (...)
   end
   # generate the static polynomial
   return StBLNBody(V.D, StaticPolynomials.Polynomial(V.c, exps), V.valN)
end

cutoff(V::StBLNBody) = cutoff(V.D)

fast(Vn::StBLNBody)  = Vn
fast(Vn::BLNBody) =  StBLNBody(Vn)
fast(Vn::BLNBody{1}) = Vn

# --------- EVALUATION OF THE BLNBody and StBLNBody potentials ---------
include("eval_blnbody.jl")


# ==================================================================
#           Generate a BLNBody Basis
# ==================================================================


"""
compute the total degree of the polynomial represented by α.
Note that `M = K-1` where `K` is the tuple length while
`M` is the number of edges.
"""
function tdegree(α)
   K = length(α)
   degs1, degs2 = tdegrees(Val(edges2bo(K-1)))
   # primary invariants
   d = sum(α[j] * degs1[j] for j = 1:K-1)
   # secondary invariants
   d += degs2[1+α[end]]
   return d
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
gen_tuples(N, deg; tuplebound = (α -> (0 < tdegree(α) <= deg))) =
   gen_tuples(Val(N), Val(bo2edges(Val(N))+1), deg, tuplebound)

function gen_tuples(vN::Val{N}, vK::Val{K}, deg, tuplebound) where {N, K}
   A = Tup{K}[]
   degs1, degs2 = tdegrees(vN)

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


"""
* `bl_basis(N, D, deg; tuplebound = ...)`
* `bl_basis(N, strans, scut, deg; tuplebound = ...)`

generates a basis set of
`N`-body functions, with dictionary `D`, maximal degree `deg`; the precise
set of basis functions constructed depends on `tuplebound` (see `?gen_tuples`)
"""
bl_basis(N::Integer, D::BLDictionary, deg; kwargs...) =
   bl_basis(gen_tuples(N, deg; kwargs...), D)

bl_basis(ts::VecTup, D::BLDictionary) = [BLNBody(t, 1.0, D) for t in ts]

bl_basis(N::Integer, strans::String, scut::String, deg; kwargs...) =
   bl_basis(N, BLDictionary(strans, scut), deg; kwargs...)



end # module

# TODO TODO

# include("poly_regularise.jl")


# # ----------------- some simplified access functions ------------------
#
# evaluate(V::NBodyFunction{2}, r::Number) = evaluate(V, SVector(r))
#
# evaluate_d(V::NBodyFunction{2}, r::Number) = evaluate_d(V, SVector(r))[1]
#
# evaluate_dd(V::NBodyFunction{2}, r::Number) =
#       ((@D V(r+1e-5)) - (@D V(r-1e-5))) / 1e-5
#
# evaluate(V::NBodyFunction{3}, r1::Number, r2::Number, r3::Number) =
#       evaluate(V, SVector(r1, r2, r3))



# # =============== Experimental:
# #   evaluate NBodyIP
#
# (V::NBodyIP)(args...) = evaluate(V, args...)
#
# evaluate(V::NBodyIP, r::Number) = evaluate(V::NBodyIP, SVector(r))
#
# evaluate(V::NBodyIP, r1::T, r2::T, r3::T) where {T <: Number} =
#       evaluate(V::NBodyIP, SVector(r1, r2, r3))
#
# function evaluate(V::NBodyIP, r::SVector{N, T}) where {N, T}
#    v = zero(T)
#    for Vn in V.orders
#       if bo2edges(bodyorder(Vn)) == N
#          v += Vn(r)
#       end
#    end
#    return v
# end

# dim(V::BLNBody{N,M}) where {N, M} = M-1
