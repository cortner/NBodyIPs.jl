
"""
`module Polys`

The exported symbols are
* `NBody`: an N-body function wrapped into a JuLIP calculator
* `blpoly_basis` : generate a basis of N-body functions
* `bapoly_basis` : generate a basis of N-body functions

## Usage

## Notation

* `N` : N-body
* `M` : number of edges, i.e., M = N (N-1)/2
* `K` : length of the tuples defining polynomials, K = M+1
"""
module Polys

import StaticPolynomials

using StaticArrays

using NBodyIPs: NBodyFunction,
                bodyorder,
                _decode_dict,
                BondLengthDesc,
                edges2bo,
                bo2edges,
                tdegrees,
                invariants


import Base:              length,
                          Dict,
                          ==
import JuLIP:             cutoff
import JuLIP.Potentials:  @pot
import NBodyIPs:          fast,
                          degree,
                          combinebasis,
                          descriptor,
                          evaluate_many!,
                          evaluate_many_d!,
                          evaluate_I,
                          evaluate_I_ed

const Tup{M} = NTuple{M, Int}
const VecTup{M} = Vector{NTuple{M, Int}}

export NBPoly,
       StNBPoly,
       blpolys


# ==================================================================
#           Polynomials of Invariants
# ==================================================================

@pot struct NBPoly{N, M, T, TD} <: NBodyFunction{N}
   t::VecTup{M}               # tuples M = #edges + 1
   c::Vector{T}               # coefficients
   D::TD                      # Descriptor
   valN::Val{N}               # encodes that this is an N-body function
end

"""
`struct NBPoly`  (N-Body Polynomial, slow implementation)

A `struct` storing the information for a (pure) N-body potential, i.e.,
containing *only* terms of a specific body-order. Several `NBPoly`s can be
combined into an interatomic potential via `NBodyIP`.

### Fields

* `t::Vector{NTuple{M,TI}}` : list of M-tuples containing basis function information
e.g., if M = 7, α = t[n] is a 7-vector then this corresponds to the basis function
```
I2[α[7]] * ∏_{j=1..6} I1[j]^α[j]
```
where `I1, I2` are the 4-body invariants.

* `c::Vector{T}`: vector of coefficients for the basis functions

* `D`: a descriptor (cf `NBodyIPs.AbstractDescriptor`)
"""
NBPoly

==(V1::NBPoly, V2::NBPoly) = ( (V1.t == V2.t) && (V1.c == V2.c) && (V1.D == V2.D) )

descriptor(V::NBPoly) = V.D

# standard constructor (N can be inferred)
NBPoly(t::VecTup{K}, c, D) where {K} = NBPoly(t, c, D, Val(edges2bo(K-1)))

# NBPoly made from a single basis function rather than a collection
NBPoly(t::Tup, c, D) = NBPoly([t], [c], D)

# collect multiple basis functions represented as NBPoly's into a single NBPoly
# (for performance reasons)
# TODO: this is not general enough!
NBPoly(B::Vector{TB}, c, D) where {TB <: NBPoly} =
      NBPoly([b.t[1] for b in B], c .* [b.c[1] for b in B], D)

# 1-body term (on-site energy)
NBPoly(c::Float64) = NBPoly([Tup{0}()], [c], nothing, Val(1))

# number of basis functions which this term is made from
length(V::NBPoly) = length(V.t)

cutoff(V::NBPoly) = cutoff(V.D)

function match_dictionary(V::NBPoly, V1::NBPoly)
   if V.D != V1.D
      if V.D.s != V1.D.s
         warn("matching two non-matching dictionaries!")
      end
   end
   return NBPoly(V.t, V.c, V1.D, V.valN)
end

combinebasis(basis::AbstractVector{TV}, coeffs) where {TV <: NBPoly} =
      NBPoly(basis, coeffs, basis[1].D)

function degree(V::NBPoly)
   if length(V) == 1
      return tdegree(V.t[1])
   end
   error("`degree` is only defined for `NBPoly` basis functions, length == 1")
end

degree(V::NBPoly{1}) = 0

function Base.info(B::Vector{T}; indent = 2) where T <: NBPoly
   ind = repeat(" ", indent)
   println(ind * "body-order = $(bodyorder(B[1]))")
   println(ind * "    length = $(length(B))")
   if bodyorder(B[1]) > 1
      println(ind * " transform = $(B[1].D.s[1])")
      println(ind * "    cutoff = $(B[1].D.s[2])")
   end
end


include("fast_monomials.jl")


# ---------------  evaluate the n-body terms ------------------

# a tuple α = (α1, …, α6, α7) means the following:
# with f[0] = 1, f[1] = I7, …, f[5] = I11 we construct the basis set
#   f[α7] * g(I1, …, I6)
# this means, that gen_tuples must generate 7-tuples instead of 6-tuples
# with the first 6 entries restricted by degree and the 7th tuple must
# be in the range 0, …, 5

evaluate_I(V::NBPoly{1}, args...) = sum(V.c)

function evaluate_I(V::NBPoly, II)
   I1, I2 = II
   E = zero(eltype(I1))
   for (α, c) in zip(V.t, V.c)
      E += c * I2[1+α[end]] * monomial(α, I1)
   end
   return E
end

function evaluate_I_ed(V::NBPoly, II)
   I1, I2, dI1, dI2 = II
   E = zero(eltype(I1))
   dM = zero(typeof(I1))
   dE = zero(typeof(I1))
   #
   for (α, c) in zip(V.t, V.c)
      m, m_d = monomial_d(α, I1)
      E += c * I2[1+α[end]] * m        # just the value of the function itself
      dM += (c * I2[1+α[end]]) * m_d   # the I2 * ∇m term without the chain rule
      dE += (c * m) * dI2[1+α[end]]    # the ∇I2 * m term
   end
   # chain rule
   for i = 1:length(dI1)   # dI1' * dM
      dE += dM[i] * dI1[i]
   end
   return E, dE
end


# -------------- Infrastructure to read/write NBPoly  --------


Dict(V::NBPoly{N}) where {N} = Dict( "__id__" => "NBPoly",
                                      "t" => V.t,
                                      "c" => V.c,
                                      "D" => Dict(V.D),
                                      "N" => N )

NBPoly(D::Dict) = NBPoly([ tuple(ti...) for ti in D["t"] ],
                           Vector{Float64}(D["c"]),
                           _decode_dict(D["D"]),
                           Val(D["N"]))

Base.convert(::Val{:NBPoly}, D::Dict) = NBPoly(D)


# ==================================================================
#    StNBPoly
# ==================================================================

@pot struct StNBPoly{N, TD, TP} <: NBodyFunction{N}
   D::TD       # Descriptor
   P::TP       # a static polynomial
   valN::Val{N}
end

"""
`struct StNBPoly`  (N-Body Bond-length Polynomial)

fast evaluation of the outer polynomial using `StaticPolynomials`
"""
StNBPoly

descriptor(V::StNBPoly) = V.D

function StNBPoly(V::NBPoly{N}) where {N}
   M = bo2edges(N)  # number of edges for body-order N
   I1, I2 = invariants(V.D, SVector(ones(M)...))  # how many invariants
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
   return StNBPoly(V.D, StaticPolynomials.Polynomial(V.c, exps), V.valN)
end

cutoff(V::StNBPoly) = cutoff(V.D)

fast(Vn::StNBPoly)  = Vn
fast(Vn::NBPoly) =  StNBPoly(Vn)
fast(Vn::NBPoly{1}) = Vn


evaluate_I(V::StNBPoly, II) =
      StaticPolynomials.evaluate(V.P, vcat(II...))


function evaluate_I_ed(V::StNBPoly, II)
   V, dV_dI = StaticPolynomials.evaluate_and_gradient(V.P, vcat(II[1], II[2]))
   return V, dot(vcat(II[3], II[4]), dV_dI)  # (dI' * dV_dI)
end


# ==================================================================
#           Generate a NBPoly Basis
# ==================================================================


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


"""
* `blpolys(N, D, deg; tuplebound = ...)`
* `blpolys(N, strans, scut, deg; tuplebound = ...)`

generates a basis set of
`N`-body functions, with dictionary `D`, maximal degree `deg`; the precise
set of basis functions constructed depends on `tuplebound` (see `?gen_tuples`)
"""
blpolys(N::Integer, desc::BondLengthDesc, deg; kwargs...) =
   blpolys(gen_tuples(desc, N, deg; kwargs...), desc)

blpolys(ts::VecTup, desc::BondLengthDesc) = [NBPoly(t, 1.0, desc) for t in ts]

blpolys(N::Integer, strans::String, scut::String, deg; kwargs...) =
   blpolys(N, BondLengthDesc(strans, scut), deg; kwargs...)



end # module
