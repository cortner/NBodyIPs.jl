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
                invariants,
                edges2bo,
                bo2edges,
                ninvariants


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

@pot struct NBPoly{N, M, T, TD} <: NBodyFunction{N, TD}
   t::VecTup{M}               # tuples M = #edges + 1
   c::Vector{T}               # coefficients
   D::TD                      # Descriptor
   valN::Val{N}               # encodes that this is an N-body function

   NBPoly(t::VecTup{M}, c::Vector{T}, D::TD, valN::Val{N}) where {N, M, T, TD} = (
      N <= 1 ? error("""NBPoly must have body-order 2 or larger;
                        use `NBodyIPs.OneBody{T}` for 1-body.""")
             : new{N, M, T, TD}(t, c, D, valN))
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

* `D`: a descriptor (cf `NBodyIPs.NBodyDescriptor`)
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

@pot struct StNBPoly{N, TD, TP} <: NBodyFunction{N, TD}
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
   nI1, nI2 = ninvariants(V.D, N)
   nI = nI1 + nI2
   nmonomials = length(V.c)
   # generate the exponents for the StaticPolynomial
   exps = zeros(Int, nI, nmonomials)
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

evaluate_I(V::StNBPoly, II) =
      StaticPolynomials.evaluate(V.P, vcat(II...))

function evaluate_I_ed(V::StNBPoly, II)
   V, dV_dI = StaticPolynomials.evaluate_and_gradient(V.P, vcat(II[1], II[2]))
   if length(dV_dI) != length(II[3]) + length(II[4])
      @show size.(II)
      @show size(dV_dI)
      @show size(vcat(II[1], II[2]))
   end
   return V, dot(vcat(II[3], II[4]), dV_dI)  # (dI' * dV_dI)
end




end # module
