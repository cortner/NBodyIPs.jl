
"""
`module Polynomials`

The main types are
* `Dictionary`: collects all the information about a basis
* `NBody`: an N-body function wrapped into a JuLIP calculator
* `NBodyIP`: a collection of N-body functions (of possibly different
body-order and cut-off) wrapped into a JuLIP calculator.

## Notation

* `N` : N-body
* `M` : number of edges, i.e., M = N (N-1)/2
* `K` : length of the tuples defining polynomials, K = M+1

## Examples

"""
module Polynomials

using Reexport

using JuLIP, NeighbourLists, StaticArrays, ForwardDiff, Calculus

import JuLIP.Potentials: @analytic
using JuLIP.Potentials: cutsw, cutsw_d, coscut, coscut_d
const cutsp = JuLIP.Potentials.fcut
const cutsp_d = JuLIP.Potentials.fcut_d

import Base: length
import JuLIP: cutoff, energy, forces
import JuLIP.Potentials: evaluate, evaluate_d

const CRg = CartesianRange
const CInd = CartesianIndex
const Tup{M} = NTuple{M, Int}
const VecTup{M} = Vector{NTuple{M, Int}}

export NBody, NBodyIP, Dictionary,
       gen_tuples, gen_basis

# somehow @reexport doesn't catch this one
export @analytic



@pot struct Dictionary{TT, TDT, TC, TDC, T}
   transform::TT              # distance transform
   transform_d::TDT            # derivative of distance transform
   fcut::TC                   # cut-off function
   fcut_d::TDC                 # cut-off function derivative
   rcut::T                    # cutoff radius
end


# generated functions for fast evaluation of monomials
include("fast_monomials.jl")

# import the raw symmetry invariant polynomials
include("invariants.jl")


# ==================================================================
#           DICTIONARY
# ==================================================================


"""
`struct Dictionary` : specifies all details about the basis functions

## Convenience Constructor

`Dictionary(ftrans, fcut[, rcut])`, where `ftrans`, `fcut` specify in one of several
ways how the dictionary is defined and `rcut` is the cutoff radius. If
the radius is not provided then `Dictionary` will try to infer it from the
`fcut` argument. E.g.,
```julia
D = Dictionary((@analytic r -> 1/r), (r -> (r-rcut)^2), rcut)
D = Dictionary( (:poly, -1), (:square, rcut) )
D = Dictionary( (:exp, 3.0), (:sw, L, rcut) )
```
Known symbols for the transformation are
```
[:poly, :exp]
```
Known symbols for the cutoff are
```
[:cos, :sw, :spline, :square]
```


## Developer Doc: Methods associated with a `D::Dictionary`:

* `invariants`, `jac_invariants`: compute invariants and jacobian
   in transformed coordinates defined by `D.transform`
* `evaluate`, `evaluate_d`: evaluate the (univariate) basis
   function associated with this dictionary; at the moment only
   standard polynomials are admitted
* `fcut, fcut_d`: evulate the cut-off function and its derivative / gradient
   when interpreted as a product of cut-off functions
"""
Dictionary

@inline invariants(D::Dictionary, r) = invariants(D.transform.(r))
@inline function invariants_d(D::Dictionary, r)
   DI1, DI2 = invariants_d(D.transform.(r))
   t_d = D.transform_d.(r)'
   return DI1 .* t_d, DI2 .* t_d
end

# @inline evaluate(D::Dictionary, α, r) = _m1(α, r)
# @inline evaluate_d(D::Dictionary, α, r) = _m1d(α, r)


@inline fcut(D::Dictionary, r::Number) = D.fcut(r)
@inline fcut_d(D::Dictionary, r::Number) = D.fcut_d(r)

@inline cutoff(D::Dictionary) = D.rcut


# -------------------- generate dictionaries -------------------

Dictionary(ftrans::AnalyticFunction, fcut::AnalyticFunction, rcut::AbstractFloat) =
      Dictionary(ftrans.f, ftrans.f_d, fcut.f, fcut.f_d, rcut)

Dictionary(ftrans::Any, fcut::Any) =
      Dictionary(ftrans_analyse(ftrans), fcut_analyse(fcut)...)

ftrans_analyse(x::AnalyticFunction) = x

function ftrans_analyse(args)
   if args isa Symbol || length(args) == 1
      sym = args
      p = nothing
   elseif length(args) == 2
      sym, p = args
   else
      sym, p = args[1], args[2:end]
   end
   if Symbol(sym) == :inv
      return @analytic r -> 1/r
   elseif Symbol(sym) == :invsqrt
      return @analytic r -> 1/sqrt(r)
   elseif Symbol(sym) == :invsquare
      return @analytic r -> (1/r)^2
   elseif Symbol(sym) == :poly
      return let p=p; @analytic r -> r^p; end
   elseif Symbol(sym) == :exp
      return let p=p; @analytic r -> exp(-p * r); end
   else
      error("Dictionary: unknown symbol $(sym) for transformation.")
   end
end

fcut_analyse(x::AnalyticFunction) = x

function fcut_analyse(args::Tuple)
   sym = args[1]::Symbol
   args = args[2:end]
   if Symbol(sym) == :cos
      rc1, rc2 = args
      return let rc1=rc1, rc2=rc2
                  AnalyticFunction( r -> coscut(r, rc1, rc2),
                                    r -> coscut_d(r, rc1, rc2),
                                    nothing ); end, rc2
   elseif Symbol(sym) == :sw
      L, rcut = args
      return let L=L, rcut=rcut
                  AnalyticFunction( r -> cutsw(r, rcut, L),
                                    r -> cutsw_d(r, rcut, L),
                                    nothing ); end, rcut
   elseif Symbol(sym) == :spline
      rc1, rc2 = args
      return let rc1=rc1, rc2=rc2
                  AnalyticFunction( r -> cutsp(r, rc1, rc2),
                                    r -> cutsp_d(r, rc1, rc2),
                                    nothing );end , rc2
   elseif Symbol(sym) == :square
      rcut = args[1]
      return let rcut=rcut; (@analytic r -> (r - rcut)^2); end, rcut
   else
      error("Dictionary: unknown symbol $(sym) for fcut.")
   end
end


# ==================================================================
#           Polynomials of Invariants
# ==================================================================


@pot struct NBody{N, M, T, TD} <: AbstractCalculator
   t::VecTup{M}               # tuples M = #edges + 1
   c::Vector{T}               # coefficients
   D::TD                      # Dictionary (or nothing)
   valN::Val{N}               # encodes that this is an N-body term
end

"""
`struct NBody`

A struct storing the information for a pure N-body potential, i.e., containing
*only* terms of a specific body-order. Several `NBody`s can be
combined into an interatomic potential via `NBodyIP`.

### Fields

* `t::Vector{NTuple{M,TI}}` : list of M-tuples containing basis function information
e.g., if M = 3, α = t[1] is a 3-vector then this corresponds to the basis function
`f[α[1]](Q[1]) * f[α[2]](Q[2]) * f[α[3]](Q[3])` where `Q` are the 3-body invariants.

* `c::Vector{T}`: vector of coefficients for the basis functions

* `D`: a `Dictionary`
"""
NBody

# standad constructor (N can be inferred)
NBody(t::VecTup{K}, c, D) where {K} =
      NBody(t, c, D, Val(edges2bo(K-1)))

# NBody made from a single basis function rather than a collection
NBody(t::Tup, c, D) =
      NBody([t], [c], D)

# collect multiple basis functions represented as NBody's into a single NBody
# (mostly for performance reasons)
NBody(B::Vector{TB}, c, D) where {TB <: NBody} =
      NBody([b.t[1] for b in B], c, D)

# 1-body term (on-site energy)
NBody(c::Float64) =
      NBody([Tup{0}()], [c], nothing, Val(1))

bodyorder(V::NBody{N}) where {N} = N
dim(V::NBody{N,M}) where {N, M} = M-1
# number of basis functions which this term is made from
length(V::NBody) = length(V.t)

# --------------------- Calculator Functionality for NBody
# TODO: site energies

cutoff(V::NBody) = cutoff(V.D)

# energy and forces of the 0-body term
energy(V::NBody{1,0,T}, at::Atoms{T}) where {T} =
      sum(V.c) * length(at)

forces(V::NBody{1,0,T}, at::Atoms{T}) where {T} =
      zeros(SVector{3, T}, length(at))


function energy(V::NBody{N, M, T}, at::Atoms{T}) where {N, M, T}
   nlist = neighbourlist(at, cutoff(V))
   Es = maptosites!(r -> V(r), zeros(length(at)), nbodies(N, nlist))
   return sum_kbn(Es)
end

function forces(V::NBody{N, M, T}, at::Atoms{T}) where {N, M, T}
   nlist = neighbourlist(at, cutoff(V))
   return scale!(maptosites_d!(r -> (@D V(r)),
                 zeros(SVector{3, T}, length(at)),
                 nbodies(N, nlist)), -1)
end

function forces(V::NBody{4, M, T}, at::Atoms{T}) where {M, T}
   nlist = neighbourlist(at, cutoff(V))
   return scale!(maptosites_d!(r -> (@D V(r)),
                 zeros(SVector{3, T}, length(at)),
                 nbodies(4, nlist)), -1)
end


# ---------------  evaluate the n-body terms ------------------

# a tuple α = (α1, …, α6, α7) means the following:
# with f[0] = 1, f[1] = I7, …, f[5] = I11 we construct the basis set
#   f[α7] * g(I1, …, I6)
# this means, that gen_tuples must generate 7-tuples instead of 6-tuples
# with the first 6 entries restricted by degree and the 7th tuple must
# be in the range 0, …, 5


function evaluate(V::NBody, r::SVector{M, T}) where {M, T}
   # @assert ((N*(N-1))÷2 == M == K-1 == dim(V) == M)
   E = zero(T)
   I1, I2 = invariants(V.D, r)         # SVector{NI, T}
   for (α, c) in zip(V.t, V.c)
      E += c * I2[1+α[end]] * monomial(α, I1)
   end
   return E * fcut(V.D, r)
end

# without AD
function evaluate_d(V::NBody, r::SVector{M, T}) where {M, T}
   E = zero(T)
   dM = zero(SVector{M, T})
   dE = zero(SVector{M, T})
   D = V.D
   I1, I2 = invariants(D, r)              # SVectors
   dI1, dI2 = invariants_d(D, r)          # SMatrices
   #
   for (α, c) in zip(V.t, V.c)
      m, m_d = monomial_d(α, I1)
      E += c * I2[1+α[end]] * m        # just the value of the function itself
      dM += (c * I2[1+α[end]]) * m_d   # the I2 * ∇m term without the chain rule
      dE += (c * m) * dI2[1+α[end],:]  # the ∇I2 * m term
   end
   # chain rule
   dE += dI1' * dM

   fc, fc_d = fcut_d(D, r)
   return dE * fc + E * fc_d
end



# ==================================================================
#           The Final Interatomic Potential
# ==================================================================
# TODO: site energies

"""
`NBodyIP` : wraps `NBody`s into a JuLIP calculator, defining
`energy`, `forces` and `cutoff`.
"""
struct NBodyIP <: AbstractCalculator
   orders::Vector{NBody}
end

cutoff(V::NBodyIP) = maximum( cutoff.(V.orders) )
energy(V::NBodyIP, at::Atoms) = sum( energy(Vn, at)  for Vn in V.orders )
forces(V::NBodyIP, at::Atoms) = sum( forces(Vn, at)  for Vn in V.orders )

function NBodyIP(basis, coeffs)
   orders = NBody[]
   bos = bodyorder.(basis)
   for N = 1:maximum(bos)
      Ibo = find(bos .== N)  # find all basis functions that have the right bodyorder
      if length(Ibo) > 0
         D = basis[Ibo[1]].D
         V_N = NBody(basis[Ibo], coeffs[Ibo], D)
         push!(orders, V_N)  # collect them
      end
   end
   return NBodyIP(orders)
end



# ==================================================================
#           Generate a Basis
# ==================================================================

nedges(::Val{N}) where {N} = (N*(N-1)) ÷ 2

"""
compute the total degree of the polynomial represented by α.
Note that `M = K-1` where `K` is the tuple length while
`M` is the number of edges.
"""
function degree(α::Tup{K}) where {K}
   degs1, degs2 = degrees(Val(edges2bo(K-1)))
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
`α -> (degree(α) <= deg)` i.e. the standard monomial degree. (note this is
the degree w.r.t. lengths, not w.r.t. invariants!)
"""
gen_tuples(N, deg; tuplebound = (α -> (0 < degree(α) <= deg))) =
   gen_tuples(Val(N), Val(nedges(Val(N))+1), deg, tuplebound)

# ------------- 2-body tuples -------------

# TODO: need to eventually generate 2-tuples
gen_tuples(vN::Val{2}, vK::Val{2}, deg, tuplebound) =
   Tup{2}[ (j, 0)   for j = 1:deg ]

# ------------- 3-body tuples -------------

function gen_tuples(vN::Val{N}, vK::Val{K}, deg, tuplebound) where {N, K}
   t = Tup{K}[]
   degs1, degs2 = degrees(vN)
   Ilo = CInd{K}(zeros(Int, K)...)
   idegs = [ceil.(Int, deg ./ [degs1...]); length(degs2)-1]
   Ihi = CInd{K}(idegs...)
   for I in CRg(Ilo, Ihi)
      if tuplebound(I.I)
         push!(t, I.I)
      end
   end
   return t
end

# ------------- 4-body tuples -------------


"""
`gen_basis(N, D, deg; tuplebound = ...)` : generates a basis set of
`N`-body functions, with dictionary `D`, maximal degree `deg`; the precise
set of basis functions constructed depends on `tuplebound` (see `?gen_tuples`)
"""
gen_basis(N, D, deg; kwargs...) = gen_basis(gen_tuples(N, deg; kwargs...), D)

gen_basis(ts::VecTup, D::Dictionary) = [NBody(t, 1.0, D) for t in ts]

end
