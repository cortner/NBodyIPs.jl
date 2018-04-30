

"""
`invariants(r::SVector{M,T}) -> SVector` : computes the invariant descriptors as a function of the
lengths in a simplex. The order is lexicographical, i.e.,

* 2-body: `r::SVector{1}`
* 3-body: `r::SVector{3}`, order is irrelevant, but `r = [r12, r13, r23]`
* 4-body: `r::SVector{6}`, order is `r = [r12, r13, r14, r23, r24, r34]`
* n-body: `r::SVector{n*(n-1)/2}`, ...

The first `n*(n-1)/2` invariants are the primary invariants; the remainind
ones are the secondary invariants.
"""
function invariants end

"""
`invariants(r::SVector{M,T}) -> SMatrix` : computes the jacobian of
`invariants`
"""
function invariants_d end

"""
`degrees(::Val{N})` where `N` is the body-order returns a
tuple of polynomials degrees corresponding to the degrees of the
individual invariants.  E.g. for 3-body, the invariants are
r1 + r2 + r3, r1 r2 + r1 r3 + r2 r3, r1 r2 r3, and the corresponding
degrees are `(1, 2, 3)`.
"""
function degrees end

"""
`bo2edges(N)` : bodyorder-to-edges
"""
bo2edges(N::Integer) = (N * (N-1)) ÷ 2

"""
`edges2bo(M)`: "edges-to-bodyorder", an internal function that translates
the number of edges in a simplex into the body-order
"""
edges2bo(M::Integer) = (M <= 0) ? 1 : round(Int, 0.5 + sqrt(0.25 + 2 * M))


# ------------------------------------------------------------------------
#             2-BODY Invariants
# ------------------------------------------------------------------------

invariants(r::SVector{1, T}) where {T} =
   copy(r), SVector{1, T}(1.0)

invariants_d(r::SVector{1, T}) where {T} =
   (@SMatrix [one(T)]), (@SMatrix [zero(T)])

degrees(::Val{2}) = (1,), (0,)

# ------------------------------------------------------------------------
#             3-BODY Invariants
# ------------------------------------------------------------------------

# the 1.0 is a "secondary invariant"
invariants(r::SVector{3, T}) where {T} =
      (@SVector T[ r[1]+r[2]+r[3],
                   r[1]*r[2] + r[1]*r[3] + r[2]*r[3],
                   r[1]*r[2]*r[3] ]),
      (@SVector T[ 1.0 ])


invariants_d(  r::SVector{3, T}) where {T} =
      (@SMatrix T[ 1.0        1.0         1.0;
                   r[2]+r[3]  r[1]+r[3]   r[1]+r[2];
                   r[2]*r[3]  r[1]*r[3]   r[1]*r[2] ]),
      (@SMatrix T[ 0.0        0.0         0.0  ])

degrees(::Val{3}) = (1, 2, 3), (0,)


# ------------------------------------------------------------------------
#             4-BODY Invariants
#
# this implementation is based on
#    Schmelzer, A., Murrell, J.N.: The general analytic expression for
#    S4-symmetry-invariant potential functions of tetra-atomic homonuclear
#    molecules. Int. J. Quantum Chem. 28, 287–295 (1985).
#    doi:10.1002/qua.560280210
#
# ------------------------------------------------------------------------
# TODO: reorder to obtain increasing degree?

degrees(::Val{4}) = (1, 2, 3, 4, 2, 3), (0, 3, 4, 5, 6, 9)

const _2 = 2.0^(-0.5)
const _3 = 3.0^(-0.5)
const _6 = 6.0^(-0.5)
const _12 = 12.0^(-0.5)

# permutation to account for the different ordering used here, vs Schmelzer et al.
const r2ρ = @SMatrix [   1 0 0 0 0 0
                         0 1 0 0 0 0
                         0 0 1 0 0 0
                         0 0 0 0 0 1
                         0 0 0 0 1 0
                         0 0 0 1 0 0 ]

const R2Q = @SMatrix [ _6     _6     _6    _6     _6     _6
                        _2      0      0   -_2      0      0
                         0     _2      0     0    -_2      0
                         0      0     _2     0      0    -_2
                         0    0.5   -0.5     0    0.5   -0.5
                        _3   -_12   -_12    _3   -_12   -_12 ]

const R2Qxr2ρ = R2Q * r2ρ

@inline invariants(r::SVector{6}) = _invariants_Q6(R2Qxr2ρ * r)

@inline function invariants_d(r::SVector{6, T}) where {T}
   J12 = ForwardDiff.jacobian( r_ -> vcat(invariants(r_)...), r )
   I1 = @SVector [1,2,3,4,5,6]
   I2 = @SVector [7,8,9,10,11,12]
   return J12[I1,:], J12[I2,:]
end

function _invariants_Q6(Q::SVector{6, T}) where {T}
   Q2 = Q .* Q
   Q2_34, Q2_24, Q2_23 = Q2[3] * Q2[4], Q2[2] * Q2[4], Q2[2] * Q2[3]
   rt3 = sqrt(3.0)
   Q_56 = Q[5] * Q[6]

   return SVector{6, T}(
      # ---------------------------- primary invariants
      # I1
      (Q[1]),
      # I2
      (Q2[2] + Q2[3] + Q2[4]),
      # I3
      (Q[2] * Q[3] * Q[4]),
      # I4
      (Q2_34 + Q2_24 + Q2_23),
      # I5
      (Q2[5] + Q2[6]),
      # I6
      (Q[6] * (Q2[6] - 3*Q2[5]))
   ),
      # ---------------------------- secondary invariants
   SVector{6, T}(
      # sneak in an additional secondary "invariant"
      1.0,
      # I7
      Q[6] * (2*Q2[2] - Q2[3] - Q2[4]) + rt3 * Q[5] * (Q2[3] - Q2[4]),
      # I8
      (Q2[6] - Q2[5]) * (2*Q2[2] - Q2[3] - Q2[4]) - 2 * rt3 * Q_56 * (Q2[3] - Q2[4]),
      # I9
      Q[6] * (2*Q2_34 - Q2_24 - Q2_23) + rt3 * Q[5] * (Q2_24 - Q2_23),
      # I10
      (Q2[6] - Q2[5])*(2*Q2_34 - Q2_24 - Q2_23) - 2 * rt3 * Q_56 * (Q2_24 - Q2_23),
      # I11
      (Q2[3] - Q2[4]) * (Q2[4] - Q2[2]) * (Q2[2] - Q2[3]) * Q[5] * (3*Q2[6] - Q2[5])
   )
end

import StaticPolynomials
using DynamicPolynomials: @polyvar

@polyvar Q1 Q2 Q3 Q4 Q5 Q6

const INV6Q = StaticPolynomials.system(
   [  Q1,
      Q2^2 + Q3^2 + Q4^2,
      Q2 * Q3 * Q4,
      Q3^2 * Q4^2 + Q2^2 * Q4^2 + Q2^2 * Q3^2,
      Q5^2 + Q6^2,
      Q6^3 - 3*Q5^2 * Q6,
      1.0,
      Q6 * (2*Q2^2 - Q3^2 - Q4^2) + √3 * Q5 * (Q3^2 - Q4^2),
      (Q6^2 - Q5^2) * (2*Q2^2 - Q3^2 - Q4^2) - 2 * √3 * Q5 * Q6 * (Q3^2 - Q4^2),
      Q6 * (2*Q3^2 * Q4^2 - Q2^2 * Q4^2 - Q2^2 * Q3^2) + √3 * Q2 * (Q2^2 * Q4^2 - Q2^2 * Q3^2),
      (Q6^2 - Q5^2)*(2*Q3^2*Q4^2 - Q2^2*Q4^2 -Q2^2*Q3^2) - 2*√3 * Q5 * Q6 * (Q2^2*Q4^2 - Q2^2*Q3^2),
      (Q3^2 - Q4^2) * (Q4^2 - Q2^2) * (Q2^2 - Q3^2) * Q5 * (3*Q6^2 - Q5^2)
   ])

@inline _invQ6_2_(Q::SVector{6}) = StaticPolynomials.evaluate(INV6Q, Q)
@inline _invQ6_2_d(Q::SVector{6}) = StaticPolynomials.jacobian(INV6Q, Q)

function sp_invariants(r::SVector{6, T}) where {T}
   I12 = _invQ6_2_(R2Qxr2ρ * r)
   I1 = @SVector [1,2,3,4,5,6]
   I2 = @SVector [7,8,9,10,11,12]
   return I12[I1], I12[I2]
end


function sp_invariants_d(r::SVector{6, T}) where {T}
   J12 = _invQ6_2_d(R2Qxr2ρ * r) * R2Qxr2ρ
   I1 = @SVector [1,2,3,4,5,6]
   I2 = @SVector [7,8,9,10,11,12]
   return J12[I1,:], J12[I2,:]
end

@polyvar x1 x2 x3 x4 x5 x6

const FINV6a = StaticPolynomials.system(
  [ x1 + x2 + x3 + x4 + x5 + x6,
    x1^2 + x2^2 + x3^2 + x4^2 + x5^2 + x6^2,
    x1*x2 + x1*x3 + x1*x4 + x1*x5 + x2*x3 + x2*x4 + x2*x6 + x3*x5 + x3*x6 + x4*x5 + x4*x6 + x5*x6,
    x1^3 + x2^3 + x3^3 + x4^3 + x5^3 + x6^3,
    x1^2*x2 + x1^2*x3 + x1^2*x4 + x1^2*x5 + x1*x2^2 + x1*x3^2 + x1*x4^2 + x1*x5^2 + x2^2*x3 + x2^2*x4 + x2^2*x6 + x2*x3^2 + x2*x4^2 + x2*x6^2 + x3^2*x5 + x3^2*x6 + x3*x5^2 + x3*x6^2 + x4^2*x5 + x4^2*x6 + x4*x5^2 + x4*x6^2 + x5^2*x6 + x5*x6^2,
    x1*x2*x3 + x1*x4*x5 + x2*x4*x6 + x3*x5*x6,
    x1^4 + x2^4 + x3^4 + x4^4 + x5^4 + x6^4,
    x1^3*x2 + x1^3*x3 + x1^3*x4 + x1^3*x5 + x1*x2^3 + x1*x3^3 + x1*x4^3 + x1*x5^3 + x2^3*x3 + x2^3*x4 + x2^3*x6 + x2*x3^3 + x2*x4^3 + x2*x6^3 + x3^3*x5 + x3^3*x6 + x3*x5^3 + x3*x6^3 + x4^3*x5 + x4^3*x6 + x4*x5^3 + x4*x6^3 + x5^3*x6 + x5*x6^3,
    x1^5 + x2^5 + x3^5 + x4^5 + x5^5 + x6^5 ])

const FINV6b = StaticPolynomials.system(
   [ x1 + x2 + x3 + x4 + x5 + x6,
   x1^2 + x2^2 + x3^2 + x4^2 + x5^2 + x6^2,
   x1*(x2+x3+x4+x5) + x2*(x3+x4+x6) + (x3+x4)*(x5+x6) + x5*x6,
   x1^3 + x2^3 + x3^3 + x4^3 + x5^3 + x6^3,
   x1^2*(x2+x3+x4+x5) + x1*(x2^2+x3^2+x4^2+x5^2) + x2^2*(x3+x4+x6) + x2*(x3^2+x4^2+x6^2) + (x3^2+x4^2)*(x5+x6) + (x3+x4)*(x5^2+x6^2) + x5^2*x6 + x5*x6^2,
   x1*x2*x3 + x1*x4*x5 + x2*x4*x6 + x3*x5*x6,
   x1^4 + x2^4 + x3^4 + x4^4 + x5^4 + x6^4,
   x1^3*(x2+x3+x4+x5) + x1*(x2^3+x3^3+x4^3+x5^3) + x2^3*(x3+x4+x6) + x2*(x3^3+x4^3+x6^3) + (x3^3+x4^3)*(x5+x6) + (x3+x4)*(x5^3+x6^3) +  x5^3*x6 + x5*x6^3,
   x1^5 + x2^5 + x3^5 + x4^5 + x5^5 + x6^5 ])

sp_fundamentals(r::SVector{6}) = StaticPolynomials.evaluate(FINV6b, r)
sp_fundamentals_d(r::SVector{6}) = StaticPolynomials.jacobian(FINV6b, r)

# import XGrad

# exF4 = quote
#    x2 = x .* x
#    x3 = x2 .* x
#    x4 = x3 .* x
#    x5 = x4 .* x
#    I1 = ((x[1]+x[2]+x[3])+x[4]+x[5])+x[6]
#    I2 = ((x2[1]+x2[2]+x2[3])+x2[4]+x2[5])+x2[6]
#    I3 = x[1]*(x[2]+x[3]+x[4]+x[5]) + x[2]*(x[3]+x[4]+x[6]) +
#          (x[3]+x[4])*(x[5]+x[6]) + x[5]*x[6]
#    I4 = ((x3[1]+x3[2]+x3[3])+x3[4]+x3[5])+x3[6]
#    I5a = x2[1]*(x[2]+x[3]+x[4]+x[5]) + x[1]*(x2[2]+x2[3]+x2[4]+x2[5])
#       I5b = I5a + x2[2]*(x[3]+x[4]+x[6]) + x[2]*(x2[3]+x2[4]+x2[6])
#       I5c = I5b + (x2[3]+x2[4])*(x[5]+x[6]) + (x[3]+x[4])*(x2[5]+x2[6])
#       I5 = I5c + x2[5]*x[6] + x[5]*x2[6]
#    I6a = x[1]*x[2]*x[3] + x[1]*x[4]*x[5]
#    I6 = I6a + x[2]*x[4]*x[6] + x[3]*x[5]*x[6]
#    I7 = ((x4[1]+x4[2]+x4[3])+x4[4]+x4[5])+x4[6]
#    I8a = x3[1]*(x[2]+x[3]+x[4]+x[5]) + x[1]*(x3[2]+x3[3]+x3[4]+x3[5])
#       I8b = I8a + x3[2]*(x[3]+x[4]+x[6]) + x[2]*(x3[3]+x3[4]+x3[6])
#       I8c = I8b + (x3[3]+x3[4])*(x[5]+x[6]) + (x[3]+x[4])*(x3[5]+x3[6])
#       I8 = I8c +  x3[5]*x[6] + x[5]*x3[6]
#    I9 = ((x5[1]+x5[2]+x5[3])+x5[4]+x5[5])+x5[6]
# end
#
# XGrad.@diffrule +(x::Real, y::Real, z::Real) x ds
# XGrad.@diffrule +(x::Real, y::Real, z::Real) y ds
# XGrad.@diffrule +(x::Real, y::Real, z::Real) z ds
# XGrad.@diffrule *(x::Real, y::Real, z::Real) x y*z*ds
# XGrad.@diffrule *(x::Real, y::Real, z::Real) y x*z*ds
# XGrad.@diffrule *(x::Real, y::Real, z::Real) z x*y*ds
#
# x_ = @SVector rand(6)
# ctx = Dict(:codegen => XGrad.VectorCodeGen())
# exF4_d = quote end
# for n = 1:9
#    exn = copy(exF4); push!(exn.args, parse("z = I$n"))
#    dexn = XGrad.xdiff(exn, ctx=ctx, x=x_.data)
#    # push!(exF4_d, dexn)
#    for subex in dexn.args
#       push!(exF4_d.args, subex)
#    end
# end
# exF4_d
#
# # joining diff expressions for the 2 components
# dex = copy(dex1)
# for subex in dex2.args
#     push!(dex.args, subex)
# end
# # auxiliary variable to return all the results
# push!(dex.args, :(result = ($(dex1.args[end].args[1]), $(dex2.args[end].args[2]))))


function fundamentals(x::SVector{6})
   x2 = x .* x
   x3 = x2 .* x
   x4 = x3 .* x
   x5 = x4 .* x
   I1 = sum(x)
   I2 = sum(x2)
   I3 = x[1]*(x[2]+x[3]+x[4]+x[5]) + x[2]*(x[3]+x[4]+x[6]) +
         (x[3]+x[4])*(x[5]+x[6]) + x[5]*x[6]
   I4 = sum(x3)
   I5 = x2[1]*(x[2]+x[3]+x[4]+x[5]) + x[1]*(x2[2]+x2[3]+x2[4]+x2[5])
      I5 += x2[2]*(x[3]+x[4]+x[6]) + x[2]*(x2[3]+x2[4]+x2[6])
      I5 += (x2[3]+x2[4])*(x[5]+x[6]) + (x[3]+x[4])*(x2[5]+x2[6])
      I5 += x2[5]*x[6] + x[5]*x2[6]
   I6 = x[1]*x[2]*x[3] + x[1]*x[4]*x[5] + x[2]*x[4]*x[6] + x[3]*x[5]*x[6]
   I7 = sum(x4)
   I8 = x3[1]*(x[2]+x[3]+x[4]+x[5]) + x[1]*(x3[2]+x3[3]+x3[4]+x3[5])
      I8 += x3[2]*(x[3]+x[4]+x[6]) + x[2]*(x3[3]+x3[4]+x3[6])
      I8 += (x3[3]+x3[4])*(x[5]+x[6]) + (x[3]+x[4])*(x3[5]+x3[6])
      I8 +=  x3[5]*x[6] + x[5]*x3[6]
   I9 = sum(x5)
   return SVector{9}(I1, I2, I3, I4, I5, I6, I7, I8, I9)
end

fundamentals_ad(x) = ForwardDiff.jacobian(fundamentals, x)


const PI43 = StaticPolynomials.Polynomial(x1*x2 + x1*x3 + x1*x4 + x1*x5 + x2*x3 + x2*x4 + x2*x6 + x3*x5 + x3*x6 + x4*x5 + x4*x6 + x5*x6)
const PI45 = StaticPolynomials.Polynomial(x1^2*(x2+x3+x4+x5) + x1*(x2^2+x3^2+x4^2+x5^2) + x2^2*(x3+x4+x6) + x2*(x3^2+x4^2+x6^2) + (x3^2+x4^2)*(x5+x6) + (x3+x4)*(x5^2+x6^2) + x5^2*x6 + x5*x6^2)
const PI46 = StaticPolynomials.Polynomial(x1*x2*x3 + x1*x4*x5 + x2*x4*x6 + x3*x5*x6)
const PI48 = StaticPolynomials.Polynomial(x1^3*(x2+x3+x4+x5) + x1*(x2^3+x3^3+x4^3+x5^3) + x2^3*(x3+x4+x6) + x2*(x3^3+x4^3+x6^3) + (x3^3+x4^3)*(x5+x6) + (x3+x4)*(x5^3+x6^3) +  x5^3*x6 + x5*x6^3)

function fundamentals_d(x::SVector{6, T}) where {T}
   x2 = x .* x
   x3 = x2 .* x
   x4 = x3 .* x
   o = SVector{6, T}(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
   ∇I3 = StaticPolynomials.gradient(PI43, x)
   ∇I5 = StaticPolynomials.gradient(PI45, x)
   ∇I6 = StaticPolynomials.gradient(PI46, x)
   ∇I8 = StaticPolynomials.gradient(PI48, x)
   return o, 2*x, ∇I3, 3*x3, ∇I5, ∇I6, 4*x3, ∇I8, 5*x4
end

using StaticArrays, BenchmarkTools
using DynamicPolynomials: @polyvar
import StaticPolynomials, XGrad


function fba(rθ::SVector{6,T}) where {T}
   Ir = SVector{3, Int}(1,2,3)
   Iθ = SVector{3, Int}(4,5,6)
   r = rθ[Ir]
   θ = rθ[Iθ]
   r2 = r .* r
   r3 = r2 .* r
   θ2 = θ .* θ
   θ3 = θ2 .* θ
   r_rev = SVector{3,T}(r[3], r[2], r[1])
   θ_rev = SVector{3,T}(θ[3], θ[2], θ[1])

   return SVector{15, T}(
      sum(r), sum(r2), sum(r3),
      sum(θ), sum(θ2), sum(θ3),
      dot(r, θ), dot(r, θ_rev),
      dot(r2, θ), dot(r2, r_rev),
      dot(r, θ2), dot(θ2, θ_rev),
       r[1]*r[2]*θ[2] + r[1]*r[3]*θ[1] + r[2]*r[3]*θ[3],
       r[1]*r2[2] + r[2]*r2[3] + r[3]*r2[1],
       θ2[1]*θ[3] + θ[1]*θ2[2] + θ[2]*θ2[3]
   )
end

function fba_d(rθ::SVector{6,T}) where {T}
   Ir = SVector{3, Int}(1,2,3)
   Iθ = SVector{3, Int}(4,5,6)
   r = rθ[Ir]
   θ = rθ[Iθ]
   r2 = r .* r
   r3 = r2 .* r
   θ2 = θ .* θ
   θ3 = θ2 .* θ
   r_rev = @SVector [r[3], r[2], r[1]]
   θ_rev = @SVector [θ[3], θ[2], θ[1]]

   return (@SMatrix [ 1.0 1.0 1.0 0.0 0.0 0.0;   # sum(r)
      2*r[1] 2*r[2] 2*r[3] 0.0 0.0 0.0;   # sum(r2)
      3*r2[1] 3*r2[2] 3*r2[3] 0.0 0.0 0.0;  # sum(r3)
      0.0 0.0 0.0 1.0 1.0 1.0;   # sum(θ)
      0.0 0.0 0.0 2*θ[1] 2*θ[2] 2*θ[3];  # sum(θ2)
      0.0 0.0 0.0 3*θ2[1] 3*θ2[2] 3*θ2[3];  # sum(θ3)
      θ[1] θ[2] θ[3] r[1] r[2] r[3];   # dot(r, θ)
      θ_rev[1] θ_rev[2] θ_rev[3] r_rev[1] r_rev[2] r_rev[3];  # dot(r, θ_rev)
      2*r[1]*θ[1] 2*r[2]*θ[2] 2*r[3]*θ[3] r2[1] r2[2] r2[3];  # dot(r2, θ)
      2*r[1]*r[3] 2*r[2]*r[2] 2*r[3]*r[1] 0.0 0.0 0.0;  #  dot(r2, r_rev)
      r2[1] r2[2] r2[3] 2*r[1]*θ[1] 2*r[2]*θ[2] 2*r[3]*θ[3];  # dot(r, θ2)
      0.0 0.0 0.0 2*θ[1]*θ[3] 2*θ[2]*θ[2] 2*θ[3]*θ[1]; # dot(θ2, θ_rev),
      r[2]*θ[2]+r[3]*θ[1] r[1]*θ[2]+r[3]*θ[3] r[1]*θ[1]+r[2]*θ[3] r[1]*r[3] r[1]*r[2] r[2]*r[3];  # r[1]*r[2]*θ[2] + r[1]*r[3]*θ[1] + r[2]*r[3]*θ[3],
      r2[2]+2*r[1]*r[3] 2*r[1]*r[2]+r2[3] 2*r[2]*r[3]+r2[1] 0.0 0.0 0.0; # r[1]*r2[2] + r[2]*r2[3] + r[3]*r2[1],
      0.0 0.0 0.0 θ2[2]+2*θ[1]*θ[3] 2*θ[1]*θ[2]+θ2[3] 2*θ[2]*θ[3]+θ2[1] ])
end



rθ = @SVector rand(6)
@btime fba($rθ)
@btime fba_d($rθ)

# fba_d = XGrad.xdiff(fba, rθ=rθ)


@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10

FIRT5a = StaticPolynomials.system([   x1 + x2 + x3 + x4, x5 + x6 + x7 + x8 + x9 + x10,
    x1^2 + x2^2 + x3^2 + x4^2,
    x1*x5 + x1*x6 + x1*x7 + x2*x5 + x2*x8 + x2*x9 + x3*x6 + x3*x8 + x3*x10 + x4*x7 + x4*x9 + x4*x10,
    x5^2 + x6^2 + x7^2 + x8^2 + x9^2 + x10^2,
    x5*x6 + x5*x7 + x5*x8 + x5*x9 + x6*x7 + x6*x8 + x6*x10 + x7*x9 + x7*x10 + x8*x9 + x8*x10 + x9*x10,
    x1^3 + x2^3 + x3^3 + x4^3,
    x1^2*x5 + x1^2*x6 + x1^2*x7 + x2^2*x5 + x2^2*x8 + x2^2*x9 + x3^2*x6 + x3^2*x8 + x3^2*x10 + x4^2*x7 + x4^2*x9 + x4^2*x10,
    x1*x2*x5 + x1*x3*x6 + x1*x4*x7 + x2*x3*x8 + x2*x4*x9 + x3*x4*x10,
    x1*x5^2 + x1*x6^2 + x1*x7^2 + x2*x5^2 + x2*x8^2 + x2*x9^2 + x3*x6^2 + x3*x8^2 + x3*x10^2 + x4*x7^2 + x4*x9^2 + x4*x10^2,
    x1*x5*x6 + x1*x5*x7 + x1*x6*x7 + x2*x5*x8 + x2*x5*x9 + x2*x8*x9 + x3*x6*x8 + x3*x6*x10 + x3*x8*x10 + x4*x7*x9 + x4*x7*x10 + x4*x9*x10,
    x5^3 + x6^3 + x7^3 + x8^3 + x9^3 + x10^3,
    (x5^2*x6 + x5^2*x7 + x5^2*x8 + x5^2*x9 + x5*x6^2 + x5*x7^2 + x5*x8^2 +
      x5*x9^2 + x6^2*x7 + x6^2*x8 + x6^2*x10 + x6*x7^2 + x6*x8^2 + x6*x10^2 +
      x7^2*x9 + x7^2*x10 + x7*x9^2 + x7*x10^2 + x8^2*x9 + x8^2*x10 + x8*x9^2 +
      x8*x10^2 + x9^2*x10 + x9*x10^2),
    x5*x6*x7 + x5*x8*x9 + x6*x8*x10 + x7*x9*x10,
    x1^4 + x2^4 + x3^4 + x4^4,
    x1^3*x5 + x1^3*x6 + x1^3*x7 + x2^3*x5 + x2^3*x8 + x2^3*x9 + x3^3*x6 + x3^3*x8 + x3^3*x10 + x4^3*x7 + x4^3*x9 + x4^3*x10,
    x1^2*x5^2 + x1^2*x6^2 + x1^2*x7^2 + x2^2*x5^2 + x2^2*x8^2 + x2^2*x9^2 + x3^2*x6^2 + x3^2*x8^2 + x3^2*x10^2 + x4^2*x7^2 + x4^2*x9^2 + x4^2*x10^2,
    x1^2*x5*x6 + x1^2*x5*x7 + x1^2*x6*x7 + x2^2*x5*x8 + x2^2*x5*x9 + x2^2*x8*x9 + x3^2*x6*x8 + x3^2*x6*x10 + x3^2*x8*x10 + x4^2*x7*x9 + x4^2*x7*x10 + x4^2*x9*x10,
    x1*x2*x5^2 + x1*x3*x6^2 + x1*x4*x7^2 + x2*x3*x8^2 + x2*x4*x9^2 + x3*x4*x10^2,
    x1*x5^3 + x1*x6^3 + x1*x7^3 + x2*x5^3 + x2*x8^3 + x2*x9^3 + x3*x6^3 + x3*x8^3 + x3*x10^3 + x4*x7^3 + x4*x9^3 + x4*x10^3,
    x1*x5^2*x6 + x1*x5^2*x7 + x1*x5*x6^2 + x1*x5*x7^2 + x1*x6^2*x7 + x1*x6*x7^2 + x2*x5^2*x8 + x2*x5^2*x9 + x2*x5*x8^2 + x2*x5*x9^2 + x2*x8^2*x9 + x2*x8*x9^2 + x3*x6^2*x8 + x3*x6^2*x10 + x3*x6*x8^2 + x3*x6*x10^2 + x3*x8^2*x10 + x3*x8*x10^2 + x4*x7^2*x9 + x4*x7^2*x10 + x4*x7*x9^2 + x4*x7*x10^2 + x4*x9^2*x10 + x4*x9*x10^2,
    x1*x5^2*x8 + x1*x5^2*x9 + x1*x6^2*x8 + x1*x6^2*x10 + x1*x7^2*x9 + x1*x7^2*x10 + x2*x5^2*x6 + x2*x5^2*x7 + x2*x6*x8^2 + x2*x7*x9^2 + x2*x8^2*x10 + x2*x9^2*x10 + x3*x5*x6^2 + x3*x5*x8^2 + x3*x6^2*x7 + x3*x7*x10^2 + x3*x8^2*x9 + x3*x9*x10^2 + x4*x5*x7^2 + x4*x5*x9^2 + x4*x6*x7^2 + x4*x6*x10^2 + x4*x8*x9^2 + x4*x8*x10^2,
    x5^4 + x6^4 + x7^4 + x8^4 + x9^4 + x10^4,
    x5^3*x6 + x5^3*x7 + x5^3*x8 + x5^3*x9 + x5*x6^3 + x5*x7^3 + x5*x8^3 + x5*x9^3 + x6^3*x7 + x6^3*x8 + x6^3*x10 + x6*x7^3 + x6*x8^3 + x6*x10^3 + x7^3*x9 + x7^3*x10 + x7*x9^3 + x7*x10^3 + x8^3*x9 + x8^3*x10 + x8*x9^3 + x8*x10^3 + x9^3*x10 + x9*x10^3,
    x1^3*x2*x5 + x1^3*x3*x6 + x1^3*x4*x7 + x1*x2^3*x5 + x1*x3^3*x6 + x1*x4^3*x7 + x2^3*x3*x8 + x2^3*x4*x9 + x2*x3^3*x8 + x2*x4^3*x9 + x3^3*x4*x10 + x3*x4^3*x10,
    x1^3*x5^2 + x1^3*x6^2 + x1^3*x7^2 + x2^3*x5^2 + x2^3*x8^2 + x2^3*x9^2 + x3^3*x6^2 + x3^3*x8^2 + x3^3*x10^2 + x4^3*x7^2 + x4^3*x9^2 + x4^3*x10^2,
    x1^2*x5^3 + x1^2*x6^3 + x1^2*x7^3 + x2^2*x5^3 + x2^2*x8^3 + x2^2*x9^3 + x3^2*x6^3 + x3^2*x8^3 + x3^2*x10^3 + x4^2*x7^3 + x4^2*x9^3 + x4^2*x10^3,
    x1*x2*x5^3 + x1*x3*x6^3 + x1*x4*x7^3 + x2*x3*x8^3 + x2*x4*x9^3 + x3*x4*x10^3,
    x1*x5^4 + x1*x6^4 + x1*x7^4 + x2*x5^4 + x2*x8^4 + x2*x9^4 + x3*x6^4 + x3*x8^4 + x3*x10^4 + x4*x7^4 + x4*x9^4 + x4*x10^4,
    x1*x5^3*x6 + x1*x5^3*x7 + x1*x5*x6^3 + x1*x5*x7^3 + x1*x6^3*x7 + x1*x6*x7^3 + x2*x5^3*x8 + x2*x5^3*x9 + x2*x5*x8^3 + x2*x5*x9^3 + x2*x8^3*x9 + x2*x8*x9^3 + x3*x6^3*x8 + x3*x6^3*x10 + x3*x6*x8^3 + x3*x6*x10^3 + x3*x8^3*x10 + x3*x8*x10^3 + x4*x7^3*x9 + x4*x7^3*x10 + x4*x7*x9^3 + x4*x7*x10^3 + x4*x9^3*x10 + x4*x9*x10^3,
    x5^5 + x6^5 + x7^5 + x8^5 + x9^5 + x10^5 ])

xx = @SVector rand(10)
@btime StaticPolynomials.evaluate($FIRT5a, $xx)
@btime StaticPolynomials.jacobian($FIRT5a, $xx)

function fba(x::SVector{10,T}) where {T}
   x2 = x .* x
   x3 = x2 .* x
   x4 = x3 .* x
   x5 = x4 .* x
   return SVector(
    x[1] + x[2] + x[3] + x[4],
    x[5] + x[6] + x[7] + x[8] + x[9] + x[10],
    x2[1] + x2[2] + x2[3] + x2[4],
    ((x[1]*x[5] + x[1]*x[6] + x[1]*x[7] + x[2]*x[5] + x[2]*x[8] + x[2]*x[9]) +
      (x[3]*x[6] + x[3]*x[8] + x[3]*x[10] + x[4]*x[7] + x[4]*x[9] + x[4]*x[10])),
    x2[5] + x2[6] + x2[7] + x2[8] + x2[9] + x2[10],
    ((x[5]*x[6] + x[5]*x[7] + x[5]*x[8] + x[5]*x[9] + x[6]*x[7] + x[6]*x[8]) +
      (x[6]*x[10] + x[7]*x[9] + x[7]*x[10] + x[8]*x[9] + x[8]*x[10] + x[9]*x[10])),
    x3[1] + x3[2] + x3[3] + x3[4],
    ((x2[1]*x[5] + x2[1]*x[6] + x2[1]*x[7] + x2[2]*x[5] + x2[2]*x[8] + x2[2]*x[9]) +
      (x2[3]*x[6] + x2[3]*x[8] + x2[3]*x[10] + x2[4]*x[7] + x2[4]*x[9] + x2[4]*x[10])),
    x[1]*x[2]*x[5] + x[1]*x[3]*x[6] + x[1]*x[4]*x[7] + x[2]*x[3]*x[8] + x[2]*x[4]*x[9] + x[3]*x[4]*x[10],
    ((x[1]*x2[5] + x[1]*x2[6] + x[1]*x2[7] + x[2]*x2[5] + x[2]*x2[8] + x[2]*x2[9]) +
      (x[3]*x2[6] + x[3]*x2[8] + x[3]*x2[10] + x[4]*x2[7] + x[4]*x2[9] + x[4]*x2[10])),
    ((x[1]*x[5]*x[6] + x[1]*x[5]*x[7] + x[1]*x[6]*x[7] + x[2]*x[5]*x[8] + x[2]*x[5]*x[9] + x[2]*x[8]*x[9]) +
     (x[3]*x[6]*x[8] + x[3]*x[6]*x[10] + x[3]*x[8]*x[10] + x[4]*x[7]*x[9] + x[4]*x[7]*x[10] + x[4]*x[9]*x[10])),
    x3[5] + x3[6] + x3[7] + x3[8] + x3[9] + x3[10],
    ((x2[5]*x[6] + x2[5]*x[7] + x2[5]*x[8] + x2[5]*x[9] + x[5]*x2[6] + x[5]*x2[7]) +
     (x[5]*x2[8] + x[5]*x2[9] + x2[6]*x[7] + x2[6]*x[8] + x2[6]*x[10] + x[6]*x2[7]) +
     (x[6]*x2[8] + x[6]*x2[10] + x2[7]*x[9] + x2[7]*x[10] + x[7]*x2[9] + x[7]*x2[10]) +
     (x2[8]*x[9] + x2[8]*x[10] + x[8]*x2[9] + x[8]*x2[10] + x2[9]*x[10] + x[9]*x2[10])),
    x[5]*x[6]*x[7] + x[5]*x[8]*x[9] + x[6]*x[8]*x[10] + x[7]*x[9]*x[10],
    x4[1] + x4[2] + x4[3] + x4[4],
    ((x3[1]*x[5] + x3[1]*x[6] + x3[1]*x[7] + x3[2]*x[5] + x3[2]*x[8] + x3[2]*x[9]) +
       (x3[3]*x[6] + x3[3]*x[8] + x3[3]*x[10] + x3[4]*x[7] + x3[4]*x[9] + x3[4]*x[10])),
    # x2[1]*x2[5] + x2[1]*x2[6] + x2[1]*x2[7] + x2[2]*x2[5] + x2[2]*x2[8] + x2[2]*x2[9] + x2[3]*x2[6] + x2[3]*x2[8] + x2[3]*x2[10] + x2[4]*x2[7] + x2[4]*x2[9] + x2[4]*x2[10],
    # x2[1]*x[5]*x[6] + x2[1]*x[5]*x[7] + x2[1]*x[6]*x[7] + x2[2]*x[5]*x[8] + x2[2]*x[5]*x[9] + x2[2]*x[8]*x[9] + x2[3]*x[6]*x[8] + x2[3]*x[6]*x[10] + x2[3]*x[8]*x[10] + x2[4]*x[7]*x[9] + x2[4]*x[7]*x[10] + x2[4]*x[9]*x[10],
    # x[1]*x[2]*x2[5] + x[1]*x[3]*x2[6] + x[1]*x[4]*x2[7] + x[2]*x[3]*x2[8] + x[2]*x[4]*x2[9] + x[3]*x[4]*x2[10],
    # x[1]*x3[5] + x[1]*x3[6] + x[1]*x3[7] + x[2]*x3[5] + x[2]*x3[8] + x[2]*x3[9] + x[3]*x3[6] + x[3]*x3[8] + x[3]*x3[10] + x[4]*x3[7] + x[4]*x3[9] + x[4]*x3[10],
    # x[1]*x2[5]*x[6] + x[1]*x2[5]*x[7] + x[1]*x[5]*x2[6] + x[1]*x[5]*x2[7] + x[1]*x2[6]*x[7] + x[1]*x[6]*x2[7] + x[2]*x2[5]*x[8] + x[2]*x2[5]*x[9] + x[2]*x[5]*x2[8] + x[2]*x[5]*x2[9] + x[2]*x2[8]*x[9] + x[2]*x[8]*x2[9] + x[3]*x2[6]*x[8] + x[3]*x2[6]*x[10] + x[3]*x[6]*x2[8] + x[3]*x[6]*x2[10] + x[3]*x2[8]*x[10] + x[3]*x[8]*x2[10] + x[4]*x2[7]*x[9] + x[4]*x2[7]*x[10] + x[4]*x[7]*x2[9] + x[4]*x[7]*x2[10] + x[4]*x2[9]*x[10] + x[4]*x[9]*x2[10],
    # x[1]*x2[5]*x[8] + x[1]*x2[5]*x[9] + x[1]*x2[6]*x[8] + x[1]*x2[6]*x[10] + x[1]*x2[7]*x[9] + x[1]*x2[7]*x[10] + x[2]*x2[5]*x[6] + x[2]*x2[5]*x[7] + x[2]*x[6]*x2[8] + x[2]*x[7]*x2[9] + x[2]*x2[8]*x[10] + x[2]*x2[9]*x[10] + x[3]*x[5]*x2[6] + x[3]*x[5]*x2[8] + x[3]*x2[6]*x[7] + x[3]*x[7]*x2[10] + x[3]*x2[8]*x[9] + x[3]*x[9]*x2[10] + x[4]*x[5]*x2[7] + x[4]*x[5]*x2[9] + x[4]*x[6]*x2[7] + x[4]*x[6]*x2[10] + x[4]*x[8]*x2[9] + x[4]*x[8]*x2[10],
    # x4[5] + x4[6] + x4[7] + x4[8] + x4[9] + x4[10],
    # x3[5]*x[6] + x3[5]*x[7] + x3[5]*x[8] + x3[5]*x[9] + x[5]*x3[6] + x[5]*x3[7] + x[5]*x3[8] + x[5]*x3[9] + x3[6]*x[7] + x3[6]*x[8] + x3[6]*x[10] + x[6]*x3[7] + x[6]*x3[8] + x[6]*x3[10] + x3[7]*x[9] + x3[7]*x[10] + x[7]*x3[9] + x[7]*x3[10] + x3[8]*x[9] + x3[8]*x[10] + x[8]*x3[9] + x[8]*x3[10] + x3[9]*x[10] + x[9]*x3[10],
    # x3[1]*x[2]*x[5] + x3[1]*x[3]*x[6] + x3[1]*x[4]*x[7] + x[1]*x3[2]*x[5] + x[1]*x3[3]*x[6] + x[1]*x3[4]*x[7] + x3[2]*x[3]*x[8] + x3[2]*x[4]*x[9] + x[2]*x3[3]*x[8] + x[2]*x3[4]*x[9] + x3[3]*x[4]*x[10] + x[3]*x3[4]*x[10],
    # x3[1]*x2[5] + x3[1]*x2[6] + x3[1]*x2[7] + x3[2]*x2[5] + x3[2]*x2[8] + x3[2]*x2[9] + x3[3]*x2[6] + x3[3]*x2[8] + x3[3]*x2[10] + x3[4]*x2[7] + x3[4]*x2[9] + x3[4]*x2[10],
    # x2[1]*x3[5] + x2[1]*x3[6] + x2[1]*x3[7] + x2[2]*x3[5] + x2[2]*x3[8] + x2[2]*x3[9] + x2[3]*x3[6] + x2[3]*x3[8] + x2[3]*x3[10] + x2[4]*x3[7] + x2[4]*x3[9] + x2[4]*x3[10],
    # x[1]*x[2]*x3[5] + x[1]*x[3]*x3[6] + x[1]*x[4]*x3[7] + x[2]*x[3]*x3[8] + x[2]*x[4]*x3[9] + x[3]*x[4]*x3[10],
    # x[1]*x4[5] + x[1]*x4[6] + x[1]*x4[7] + x[2]*x4[5] + x[2]*x4[8] + x[2]*x4[9] + x[3]*x4[6] + x[3]*x4[8] + x[3]*x4[10] + x[4]*x4[7] + x[4]*x4[9] + x[4]*x4[10],
    # x[1]*x3[5]*x[6] + x[1]*x3[5]*x[7] + x[1]*x[5]*x3[6] + x[1]*x[5]*x3[7] + x[1]*x3[6]*x[7] + x[1]*x[6]*x3[7] + x[2]*x3[5]*x[8] + x[2]*x3[5]*x[9] + x[2]*x[5]*x3[8] + x[2]*x[5]*x3[9] + x[2]*x3[8]*x[9] + x[2]*x[8]*x3[9] + x[3]*x3[6]*x[8] + x[3]*x3[6]*x[10] + x[3]*x[6]*x3[8] + x[3]*x[6]*x3[10] + x[3]*x3[8]*x[10] + x[3]*x[8]*x3[10] + x[4]*x3[7]*x[9] + x[4]*x3[7]*x[10] + x[4]*x[7]*x3[9] + x[4]*x[7]*x3[10] + x[4]*x3[9]*x[10] + x[4]*x[9]*x3[10],
    # x5[5] + x5[6] + x5[7] + x5[8] + x5[9] + x5[10]
    )
end

@btime fba($r10)
r10 = @SVector rand(10)
@btime fba($r10)

# ------------------------------------------------------------------------
#             5-BODY Invariants
# ------------------------------------------------------------------------

# TODO

# COPIED FROM SLACK:
# For the 5-body, if I am right, the number of primary invariants is 10 with
# degrees : [ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6 ]. The number of secondary for each
# degree (starting from 0) is [ 1, 0, 0, 2, 5, 8, 15, 23, 33, 46 ] which means
# in total 133 secondary of degree less than 9. That’s quite a lot. Some of
# them have very long expressions.

degrees(::Val{5}) = (1, 2, 2, 3, 3, 4, 4, 5, 5, 6),
   (0,     # 1 x degree 0
    3, 3,  # 2 x degree 3
    4, 4, 4, 4, 4,   # 5 x degree 4
    5, 5, 5, 5, 5, 5, 5, 5,  # 8 x degree 5
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)   # 15 x degree 6
