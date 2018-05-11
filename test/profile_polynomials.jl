#
# This short script demonstrates that we should really be using
# StaticPolynomials to evaluate the final NBodyIP
#

import StaticPolynomials

using StaticArrays, BenchmarkTools
using NBodyIPs

using NBodyIPs.Polys: monomial, monomial_d

function eval_poly(A::Vector{NTuple{7,Int}},
                   C::Vector{Float64},
                   x::SVector{12, Float64})
   E = 0.0
   x6 = x[SVector(1,2,3,4,5,6)]
   for (α, c) in zip(A, C)
      E += c * x[7+α[end]] * monomial(α, x6)
   end
   return E
end

function eval_poly_d(A::Vector{NTuple{7,Int}},
                     C::Vector{Float64},
                     x::SVector{12, Float64})
   E = zero(Float64)
   dM = zero(SVector{6, Float64})
   dE = zero(SVector{6, Float64})
   x6 = x[SVector(1,2,3,4,5,6)]
   #
   for (α, c) in zip(A, C)
      m, m_d = monomial_d(α, x6)
      E += c * x[7+α[end]] * m        # just the value of the function itself
      dM += (c * x[7+α[end]]) * m_d   # the I2 * ∇m term without the chain rule
      dE += (c * m) * x6              # a made-up line
   end
   # chain rule
   dE += dM
   return dE
end

function sp_exps(A)
   exps = zeros(Int, (length(A), 12))
   for (i, a) in enumerate(A)
      exps[i,1:6] .= a[1:6]
      exps[i,7+a[end]] = 1
   end
   return exps
end

A = NBodyIPs.gen_tuples(4, 10)
C = rand(length(A))
SP = StaticPolynomials.Polynomial(C, sp_exps(A)')

x = @SVector rand(12)

print("   hand-coded:"); @btime eval_poly($A, $C, $x)
print("  StaticPolys:"); @btime StaticPolynomials.evaluate($SP, $x)
print(" ∇-hand-coded:"); @btime eval_poly_d($A, $C, $x)
print("∇-StaticPolys:"); @btime StaticPolynomials.gradient($SP, $x)
