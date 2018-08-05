

import NBodyIPs: evaluate_many!,
                 evaluate_many_d!,
                 invariants,
                 invariants_d,
                 invariants_ed

import JuLIP.Potentials: evaluate,
                         evaluate_d


include("fast_monomials.jl")


# ---------------  evaluate the n-body terms ------------------

# a tuple α = (α1, …, α6, α7) means the following:
# with f[0] = 1, f[1] = I7, …, f[5] = I11 we construct the basis set
#   f[α7] * g(I1, …, I6)
# this means, that gen_tuples must generate 7-tuples instead of 6-tuples
# with the first 6 entries restricted by degree and the 7th tuple must
# be in the range 0, …, 5

evaluate_I(V::NBPoly{1}, args...) = sum(V.c)


evaluate_I(V::NBPoly, I1, I2, fc) = sum( c * I2[1+α[end]] * monomial(α, I1)
                                         for  (α, c) in zip(V.t, V.c) ) * fc


function evaluate_I_d(V::NBPoly, I1::SVector{M,T}, I2, I1_d, I2_d, fc, fc_d
                      ) where {M, T}
   E = zero(T)
   dM = zero(SVector{M, T})
   dE = zero(SVector{M, T})
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
   fc, fc_d = fcut_d(D, r)
   return dE * fc + E * fc_d
end


evaluate_I(V::StNBPoly, I1, I2, fc) =
      StaticPolynomials.evaluate(V.P, vcat(I1, I2)) * fc


function evaluate_I_d(V::StNBPoly, I1::SVector{M, T}, I2, I1_d, I2_d,
                                   fc, fc_dr) where {M, T}
   V, dV_dI = StaticPolynomials.evaluate_and_gradient(V.P, vcat(I1, I2))
   return V * fc_d + fc * dot(vcat(dI1, dI2), dV_dI)  # (dI' * dV_dI)
end




# function evaluate_many!(temp, B::Vector{TB}, r::SVector{M, T}
#                ) where {TB <: NBPoly{N}} where {N, M, T}
#    E = temp # zeros(T, length(B))
#    D = B[1].D
#    # it is assumed implicitly that all basis functions use the same dictionary!
#    # and that each basis function NBPoly contains just a single c and a single t
#    I1, I2 = invariants(D, r)
#    for ib = 1:length(B)
#       E[ib] = B[ib].c[1] * I2[1+B[ib].t[1][end]] * monomial(B[ib].t[1], I1)
#    end
#    fc = fcut(D, r)
#    scale!(E, fc)
#    return E
# end
#
#
# function evaluate_many_d!(temp, B::Vector{TB}, r::SVector{M, T}
#                ) where {TB <: NBPoly{N}} where {N, M, T}
#    E, dM, dE, dEfinal = temp
#    fill!(E, 0.0)     # size (nB,)
#    fill!(dE, 0.0)    # size (M, nB)
#    fill!(dM, 0.0)    # size (M, nB)
#    nB = length(B)
#
#    D = B[1].D
#    # it is assumed implicitly that all basis functions use the same dictionary!
#    # and that each basis function NBPoly contains just a single c and a single t
#    I1, I2, dI1, dI2 = invariants_ed(D, r)
#    for ib = 1:nB    # (ib, b) in enumerate(B)
#       α = B[ib].t[1]
#       c = B[ib].c[1]
#       m, m_d = monomial_d(α, I1)
#       E[ib] += c * I2[1+α[end]] * m        # just the value of the function itself
#       for t = 1:M
#          dM[t,ib] += (c * I2[1+α[end]]) * m_d[t]   # the I2 * ∇m term without the chain rule
#          dE[t,ib] += (c * m) * dI2[1+α[end]][t]  # the ∇I2 * m term
#       end
#    end
#    # chain rule
#    fc, fc_d = fcut_d(D, r)
#    for ib = 1:nB
#       for i = 1:M, t = 1:M  # dE[ib] += dI1' * dM[ib]
#          dE[t,ib] += dI1[i][t] * dM[i,ib]
#       end
#       for t = 1:M
#          dE[t,ib] = dE[t,ib] * fc + E[ib] * fc_d[t]
#       end
#    end
#    # write into an output vector
#    for i = 1:nB, t = 1:M
#       dEfinal[i][t] = dE[t, i]
#    end
#    return dEfinal
# end
