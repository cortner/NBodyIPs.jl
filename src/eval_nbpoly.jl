

import NBodyIPs: evaluate_many!,
                 evaluate_many_d!,
                 evaluate_I,
                 evaluate_I_d


include("fast_monomials.jl")


# ---------------  evaluate the n-body terms ------------------

# a tuple α = (α1, …, α6, α7) means the following:
# with f[0] = 1, f[1] = I7, …, f[5] = I11 we construct the basis set
#   f[α7] * g(I1, …, I6)
# this means, that gen_tuples must generate 7-tuples instead of 6-tuples
# with the first 6 entries restricted by degree and the 7th tuple must
# be in the range 0, …, 5

evaluate_I(V::NBPoly{1}, args...) = sum(V.c)


function evaluate_I(V::NBPoly, I1, I2, fc)
   E = zero(eltype(I1))
   for (α, c) in zip(V.t, V.c)
      E += c * I2[1+α[end]] * monomial(α, I1)
   end
   return E * fc
end

function evaluate_I_d(V::NBPoly, I1::SVector{M,T}, I2, dI1, dI2, fc, fc_d
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
   return dE * fc + E * fc_d
end


evaluate_I(V::StNBPoly, I1, I2, fc) =
      StaticPolynomials.evaluate(V.P, vcat(I1, I2)) * fc


function evaluate_I_d(V::StNBPoly, I1::SVector{M, T}, I2, dI1, dI2,
                                   fc, fc_d) where {M, T}
   V, dV_dI = StaticPolynomials.evaluate_and_gradient(V.P, vcat(I1, I2))
   return V * fc_d + fc * dot(vcat(dI1, dI2), dV_dI)  # (dI' * dV_dI)
end
