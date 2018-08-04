

using NBodyIPs.FastPolys: fpoly,
                          fpoly_d

import NBodyIPs: evaluate_many!,
                 evaluate_many_d!

import JuLIP.Potentials: evaluate,
                         evaluate_d

import NBodyIPs.BLPolys.BLInvariants: invariants,
                                      invariants_d,
                                      invariants_ed


include("fast_monomials.jl")

@generated function _sdot(a::T, B::SVector{N, T}) where {N, T}
   code = "@SVector $T["
   for n = 1:N
      code *= "a .* B[$n],"
   end
   code = code[1:end-1] * "]"
   ex = parse(code)
   quote
      $ex
   end
end

@inline invariants(D::BLDictionary, r) = invariants(transform.(D, r))

@inline function invariants_d(D::BLDictionary, r)
   DI1, DI2 = invariants_d(transform.(D, r))
   t_d = transform_d.(D, r)
   return _sdot(t_d, DI1), _sdot(t_d, DI2)
end

@inline function invariants_ed(D::BLDictionary, r)
   t = transform.(D, r)
   I1, I2, DI1, DI2 = invariants_ed(t)
   t_d = transform_d.(D, r)
   return I1, I2, _sdot(t_d, DI1), _sdot(t_d, DI2)
end


# ---------------  evaluate the n-body terms ------------------

# a tuple α = (α1, …, α6, α7) means the following:
# with f[0] = 1, f[1] = I7, …, f[5] = I11 we construct the basis set
#   f[α7] * g(I1, …, I6)
# this means, that gen_tuples must generate 7-tuples instead of 6-tuples
# with the first 6 entries restricted by degree and the 7th tuple must
# be in the range 0, …, 5

evaluate(V::BLNBody{1}) = sum(V.c)

function evaluate(V::BLNBody, r::SVector{M, T}) where {M, T}
   # @assert ((N*(N-1))÷2 == M == K-1 == dim(V) == M)
   E = zero(T)
   I1, I2 = invariants(V.D, r)         # SVector{NI, T}
   for (α, c) in zip(V.t, V.c)
      E += c * I2[1+α[end]] * monomial(α, I1)
   end
   return E * fcut(V.D, r)
end


function evaluate_d(V::BLNBody, r::SVector{M, T}) where {M, T}
   E = zero(T)
   dM = zero(SVector{M, T})
   dE = zero(SVector{M, T})
   D = V.D
   I1, I2, dI1, dI2 = invariants_ed(D, r)
   #
   for (α, c) in zip(V.t, V.c)
      m, m_d = monomial_d(α, I1)
      E += c * I2[1+α[end]] * m        # just the value of the function itself
      dM += (c * I2[1+α[end]]) * m_d   # the I2 * ∇m term without the chain rule
      dE += (c * m) * dI2[1+α[end]]  # the ∇I2 * m term
   end
   # chain rule
   for i = 1:length(dI1)   # dI1' * dM
      dE += dM[i] * dI1[i]
   end
   fc, fc_d = fcut_d(D, r)
   return dE * fc + E * fc_d
end


function evaluate(V::StBLNBody, r::SVector{M, T}) where {M, T}
   I = vcat(invariants(V.D, r)...)
   E = StaticPolynomials.evaluate(V.P, I)
   return E * fcut(V.D, r)
end

function evaluate_d(V::StBLNBody, r::SVector{M, T}) where {M, T}
   D = V.D
   # I = vcat(invariants(D, r)...)   # TODO: combine into a single evaluation
   # dI = vcat(invariants_d(D, r)...)
   I1, I2, dI1, dI2 = invariants_ed(D, r)
   I = vcat(I1, I2)
   dI = vcat(dI1, dI2)
   V, dV_dI = StaticPolynomials.evaluate_and_gradient(V.P, I)
   fc, fc_d = fcut_d(D, r)
   return V * fc_d + fc * dot(dI, dV_dI)  # (dI' * dV_dI)
end




function evaluate_many!(temp, B::Vector{TB}, r::SVector{M, T}
               ) where {TB <: BLNBody{N}} where {N, M, T}
   E = temp # zeros(T, length(B))
   D = B[1].D
   # it is assumed implicitly that all basis functions use the same dictionary!
   # and that each basis function BLNBody contains just a single c and a single t
   I1, I2 = invariants(D, r)
   for ib = 1:length(B)
      E[ib] = B[ib].c[1] * I2[1+B[ib].t[1][end]] * monomial(B[ib].t[1], I1)
   end
   fc = fcut(D, r)
   scale!(E, fc)
   return E
end


function evaluate_many_d!(temp, B::Vector{TB}, r::SVector{M, T}
               ) where {TB <: BLNBody{N}} where {N, M, T}
   E, dM, dE, dEfinal = temp
   fill!(E, 0.0)     # size (nB,)
   fill!(dE, 0.0)    # size (M, nB)
   fill!(dM, 0.0)    # size (M, nB)
   nB = length(B)

   D = B[1].D
   # it is assumed implicitly that all basis functions use the same dictionary!
   # and that each basis function BLNBody contains just a single c and a single t
   I1, I2, dI1, dI2 = invariants_ed(D, r)
   for ib = 1:nB    # (ib, b) in enumerate(B)
      α = B[ib].t[1]
      c = B[ib].c[1]
      m, m_d = monomial_d(α, I1)
      E[ib] += c * I2[1+α[end]] * m        # just the value of the function itself
      for t = 1:M
         dM[t,ib] += (c * I2[1+α[end]]) * m_d[t]   # the I2 * ∇m term without the chain rule
         dE[t,ib] += (c * m) * dI2[1+α[end]][t]  # the ∇I2 * m term
      end
   end
   # chain rule
   fc, fc_d = fcut_d(D, r)
   for ib = 1:nB
      for i = 1:M, t = 1:M  # dE[ib] += dI1' * dM[ib]
         dE[t,ib] += dI1[i][t] * dM[i,ib]
      end
      for t = 1:M
         dE[t,ib] = dE[t,ib] * fc + E[ib] * fc_d[t]
      end
   end
   # write into an output vector
   for i = 1:nB, t = 1:M
      dEfinal[i][t] = dE[t, i]
   end
   return dEfinal
end
