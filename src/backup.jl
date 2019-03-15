# ------- distribution functions on the invariants --------------

function rdf(at::Atoms, rcut, transform=identity)
   if transform != identity
      return idf(2, at, rcut, transform)[1][1]
   end
   nlist = neighbourlist(at, rcut)
   return nlist.r
end

"""
`idf(N::Integer, at, rcut, transform)`

invariants distribution function => accumulates all the invariants values
arising during an N-body assembly over at.
"""
idf(N::Integer, at, rcut, transform) =
      _idf(Val(N), Val(bo2edges(N)), at, rcut, transform)

function _idf(valN::Val{N}, valM::Val{M}, at::Atoms{T}, rcut::T,
              transform) where {N, M, T}
   # compute invariants vectors to learn how many there are
   x = rand(SVector{M,T})
   I1, I2 = invariants(x)

   I1acc = [ T[] for n = 1:length(I1) ]
   I2acc = [ T[] for n = 1:length(I2) ]
   Iacc = (I1acc, I2acc)

   for (i, j, r, R) in sites(at, rcut)
      eval_site_nbody!(valN, R, rcut,
                  (Iacc, s, _1, _2, _3) -> idf_accumulator(Iacc, transform.(s)),
                  Iacc, nothing)
   end

   return I1acc, I2acc
end

function idf_accumulator(Iacc, x)
   I1, I2 = invariants(x)
   for n = 1:length(I1)
      push!(Iacc[1][n], I1[n])
   end
   for n = 1:length(I2)
      push!(Iacc[2][n], I2[n])
   end
   return Iacc
end




# ---------------------- Regularisation of EnvIPs --------------------------

function regularise(N::Integer, P::Integer, B::Vector, r0, r1; kwargs...)
   I = find(isa.(B, EnvBL{N, P}))
   Br = [b.Vr for b in B[I]]
   Ψr = Polys.regularise(N, Br, r0, r1; kwargs...)
   Ψ = zeros(size(Ψr, 1), length(B))
   Ψ[:, I] = Ψr
   return Ψ
end

function regularise(N::Integer, B::Vector, r0, r1; kwargs...)
   # find all EnvBL basis functions with body-order N
   I = find(isa.(B, EnvBL{N}))
   max_P = maximum( b.t for b in B[I] )
   # each polynomial degree p => separate V_Np, hence
   # compute multiple regularisations
   Ψ = zeros(0, length(B))
   for p = 0:max_P
      Ψp = regularise(N, p, B, r0, r1; kwargs...)
      Ψ = [Ψ; Ψp]
   end
   return Ψ
end


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
#    for Vn in V.components
#       if bo2edges(bodyorder(Vn)) == N
#          v += Vn(r)
#       end
#    end
#    return v
# end

# dim(V::NBPoly{N,M}) where {N, M} = M-1
