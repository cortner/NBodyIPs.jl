struct BondLengthDesc{TT, TC} <: NBSiteDescriptor
   transform::TT
   cutoff::TC
end

BondLengthDesc(transform::String, cutoff::Union{String, Tuple}) =
         BondLengthDesc(SpaceTransform(transform), Cutoff(cutoff))

# -------------- IO -------------------
Dict(D::BondLengthDesc) = Dict( "__id__"    =>  "BondLengthDesc",
                                "transform" =>  Dict(D.transform),
                                "cutoff"    =>  Dict(D.cutoff) )

BondLengthDesc(D::Dict) = BondLengthDesc( SpaceTransform(D["transform"]),
                                          Cutoff(D["cutoff"]) )

==(D1::BondLengthDesc, D2::BondLengthDesc) =
      ( (D1.transform == D2.transform) && (D1.cutoff == D2.cutoff) )

Base.convert(::Val{:BondLengthDesc}, D::Dict) = BondLengthDesc(D)

# ------------- Interface Code ---------------

tdegrees(::BondLengthDesc, vN::Val{N}) where {N} = BLInvariants.tdegrees(vN)

@inline ricoords(D::BondLengthDesc, Rs, J) = edge_lengths(Rs, J)

@inline gradri2gradR!(desc::BondLengthDesc, dVsite, dV_dr, Rs, J, r) =
   _grad_len2pos!(dVsite, dV_dr, Rs, J, r)

@inline fcut(D::BondLengthDesc, r) = fcut(D.cutoff, r)
@inline fcut_d(D::BondLengthDesc, r) = fcut_d(D.cutoff, r)

@inline skip_simplex(D::BondLengthDesc, r) = (maximum(r) > cutoff(D.cutoff))

@inline invariants(D::BondLengthDesc, r) = BLInvariants.invariants(transform.(D, r))

@inline function invariants_ed(D::BondLengthDesc, r)
   x = transform.(D, r)
   I1, I2, DI1, DI2 = BLInvariants.invariants_ed(x)
   x_d = transform_d.(D, r)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end

# ------------- main evaluation code -----------

function evaluate(V::NBodyFunction{N},
                  desc::NBSiteDescriptor,
                  Rs::AbstractVector{JVec{T}},
                  J::SVector{K, Int}) where {N, T, K}

   # get the physical descriptor: bond-lengths (+ bond-angles)
   rθ = ricoords(desc, Rs, J)
   # check whether to skip this N-body term?
   skip_simplex(desc, rθ) && zero(T)
   # compute the cutoff (and skip this site if the cutoff is zero)
   fc = fcut(desc, rθ)
   fc == 0 && return zero(T)
   # compute the invariants (this also applies the transform)
   II = invariants(desc, rθ)
   # evaluate the inner potential function (e.g. polynomial)
   return evaluate_I(V, II) * fc
end

# (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(V, s)/N, J, S),
function evaluate_d!(dVsite,
                     V::NBodyFunction{N},
                     desc::BondLengthDesc,
                     Rs,
                     J) where {N}
   # get the physical descriptor: bond-lengths (+ bond-angles)
   rθ = ricoords(desc, Rs, J)
   # check whether to skip this N-body term?
   skip_simplex(desc, rθ) && dVsite
   # compute the cutoff (and skip this site if the cutoff is zero)
   fc, fc_d = fcut_d(desc, rθ)
   fc == 0 && return dVsite
   # get the invariants
   II = invariants_ed(desc, rθ)
   # evaluate the inner potential function (e.g. polynomial)
   V, dV_drθ = evaluate_I_ed(V, II)
   # convert to gradient w.r.t. rθ
   dV_drθ = dV_drθ * fc + V * fc_d
   # convert to gradient w.r.t. dR
   return gradri2gradR!(desc, dVsite, dV_drθ, Rs, J, rθ)
end


function evaluate_many!(Es,
                        B::AbstractVector{TB},
                        desc::BondLengthDesc,
                        Rs, J)  where {TB <: NBodyFunction{N}} where {N}

   rθ = ricoords(desc, Rs, J)
   skip_simplex(desc, rθ) && return Es
   fc = fcut(desc, rθ)
   fc == 0 && return Es
   II = invariants(desc, rθ)

   # evaluate the inner potential function (e.g. polynomial)
   for n = 1:length(B)
      Es[n] += evaluate_I(B[n], II) * fc
   end
   return Es
end


function evaluate_many_d!(dVsite::AbstractVector,
                          B::AbstractVector{TB},
                          desc::BondLengthDesc,
                          Rs,
                          J)  where {TB <: NBodyFunction{N}} where {N}
   rθ = ricoords(desc, Rs, J)
   skip_simplex(desc, rθ) && return dVsite
   fc, fc_d = fcut_d(desc, rθ)
   fc == 0 && return dVsite
   II = invariants_ed(desc, rθ)

   # evaluate the inner potential function (e.g. polynomial)
   for n = 1:length(B)
      V, dV_drθ = evaluate_I_ed(B[n], II)
      dV_drθ = dV_drθ * fc + V * fc_d
      gradri2gradR!(desc, dVsite[n], dV_drθ, Rs, J, rθ)
   end
   return dVsite
end





# ---------------- inner kernels for the edge length descriptor --------------

"""
`edge_lengths(Rs::AbstractVector, J::SVector{K, Int})`

Take collection `Rs::AbstractVector{JVec}` representing a simplex
and convert into an ordered vector of bondlengths, e.g.,
```
[ r01, r02, r03, r12, r13, r23 ]
[ r01, r02, r03, r04, r12, r13, r14, r23, r24, r34 ]
```

* `Rs` : collection of all R-vectors within a site energy computation
* `J` : indices of the vectors involved in the current nbody computation
"""
@generated function edge_lengths(Rs::AbstractVector, J::SVector{K, Int}) where K
   # note K = N-1 ]
   code = Expr[]
   idx = 0
   for n = 1:K
      idx += 1
      push_str!(code, "s_$idx = norm(Rs[J[$n]])")  # r_0n
   end
   for n = 1:K-1, m = (n+1):K
      idx += 1
      push_str!(code, "s_$idx = norm(Rs[J[$n]] - Rs[J[$m]])")   # r_nm
   end
   s_str = "s = @SVector [s_1"
   if idx > 1; s_str *= prod(", s_$i" for i = 2:idx); end
   push_str!(code, s_str * "]")
   quote
      $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, code...))
      return s
   end
end



"""
convert ∇V (where ∇ is the gradient w.r.t. bond-lengths) into forces, i.e., into
∇Vsite (where ∇ is the gradient w.r.t. positions)

J : neighbour (sub-) indices
"""
@generated function _grad_len2pos!(dVsite, dV, Rs, J::SVector{K, Int}, r) where {K}
   # K is the number of neighbours, i.e. N = K+1 counting also the center atom
   # length(dV) == length(s) == length(S) == K * (K+1)/2
   # ------
   code = Expr[]
   idx = 0
   # the first K entries of dV, s, S are the |Ri - 0|
   for k = 1:K
      idx += 1   # idx == k of course
      push!(code, :( dVsite[J[$k]] += dV[$idx] * (Rs[J[$k]]/r[$idx]) ))
   end
   # the remaining ones are |R_i - R_j|
   for n = 1:K-1, m = (n+1):K
      idx += 1
      push!(code, :( S = (Rs[J[$n]]-Rs[J[$m]]) / r[$idx] ))
      push!(code, :( dVsite[J[$n]] += dV[$idx] * S ))
      push!(code, :( dVsite[J[$m]] -= dV[$idx] * S ))
   end
   quote
      $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, code...))
      return dVsite
   end
end
