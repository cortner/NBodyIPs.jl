struct BondLengthDesc{TT, TC} <: NBSiteDescriptor
   transform::TT
   cutoff::TC
end

BondLengthDesc(transform::String, cutoff::Union{String, Tuple}) =
         BondLengthDesc(SpaceTransform(transform), Cutoff(cutoff))

Dict(D::BondLengthDesc) = Dict( "__id__"    =>  "BondLengthDesc",
                                "transform" =>  Dict(D.transform),
                                "cutoff"    =>  Dict(D.cutoff) )

BondLengthDesc(D::Dict) = BondLengthDesc( SpaceTransform(D["transform"]),
                                          Cutoff(D["cutoff"]) )

==(D1::BondLengthDesc, D2::BondLengthDesc) =
      ( (D1.transform == D2.transform) && (D1.cutoff == D2.cutoff) )

Base.convert(::Val{:BondLengthDesc}, D::Dict) = BondLengthDesc(D)


tdegrees(::BondLengthDesc, vN::Val{N}) where {N} = BLInvariants.tdegrees(vN)

# @inline ricoords(D::BondLengthDesc, Rs, J) = 

@inline invariants(D::BondLengthDesc, r) = BLInvariants.invariants(transform.(D, r))

@inline function invariants_d(D::BondLengthDesc, r)
   DI1, DI2 = BLInvariants.invariants_d(transform.(D, r))
   t_d = transform_d.(D, r)
   return _sdot(t_d, DI1), _sdot(t_d, DI2)
end

@inline function invariants_ed(D::BondLengthDesc, r)
   x = transform.(D, r)
   I1, I2, DI1, DI2 = BLInvariants.invariants_ed(x)
   x_d = transform_d.(D, r)
   return I1, I2, _sdot(x_d, DI1), _sdot(x_d, DI2)
end

"""
`_simplex_edges(Rs::AbstractVector, J::SVector{K, Int})`

Take collection `Rs::AbstractVector{JVec}` representing a simplex
and convert into an ordered vector of bondlengths, e.g.,
```
[ r01, r02, r03, r12, r13, r23 ]
[ r01, r02, r03, r04, r12, r13, r14, r23, r24, r34 ]
```

* `Rs` : collection of all R-vectors within a site energy computation
* `J` : indices of the vectors involved in the current nbody computation
"""
@generated function _simplex_edges(Rs::AbstractVector, J::SVector{K, Int}) where K
   # note K = N-1 ]
   code = Expr[]
   idx = 0
   for n = 1:K
      idx += 1
      push_str!(code, "s_$idx = norm(Rs[J[$n]])")  # r_0n
      push_str!(code, "S_$idx = Rs[J[$n]] / s_$idx")
   end
   for n = 1:K-1, m = (n+1):K
      idx += 1
      push_str!(code, "s_$idx = norm(Rs[J[$n]] - Rs[J[$m]])")   # r_nm
      push_str!(code, "S_$idx = (Rs[J[$n]] - Rs[J[$m]])/s_$idx")   # S_nm ∝ Rn - Rm
   end
   str_s = "s = @SVector [s_1"
   str_S = "S = @SVector [S_1"
   for i = 2:idx
      str_s *= ", s_$i"
      str_S *= ", S_$i"
   end
   push_str!(code, str_s * "]")
   push_str!(code, str_S * "]")
   quote
      $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, code...))
      return s, S
   end
end

"""
convert ∇V (where ∇ is the gradient w.r.t. bond-lengths) into forces, i.e., into
∇Vsite (where ∇ is the gradient w.r.t. positions)

J : neighbour (sub-) indices
"""
@generated function _grad_len2pos!(dVsite, dV, J::SVector{K, Int}, S) where {K}
   # K is the number of neighbours, i.e. N = K+1 counting also the center atom
   # length(dV) == length(s) == length(S) == K * (K+1)/2
   # ------
   code = Expr[]
   idx = 0
   # the first K entries of dV, s, S are the |Ri - 0|
   for k = 1:K
      idx += 1   # idx == k of course
      push!(code, :( dVsite[J[$k]] += dV[$idx] * S[$idx] ))
   end
   # the remaining ones are |R_i - R_j|
   for n = 1:K-1, m = (n+1):K
      idx += 1
      push!(code, :( dVsite[J[$n]] += dV[$idx] * S[$idx] ))
      push!(code, :( dVsite[J[$m]] -= dV[$idx] * S[$idx] ))
   end
   quote
      $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, code...))
      return dVsite
   end
end



function evaluate(V::NBodyFunction{N},
                  desc::BondLengthDesc,
                  Rs::AbstractVector{JVec{T}},
                  J::SVector{K, Int}) where {N, T, K}

   # get bond-lengths (and bond directions => throw those away)
   rcut = cutoff(desc.cutoff)
   r, _ = _simplex_edges(Rs, J)
   if maximum(r) > rcut
      return zero(T)
   end
   # compute the cut-off
   fc = fcut(desc.cutoff, r)
   if fc == 0
      return zero(T)
   end
   # get the invariants
   I1, I2 = invariants(desc, r)
   # evaluate the inner potential function (e.g. polynomial)
   return evaluate_I(V, I1, I2, fc)
end

# (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(V, s)/N, J, S),
function evaluate_d!(dVsite,
                     V::NBodyFunction{N},
                     desc::BondLengthDesc,
                     Rs,
                     J) where {N}
   # get bond-lengths (and bond directions => throw those away)
   rcut = cutoff(desc.cutoff)
   r, S = _simplex_edges(Rs, J)
   if maximum(r) > rcut
      return dVsite
   end
   # compute the cut-off
   fc, fc_d = fcut_d(desc.cutoff, r)
   if fc == 0
      return dVsite
   end
   # get the invariants
   I1, I2, I1_d, I2_d = invariants_ed(desc, r)
   # evaluate the inner potential function (e.g. polynomial)
   dV_dr = evaluate_I_d(V, I1, I2, I1_d, I2_d, fc, fc_d)
   # convert to gradient w.r.t. (relative) position vectors and write into dVsite
   return _grad_len2pos!(dVsite, dV_dr, J, S)
end


function evaluate_many!(Es,
                        B::AbstractVector{TB},
                        desc::BondLengthDesc,
                        Rs, J)  where {TB <: NBodyFunction{N}} where {N}

   rcut = cutoff(desc.cutoff)
   r, _ = _simplex_edges(Rs, J)
   maximum(r) > rcut && return Es
   fc = fcut(desc.cutoff, r)
   fc == 0 && return Es

   # get the invariants
   I1, I2 = invariants(desc, r)
   # evaluate the inner potential function (e.g. polynomial)
   for n = 1:length(B)
      Es[n] += evaluate_I(B[n], I1, I2, fc)
   end
   return Es
end


function evaluate_many_d!(dVsite::AbstractVector,
                          B::AbstractVector{TB},
                          desc::BondLengthDesc,
                          Rs,
                          J)  where {TB <: NBodyFunction{N}} where {N}
   rcut = cutoff(desc.cutoff)
   r, S = _simplex_edges(Rs, J)
   maximum(r) > rcut && return dVsite
   fc, fc_d = fcut_d(desc.cutoff, r)
   fc == 0 && return dVsite

   # get the invariants
   I1, I2, I1_d, I2_d = invariants_ed(desc, r)
   # evaluate the inner potential function (e.g. polynomial)
   for n = 1:length(B)
      dV_dr = evaluate_I_d(B[n], I1, I2, I1_d, I2_d, fc, fc_d)
      _grad_len2pos!(dVsite[n], dV_dr, J, S)
   end
   return dVsite
end
