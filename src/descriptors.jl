
# ============= General Utility functions ==================

ninvariants(D::NBodyDescriptor, N::Integer) = ninvariants(D, Val(N))
ninvariants(D::NBodyDescriptor, vN::Val{N}) where {N} = length.(tdegrees(D, vN))

@inline transform(D::NBodyDescriptor, r::Number) = D.transform.f(r)
@inline transform_d(D::NBodyDescriptor, r::Number) = D.transform.f_d(r)
@inline cutoff(D::NBodyDescriptor) = D.cutoff.rcut

combiscriptor(D::NBodyDescriptor) =
      (combiscriptor(D.transform), combiscriptor(D.cutoff))

evaluate(V::NBodyFunction, Rs, i, J) =
      evaluate(V, descriptor(V), Rs, i, J)

evaluate_d!(dVsite, V::NBodyFunction, Rs, i, J) =
      evaluate_d!(dVsite, V, descriptor(V), Rs, i, J)


evaluate_many!(out, B, Rs, i, J) =
      evaluate_many!(out, B, descriptor(B[1]), Rs, i, J)

evaluate_many_d!(out, B, Rs, i, J) =
      evaluate_many_d!(out, B, descriptor(B[1]), Rs, i, J)


# TODO: switch to tuples??
struct TempEdV{T1, T2}
   Etemp::T1
   dVsite::T2
end

evaluate_many_ed!(out::TempEdV, B, Rs, i, J) =
      evaluate_many_ed!(out, B, descriptor(B[1]), Rs, i, J)


"""
`_sdot(a::T, B::SVector{N, T})`: efficiently compute `[ a .* b for b in B ]`
"""
@generated function _sdot(a::T, B::SVector{N, T}) where {N, T}
   code = "@SVector $T["
   for n = 1:N
      code *= "a .* B[$n],"
   end
   code = code[1:end-1] * "]"
   ex = Meta.parse(code)
   quote
      $ex
   end
end


# ------------- main evaluation code -----------


function evaluate(V::NBodyFunction{N},
                  desc::NBSiteDescriptor,
                  Rs::AbstractVector{JVec{T}},
                  i::Int,
                  J::SVector{K, Int}) where {N, T, K}
   # get the physical descriptor: bond-lengths (+ bond-angles)
   rθ = ricoords(desc, Rs, J)
   # check whether to skip this N-body term?
   skip_simplex(desc, rθ) && return zero(T)
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
                     desc::NBSiteDescriptor,
                     Rs,
                     i::Int,
                     J) where {N}
   # get the physical descriptor: bond-lengths (+ bond-angles)
   rθ = ricoords(desc, Rs, J)
   # check whether to skip this N-body term?
   skip_simplex(desc, rθ) && return dVsite
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
                        desc::NBSiteDescriptor,
                        Rs, i, J)  where {TB <: NBodyFunction{N}} where {N}
   rθ = ricoords(desc, Rs, J)
   skip_simplex(desc, rθ) && return Es
   return _evaluate_many_ricoords!(Es, B, desc, rθ)
end

evaluate_many_ricoords!(Es,
               B::AbstractVector{TB},
               rθ) where {TB <: NBodyFunction{N}} where {N} =
   _evaluate_many_ricoords!(Es, B, descriptor(B[1]), rθ)

function _evaluate_many_ricoords!(Es,
               B::AbstractVector{TB}, desc,
               rθ) where {TB <: NBodyFunction{N}} where {N}
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
                          desc::NBSiteDescriptor,
                          Rs,
                          i, J)  where {TB <: NBodyFunction{N}} where {N}
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


function evaluate_many_ed!(out::TempEdV,
                           B::AbstractVector{TB},
                           desc::NBSiteDescriptor,
                           Rs,
                           i, J)  where {TB <: NBodyFunction{N}} where {N}
   rθ = ricoords(desc, Rs, J)
   skip_simplex(desc, rθ) && return out
   fc, fc_d = fcut_d(desc, rθ)
   fc == 0 && return out
   II = invariants_ed(desc, rθ)

   # evaluate the inner potential function (e.g. polynomial)
   for n = 1:length(B)
      V, dV_drθ = evaluate_I_ed(B[n], II)
      out.Etemp[n] += V * fc
      dV_drθ = dV_drθ * fc + V * fc_d
      gradri2gradR!(desc, out.dVsite[n], dV_drθ, Rs, J, rθ)
   end
   return out
end




# bond length descriptor
include("bldescriptor.jl")

# bond angle descriptor
include("badescriptor.jl")
