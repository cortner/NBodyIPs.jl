
# ============= General Utility functions ==================

@inline transform(D::AbstractDescriptor, r::Number) = D.transform.f(r)
@inline transform_d(D::AbstractDescriptor, r::Number) = D.transform.f_d(r)
@inline cutoff(D::AbstractDescriptor) = D.cutoff.rcut


evaluate(V::NBodyFunction, Rs, J) =
      evaluate(V, descriptor(V), Rs, J)

evaluate_d!(dVsite, V::NBodyFunction, Rs, J) =
      evaluate_d!(dVsite, V, descriptor(V), Rs, J)


evaluate_many!(out, B, Rs, J) =
      evaluate_many!(out, B, descriptor(B[1]), Rs, J)

evaluate_many_d!(out, B, Rs, J) =
      evaluate_many_d!(out, B, descriptor(B[1]), Rs, J)


"""
`_sdot(a::T, B::SVector{N, T})`: efficiently compute `[ a .* b for b in B ]`
"""
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

# bond length descriptor
include("bldescriptor.jl")

# # bond angle descriptor
# include("badescriptor.jl")
