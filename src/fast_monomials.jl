

using StaticArrays

# univariate monomial
@inline _m1(α::Number, x::T) where {T <: Number} =
      (@fastmath x^α)

# univariate monomial derivative
@inline _m1d(α::Number, x::T) where {T <: Number} =
      α == 0 ? zero(T) : (@fastmath α * x^(α-1))


function monomial(α, x::SVector{K, T}) where {K, T}
   m = one(T)
   for i = 1:K
      @inbounds m *= _m1(α[i], x[i])
   end
   return m
end


@generated function monomial_d(α, x::SVector{K, T}) where {K, T}
   # @assert length(α) >= K
   # evaluate the scalar monomials
   #  f = SVector{...}( x[1]^α[1], ...)
   # df = SVector{...}( α[1] * x[1]^(α[1]-1), ...)
   ex_f = "f = SVector{$K, $T}("
   ex_df = "df = SVector{$K, $T}("
   for i = 1:K
      ex_f  *= " _m1(α[$i], x[$i]), "
      ex_df *= "_m1d(α[$i], x[$i]), "
   end
   ex_f  =  ex_f[1:end-2] * ")"
   ex_df = ex_df[1:end-2] * ")"

   # evaluate the multi-variate monomial
   ex_m = "m = "
   for i = 1:K
      ex_m *= "f[$i] * "
   end
   ex_m = ex_m[1:end-3]

   # evaluate the derivative
   ex_dm = "m_d = SVector{$K, $T}("
   for i = 1:K
      dmj = ""
      for j = 1:K
         if i == j
            dmj *= " df[$j] *"
         else
            dmj *= "  f[$j] *"
         end
      end
      ex_dm *= dmj[1:end-2] * ", "
   end
   ex_dm = ex_dm[1:end-2] * ")"

   quote
      $(parse(ex_f))
      $(parse(ex_df))
      $(parse(ex_m))
      $(parse(ex_dm))
      return m, m_d
   end
end



# ---------------------------------------------------------------
# NOT TECHNICALLY MONOMIALS, BUT VERY SIMILAR OBJECT:
#   TODO: combine fcut_d and monomial_d into a single function 
# ---------------------------------------------------------------

function fcut(D::Dictionary, r::SVector{M, T}) where {M, T}
   fc = one(T)
   for i = 1:M
      @fastmath fc *= fcut(D, r[i])
   end
   return fc
end


@generated function fcut_d(D::Dictionary, r::SVector{M, T}) where {M, T}

   # evaluate the scalar monomials
   #  f = SVector{...}( x[1]^α[1], ...)
   # df = SVector{...}( α[1] * x[1]^(α[1]-1), ...)
   ex_f = "f = SVector{$M, $T}("
   ex_df = "df = SVector{$M, $T}("
   for i = 1:M
      ex_f  *= "  fcut(D, r[$i]), "
      ex_df *= "fcut_d(D, r[$i]), "
   end
   ex_f  =  ex_f[1:end-2] * ")"
   ex_df = ex_df[1:end-2] * ")"

   # evaluate the multi-variate monomial
   ex_m = "fc = "
   for i = 1:M
      ex_m *= "f[$i] * "
   end
   ex_m = ex_m[1:end-3]

   # evaluate the derivative
   ex_dm = "fc_d = SVector{$M, $T}("
   for i = 1:M
      dmj = ""
      for j = 1:M
         if i == j
            dmj *= " df[$j] *"
         else
            dmj *= "  f[$j] *"
         end
      end
      ex_dm *= dmj[1:end-2] * ", "
   end
   ex_dm = ex_dm[1:end-2] * ")"

   quote
      $(parse(ex_f))
      $(parse(ex_df))
      $(parse(ex_m))
      $(parse(ex_dm))
      return fc, fc_d
   end
end
