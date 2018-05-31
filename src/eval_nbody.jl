

@generated function simplex_edges(Rs::AbstractVector, J::SVector{K, Int}) where K
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


function _get_loop_ex(N)
   # generate the expression for the multiple for loops, e.g. for 4B:
   # for i_1 = 1:(nR-2), i_2 = (i_1+1):(nR-1), i_3 = (i_2+1):nR
   str_loop = "for i_1 = 1:(nR-$(N-2))"
   for n = 2:N-1
      str_loop *= ", i_$n = (i_$(n-1)+1):(nR-$(N-1-n))"
   end
   str_loop *= "\n end"
   return parse(str_loop)
end


function _get_Jvec_ex(N)
   # inside these N-1 loops we collect the loop indices into an SVector, e.g.
   # J = @SVector [i_1, i_2, i_3]
   str = "J = SVector(i_1"
   for n = 2:N-1
      str *= ", i_$n"
   end
   return parse(str * ")")
end


@generated function eval_site(V::NBodyFunction{N},
                              Rs::AbstractVector{JVec{T}}) where {N, T}
   code = Expr[]
   # initialise the output
   push!(code, :( E = zero(T)      ))
   push!(code, :( nR = length(Rs)  ))
   push!(code, :( rcut = cutoff(V) ))

   # generate the multiple for-loop
   ex_loop = _get_loop_ex(N)

   # inside the loop
   code_inner = Expr[]
   # collect the indices into a vector
   push!(code_inner,      _get_Jvec_ex(N) )
   # collect the edge lengths and edge directions (lexicographical ordering)
   push!(code_inner, :(   (s, S) = simplex_edges(Rs, J) ))
   push!(code_inner, :(   if maximum(s) > rcut; continue; end ))
   # now call `V` with the simplex-lengths and add this to the site energy
   # the normalisation is due to the fact that this term actually appears
   # in N site energies. (once for each corner of the simplex)
   push!(code_inner, :(   E += V(s) / $N ))

   # put code_inner into the loop expression
   ex_loop.args[2] = Expr(:block, code_inner...)

   # now append the loop to the main code
   push!(code, ex_loop)

   quote
      # $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, code...))
      return E
   end
end


"""
convert ∇V (where ∇ is the gradient w.r.t. bond-lengths) into forces, i.e., into
∇Vsite (where ∇ is the gradient w.r.t. positions)

J : neighbour (sub-) indices
"""
@generated function _grad_len2pos!(dVsite, dV, J::SVector{K, Int}, s, S) where {K}
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
      # $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, code...))
      return nothing
   end
end



@generated function eval_site_d!(dVsite::AbstractVector{JVec{T}},
                                 V::NBodyFunction{N},
                                 Rs::AbstractVector{JVec{T}}) where {N, T}
   code = Expr[]
   # initialise the output
   push!(code, :( nR = length(Rs) ))
   push!(code, :( for n = 1:nR; dVsite[n] *= 0.0; end ))

   # generate the multiple for-loop
   ex_loop = _get_loop_ex(N)

   # inside the loop
   code_inner = Expr[]
   # collect the indices into a vector
   push!(code_inner,     _get_Jvec_ex(N))
   # collect the edge lengths and edge directions (lexicographical ordering)
   push!(code_inner, :(  (s, S) = simplex_edges(Rs, J) ))
   push!(code_inner, :(  if maximum(s) > rcut; continue; end ))
   # now call `V` with the simplex-lengths and add this to the site energy
   # the normalisation is due to the fact that this term actually appears
   # in N site energies. (once for each corner of the simplex)
   push!(code_inner, :(  dV = evaluate_d(V, s) / $N ))
   # still inside the loop: convert ∇V (where ∇ is the gradient w.r.t.
   # bond-lengths) into ∇Vsite (where ∇ is the gradient w.r.t. positions)
   push!(code_inner, :(  _grad_len2pos!(dVsite, dV, J, s, S) ))

   # put code_inner into the loop expression
   ex_loop.args[2] = Expr(:block, code_inner...)

   # now append the loop to the main code
   push!(code, ex_loop)

   quote
      # $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, code...))
      return nothing
   end
end
