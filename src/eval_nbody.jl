

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




@generated function eval_site_nbody!( ::Val{N},
                                      Rs::AbstractVector{JVec{T}},
                                      rcut::T,
                                      accumfun,
                                      out,
                                      temp ) where {N, T}
   code = Expr[]
   # initialise the output
   push!(code, :( nR = length(Rs)  ))

   # generate the multiple for-loop
   ex_loop = _get_loop_ex(N)

   # inside the loop
   code_inner = Expr[]
   # collect the indices into a vector
   push!(code_inner,      _get_Jvec_ex(N) )
   # collect the edge lengths and edge directions (lexicographical ordering)
   push!(code_inner, :(   (s, S) = _simplex_edges(Rs, J) ))
   push!(code_inner, :(   if maximum(s) > rcut; continue; end ))
   # now call `V` with the simplex-lengths and add this to the site energy
   # the normalisation is due to the fact that this term actually appears
   # in N site energies. (once for each corner of the simplex)
   push!(code_inner, :(   out = accumfun(out, s, S, J, temp) ))

   # put code_inner into the loop expression
   ex_loop.args[2] = Expr(:block, code_inner...)

   # now append the loop to the main code
   push!(code, ex_loop)

   quote
      @inbounds $(Expr(:block, code...))
      return out
   end
end




function site_energies(V::NBodyFunction{N}, at::Atoms{T}) where {N, T}
   Es = zeros(T, length(at))
   for (i, j, r, R) in sites(at, cutoff(V))
      Es[i] = eval_site_nbody!(Val(N), R, cutoff(V),
                               ((out, s, _,_1,_2) -> out + evaluate(V, s)),
                               zero(T), nothing)
   end
   return Es/N
end



# this is probably already in JuLIP??? if not, it should be moved to JuLIP??
energy(V::NBodyFunction, at::Atoms) =
      sum_kbn(site_energies(V, at))

# this appears to be a nice generic implementation of forces with a
# temporary array => move this to JuLIP!
function forces(V::NBodyFunction{N}, at::Atoms{T}) where {N, T}
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   for (i, j, r, R) in sites(nlist)
      dVsite .*= 0.0
      eval_site_nbody!(
            Val(N), R, cutoff(V),
            (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(V, s)/N, J, S),
            dVsite, nothing )
      # write site energy gradient into forces
      for n = 1:length(j)
         F[j[n]] -= dVsite[n]
         F[i] += dVsite[n]
      end
   end
   return F
end



function virial(V::NBodyFunction{N}, at::Atoms{T}) where {N, T}
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   S = @SMatrix zeros(3,3)
   dVsite = zeros(JVec{T}, maxneigs)
   for (i, j, r, R) in sites(nlist)
      dVsite .*= 0.0
      # eval_site_d!(dVsite, V, R)
      eval_site_nbody!(
            Val(N), R, cutoff(V),
            (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(V, s)/N, J, S),
            dVsite, nothing )
      S += JuLIP.Potentials.site_virial(dVsite, R)
   end
   return S
end


# ------ special treatment of 1-body functions

site_energies(V::NBodyFunction{1}, at::Atoms) =
      fill(V(), length(at))

forces(V::NBodyFunction{1}, at::Atoms{T}) where {T} =
      zeros(SVector{3, T}, length(at))

virial(V::NBodyFunction{1}, at::Atoms{T}) where {T} =
      zero(SMatrix{3, 3, T})



# ========================= assembly support for LSQ system ====================

# For assembling the LSQ system efficiently we need a way to evaluate all basis
# functions of the same body-order at the same time. Otherwise we would be
# re-computing the invariants many many times, which is very expensive.
# To achieve this we just wrap all basis functions of a body-order into
# a new type `NBBasis` which evaluates to a long vector
#
# at the moment, it seems we need to hard-code this to the Polys
# sub-module, but it would be good if this can be fixed, so we keep this
# "interface" here.

function evaluate_many! end
function evaluate_many_d! end


energy(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ energy(b, at) for b in B ]

function energy(B::AbstractVector{TB}, at::Atoms{T}
                ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   temp = zeros(T, length(B))
   E = zeros(T, length(B))
   for (i, j, r, R) in sites(nlist)
      # evaluate all the site energies at the same time
      # for each simples, write the nB energies into temp
      # then add them to E, which is just passed through all the
      # various loops, so no need to update it here again
      eval_site_nbody!(Val(N), R, rcut,
                       (out, s, S, _1, temp) -> (out .+= evaluate_many!(temp, B, s)),
                       E, temp)
   end
   # rescale to account for permutations
   E ./= N
   return E
end



function _acc_manyfrcs(B, dVsite, s, S, J, temp)
   dV = evaluate_many_d!(temp, B, s)
   for ib = 1:length(dVsite)
      _grad_len2pos!(dVsite[ib], dV[ib], J, S)
   end
   return dVsite
end

function forces(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nedges = (N*(N-1))÷2
   nB = length(B)
   # forces
   F =      [ zeros(JVec{T}, length(at)) for n = 1:nB ]
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:nB ]
   # n-body gradients
   dV =     [ zeros(T, nedges)      for n = 1:nB ]
   # temporary arrays to compute the site gradients
   temp = ( zeros(T, nB),
            zeros(T, nedges, nB),
            zeros(T, nedges, nB),
            dV )
   accum_fun = let B=B
      (out, s, S, J, temp) -> _acc_manyfrcs(B, out, s, S, J, temp)
   end

   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:nB
         dVsite[n] .*= 0.0
      end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut, accum_fun, dVsite, temp)
      # write it into the force vectors
      for ib = 1:nB, n = 1:length(j)
         F[ib][j[n]] -= dVsite[ib][n]/N
         F[ib][i] += dVsite[ib][n]/N
      end
   end
   return F
end

forces(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ forces(b, at) for b in B ]


function virial(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nedges = (N*(N-1))÷2
   # virials (main output)
   S = fill((@SMatrix zeros(3,3)), length(B))
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:length(B) ]
   # n-body gradients
   dV =     [ zeros(T, nedges)      for n = 1:length(B) ]
   # temporary arrays to compute the site gradients
   temp = ( zeros(T, length(B)),
            zeros(T, nedges, length(B)),
            zeros(T, nedges, length(B)),
            dV )
   accum_fun = let B=B
      (out, s, S, J, temp) -> _acc_manyfrcs(B, out, s, S, J, temp)
   end

   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:length(dVsite)
         dVsite[n] .*= 0.0
      end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut, accum_fun, dVsite, temp)
      # update the virials
      for iB = 1:length(B)
         S[iB] += JuLIP.Potentials.site_virial(dVsite[iB], R) / N
      end
   end
   return S
end


virial(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ virial(b, at) for b in B ]
