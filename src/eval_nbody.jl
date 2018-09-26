using JuLIP: neighbourlist,
             @D,
             JVec

import JuLIP.Potentials: evaluate,
                         evaluate_d,
                         evaluate_dd,
                         evaluate_d!

using NeighbourLists: nbodies,
                      maptosites!,
                      maptosites_d!,
                      virial!,
                      max_neigs,
                      sites


function _get_loop_ex(N)
   # generate the expression for the multiple for loops, e.g. for 4B:
   # for i_1 = 1:(nR-2), i_2 = (i_1+1):(nR-1), i_3 = (i_2+1):nR
   str_loop = "for i_1 = 1:(nR-$(N-2))"
   for n = 2:N-1
      str_loop *= ", i_$n = (i_$(n-1)+1):(nR-$(N-1-n))"
   end
   str_loop *= "\n end"
   return Meta.parse(str_loop)
end


function _get_Jvec_ex(N)
   # inside these N-1 loops we collect the loop indices into an SVector, e.g.
   # J = @SVector [i_1, i_2, i_3]
   str = "J = SVector(i_1"
   for n = 2:N-1
      str *= ", i_$n"
   end
   return Meta.parse(str * ")")
end


skipunorderedsimplices(V::NBodyFunction) =
      skipunorderedsimplices(descriptor(V))

skipunorderedsimplices(D::NBSiteDescriptor) = false

skipunorderedsimplices(D::NBClusterDescriptor) = true

function unordered(i, j, J)
   for jj in J
      i > j[jj] && return true
   end
   return false
end

@generated function eval_site_nbody!( ::Val{N},
                                      ii::Int, jj::AbstractVector{TI},
                                      Rs::AbstractVector{JVec{T}},
                                      rcut::T,
                                      skipunordered::Bool,
                                      reducefun,
                                      out,
                                      temp ) where {N, T, TI <: Integer}
   code = Expr[]
   # initialise the output
   push!(code, :( nR = length(Rs)  ))
   push!(code, :( skip = false ))

   # generate the multi-for-loop
   ex_loop = _get_loop_ex(N)

   # inside the loop
   # ---------------
   code_inner = Expr[]
   # collect the indices into a vector
   push!(code_inner,      _get_Jvec_ex(N) )
   push!(code_inner, :(   if skipunordered; skip = unordered(ii, jj, J); end   ))
   # now call `V` with the simplex-corner vectors and "add" this to the site energy
   push!(code_inner, :(   if !skip; out = reducefun(out, Rs, ii, J, temp); end ))

   # put code_inner into the loop expression
   ex_loop.args[2] = Expr(:block, code_inner...)

   # now append the loop to the main code
   push!(code, ex_loop)

   quote
      @inbounds $(Expr(:block, code...))
      return out
   end
end



function site_energies(V::NBodyFunction{N, DT}, at::Atoms{T}
                  ) where {N, T, DT <: NBSiteDescriptor}
   Es = zeros(T, length(at))
   for (i, j, r, R) in sites(at, cutoff(V))
      Es[i] = eval_site_nbody!(Val(N), i, j, R, cutoff(V),
                               skipunorderedsimplices(V),
                               ((out, R, ii, J, temp) -> out + evaluate(V, R, ii, J)),
                               zero(T), nothing)
   end
   return Es
end


# this is probably already in JuLIP??? if not, it should be moved to JuLIP??
energy(V::NBodyFunction, at::Atoms) = sum_kbn(site_energies(V, at))

# this appears to be a nice generic implementation of forces with a
# temporary array => move this to JuLIP!
function forces(V::NBodyFunction{N}, at::Atoms{T}) where {N, T}
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   for (i, j, r, R) in sites(nlist)
      fill!(dVsite, zero(JVec{T}))
      eval_site_nbody!(
            Val(N), i, j, R, cutoff(V), skipunorderedsimplices(V),
            (out, R, ii, J, temp) -> evaluate_d!(out, V, R, ii, J),
            dVsite, nothing )   # dVsite == out, nothing == temp
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
            Val(N), i, j, R, cutoff(V), skipunorderedsimplices(V),
            (out, R, ii, J, temp) -> evaluate_d!(out, V, R, ii, J),
            dVsite, nothing )
      S += JuLIP.Potentials.site_virial(dVsite, R)
   end
   return S
end



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



function energy(B::AbstractVector{TB}, at::Atoms{T}
                ) where {TB <: NBodyFunction{N}, T} where {N}
   # TODO: assert that all B[j] have the same invariants
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   E = zeros(T, length(B))
   for (i, j, r, R) in sites(nlist)
      # evaluate all the site energies at the same time
      # for each simplex, write the nB energies into temp
      # then add them to E, which is just passed through all the
      # various loops, so no need to update it here again
      eval_site_nbody!(Val(N), i, j, R, rcut, skipunorderedsimplices(B[1]),
                       (out, R, ii, J, temp) -> evaluate_many!(out, B, R, ii, J),
                       E, nothing)
   end
   return E
end


function forces(B::AbstractVector{TB}, at::Atoms{T}
              ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nB = length(B)
   # forces
   F =      [ zeros(JVec{T}, length(at)) for n = 1:nB ]
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:nB ]

   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:nB; fill!(dVsite[n], zero(JVec{T})); end
      # fill dVsite
      eval_site_nbody!(Val(N), i, j, R, rcut, skipunorderedsimplices(B[1]),
                       (out, R, ii, J, temp) -> evaluate_many_d!(out, B, R, ii, J),
                       dVsite, nothing)
      # write it into the force vectors
      for ib = 1:nB, n = 1:length(j)
         F[ib][j[n]] -= dVsite[ib][n]
         F[ib][i] += dVsite[ib][n]
      end
   end
   return F
end


function virial(B::AbstractVector{TB}, at::Atoms{T}
                ) where {TB <: NBodyFunction{N}, T} where {N}
   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nB = length(B)
   # virials (main output)
   S = fill((@SMatrix zeros(3,3)), nB)
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:nB ]

   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:nB; dVsite[n] .*= 0.0; end
      # fill dVsite
      eval_site_nbody!(Val(N), i, j, R, rcut, skipunorderedsimplices(B[1]),
                       (out, R, ii, J, temp) -> evaluate_many_d!(out, B, R, ii, J),
                       dVsite, nothing)
      # update the virials
      for iB = 1:nB
         S[iB] += JuLIP.Potentials.site_virial(dVsite[iB], R)
      end
   end
   return S
end
