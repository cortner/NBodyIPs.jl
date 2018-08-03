



function cutoff(V::EnvBLFunction)
   @assert cutoff(Vn(V)) <= cutoff(Vr(V))
   return cutoff(Vr(V::EnvBL))
end

n_fun(V::EnvBL, n) = (V.t == 0) ? 1.0 : n^V.t

function n_fun_d(V::EnvBL, n)
   if V.t == 0
      return 0.0
   elseif V.t == 1
      return 1.0
   else
      return V.t * n^(V.t-1)
   end
end

site_ns(V::EnvBL, at) = n_fun.(V, site_energies(Vn(V), at))

function site_ns_ed(V::EnvBL, at)
   Vns = site_energies(Vn(V), at)
   return n_fun.(V, Vns), n_fun_d.(V, Vns)
end

function site_n_d!(dVn, V::EnvBL, r, R, Ni, dNi)
   for n = 1:length(r)
      dVn[n] = 0.5 * dNi * evaluate_d(Vn(V), r[n]) * R[n] / r[n]
   end
   return dVn
end

function energy(V::EnvBLFunction, at::Atoms)
   # compute the n values
   Ns = site_ns(V, at)
   # compute the inner v values
   Vs = site_energies(Vr(V), at)
   # the energy is just the dot-product now
   return dot(Ns, Vs)
end


function forces(V::EnvBLFunction{N}, at::Atoms{T}) where {N, T}

   # compute the n values
   Ns, dNs = site_ns_ed(V, at)
   # compute the inner v values
   Vs = site_energies(Vr(V), at)

   # now assemble site forces and use those to create
   # total forces by mixing N and V
   cutoff_n = cutoff(Vn(V))
   nlist = neighbourlist(at, cutoff(V))  # this checks that cutoff(Vn) <= cutoff(Vr)
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   dVn = zeros(JVec{T}, maxneigs)

   for (i, j, r, R) in sites(nlist)
      dVsite .*= 0.0
      eval_site_nbody!(
            Val(N), R, cutoff(Vr(V)),
            (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(Vr(V), s)/N, J, S),
            dVsite, nothing )

      site_n_d!(dVn, V, r, R, Ns[i], dNs[i])

      # write site energy gradient into forces
      for n = 1:length(j)
         f = Ns[i] * dVsite[n] + Vs[i] * dVn[n]
         F[j[n]] -= f
         F[i] += f
      end
   end
   return F
end



function virial(V::EnvBLFunction{N}, at::Atoms{T}) where {N, T}
   # compute the n values
   Ns, dNs = site_ns_ed(V, at)
   # compute the inner v values
   Vs = site_energies(Vr(V), at)

   # now assemble site forces and use those to create
   # total forces by mixing N and V
   cutoff_n = cutoff(Vn(V))
   nlist = neighbourlist(at, cutoff(V))  # this checks that cutoff(Vn) <= cutoff(Vr)
   maxneigs = max_neigs(nlist)
   dVsite = zeros(JVec{T}, maxneigs)
   dVn = zeros(JVec{T}, maxneigs)

   S = @SMatrix zeros(3,3)

   for (i, j, r, R) in sites(nlist)
      dVsite .*= 0.0
      dVn .*= 0.0
      eval_site_nbody!(
            Val(N), R, cutoff(V),
            (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(Vr(V), s)/N, J, S),
            dVsite, nothing )
      site_n_d!(dVn, V, r, R, Ns[i], dNs[i])
      for n = 1:length(j)
         dVsite[n] = Ns[i] * dVsite[n] + Vs[i] * dVn[n]
      end

      S += JuLIP.Potentials.site_virial(dVsite, R)
   end
   return S
end



function energy(B::Vector{TB}, at::Atoms{T}, typewarn=true
                ) where {TB <: EnvBL{N}, T} where {N}
   if typewarn
      !isleaftype(TB) && warn("TB is not a leaf type")
   end

   Br = [Vr(b) for b in B]

   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   temp = zeros(T, length(B))
   E = zeros(T, length(B))
   Etemp = zeros(T, length(B))
   # all the site energies => should be trivial in terms of cost
   Ns = [ site_ns(V, at) for V in B ]

   for (i, j, r, R) in sites(nlist)
      # evaluate all the site energies at the same time
      # for each simples, write the nB energies into temp
      # then add them to E, which is just passed through all the
      # various loops, so no need to update it here again
      Etemp .*= 0.0
      eval_site_nbody!(Val(N), R, rcut,
                       (out, s, S, _1, temp) -> (out .+= evaluate_many!(temp, Br, s)),
                       Etemp, temp)
      #
      for nb = 1:length(B)
         E[nb] += Etemp[nb] * Ns[nb][i]
      end
   end
   # rescale to account for permutations
   E ./= N
   return E
end



function forces(B::AbstractVector{TB}, at::Atoms{T}, typewarn=true
              )where {TB <: EnvBL{N}, T} where {N}
   if typewarn
      !isleaftype(TB) && warn("TB is not a leaf type")
   end

   Br = [Vr(b) for b in B]

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
      (out, s, S, J, temp) -> _acc_manyfrcs(Br, out, s, S, J, temp)
   end

   # extras dfor Env
   dVn = zeros(JVec{T}, maxneigs)
   Etemp = zeros(T, length(B))
   temp2 = zeros(T, length(B))

   # compute the N-components
   Ns = [ site_ns(V, at) for V in B ]
   dNs = [ site_ns_ed(V, at)[2] for V in B ]

   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:nB
         dVsite[n] .*= 0.0
      end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut, accum_fun, dVsite, temp)
      # compute the site energies as well
      Etemp .*= 0.0
      eval_site_nbody!(Val(N), R, rcut,
                       (out, s, S, _1, temp) -> (out .+= evaluate_many!(temp, Br, s)),
                       Etemp, temp2)

      # write it into the force vectors
      for ib = 1:nB, n = 1:length(j)
         site_n_d!(dVn, B[ib], r, R, Ns[ib][i], dNs[ib][i])
         # @show Ns[ib][i], dVsite[ib][n], dVn[n], Etemp[ib]
         f = (Ns[ib][i] * dVsite[ib][n] + dVn[n] * Etemp[ib]) / N
         F[ib][j[n]] -= f
         F[ib][i] += f
      end
   end
   return F
end


function virial(B::AbstractVector{TB}, at::Atoms{T}, typewarn=true
              )where {TB <: EnvBL{N}, T} where {N}
   if typewarn
      !isleaftype(TB) && warn("TB is not a leaf type")
   end

   Br = [Vr(b) for b in B]

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
      (out, s, S, J, temp) -> _acc_manyfrcs(Br, out, s, S, J, temp)
   end

   # extras dfor Env
   dVn = zeros(JVec{T}, maxneigs)
   Etemp = zeros(T, length(B))
   temp2 = zeros(T, length(B))

   # compute the N-components
   Ns = [ site_ns(V, at) for V in B ]
   dNs = [ site_ns_ed(V, at)[2] for V in B ]


   for (i, j, r, R) in sites(nlist)
      # clear dVsite
      for n = 1:length(dVsite)
         dVsite[n] .*= 0.0
      end
      # fill dVsite
      eval_site_nbody!(Val(N), R, rcut, accum_fun, dVsite, temp)

      # compute the site energies as well
      Etemp .*= 0.0
      eval_site_nbody!(Val(N), R, rcut,
                       (out, s, S, _1, temp) -> (out .+= evaluate_many!(temp, Br, s)),
                       Etemp, temp2)

      # use Etemp (site energies), the Ns and dVn to convert dVsite in
      # proper site energies
      for ib = 1:length(B)
         site_n_d!(dVn, B[ib], r, R, Ns[ib][i], dNs[ib][i])
         for n = 1:length(j)
            dVsite[ib][n] = Ns[ib][i] * dVsite[ib][n] + Etemp[ib] * dVn[n]
         end
      end

      # update the virials
      for iB = 1:length(B)
         S[iB] += JuLIP.Potentials.site_virial(dVsite[iB], R) / N
      end
   end
   return S
end
