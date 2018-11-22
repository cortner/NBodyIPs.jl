
using JuLIP: Atoms,
             neighbourlist,
             site_energies,
             sites,
             JVec

using JuLIP.Potentials: site_virial,
                        evaluate,
                        evaluate_d,
                        evaluate_d!

using NeighbourLists: max_neigs

import JuLIP: site_energies,
              energy,
              forces,
              virial,
              cutoff

import NBodyIPs: evaluate_many!,
                 evaluate_many_d!,
                 evaluate_many_ed!,
                 eval_site_nbody!,
                 _grad_len2pos!,
                 TempEdV

function cutoff(V::AbstractEnvIP)
   @assert cutoff(Vn(V)) <= cutoff(Vr(V))
   return cutoff(Vr(V::EnvIP))
end

n_fun(V::EnvIP, n) = (V.t == 0) ? 1.0 : n^V.t

function n_fun_d(V::EnvIP, n)
   if V.t == 0
      return 0.0
   elseif V.t == 1
      return 1.0
   else
      return V.t * n^(V.t-1)
   end
end

site_ns(V::EnvIP, at) = n_fun.(V, site_energies(Vn(V), at))

function site_ns_ed(V::EnvIP, at)
   Vns = site_energies(Vn(V), at)
   return n_fun.(V, Vns), n_fun_d.(V, Vns)
end

function site_n_d!(dVn, V::EnvIP, r, R, Ni, dNi)
   for n = 1:length(r)
      dVn[n] = 0.5 * dNi * evaluate_d(Vn(V), r[n]) * R[n] / r[n]
   end
   return dVn
end

site_energies(V::AbstractEnvIP, at::Atoms) =
      site_ns(V, at) .* site_energies(Vr(V), at)


energy(V::AbstractEnvIP, at::Atoms) = sum_kbn(site_energies(V, at))


function forces(V::AbstractEnvIP{N}, at::Atoms{T}) where {N, T}

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
      # compute the site energy gradients
      fill!(dVsite, zero(JVec{T}))
      eval_site_nbody!(
            Val(N), i, j, R, cutoff(V), false,
            (out, R, i, J, temp) -> evaluate_d!(out, Vr(V), R, i, J),
            dVsite, nothing )   # dVsite == out, nothing == temp
      # compute the neighbour count gradients
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



function virial(V::AbstractEnvIP{N}, at::Atoms{T}) where {N, T}
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
      # compute the site energy gradients
      fill!(dVsite, zero(JVec{T}))
      eval_site_nbody!(
            Val(N), R, cutoff(V),
            (out, R, J, temp) -> evaluate_d!(out, Vr(V), R, J),
            dVsite, nothing )   # dVsite == out, nothing == temp
      # compute the neighbour count gradients
      site_n_d!(dVn, V, r, R, Ns[i], dNs[i])

      # convert the two into a single site energy gradient
      # so that we can call site_virial on it
      for n = 1:length(j)
         dVsite[n] = Ns[i] * dVsite[n] + Vs[i] * dVn[n]
      end

      S += site_virial(dVsite, R)
   end
   return S
end



function energy(B::Vector{TB}, at::Atoms{T}, typewarn=true
                ) where {TB <: EnvIP{N}, T} where {N}
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
      fill!(Etemp, zero(T))
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp) -> evaluate_many!(out, Br, R, J),
                       Etemp, nothing)
      #
      for nb = 1:length(B)
         E[nb] += Etemp[nb] * Ns[nb][i]
      end
   end
   return E
end



function forces(B::AbstractVector{TB}, at::Atoms{T}, typewarn=true
              )where {TB <: EnvIP{N}, T} where {N}
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

   # extras dfor Env
   dVn = zeros(JVec{T}, maxneigs)
   Etemp = zeros(T, length(B))

   # compute the N-components
   Ns = [ site_ns(V, at) for V in B ]
   dNs = [ site_ns_ed(V, at)[2] for V in B ]

   for (i, j, r, R) in sites(nlist)
      # clear dVsite and Etemp
      for n = 1:nB; fill!(dVsite[n], zero(JVec{T})); end
      fill!(Etemp, zero(T))
      # fill site energy and dVsite
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp) -> evaluate_many!(out, Br, R, J),
                       Etemp, nothing)
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp) -> evaluate_many_d!(out, Br, R, J),
                       dVsite, nothing)

      # write it into the force vectors
      for ib = 1:nB, n = 1:length(j)
         site_n_d!(dVn, B[ib], r, R, Ns[ib][i], dNs[ib][i])
         f = (Ns[ib][i] * dVsite[ib][n] + dVn[n] * Etemp[ib])
         F[ib][j[n]] -= f
         F[ib][i] += f
      end
   end
   return F
end


function virial(B::AbstractVector{TB}, at::Atoms{T}, typewarn=true
              )where {TB <: EnvIP{N}, T} where {N}
   if typewarn
      !isleaftype(TB) && warn("TB is not a leaf type")
   end

   Br = [Vr(b) for b in B]

   rcut = cutoff(B[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nedges = (N*(N-1))÷2
   nB = length(B)
   # virials (main output)
   S = fill((@SMatrix zeros(3,3)), length(B))
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs)   for n = 1:length(B) ]

   # extras dfor Env
   dVn = zeros(JVec{T}, maxneigs)
   Etemp = zeros(T, length(B))

   # compute the N-components
   Ns = [ site_ns(V, at) for V in B ]
   dNs = [ site_ns_ed(V, at)[2] for V in B ]

   for (i, j, r, R) in sites(nlist)
      # clear dVsite and Etemp
      for n = 1:nB; fill!(dVsite[n], zero(JVec{T})); end
      fill!(Etemp, zero(T))
      # fill site energy and dVsite
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp) -> evaluate_many!(out, Br, R, J),
                       Etemp, nothing)
      eval_site_nbody!(Val(N), R, rcut,
                       (out, R, J, temp) -> evaluate_many_d!(out, Br, R, J),
                       dVsite, nothing)

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
         S[iB] += site_virial(dVsite[iB], R)
      end
   end
   return S
end



# ============================================================
#     Re-Implement from Scratch for Faster `EnvPolys`
# ============================================================



function site_energies(V::EnvPoly{N}, at::Atoms{T}) where {N, T <: AbstractFloat}
   Es = zeros(T, length(at))

   # all the site energies => should be trivial in terms of cost
   Ns = site_energies(V.Vn, at)

   # storage
   Etemp = zeros(T, length(V.Vr))

   rcut = cutoff(V.Vr[1])
   for (i, j, r, R) in sites(at, rcut)
      fill!(Etemp, zero(T))
      eval_site_nbody!(Val(N), i, j, R, rcut, false,
                       (out, R, ii, J, temp) -> evaluate_many!(out, V.Vr, R, ii, J),
                       Etemp, nothing)
      for P = 1:length(V.Vr)
         Es[i] += Ns[i]^(P-1) * Etemp[P]
      end
   end

   return Es
end


function forces(V::EnvPoly{N}, at::Atoms{T}) where {N, T}
   rcut = cutoff(V.Vr[1])
   nlist = neighbourlist(at, rcut)
   maxneigs = max_neigs(nlist)
   nedges = (N*(N-1))÷2
   nVr = length(V.Vr)

   # forces
   F = zeros(JVec{T}, length(at))
   # site gradient
   dVsite = [ zeros(JVec{T}, maxneigs) for _ = 1:nVr ]

   # extras for Env
   Etemp = zeros(T, nVr)
   out = TempEdV(Etemp, dVsite)

   # compute the N-components
   Ns = site_energies(V.Vn, at)

   for (i, j, r, R) in sites(nlist)
      # clear dVsite and Etemp
      for n = 1:nVr; fill!(dVsite[n], zero(JVec{T})); end
      fill!(Etemp, zero(T))
      # fill site energy and dVsite
      eval_site_nbody!(Val(N), i, j, R, rcut, false,
                       (out, R, ii, J, temp) -> evaluate_many_ed!(out, V.Vr, R, ii, J),
                       out, nothing)

      # write it into the force vectors
      # first for P = 0
      for n = 1:length(j)
         F[j[n]] -= dVsite[1][n]
         F[i] += dVsite[1][n]
      end

      for P = 2:nVr, n = 1:length(j)
         # energy contribution:
         #   Es[i] += Ns[i]^(P-1) * Etemp[P]
         dVn = 0.5 * evaluate_d(V.Vn, r[n]) * R[n] / r[n]
         f = (  Ns[i]^(P-1) * dVsite[P][n]
              + (P-1) * Ns[i]^(P-2) * dVn * Etemp[P] )
         F[j[n]] -= f
         F[i] += f
      end
   end
   return F
end
