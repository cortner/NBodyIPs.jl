
module EnvBLs

using JuLIP, StaticArrays
using JuLIP.Potentials: evaluate, evaluate_d

using NBodyIPs: max_neigs, eval_site_nbody!, _grad_len2pos!
using NBodyIPs.Polys: Dictionary, NBody

abstract type EnvBLFunction{N} <: AbstractCalculator end

import JuLIP: cutoff, energy, forces, virial

import NBodyIPs: NBodyIP, bodyorder, fast,
                 evaluate_many!, evaluate_many_d!,
                 saveas, loadas,
                 saveas_json, loadas_json


@pot struct EnvBL{N, T, TVR, TVN} <: EnvBLFunction{N}
   t::Vector{Int}
   c::Vector{T}
   Vr::TVR
   Vn::TVN
   valN::Val{N}
   #
   # function EnvBL(t, c::Vector{T}, Vr::TVR, Vn::TVN, valN::Val{N}
   #               ) where {N, T, TVR, TVN}
   #    if N != bodyorder(Vn)
   #       error("EnvBL : body-orders must match")
   #    end
   #    return new(t, c, Vr, Vn, valN)
   # end
end

"""
An EnvBL is a product of a polynomial in `n` and an
`NBody`
"""
EnvBL


EnvBL(t, c, Vr, Vn) = EnvBL(t, c, Vr, Vn, Val(bodyorder(Vr)))

Vn(V::EnvBL) = V.Vn
Vr(V::EnvBL) = V.Vr

bodyorder(V::EnvBLFunction) = bodyorder(Vr(V))

function cutoff(V::EnvBLFunction)
   @assert cutoff(Vn(V)) <= cutoff(Vr(V))
   return cutoff(Vr(V::EnvBL))
end

# evaluate(V::EnvBL, n, r) =
#    evaluate(V.Vn, n) * evaluate(V.Vr, r)
#
# evaluate_d(V::EnvBL, n, r) =
#    evaluate_d(V, n, r) * evaluate(V.Vr, r),
#    evaluate(V.Vn, n) * evaluate_d(V.Vr, r)

n_fun(V::EnvBL, n) = sum( c * n^t  for (c,t) in zip(V.c, V.t) )

n_fun_d(V::EnvBL, n) = sum( t * c * n^(t-1)  for (c,t) in zip(V.c, V.t) )

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


end


# TODO TESTIN
#  * energy = force consistency
#  * This is actually simply a product of two site potentials;
#    can we not just move this into JuLIP?
