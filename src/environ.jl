
module EnvBLs

using JuLIP

using NBodyIPs.Polys: Dictionary, NBody

abstract type EnvBLFunction{N} <: AbstractCalculator


"""
An EnvBL is a product of a polynomial in `n` and an
`NBody`
"""
struct EnvBL{N, T, TVR <: NBody{N}, TVN}
   t::Vector{Int}
   c::Vector{T}
   Vr::TVR
   Vn::TVN
   valN::Val{N}
end

Vn(V::EnvBL) = V.Vn
Vr(V::EnvBL) = V.Vr

# turn it into a potential
@pot EnvBL

# evaluate(V::EnvBL, n, r) =
#    evaluate(V.Vn, n) * evaluate(V.Vr, r)
#
# evaluate_d(V::EnvBL, n, r) =
#    evaluate_d(V, n, r) * evaluate(V.Vr, r),
#    evaluate(V.Vn, n) * evaluate_d(V.Vr, r)

n_fun(V::EnvBL, n) = sum( c * n^t  for (c,t) in zip(V.c, V.t) )

n_fun_d(V::EnvBL, n) = sum( t * c * n^(t-1)  for (c,t) in zip(V.c, V.t) )

site_ns(V::EnvBLFunction{N}, at) =
   nfun.(V, site_energies(Vn(V), at))

function site_n_d!(dVn, V::EnvBLFunction{N}, r, R, Ni)
   dNi = n_fun_d(Ni)
   for n = 1:length(r)
      dVn[n] = dNi * grad(Vn(V), r[n], R[n])
   end
   return dVn
end

function energy(V::EnvBLFunction{N}, at::Atoms)
   # compute the n values
   Ns = site_ns(V, at)
   # compute the inner v values
   Vs = site_energies(Vr(V), at)
   # the energy is just the dot-product now
   return dot(Ns, Vs)
end


function forces(V::EnvBLFunction{N}, at::Atoms)
   @assert cutoff(V.Vn) <= cutoff(V.Vr)

   # compute the n values
   Ns = site_ns(V, at)
   # compute the inner v values
   Vs = site_energies(Vr(V), at)

   # now assemble site forces and use those to create
   # total forces by mixing N and V
   nlist = neighbourlist(at, cutoff(V))
   maxneigs = max_neigs(nlist)
   F = zeros(JVec{T}, length(at))
   dVsite = zeros(JVec{T}, maxneigs)
   dVn = zeros(JVec{T}, maxneigs)

   for (i, j, r, R) in sites(nlist)
      dVsite .*= 0.0
      eval_site_nbody!(
            Val(N), R, cutoff(V),
            (out, s, S, J, _) -> _grad_len2pos!(out, evaluate_d(V, s)/N, J, S),
            dVsite, nothing )

      site_n_d!(dVn, V, r, R, Ns[i])

      # write site energy gradient into forces
      for n = 1:length(j)
         f = Ns[i] * dVsite[n] + Vs[i] * dVn[n]
         F[j[n]] -= f
         F[i] += f
      end
   end
   return F
end


end


# TODO TESTIN
#  * energy = force consistency
#  * This is actually simply a product of two site potentials;
#    can we not just move this into JuLIP?
