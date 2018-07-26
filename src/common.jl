
using StaticArrays

using JuLIP: AbstractCalculator, Atoms, neighbourlist, @D, JVec
using NeighbourLists: nbodies, maptosites!, maptosites_d!, virial!, max_neigs

import Base:Dict
import JuLIP.Potentials: evaluate, evaluate_d, evaluate_dd
import JuLIP: cutoff, energy, forces, site_energies, virial, stress

export NBodyIP,
       bodyorder,
       fast,
       dictionary,
       match_dictionary,
       rdf, idf,
       recover_basis,
       degree

"""
`NBodyFunction` : abstract supertype of all "pure" N-body functions.
concrete subtypes must implement

* `bodyorder`
* `evaluate`
* `evaluate_d`
"""
abstract type NBodyFunction{N} <: AbstractCalculator end

# prototypes of function defined on `NBodyFunction`
function bodyorder end

"""
some measure of degree - need not be polymomial degree?
TODO: figure out exactly how to use this
"""
degree(::NBodyFunction) = 0

"""
return the object attached to an NBodyFunction that describes
the underlying basis set
"""
function dictionary end

"""
this is to allow fix of a technical issue: it can be necessary to
replace a dictionary in an NBodyFunction with another dictionary that
is effectively the same but technically a different Julia type.
This occurs e.g. when loading basis sets from a file.
"""
function match_dictionary end

"""
recover a basis from an NBodyFunction
"""
function recover_basis end

"""
combine several basis functions into a single one
"""
function combine_basis end

include("eval_nbody.jl")


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

"""
`NBodyIP` : wraps `NBodyFunction`s into a JuLIP calculator, defining
`energy`, `forces`, `virial` and `cutoff`.

TODO: `site_energies`, etc.
"""
struct NBodyIP{TV} <: AbstractCalculator
   orders::Vector{TV}
end

cutoff(V::NBodyIP) = maximum( cutoff.(V.orders) )
energy(V::NBodyIP, at::Atoms) = sum( energy(Vn, at)  for Vn in V.orders )
forces(V::NBodyIP, at::Atoms) = sum( forces(Vn, at)  for Vn in V.orders )
virial(V::NBodyIP, at::Atoms) = sum( virial(Vn, at)  for Vn in V.orders )

"""
turn a potentially slow representation of an IP into a fast one,
by switching to a different representation.
"""
fast(IP::NBodyIP) = NBodyIP( fast.(IP.orders) )


# ----------------- some simplified access functions ------------------

evaluate(V::NBodyFunction{2}, r::Number) = evaluate(V, SVector(r))

evaluate_d(V::NBodyFunction{2}, r::Number) = evaluate_d(V, SVector(r))[1]

evaluate_dd(V::NBodyFunction{2}, r::Number) =
      ((@D V(r+1e-5)) - (@D V(r-1e-5))) / 1e-5

evaluate(V::NBodyFunction{3}, r1::Number, r2::Number, r3::Number) =
      evaluate(V, SVector(r1, r2, r3))




# ==================================================================
#    construct an NBodyIP from a basis
# ==================================================================


function NBodyIP(basis, coeffs)
   components = []
   tps = typeof.(basis)
   for tp in unique(tps)
      # find all basis functions that have the same type, which in particular
      # incorporated the body-order
      Itp = find(tps .== tp)
      # construct a new basis function by combining all of them in one
      # (this assumes that we have NBody types)
      V_N = combine_basis([basis[Itp]...], coeffs[Itp])
      push!(components, V_N)
   end
   return NBodyIP([components...])
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



# ------- distribution functions on the invariants --------------

function rdf(at::Atoms, rcut, transform=identity)
   if transform != identity
      return idf(2, at, rcut, transform)[1][1]
   end
   nlist = neighbourlist(at, rcut)
   return nlist.r
end

"""
`idf(N::Integer, at, rcut, transform)`

invariants distribution function => accumulates all the invariants values
arising during an N-body assembly over at.
"""
idf(N::Integer, at, rcut, transform) =
      _idf(Val(N), Val(bo2edges(N)), at, rcut, transform)

function _idf(valN::Val{N}, valM::Val{M}, at::Atoms{T}, rcut::T,
              transform) where {N, M, T}
   # compute invariants vectors to learn how many there are
   x = rand(SVector{M,T})
   I1, I2 = invariants(x)

   I1acc = [ T[] for n = 1:length(I1) ]
   I2acc = [ T[] for n = 1:length(I2) ]
   Iacc = (I1acc, I2acc)

   for (i, j, r, R) in sites(at, rcut)
      eval_site_nbody!(valN, R, rcut,
                  (Iacc, s, _1, _2, _3) -> idf_accumulator(Iacc, transform.(s)),
                  Iacc, nothing)
   end

   return I1acc, I2acc
end

function idf_accumulator(Iacc, x)
   I1, I2 = invariants(x)
   for n = 1:length(I1)
      push!(Iacc[1][n], I1[n])
   end
   for n = 1:length(I2)
      push!(Iacc[2][n], I2[n])
   end
   return Iacc
end




# =============== Experimental:
#   evaluate NBodyIP

(V::NBodyIP)(args...) = evaluate(V, args...)

evaluate(V::NBodyIP, r::Number) = evaluate(V::NBodyIP, SVector(r))

evaluate(V::NBodyIP, r1::T, r2::T, r3::T) where {T <: Number} =
      evaluate(V::NBodyIP, SVector(r1, r2, r3))

function evaluate(V::NBodyIP, r::SVector{N, T}) where {N, T}
   v = zero(T)
   for Vn in V.orders
      if bo2edges(bodyorder(Vn)) == N
         v += Vn(r)
      end
   end
   return v
end


# ------------------ IO of NBodyIPs ---------------
# -------------- Serialising and deserialising ---------------
# TODO: Move to IO submodule

export save_ip, load_ip

using FileIO: load
import FileIO: save
using JSON

struct XJld2 end
struct XJson end
struct XJld end

Dict(IP::NBodyIP) = Dict("id" => "NBodyIP",
                         "orders" => Dict.(IP.orders))

NBodyIP(D::Dict) = NBodyIP(_decode_dict.(D["orders"]))
Base.convert(::Val{:NBodyIP}, D::Dict) = NBodyIP(D)

function _checkextension(fname)
   if fname[end-3:end] == "jld2"
      return XJld2()
   elseif fname[end-2:end] == "jld"
      return XJld()
   elseif fname[end-3:end] == "json"
      return XJson()
   end
   error("the filename should end in `jld2`")
end

save(fname::AbstractString, IP::NBodyIP) =
   save_ip(_checkextension(fname), fname, IP)

load_ip(fname::AbstractString) =
   load_ip(_checkextension(fname), fname)

save_ip(::Union{XJld2, XJld}, fname, IP) = save(fname, "IP", Dict(IP))

function load_ip(::Union{XJld2, XJld}, fname)
   IPs = nothing
   try
      IPs = load(fname, "IP")
   catch
      IPs = load(fname, "IP")
   end
   return NBodyIP(IPs)
end

function save_ip(::XJson, fname, IP)
   f = open(fname, "w")
   print(f, JSON.json(Dict(IP)))
   close(f)
end

function load_ip(::XJson, fname)
   IPj = JSON.parsefile(fname)
   @assert IPj["id"] == "NBodyIPs.NBodyIP"
   return NBodyIP(IPj)
end
