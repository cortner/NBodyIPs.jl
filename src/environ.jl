
module EnvBLs

using JuLIP, StaticArrays

using JuLIP.Potentials: evaluate, evaluate_d, Shift, @analytic

using NeighbourLists: max_neigs

using NBodyIPs: eval_site_nbody!, _grad_len2pos!, _acc_manyfrcs

using NBodyIPs.BLPolys: BLDictionary, bl_basis

abstract type EnvBLFunction{N} <: AbstractCalculator end

import Base: Dict, ==, convert

import JuLIP: cutoff, energy, forces, virial

import NBodyIPs: NBodyIP, bodyorder, fast,
                 evaluate_many!, evaluate_many_d!,
                 combinebasis

export envbl_basis

@pot struct EnvBL{N, P, TVR, TVN} <: EnvBLFunction{N}
   t::Int
   Vr::TVR
   Vn::TVN
   str_Vn::String
   valN::Val{N}
   valP::Val{P}
end

==(V1::EnvBL, V2::EnvBL) = ( (V1.t == V2.t) &&
                             (V1.Vr == V2.Vr) &&
                             (V1.str_Vn == V2.str_Vn) &&
                             (V1.valN == V2.valN) )

Dict(V::EnvBL) = Dict( "__id__" => "EnvBL",
                       "t" => V.t,
                       "Vr" => Dict(V.Vr),
                       "str_Vn" => V.str_Vn,
                       "cutoff_Vn" => cutoff(V.Vn)  )

EnvBL(D::Dict) = EnvBL( D["t"],
                        _decode_dict(D["Vr"]),
                        analyse_Vn(D["str_Vn"], D["cutoff_Vn"]) )

convert(::Val{:EnvBL}, D::Dict) = EnvBL(D)


EnvBL(t::Int, Vr, Vn, str_Vn::String) =
      EnvBL(t, Vr, Vn, str_Vn, Val(bodyorder(Vr)), Val(t))

function analyse_Vn(str_Vn, cutoff_Vn)
   Vn1 = eval(parse("@analytic r -> " * str_Vn))
   Vn = Shift(Val(1), Vn1, cutoff_Vn,
              Base.invokelatest(Vn1.f,   cutoff_Vn),
              Base.invokelatest(Vn1.f_d, cutoff_Vn),
              0.0)
   return Vn
end


function EnvBL(t, Vr, str_Vn::String, cutoff_Vn::AbstractFloat)
   Vn = analyse_Vn(str_Vn, cutoff_Vn)
   return EnvBL(t, Vr, Vn, str_Vn)
end

Vn(V::EnvBL) = V.Vn
Vr(V::EnvBL) = V.Vr

bodyorder(V::EnvBLFunction) = bodyorder(Vr(V))



# ----------------- generate basis / IP / convert ----------------

function envbl_basis(N::Integer, D::BLDictionary, deg_poly::Integer,
                     Vn_descr, deg_n::Integer; kwargs...)
   B_poly = bl_basis(N, D, deg_poly; kwargs...)
   B = EnvBL[]
   str_Vn = Vn_descr[1]
   Vn = analyse_Vn(Vn_descr...)
   for deg = 0:deg_n
      append!(B, [EnvBL(deg, Vr, Vn, str_Vn) for Vr in B_poly])
   end
   return [b for b in B]
end


function combinebasis(basis::AbstractVector{TV}, coeffs) where {TV <: EnvBL}
   @assert isleaftype(TV)
   @assert all( b.t == basis[1].t for b in basis )
   # combine the Vr components of the basis functions
   # (we get to do this because all t (=P) are the same
   Vr = combinebasis( [b.Vr for b in basis], coeffs )
   return EnvBL(basis[1].t, Vr, basis[1].Vn, basis[1].str_Vn)
end


function Base.info(B::Vector{T}; indent = 2) where T <: EnvBL
   ind = repeat(" ", indent)
   println(ind * "EnvBL with P = $(B[1].t)")
   println(ind * "           Vn : $(B[1].str_Vn)")
   println(ind * "           Vr : see below...")

   info([ b.Vr for b in B ], indent = indent+2)
end

fast(V::EnvBL) = EnvBL(V.t, fast(V.Vr), V.Vn, V.str_Vn)



# all evaluation and assembly in here:
include("eval_environ.jl")

end
