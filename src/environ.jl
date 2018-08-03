
module EnvIPs

using StaticArrays
using JuLIP:              AbstractCalculator
using JuLIP.Potentials:   Shift,
                          @analytic,
                          @pot
using NBodyIPs.BLPolys:   BLDictionary,
                          bl_basis

import Base:              Dict,
                          ==,
                          convert
import NBodyIPs:          NBodyIP,
                          bodyorder,
                          fast,
                          combinebasis,
                          _decode_dict


export envbl_basis

abstract type AbstractEnvIP{N} <: AbstractCalculator end

@pot struct EnvIP{N, P, TVR, TVN} <: AbstractEnvIP{N}
   t::Int
   Vr::TVR
   Vn::TVN
   str_Vn::String
   valN::Val{N}
   valP::Val{P}
end

==(V1::EnvIP, V2::EnvIP) = ( (V1.t == V2.t) &&
                             (V1.Vr == V2.Vr) &&
                             (V1.str_Vn == V2.str_Vn) &&
                             (V1.valN == V2.valN) )

Dict(V::EnvIP) = Dict( "__id__" => "EnvIP",
                       "t" => V.t,
                       "Vr" => Dict(V.Vr),
                       "str_Vn" => V.str_Vn,
                       "cutoff_Vn" => cutoff(V.Vn)  )

EnvIP(D::Dict) = EnvIP( D["t"],
                        _decode_dict(D["Vr"]),
                        analyse_Vn(D["str_Vn"], D["cutoff_Vn"]),
                        D["str_Vn"] )

convert(::Val{:EnvIP}, D::Dict) = EnvIP(D)


EnvIP(t::Int, Vr, Vn, str_Vn::String) =
      EnvIP(t, Vr, Vn, str_Vn, Val(bodyorder(Vr)), Val(t))

function analyse_Vn(str_Vn, cutoff_Vn)
   Vn1 = eval(parse("@analytic r -> " * str_Vn))
   Vn = Shift(Val(1), Vn1, cutoff_Vn,
              Base.invokelatest(Vn1.f,   cutoff_Vn),
              Base.invokelatest(Vn1.f_d, cutoff_Vn),
              0.0)
   return Vn
end


function EnvIP(t, Vr, str_Vn::String, cutoff_Vn::AbstractFloat)
   Vn = analyse_Vn(str_Vn, cutoff_Vn)
   return EnvIP(t, Vr, Vn, str_Vn)
end

Vn(V::EnvIP) = V.Vn
Vr(V::EnvIP) = V.Vr

bodyorder(V::AbstractEnvIP) = bodyorder(Vr(V))


# ----------------- generate basis / IP / convert ----------------

function envbl_basis(N::Integer, D::BLDictionary, deg_poly::Integer,
                     Vn_descr, deg_n::Integer; kwargs...)
   B_poly = bl_basis(N, D, deg_poly; kwargs...)
   B = EnvIP[]
   str_Vn = Vn_descr[1]
   Vn = analyse_Vn(Vn_descr...)
   for deg = 0:deg_n
      append!(B, [EnvIP(deg, Vr, Vn, str_Vn) for Vr in B_poly])
   end
   return [b for b in B]
end


function combinebasis(basis::AbstractVector{TV}, coeffs) where {TV <: EnvIP}
   @assert isleaftype(TV)
   @assert all( b.t == basis[1].t for b in basis )
   # combine the Vr components of the basis functions
   # (we get to do this because all t (=P) are the same
   Vr = combinebasis( [b.Vr for b in basis], coeffs )
   return EnvIP(basis[1].t, Vr, basis[1].Vn, basis[1].str_Vn)
end


function Base.info(B::Vector{T}; indent = 2) where T <: EnvIP
   ind = repeat(" ", indent)
   println(ind * "EnvIP with P = $(B[1].t)")
   println(ind * "           Vn : $(B[1].str_Vn)")
   println(ind * "           Vr : ...")
   info([ b.Vr for b in B ], indent = indent+5)
end

fast(V::EnvIP) = EnvIP(V.t, fast(V.Vr), V.Vn, V.str_Vn)

# all evaluation and assembly in here:
include("eval_environ.jl")

end
