# TODO:
#   - allow Vn to be an arbitrary pair potential
#   - replace the polynomial with an arbitrary family of pair potentials
#   - (or even nbody?)

module EnvIPs

using StaticArrays
using JuLIP:              AbstractCalculator
using JuLIP.Potentials:   Shift,
                          @analytic,
                          @pot
using NBodyIPs:           NBodyDescriptor
using NBodyIPs.PolyBasis: nbpolys

import Base:              Dict,
                          ==,
                          convert,
                          hash
import NBodyIPs:          NBodyIP,
                          bodyorder,
                          fast,
                          combinebasis,
                          _decode_dict,
                          descriptor,
                          degree,
                          basisname,
                          BASIS


export envpolys, info

abstract type AbstractEnvIP{N} <: AbstractCalculator end

struct EnvIP{N, P, TVR, TVN} <: AbstractEnvIP{N}
   t::Int
   Vr::TVR     # N-body potential
   Vn::TVN     # neighbour counter
   str_Vn::String  # string describing the neighbour counter
   valN::Val{N}
   valP::Val{P}
end

@pot EnvIP

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
   Vn1 = eval(Meta.parse("@analytic r -> " * str_Vn))
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

bodyorder(V::AbstractEnvIP{N}) where {N} = Int(N)

descriptor(V::EnvIP) = descriptor(V.Vr)

hash(::BASIS, V::EnvIP) = hash((EnvIP,
                                hash(BASIS(), V.Vr),
                                hash(V.str_Vn),
                                hash(V.t)))

function degree(V::EnvIP)
   if length(V.Vr) == 1
      return ( degree(V.Vr), V.t )
   end
   error("`degree` is only defined for `EnvIP` basis functions, length == 1")
end

basisname(::EnvIP) = "EnvIP"

# ----------------- generate basis / IP / convert ----------------

function envpolys(N::Integer, D::NBodyDescriptor, deg_poly::Integer,
                  Vn_descr, deg_n::Integer; kwargs...)
   B_poly = nbpolys(N, D, deg_poly; kwargs...)
   B = EnvIP[]
   str_Vn = Vn_descr[1]
   Vn = analyse_Vn(Vn_descr...)
   for deg = 0:deg_n
      append!(B, [EnvIP(deg, Vr, Vn, str_Vn) for Vr in B_poly])
   end
   return [b for b in B]
end

function combinebasis(basis::AbstractVector{TV}, coeffs) where {TV <: EnvIP}
   # @assert isleaftype(TV)
   @assert all( b.t == basis[1].t for b in basis )
   # combine the Vr components of the basis functions
   # (we get to do this because all t (=P) are the same
   Vr = combinebasis( [b.Vr for b in basis], coeffs )
   return EnvIP(basis[1].t, Vr, basis[1].Vn, basis[1].str_Vn)
end


function info(B::Vector{T}; indent = 2) where T <: EnvIP
   ind = repeat(" ", indent)
   println(ind * "EnvIP with P = $(B[1].t)")
   println(ind * "           Vn : $(B[1].str_Vn)")
   println(ind * "           Vr : ...")
   info([ b.Vr for b in B ], indent = indent+5)
end

fast(V::EnvIP) = EnvIP(V.t, fast(V.Vr), V.Vn, V.str_Vn)



# ========================================================
#    Re-Implement Fast Version of EnvIPs
# ========================================================

struct EnvPoly{N, TVR, TVN} <: AbstractEnvIP{N}
   Vr::Vector{TVR}    # N-body potentials multiplied by n^j, j = 0, 1, ...
   Vn::TVN            # neighbour counter
   str_Vn::String     # string describing the neighbour counter
   valN::Val{N}
end

@pot EnvPoly

envdegree(V::EnvPoly) = length(V.P) - 1

==(V1::EnvPoly, V2::EnvPoly) = ( (V1.Vr == V2.Vr) &&
                                 (V1.str_Vn == V2.str_Vn) &&
                                 (V1.valN == V2.valN) )

# Dict(V::EnvPoly) = Dict( "__id__" => "EnvPoly",
#                          "Vr" => Dict.(V.Vr),
#                          "str_Vn" => V.str_Vn,
#                          "cutoff_Vn" => cutoff(V.Vn)  )
#
# EnvPoly(D::Dict) = EnvPoly( _decode_dict.(D["Vr"]),
#                             analyse_Vn(D["str_Vn"], D["cutoff_Vn"]),
#                             D["str_Vn"] )

convert(::Val{:EnvPoly}, D::Dict) = EnvPoly(D)

EnvPoly(t, Vr, str_Vn::String, cutoff_Vn::AbstractFloat) =
   EnvPoly(t, Vr, analyse_Vn(str_Vn, cutoff_Vn), str_Vn)

descriptor(V::EnvPoly) = descriptor(V.Vr[1])

hash(::BASIS, V::EnvPoly) = hash(( hash(EnvPoly),
                                   hash.(Ref(BASIS()), V.Vr),
                                   hash(V.str_Vn) ))

Vn(V::EnvPoly) = V.Vn
Vr(V::EnvPoly) = V.Vr

function EnvPoly(Vr::AbstractVector, Vn, str_Vn::String)
   # for now we put the burden on the user to create a Vr vector where
   # the types match But we can still collect them here for convenience
   Vr_coll = [ v for v in Vr ]
   if !isconcretetype(typeof(Vr_coll))
      error("""For now, `EnvPoly` can only be constructed from an array of
               potentials that all share the same concrete type.""")
   end
   return EnvPoly(Vr_coll, Vn, str_Vn, Val(bodyorder(Vr[1])))
end


# function degree(V::EnvIP)
#    if length(V.Vr) == 1
#       return ( degree(V.Vr), V.t )
#    end
#    error("`degree` is only defined for `EnvIP` basis functions, length == 1")
# end

basisname(::EnvPoly) = "EnvPoly"




# all evaluation and assembly in here:
include("eval_environ.jl")


end
