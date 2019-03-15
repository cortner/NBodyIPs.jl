
module Cutoffs

using JuLIP.Potentials: cutsw,
                        cutsw_d,
                        coscut,
                        coscut_d

import JuLIP
using Roots, StaticArrays

const cutsp = JuLIP.Potentials.fcut
const cutsp_d = JuLIP.Potentials.fcut_d

import JuLIP: decode_dict
import JuLIP.Potentials: cutoff
import NBodyIPs: fcut, fcut_d, NBCutoff
import Base: Dict, convert

export CosCut, CosCut2s, PolyCut, PolyCutSym, PolyCut2sA


"""
`struct CosCut{T} <: NBCutoff`

### Constructor
```
CosCut(rc1, rc2)
```

### Cut-off function
```
fcut(r) = (r < rc1) + (rc1 <= r < rc2) * (1 + cos(...))/2
```
where `cos(...)` is a rescaled cos to blend from one to zero.
"""
struct CosCut{T} <: NBCutoff
   rc1::T
   rc2::T
end
cutoff(C::CosCut) = C.rc2
fcut(C::CosCut, r::Number) = coscut(r, C.rc1, C.rc2)
fcut_d(C::CosCut, r::Number) = coscut_d(r, C.rc1, C.rc2)
Dict(C::CosCut) = Dict("__id__" => "NBodyIPs_CosCut",
                       "rc1" => C.rc1, "rc2" => C.rc2)
CosCut(D::Dict) = CosCut(D["rc1"], D["rc2"])
convert(::Val{:NBodyIPs_CosCut}, D::Dict) = CosCut(D)


"""
`struct CosCut2s{T} <: NBCutoff`

### Constructor
```
CosCut2s(ri1, ri2, ro1, ro2)
```

### Cut-off function
```
fcut(r) = (1-coscut(r, ri1, ri2)) * coscut(r, ro1, ro2)
```
"""
struct CosCut2s{T} <: NBCutoff
   ri1::T
   ri2::T
   ro1::T
   ro2::T
end
cutoff(C::CosCut2s) = C.ro2
fcut(C::CosCut2s, r::Number) = (1-coscut(r, C.ri1, C.ri2)) * coscut(r, C.ro1, C.ro2)
fcut_d(C::CosCut2s, r::Number) = (- coscut_d(r, C.ri1, C.ri2) * coscut(r, C.ro1, C.ro2)
   + (1-coscut(r, C.ri1, C.ri2)) * coscut_d(r, C.ro1, C.ro2))
Dict(C::CosCut2s) = Dict( "__id__" => "NBodyIPs_CosCut2s",
                        "ri1" => C.ri1, "ri2" => C.ri2,
                        "ro1" => C.ro1, "ro2" => C.ro2 )
CosCut2s(D::Dict) = CosCut2s(D["ri1"], D["ri2"], D["ro1"], D["ro2"])
convert(::Val{:NBodyIPs_CosCut2s}, D::Dict) = CosCut2s(D)


"""
`struct PolyCut{TI, T} <: NBCutoff`

### Constructor
```
PolyCut(p, rc)
```

### Cut-off function
```
fcut(r) = (r/rc-1)^p
```
"""
struct PolyCut{TI, T} <: NBCutoff
   p::TI
   rc::T
end
cutoff(C::PolyCut) = C.rc
fcut(C::PolyCut, r::Number) = @fastmath((r/C.rc - 1)^C.p * (r < C.rc))
fcut_d(C::PolyCut, r::Number) = @fastmath(C.p/C.rc * (r/C.rc - 1)^(C.p-1) * (r < C.rc))
Dict(C::PolyCut) = Dict( "__id__" => "NBodyIPs_PolyCut",
                        "p" => C.p, "rc" => C.rc)
PolyCut(D::Dict) = PolyCut(D["p"], D["rc"])
convert(::Val{:NBodyIPs_PolyCut}, D::Dict) = PolyCut(D)


"""
`struct PolyCutSym{TI, T} <: NBCutoff`

### Constructor
```
PolyCutSym(p, rc)
```

### Cut-off function
```
fcut(r) = (r/rc-1)^p (r/rc+1)^p
```
"""
struct PolyCutSym{TI, T} <: NBCutoff
   p::TI
   rc::T
end
cutoff(C::PolyCutSym) = C.rc
fcut(C::PolyCutSym, r::Number) = @fastmath(((r/C.rc)^2 - 1)^C.p * (r<C.rc))
fcut_d(C::PolyCutSym, r::Number) = @fastmath(2*C.p*r/C.rc^2 * ((r/C.rc)^2 - 1)^(C.p-1) * (r<C.rc))
Dict(C::PolyCutSym) = Dict( "__id__" => "NBodyIPs_PolyCutSym",
                        "p" => C.p, "rc" => C.rc)
PolyCutSym(D::Dict) = PolyCutSym(D["p"], D["rc"])
convert(::Val{:NBodyIPs_PolyCutSym}, D::Dict) = PolyCutSym(D)


"""
`struct PolyCut2sA{T} <: NBCutoff`

experimental...
"""
struct PolyCut2sA{T} <: NBCutoff
   r0::T
   rnn::T
   rc::T
   C::T
   λ::T
end

cutoff(C::PolyCut2sA) = C.rc

function PolyCut2sA(r0, rnn, rc)
   gg = λ -> exp( λ*(rc/rnn - 1) ) + exp( λ*(r0/rnn - 1) ) - 2
   λopt = find_zero(gg, -2.5)
   if abs(λopt) < 1e-7
      @error("`Roots` found the *bad* λ parameter.")
   end
   C = 1 / (exp( λopt*(rc/rnn - 1) ) - 1)
   return PolyCut2sA(r0, rnn, rc, C, λopt)
end

function fcut(C::PolyCut2sA, r::Number)
   x = @fastmath(C.C * (exp( C.λ*(r/C.rnn - 1) ) - 1))
   return (x^2-1)^2 * (C.r0 < r < C.rc)
end

function fcut_d(C::PolyCut2sA, r::Number)
   e =  @fastmath(exp( C.λ*(r/C.rnn - 1) ))
   x = C.C * (e - 1)
   dx = C.C * e * C.λ / C.rnn
   return 4*x*(x^2-1) * dx * (C.r0 < r < C.rc)
end

Dict(C::PolyCut2sA) = Dict( "__id__" => "NBodyIPs_PolyCut2sA",
                            "r0" => C.r0, "rnn" => C.rnn, "rc" => C.rc,
                            "C" => C.C, "lam" => C.λ )

PolyCut2sA(D::Dict) = PolyCut2sA(D["r0"], D["rnn"], D["rc"], D["C"], D["lam"])

convert(::Val{:NBodyIPs_PolyCut2sA}, D::Dict) = PolyCut2sA(D)



# -------------- Some General cutoff mechanics

fcut(C::NBCutoff, r::AbstractVector) = prod(fcut.(Ref(C), r))

function fcut_d(C::NBCutoff, r::AbstractVector)
   fc = fcut.(Ref(C), r)
   fctot = prod(fc)
   if fctot == 0
      return fctot, zero(typeof(r))
   end
   fc_d = fcut_d.(Ref(C), r)
   return fctot, (fctot * fc_d) ./ fc
end

end



# struct SplineCut{T} <: NBCutoff
# end
#
# struct SplineCut2s{T} <: NBCutoff
# end
#
# struct CosEnvCut{T} <: NBCutoff
# end
#
# struct CosEnvCut2s{T} <: NBCutoff
# end
