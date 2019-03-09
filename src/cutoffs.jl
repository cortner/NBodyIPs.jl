

using JuLIP.Potentials: coscut, coscut_d

import JuLIP: decode_dict

abstract type NBCutoff end

# most of the time the type itself can server a descriptor for
# combining basis functions
combiscriptor(C::NBCutoff) = C


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
fcut(C::CosCut, r) = coscut(r, C.rc1, C.rc2)
fcut_d(C::CosCut, r) = coscut_d(r, C.rc1, C.rc2)
Dict(C::CosCut) = Dict("__id__" => "NBodyIPs_CosCut",
                       "rc1" => C.rc1, "rc2" => C.rc2)
CosCut(D::Dict) = CosCut(D["rc1"], D["rc2"])
decode_dict(::Val{:NBodyIPs_CosCut}, D::Dict) = CosCut(D)


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
fcut(C::CosCut2s, r) = (1-coscut(r, C.ri1, C.ri2)) * coscut(r, C.ro1, C.ro2)
fcut_d(C::CosCut2s, r) = (- coscut_d(r, ri1, ri2) * coscut(r, ro1, ro2)
                       + (1-coscut(r, ri1, ri2)) * coscut_d(r, ro1, ro2))
Dict(C::CosCut) = Dict( "__id__" => "NBodyIPs_CosCut2s",
                        "ri1" => C.ri1, "ri2" => C.ri2,
                        "ro1" => C.ro1, "ro2" => C.ro2 )
CosCut2s(D::Dict) = CosCut2s(D["ri1"], D["ri2"], D["ro1"], D["ro2"])
decode_dict(::Val{:NBodyIPs_CosCut2s}, D::Dict) = CosCut2s(D)


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
fcut(C::PolyCut, r) = @fastmath((r/rc - 1)^p * (r < rc))
fcut_d(C::PolyCut, r) = @fastmath(p/rc * (r/rc - 1)^(p-1) * (r-rc))
Dict(C::PolyCut) = Dict( "__id__" => "NBodyIPs_PolyCut",
                        "p" => C.p, "rc" => C.rc)
PolyCut(D::Dict) = PolyCut(D["p"], D["rc"])
decode_dict(::Val{:NBodyIPs_PolyCut}, D::Dict) = PolyCut(D)


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
fcut(C::PolyCutSym, r) = @fastmath(((r/rc)^2 - 1)^p * (r<rc))
fcut_d(C::PolyCutSym, r) = @fastmath(2*p*r/rc^2 * ((r/rc)^2 - 1)^(p-1) * (r-rc))
Dict(C::PolyCutSym) = Dict( "__id__" => "NBodyIPs_PolyCutSym",
                        "p" => C.p, "rc" => C.rc)
PolyCutSym(D::Dict) = PolyCutSym(D["p"], D["rc"])
decode_dict(::Val{:NBodyIPs_PolyCutSym}, D::Dict) = PolyCutSym(D)


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
   @show λopt = find_zero(gg, -2.5)
   C = 1 / (exp( λopt*(rc/rnn - 1) ) - 1)
   return PolyCut2sA(r0, rnn, rc, C, λopt)
end

function fcut(C::PolyCut2sA, r)
   x = @fastmath(C.C * (exp( C.λ*(r/C.rnn - 1) ) - 1))
   return (x^2-1)^2 * (r0 < r < rc)
end

function fcut_d(C::PolyCut2sA, r)
   e =  @fastmath(exp( C.λ*(r/C.rnn - 1) ))
   x = C.C * (e - 1)
   dx = C.C * e * C.λ / C.rnn
   return 4*x*(x^2-1) * dx * (r0 < r < rc)
end

Dict(C::PolyCut2sA) = Dict( "__id__" => "NBodyIPs_PolyCut2sA",
                            "r0" => C.r0, "rnn" => C.rnn, "rc" => C.rc,
                            "C" => C.C, "lam" => C.λ )

PolyCut2sA(D::Dict) = PolyCut2sA(D["r0"], D["rnn"], D["rc"], D["C"], D["lam"])

decode_dict(::Val{:NBodyIPs_PolyCut2sA}, D::Dict) = PolyCut2sA(D)



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
