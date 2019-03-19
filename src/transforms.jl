
using Base:@pure

export AnalyticTransform, ExpTransform, PolyTransform,
       IdTransform, CoulombTransform, MorseTransform
import Base: convert, Dict, ==

# TODO: rename this AnalyticTransform
struct AnalyticTransform{FT, FDT, VT} <: SpaceTransform
   id::String
   f::FT
   f_d::FDT
   valid::VT
end

@pure transform(t::AnalyticTransform, r::Number) = t.f(r)
@pure transform_d(t::AnalyticTransform, r::Number) = t.f_d(r)

AnalyticTransform(id, f, f_d) = AnalyticTransform(id, f, f_d, Val(Symbol(id)))

==(T1::AnalyticTransform, T2::AnalyticTransform) = (T1.id == T2.id)

function AnalyticTransform(strans::String; fwrap = true)
   strans0 = strans
   # if @analytic is a substring then we don't do anything
   if !occursin(r"@analytic", strans)
      # but if not, then we next look for ->
      if !occursin(r"->", strans)
         # if -> is not a substring then we assume that strans is of the form
         # "(r0/r)^4" or similar i.e. explicitly uses r as the variable.
         strans = "@analytic r -> " * strans
      else
         # @analytic is not a substring but -> is a substring. e.g.
         # r -> (r0/r)^3 or s -> (r0/s)^3. we add @analytic to compile it
         strans = "@analytic " * strans
      end
   end
   ftrans = eval(Meta.parse(strans))
   if fwrap
      return AnalyticTransform(strans0, F64fun(ftrans.f), F64fun(ftrans.f_d))
   end
   return AnalyticTransform(strans0, ftrans.f, ftrans.f_d)
end

Dict(t::AnalyticTransform) = Dict( "__id__" => "NBodyIPs_AnalyticTransform",
                                "defn" => t.id )
AnalyticTransform(D::Dict) = AnalyticTransform(D["defn"])
convert(::Val{:NBodyIPs_AnalyticTransform}, D::Dict) = AnalyticTransform(D)
hash(::BASIS, t::AnalyticTransform) = hash((AnalyticTransform, t.id))



"""
Implements the space transform
```
r -> exp( - A * (r/r0 - 1) )
```

Constructor:
```
ExpTransform(A, r0)
```
"""
struct ExpTransform{TA, T} <: SpaceTransform
   A::TA
   r0::T
end

@pure transform(t::ExpTransform, r::Number) = exp( - t.A * (r/t.r0 - 1))
@pure transform_d(t::ExpTransform, r::Number) = (-t.A/t.r0) * exp( - t.A * (r/t.r0 - 1) )

# e^{-A (r/r0 - 1)) = x
# r/r0 - 1 = - (log x) / A
# r = r0 * (1 - (log x) / A)
# TODO: write a test that checks this inverse
inv_transform(t::ExpTransform, x::Number) = t.r0 * (1 - log(x) / t.A)

Dict(t::ExpTransform) =
   Dict("__id__" => "NBodyIPs_ExpTransform", "A" => t.A, "r0" => t.r0)
ExpTransform(D::Dict) = ExpTransform(D["A"], D["r0"])
convert(::Val{:NBodyIPs_ExpTransform}, D::Dict) = ExpTransform(D)
hash(::BASIS, t::ExpTransform) = hash(t)

"""
Implements the space transform
```
r -> (r0/r)^p
```

Constructor:
```
PolyTransform(p, r0)
```
"""
struct PolyTransform{TP, T} <: SpaceTransform
   p::TP
   r0::T
end

@pure transform(t::PolyTransform, r::Number) = @fastmath((t.r0/r)^t.p)
@pure transform_d(t::PolyTransform, r::Number) = @fastmath((-t.p/t.r0) * (t.r0/r)^(t.p+1))

# x = (r0/r)^p
# r x^{1/p} = r0
inv_transform(t::PolyTransform, x::Number) = t.r0 / x^(1/t.p)

Dict(t::PolyTransform) =
   Dict("__id__" => "NBodyIPs_PolyTransform",
        "p" => t.p,
        "r0" => t.r0)
PolyTransform(D::Dict) = PolyTransform(D["p"], D["r0"])
convert(::Val{:NBodyIPs_PolyTransform}, D::Dict) = PolyTransform(D)
hash(::BASIS, t::PolyTransform) = hash(t)


"""
Implements the space transform
```
r -> 1/r
```
Constructor: `CoulombTransform()`
"""
struct CoulombTransform <: SpaceTransform
end
@pure transform(t::CoulombTransform, r::Number) = @fastmath(1/r)
@pure transform_d(t::CoulombTransform, r::Number) = @fastmath(-1/(r*r))
inv_transform(t::CoulombTransform, x::Number) = 1/x

Dict(t::CoulombTransform) = Dict("__id__" => "NBodyIPs_CoulombTransform")
CoulombTransform(D::Dict) = CoulombTransform()
convert(::Val{:NBodyIPs_CoulombTransform}, D::Dict) = CoulombTransform(D)
hash(::BASIS, t::CoulombTransform) = hash(t)


"""
Implements the space transform
```
r -> r
```
Constructor: `IdTransform()`
"""
struct IdTransform <: SpaceTransform
end

@pure transform(t::IdTransform, r::Number) = r
@pure transform_d(t::IdTransform, r::Number) = one(typeof(r))
inv_transform(t::IdTransform, x::Number) = x

Dict(t::IdTransform) = Dict("__id__" => "NBodyIPs_IdTransform")
IdTransform(D::Dict) = IdTransform()
convert(::Val{:NBodyIPs_IdTransform}, D::Dict) = IdTransform(D)
hash(::BASIS, t::IdTransform) = hash(t)
