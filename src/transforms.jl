
using Base:@pure

export SpaceTransform, ExpTransform, PolyTransform


struct SpaceTransform{FT, FDT, VT} <: AbstractTransform
   id::String
   f::FT
   f_d::FDT
   valid::VT
end

SpaceTransform(id, f, f_d) = SpaceTransform(id, f, f_d, Val(Symbol(id)))

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
struct ExpTransform{TA, T} <: AbstractTransform
   A::TA
   r0::T
end

@pure transform(t::ExpTransform, r::Number) = exp( - t.A * (r/t.r0 - 1) )
@pure transform_d(t::ExpTransform, r::Number) = (-t.A/t.r0) * exp( - t.A * (r/t.r0 - 1) )



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
struct PolyTransform{TP, T} <: AbstractTransform
   p::TP
   r0::T
end

@pure transform(t::PolyTransform, r::Number) = @fastmath((t.r0/r)^t.p)
@pure transform_d(t::PolyTransform, r::Number) = @fastmath((-t.p/t.r0) * (t.r0/r)^(t.p+1))
