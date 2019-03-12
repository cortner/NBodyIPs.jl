import Base: convert, ==
using JuLIP.Potentials: @pot
export OneBody

descriptor(::OneBody) = NullDesc()

hash(::BASIS, ::OneBody) = hash(OneBody)

evaluate(V::OneBody) = V.E0
site_energies(V::NBodyFunction{1}, at::Atoms) = fill(V(), length(at))
forces(V::NBodyFunction{1}, at::Atoms{T}) where {T} =
      zeros(SVector{3, T}, length(at))
virial(V::NBodyFunction{1}, at::Atoms{T}) where {T} =
      zero(SMatrix{3, 3, T})
energy(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ energy(b, at) for b in B ]
forces(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ forces(b, at) for b in B ]
virial(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: NBodyFunction{1}, T} =
   [ virial(b, at) for b in B ]
combinebasis(Vs::AbstractVector{<: OneBody},
             Cs::AbstractVector{<: Real}) =
             OneBody( sum( V.E0 * c for (V, c) in zip(Vs, Cs) ) )
fast(V::OneBody) = V

Dict(V::OneBody) = Dict("__id__" => "OneBody", "E0" => V.E0)
OneBody(D::Dict) = OneBody(D["E0"])
convert(::Val{:OneBody}, D::Dict) = OneBody(D)

==(V1::OneBody, V2::OneBody) = (V1.E0 == V2.E0)
