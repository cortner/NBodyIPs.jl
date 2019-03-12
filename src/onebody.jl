export OneBody

using JuLIP: Atoms
import JuLIP.Potentials: OneBody
import JuLIP: energy, forces, virial
import Base: hash
import NBodyIPs: fast, combinebasis

combinebasis(Vs::AbstractVector{<: OneBody}, Cs::AbstractVector{<: Real}) =
             OneBody( sum( V.E0 * c for (V, c) in zip(Vs, Cs) ) )

fast(V::OneBody) = V

hash(::BASIS, ::OneBody) = hash(OneBody)

bodyorder(V::OneBody) = 1 

energy(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: OneBody, T} =
   [ energy(b, at) for b in B ]

forces(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: OneBody, T} =
   [ forces(b, at) for b in B ]

virial(B::AbstractVector{TB}, at::Atoms{T}) where {TB <: OneBody, T} =
   [ virial(b, at) for b in B ]

# TODO: maybe the above should be moved to TB <: AbstractCalculator?
