

"""
`invariants(r::SVector{M,T})` : computes the invariant descriptors as a function of the
lengths in a simplex. The order is lexicographical, i.e.,

* 2-body: `r::SVector{1}`
* 3-body: `r::SVector{3}`, order is irrelevant, but `r = [r12, r13, r23]`
* 4-body: `r::SVector{6}`, order is `r = [r12, r13, r14, r23, r24, r34]`
* n-body: `r::SVector{n*(n-1)/2}`, ...
"""
function invariants end


"""
Use polynomials in r as the invariants
"""
struct PolyInvariants
end

invariants(::PolyInvariants, r::SVector{3, T}) where {T} = SVector{3, T}(
   r[1]+r[2]+r[3],
   r[1]*r[2] + r[1]*r[3] + r[2]*r[3],
   r[1]*r[2]*r[3] )

grad_invariants(::PolyInvariants, r::SVector{3, T}) where {T} = SVector{3, SVector{3,T}}(
      SVector{3,T}(1,1,1),
      SVector{3,T}(r[2]+r[3], r[1]+r[3], r[1]+r[2]),
      SVector{3,T}(r[2]*r[3], r[1]*r[3], r[1]*r[2])    )


"""
Use polynomials in r^{-1} as the invariants
"""
struct InvInvariants
end

invariants(::InvInvariants, r::SVector{M, T}) where {M, T} = invariants_inv(1 ./ r)
grad_invariants(::InvInvariants, r::SVector{M, T}) where {M, T} = grad_invariants_inv(1 ./ r)

invariants_inv(s::SVector{3, T}) where {T} = SVector{3, T}(
   s[1]+s[2]+s[3],
   s[1]*s[2] + s[1]*s[3] + s[2]*s[3],
   s[1]*s[2]*s[3] )

function grad_invariants_inv(s::SVector{3, T}) where {T}
   t = - s.^2
   return SVector{3, SVector{3,T}}(
      SVector{3,T}(t[1],t[2],t[3]),
      SVector{3,T}(t[1]*(s[2]+s[3]), t[2]*(s[1]+s[3]), t[3]*(s[1]+s[2])),
      SVector{3,T}(t[1]*s[2]*s[3], t[2]*s[1]*s[3], t[3]*s[1]*s[2])    )
end

"""
`inv_degrees(::Val{N})` where `N` is the body-order returns a
tuple of polynomials degrees corresponding to the degrees of the
individual invariants.  E.g. for 3-body, the invariants are
r1 + r2 + r3, r1 r2 + r1 r3 + r2 r3, r1 r2 r3, and the corresponding
degrees are `(1, 2, 3)`.
"""
inv_degrees(::Val{3}) = (1, 2, 3)

# TODO
# struct ExpInvariants
# end

# -------------- 4-body invariants ----------------

# TODO: reorder to obtain increasing degree?
inv_degrees(::Val{4}) = (1, 2, 3, 4, 2, 3, 3, 4, 5, 6, 9)

const _2 = 2.0^(-0.5)
const _3 = 3.0^(-0.5)
const _6 = 6.0^(-0.5)
const _12 = 12.0^(-0.5)
#                     ρ1=r1  ρ2=r2  ρ3=r3  ρ4=    ρ5=r5  ρ6=
#
const r2ρ = @SMatrix [   1 0 0 0 0 0
                         0 1 0 0 0 0
                         0 0 1 0 0 0
                         0 0 0 0 0 1
                         0 0 0 0 1 0
                         0 0 0 1 0 0 ]  # to account for the different ordering

const R2Q = @SMatrix [ _6     _6     _6    _6     _6     _6
                        _2      0      0   -_2      0      0
                         0     _2      0     0    -_2      0
                         0      0     _2     0      0    -_2
                         0    0.5   -0.5     0    0.5   -0.5
                        _3   -_12   -_12    _3   -_12   -_12 ]

const R2Qxr2ρ = R2Q * r2ρ

function invariants_inv(s::SVector{6, T}) where {T}
   Q = R2Qxr2ρ * s
   Q2 = Q .* Q
   Q2_34, Q2_24, Q2_23 = Q2[3] * Q2[4], Q2[2] * Q2[4], Q2[2] * Q2[3]
   rt3 = sqrt(3.0)
   Q_56 = Q[5] * Q[6]

   return SVector{11, T}(
      # I1
      (Q[1]),
      # I2
      (Q2[2] + Q2[3] + Q2[4]),
      # I3
      (Q[2] * Q[3] * Q[4]),
      # I4
      (Q2_34 + Q2_24 + Q2_23),
      # I5
      (Q2[5] + Q2[6]),
      # I6
      (Q[6] * (Q2[6] - 3*Q2[5])),
   # ),
   # SVector{5, T}(
      # I7
      (Q[6] * (2*Q2[2] - Q2[3] - Q2[4]) + rt3 * Q[5] * (Q2[3] - Q2[4])),
      # I8
      (( (Q2[6] - Q2[5]) * (2*Q2[2] - Q2[3] - Q2[4])
            - 2 * rt3 * Q_56 * (Q2[3] - Q2[4]) )),
      # I9
      (( Q[6] * (2*Q2_34 - Q2_24 - Q2_23) + rt3 * Q[5] * (Q2_24 - Q2_23) )),
      # I10
      (( (Q2[6] - Q2[5])*(2*Q2_34 - Q2_24 - Q2_23)
                   - 2 * rt3 * Q_56 * (Q2_24 - Q2_23) )),
      # I11
      (( (Q2[3] - Q2[4]) * (Q2[4] - Q2[2]) * (Q2[2] - Q2[3]) *
            Q[5] * (3*Q2[6] - Q2[5]) ))
   )
end
