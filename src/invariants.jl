

"""
`invariants(r::SVector{M,T}) -> SVector` : computes the invariant descriptors as a function of the
lengths in a simplex. The order is lexicographical, i.e.,

* 2-body: `r::SVector{1}`
* 3-body: `r::SVector{3}`, order is irrelevant, but `r = [r12, r13, r23]`
* 4-body: `r::SVector{6}`, order is `r = [r12, r13, r14, r23, r24, r34]`
* n-body: `r::SVector{n*(n-1)/2}`, ...

The first `n*(n-1)/2` invariants are the primary invariants; the remainind
ones are the secondary invariants.
"""
function invariants end

"""
`invariants(r::SVector{M,T}) -> SMatrix` : computes the jacobian of
`invariants`
"""
function invariants_d end

"""
`degrees(::Val{N})` where `N` is the body-order returns a
tuple of polynomials degrees corresponding to the degrees of the
individual invariants.  E.g. for 3-body, the invariants are
r1 + r2 + r3, r1 r2 + r1 r3 + r2 r3, r1 r2 r3, and the corresponding
degrees are `(1, 2, 3)`.
"""
function degrees end

"""
`bo2edges(N)` : bodyorder-to-edges
"""
bo2edges(N::Integer) = (N * (N-1)) ÷ 2

"""
`edges2bo(M)`: "edges-to-bodyorder", an internal function that translates
the number of edges in a simplex into the body-order
"""
edges2bo(M::Integer) = round(Int, 0.5 + sqrt(0.25 + 2 * M))


# ------------------------------------------------------------------------
#             2-BODY Invariants
# ------------------------------------------------------------------------

invariants(r::SVector{1, T}) where {T} = copy(r)

invariants_d(r::SVector{1, T}) where {T} = @SMatrix [one(T)]

degrees(::Val{2}) = (1,)

# ------------------------------------------------------------------------
#             3-BODY Invariants
# ------------------------------------------------------------------------

# the 1.0 is a "secondary invariant"
invariants(r::SVector{3, T}) where {T} =
      @SVector T[ r[1]+r[2]+r[3],
                  r[1]*r[2] + r[1]*r[3] + r[2]*r[3],
                  r[1]*r[2]*r[3],
                  1.0 ]

invariants_d(  r::SVector{3, T}) where {T} =
      @SMatrix T[ 1.0        1.0         1.0;
                  r[2]+r[3]  r[1]+r[3]   r[1]+r[2];
                  r[2]*r[3]  r[1]*r[3]   r[1]*r[2];
                  0.0        0.0         0.0  ]

degrees(::Val{3}) = (1, 2, 3, 0)


# ------------------------------------------------------------------------
#             4-BODY Invariants
#
# this implementation is based on
#    Schmelzer, A., Murrell, J.N.: The general analytic expression for
#    S4-symmetry-invariant potential functions of tetra-atomic homonuclear
#    molecules. Int. J. Quantum Chem. 28, 287–295 (1985).
#    doi:10.1002/qua.560280210
#
# ------------------------------------------------------------------------
# TODO: reorder to obtain increasing degree?

degrees(::Val{4}) = (1, 2, 3, 4, 2, 3, 0, 3, 4, 5, 6, 9)

const _2 = 2.0^(-0.5)
const _3 = 3.0^(-0.5)
const _6 = 6.0^(-0.5)
const _12 = 12.0^(-0.5)

# permutation to account for the different ordering used here, vs Schmelzer et al.
const r2ρ = @SMatrix [   1 0 0 0 0 0
                         0 1 0 0 0 0
                         0 0 1 0 0 0
                         0 0 0 0 0 1
                         0 0 0 0 1 0
                         0 0 0 1 0 0 ]

const R2Q = @SMatrix [ _6     _6     _6    _6     _6     _6
                        _2      0      0   -_2      0      0
                         0     _2      0     0    -_2      0
                         0      0     _2     0      0    -_2
                         0    0.5   -0.5     0    0.5   -0.5
                        _3   -_12   -_12    _3   -_12   -_12 ]

const R2Qxr2ρ = R2Q * r2ρ

@inline invariants(r::SVector{6, T}) where {T} = _invariants_Q6(R2Qxr2ρ * r)

@inline invariants_d(r::SVector{6, T}) where {T} =
            _invariants_Q6_d(R2Qxr2ρ * r) * R2Qxr2ρ

function _invariants_Q6(Q::SVector{6, T}) where {T}
   Q2 = Q .* Q
   Q2_34, Q2_24, Q2_23 = Q2[3] * Q2[4], Q2[2] * Q2[4], Q2[2] * Q2[3]
   rt3 = sqrt(3.0)
   Q_56 = Q[5] * Q[6]

   return SVector{12, T}(
      # ---------------------------- primary invariants
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
      # ---------------------------- secondary invariants
      # sneak in an additional secondary "invariant"
      1.0,
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


_invariants_Q6_d(Q) = ForwardDiff.jacobian(_invariants_Q6, Q)


# ------------------------------------------------------------------------
#             5-BODY Invariants
# ------------------------------------------------------------------------

# TODO
