
using NBodyIPs.FastPolys


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
TODO
"""
function invariants_ed end

"""
`tdegrees(::Val{N})` where `N` is the body-order returns a
tuple of polynomial <vector degrees> corresponding to the degrees of the
individual invariants.  E.g. for 3-body, the primary invariants are
r1 + r2 + r3, r1 r2 + r1 r3 + r2 r3, r1 r2 r3, while there is only a single
secondary invariant 1.0. The corresponding
degrees are `( (1,0,0), (1,1,0), (1,1,1) ), ( (0,0,0), )`. Note that only the
lexicographically leading term is included. The remaining terms are
obtained from permutation invariance.
"""
function degrees end

"""
need documentation here
"""
function corners end

"""
`tdegrees(::Val{N})` where `N` is the body-order returns a
tuple of polynomial <total degrees> corresponding to the degrees of the
individual invariants.  E.g. for 3-body, the invariants are
r1 + r2 + r3, r1 r2 + r1 r3 + r2 r3, r1 r2 r3, and the corresponding
degrees are `(1, 2, 3)`.
"""
function tdegrees(v)
   deg1, deg2 = degrees(v)
   return sum.(deg1), sum.(deg2)
end

"""
`bo2edges(N)` : bodyorder-to-edges
"""
bo2edges(N::Integer) = (N * (N-1)) ÷ 2

"""
`edges2bo(M)`: "edges-to-bodyorder", an internal function that translates
the number of edges in a simplex into the body-order
"""
edges2bo(M::Integer) = (M <= 0) ? 1 : round(Int, 0.5 + sqrt(0.25 + 2 * M))


# ------------------------------------------------------------------------
#             2-BODY Invariants
# ------------------------------------------------------------------------

# r = (r12,)

invariants(r::SVector{1, T}) where {T} =
      copy(r),
      SVector{1, T}(1.0)

invariants_d(r::SVector{1, T}) where {T} =
      (@SVector [ SVector(one(T))  ]),
      (@SVector [ SVector(zero(T)) ])

invariants_ed(r::SVector{1,T}) where {T} =
      copy(r),
      SVector{1, T}(1.0),
      (@SVector [ SVector(one(T))  ]),
      (@SVector [ SVector(zero(T)) ])


# tdegrees(::Val{2}) = (1,), (0,)

degrees(::Val{2}) = (SVector(1,),), (SVector(0,),)

corners(::Val{2}) = ( SVector(1,2), )

# ------------------------------------------------------------------------
#             3-BODY Invariants
# ------------------------------------------------------------------------

# r = (r12, r13, r23)

# the 1.0 is a "secondary invariant"
invariants(r::SVector{3, T}) where {T} =
      (@SVector T[ r[1]+r[2]+r[3],
                   r[1]*r[2] + r[1]*r[3] + r[2]*r[3],
                   r[1]*r[2]*r[3] ]),
      (@SVector T[ 1.0 ])


invariants_d(  r::SVector{3, T}) where {T} =
      (@SVector [ (@SVector [1.0, 1.0, 1.0]),
                  (@SVector [r[2]+r[3], r[1]+r[3], r[1]+r[2]]),
                  (@SVector [r[2] * r[3], r[1] * r[3], r[1] * r[2]]) ]),
      (@SVector [ (@SVector [0.0, 0.0, 0.0]) ])

function invariants_ed(r::SVector{3, T}) where {T}
   r1 = r[1]
   r2 = r[2]
   r3 = r[3]
   r12 = r[1]*r[2]
   r13 = r[1]*r[3]
   r23 = r[2]*r[3]
   r123 = r12*r3
   return (@SVector T[ r1+r2+r3, r12 + r13 + r23, r123 ]),
          (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector [1.0, 1.0, 1.0]),
                  (@SVector [r2+r3, r1+r3, r1+r2]),
                  (@SVector [r23, r13, r12]) ]),
      (@SVector [ (@SVector [0.0, 0.0, 0.0]) ] )
end

# tdegrees(::Val{3}) = (1, 2, 3), (0,)

degrees(::Val{3}) = ( SVector(1,0,0), SVector(1,1,0), SVector(1,1,1) ),
                    ( SVector(0,0,0), )

corners(::Val{3}) = ( SVector(1,2), SVector(1,3), SVector(2,3) )

# ------------------------------------------------------------------------
#             4-BODY Invariants
# ------------------------------------------------------------------------

# tdegrees(::Val{4}) = (1, 2, 2, 3, 3, 4), (0, 3, 4, 5, 6, 9)

degrees(::Val{4}) = ( SVector(1,0,0,0,0,0),
                      SVector(1,0,0,0,0,1),
                      SVector(2,0,0,0,0,0),
                      SVector(1,1,1,0,0,0),
                      SVector(3,0,0,0,0,0),
                      SVector(4,0,0,0,0,0) ),
                    ( SVector(0,0,0,0,0,0),
                      SVector(2,1,0,0,0,0),
                      SVector(3,1,0,0,0,0),
                      SVector(5,0,0,0,0,0),
                      SVector(4,2,0,0,0,0),
                      SVector(8,1,0,0,0,0) )

corners(::Val{4}) = ( SVector(1,2), SVector(1,3), SVector(1,4), SVector(2,3), SVector(2,4), SVector(3,4) )

const A = @SMatrix [0 1 1 1 1 0
                    1 0 1 1 0 1
                    1 1 0 0 1 1
                    1 1 0 0 1 1
                    1 0 1 1 0 1
                    0 1 1 1 1 0]

const P42 = Val(( (1,2,3), (6,5,4) ))

function invariants(x::SVector{6, T}) where {T}
   x2 = x .*x
   x3 = x2.*x
   x4 = x3.*x

   I1 = sum(x)
   I2 = fpoly( (x,x), P42 )  # x[1]*x[6] + x[2]*x[5] + x[3]*x[4]
   I3 = sum(x2)
   I4 = x[1]*x[2]*x[3] + x[1]*x[4]*x[5] + x[2]*x[4]*x[6] + x[3]*x[5]*x[6]
   I5 = sum(x3)
   I6 = sum(x4)

   Ax = A*x
   PV1 = dot(x2, Ax)
   PV2 = dot(x3, Ax)
   PV3 = dot(x4, x)
   I11 = PV1 * PV1
   I12 = PV2 * PV3

   return SVector(I1, I2, I3, I4, I5, I6),
          SVector(one(T), PV1, PV2, PV3, I11, I12)
end

function invariants_d(x::SVector{6, T}) where {T}
   x2 = x.*x
   x3 = x2.*x
   x4 = x3.*x
   o = @SVector ones(T, 6)
   z = @SVector zeros(T, 6)
   # ∇I2 = @SVector [x[6], x[5], x[4], x[3], x[2], x[1]]
   ∇I2 = fpoly_d( (x, x), (o, o), P42 )
   ∇I4 = @SVector [ x[2]*x[3]+x[4]*x[5],
                    x[1]*x[3]+x[4]*x[6],
                    x[1]*x[2]+x[5]*x[6],
                    x[1]*x[5]+x[2]*x[6],
                    x[1]*x[4]+x[3]*x[6],
                    x[2]*x[4]+x[3]*x[5] ]
   Ax = A*x
   PV1 = dot(x2, Ax)
   PV2 = dot(x3, Ax)
   PV3 = dot(x4, x)
   ∇PV1 = A * x2 + 2 * (x .* Ax)
   ∇PV2 = A * x3 + 3 * (x2 .* Ax)
   ∇PV3 = 5 * x4

   return SVector(o, ∇I2, 2*x, ∇I4, 3*x2, 4*x3),
          SVector(z, ∇PV1, ∇PV2, ∇PV3, 2*PV1*∇PV1, PV3*∇PV2 + PV2*∇PV3)
end


function invariants_ed(x::SVector{6, T}) where {T}
   x2 = x.*x
   x3 = x2.*x
   x4 = x3.*x
   o = @SVector ones(T, 6)
   z = @SVector zeros(T, 6)

   # primary invariants
   I1 = sum(x)
   I2 = x[1]*x[6] + x[2]*x[5] + x[3]*x[4]
   I3 = sum(x2)
   x23 = x[2]*x[3]
   x45 = x[4]*x[5]
   x24 = x[2]*x[4]
   x35 = x[3]*x[5]
   I4 = x[1]*(x23 + x45) + (x24+x35)*x[6]
   I5 = sum(x3)
   I6 = sum(x4)

   # secondary invariants
   Ax = A*x
   PV1 = dot(x2, Ax)
   PV2 = dot(x3, Ax)
   PV3 = dot(x4, x)
   I11 = PV1 * PV1
   I12 = PV2 * PV3

   # gradient precomputations
   ∇I2 = @SVector [x[6], x[5], x[4], x[3], x[2], x[1]]
   ∇I4 = @SVector [ x23+x45,
                    x[1]*x[3]+x[4]*x[6],
                    x[1]*x[2]+x[5]*x[6],
                    x[1]*x[5]+x[2]*x[6],
                    x[1]*x[4]+x[3]*x[6],
                    x24+x35 ]
   ∇PV1 = A * x2 + 2 * (x .* Ax)
   ∇PV2 = A * x3 + 3 * (x2 .* Ax)
   ∇PV3 = 5 * x4

   return SVector(I1, I2, I3, I4, I5, I6),
         SVector(one(T), PV1, PV2, PV3, I11, I12),
         SVector(o, ∇I2, 2*x, ∇I4, 3*x2, 4*x3),
         SVector(z, ∇PV1, ∇PV2, ∇PV3, 2*PV1*∇PV1, PV3*∇PV2 + PV2*∇PV3)
end


# ------------------------------------------------------------------------
#             5-BODY Invariants
# ------------------------------------------------------------------------

#Invariants up to degree 7 only.
# COPIED FROM SLACK:
# For the 5-body, if I am right, the number of primary invariants is 10 with
# degrees : [ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6 ]. The number of secondary for each
# degree (starting from 0) is [ 1, 0, 0, 2, 5, 8, 15, 23, 33, 46 ] which means
# in total 133 secondary of degree less than 9. That’s quite a lot. Some of
# them have very long expressions.

# the tdegrees are now directly in the invariant file
# tdegrees(::Val{5}) =
#       (@SVector [ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6 ]),
#       (@SVector [ 0,
#                   3, 3,
#                   4, 4, 4, 4, 4,
#                   5, 5, 5, 5, 5, 5, 5, 5,
#                   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
#                   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 ])

# tdegrees(::Val{5}) =
#           (@SVector [1, 2, 2, 3, 3, 4, 4, 5, 5, 6, ]),
#     (@SVector [0, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, ])

# tdegrees(::Val{5}) =
#           (@SVector [1, 2, 2, 3, 3, 4, 4, 5, 5, 6, ]),
#     (@SVector [0, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, ])

tdegrees(::Val{5}) =
          (@SVector [1, 2, 2, 3, 3, 4, 4, 5, 5, 6, ]),
    (@SVector [0, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, ])

include("../magma/data/NB_5_deg_12/NB_5_deg_12_invariants.jl")


@inline invariants(x::SVector{10}) = NB5I.invariants_gen(x)
@inline invariants_d(x::SVector{10}) = NB5I.invariants_d_gen(x)
@inline invariants_ed(x::SVector{10}) = NB5I.invariants_ed_gen(x)


# ------------------------------------------------------------------------
#      Utility Functions
# ------------------------------------------------------------------------


# function ispure(vN::Val{N}, t::NTuple{K, Int}) where {N, K}
#    deg1, deg2 = degrees(vN)
#    d = sum( deg1[n] * t[n] for n = 1:K-1 ) + deg2[t[end]+1]
#    contains_x = zeros(Int, N)
#    c = corners(vN)
#    for n = 1:length(d)
#       if d[n] != 0
#          contains_x[c[n]] += 1
#       end
#    end
#    if length(find(contains_x)) == N
#       return true
#    end
#    return false
# end
