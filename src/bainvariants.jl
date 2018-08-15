
module BAInvariants

using StaticArrays



# ------------------------------------------------------------------------
#             2-BODY Invariants
#             fully equivalent to BL invariants
# ------------------------------------------------------------------------

# x = (r12,)

invariants(x::SVector{1, T}) where {T} =
      copy(x),
      SVector{1, T}(1.0)

invariants_d(x::SVector{1, T}) where {T} =
      (@SVector [ SVector(one(T))  ]),
      (@SVector [ SVector(zero(T)) ])

invariants_ed(x::SVector{1,T}) where {T} =
      copy(x),
      SVector{1, T}(1.0),
      (@SVector [ SVector(one(T))  ]),
      (@SVector [ SVector(zero(T)) ])



# ------------------------------------------------------------------------
#             3-BODY Invariants
# ------------------------------------------------------------------------

# x = (r1, r2, θ12)

# the 1.0 is a "secondary invariant"
invariants(x::SVector{3, T}) where {T} =
      (@SVector T[ x[1] + x[2], x[1] * x[2], x[3] ]),
      (@SVector T[ 1.0 ])


invariants_d(x::SVector{3, T}) where {T} =
      (@SVector [ (@SVector T[1.0, 1.0, 0.0]),
                  (@SVector T[x[2], x[1], 0.0]),
                  (@SVector T[0.0, 0.0, 1.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])

invariants_ed(x::SVector{3, T}) where {T} =
      (@SVector T[ x[1] + x[2], x[1] * x[2], x[3] ]),
      (@SVector T[ 1.0 ]),
      (@SVector [ (@SVector T[1.0, 1.0, 0.0]),
                  (@SVector T[x[2], x[1], 0.0]),
                  (@SVector T[0.0, 0.0, 1.0]) ]),
      (@SVector [ (@SVector T[0.0, 0.0, 0.0]) ])



# ------------------------------------------------------------------------
#             4-BODY Invariants
# ------------------------------------------------------------------------

include("generated/BA_4B/BA_4B_invariants_hack.jl")

# ------------------------------------------------------------------------
#             5-BODY Invariants
# ------------------------------------------------------------------------

include("generated/BA_5B/BA_5B_invariants_hack.jl")

end
