
using StaticArrays, BenchmarkTools, Combinatorics

include("../data/NB_5_deg_6_non_efficient_invariants.jl")
include("../data/NB_5_deg_6_invariants.jl")

using FastPolys

x = @SVector rand(10)

(Primary_slow, Sec_slow, Irr_sec_slow) = invariants_Q10_check(x)
(Primary_fast,Irr_sec_fast) = invariants(x)

# ------------------
# Invariant check vs slow version
# ------------------

# Primary comparison
Primary - Primary_inv
maximum(abs.(SVector(Primary...) - Primary_inv))
#dont match yet...dont know why.

# Secondary comparison
SVector(Sec_Inv...) - Sec
maximum(abs.(SVector(Sec_Inv...) - Sec))

# ------------------
# Timings
# ------------------

@btime invariants_Q10_check($x)
@btime invariants_Q10($x)
