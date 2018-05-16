
using StaticArrays, BenchmarkTools, Combinatorics

include("../data/NB_5_deg_6_non_efficient_invariants.jl")
include("../data/NB_5_deg_6_invariants.jl")
# include("invariants_co_new.jl")

using FastPolys

x = @SVector rand(10)

(Primary_slow, Sec_slow, Irr_sec_slow) = invariants_Q10_check(x)
(Primary_fast,Irr_sec_fast) = invariants(x)

# ------------------
# Invariant check vs slow version
# ------------------

# Primary comparison
SVector(Primary_slow...) - Primary_fast
maximum(abs.(SVector(Primary_slow...) - Primary_fast))
#dont match yet...dont know why.

# Irreducible secondary comparison
SVector(Irr_sec_slow...) - Irr_sec_fast
maximum(abs.(SVector(Irr_sec_slow...) - Irr_sec_fast))

# # Secondary comparison
# SVector(Sec_Inv...) - Sec
# maximum(abs.(SVector(Sec_Inv...) - Sec))

# ------------------
# Timings
# ------------------

@btime invariants_Q10_check($x)
@btime invariants($x)
