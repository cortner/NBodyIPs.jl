
using StaticArrays, BenchmarkTools, Combinatorics

include("../data/NB_5_deg_6_non_efficient_invariants.jl")
include("../data/NB_5_deg_6_invariants.jl")
# include("invariants_co_new.jl")

include("fastpolys.jl")
using FastPolys

x = @SVector rand(10)

(Primary_slow, Sec_slow, Irr_sec_slow) = invariants_Q10_check(x)
(Primary_fast,Sec_fast) = invariants(x)

# ------------------
# Invariant check vs slow version
# ------------------

# Primary comparison
SVector(Primary_slow...) - Primary_fast
display(maximum(abs.(SVector(Primary_slow...) - Primary_fast)))
#dont match yet...dont know why.

# Irreducible secondary comparison
SVector(Sec_slow...) - Sec_fast
display(maximum(abs.(SVector(Sec_slow...) - Sec_fast)))

# # Secondary comparison
# SVector(Sec_Inv...) - Sec
# maximum(abs.(SVector(Sec_Inv...) - Sec))

# ------------------
# Timings
# ------------------

@btime invariants_Q10_check($x)
@btime invariants($x)
@btime invariants_d($x)
