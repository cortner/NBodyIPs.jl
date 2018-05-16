
using StaticArrays, BenchmarkTools, Combinatorics, Base.Test
using ForwardDiff

include("../data/NB_5_deg_6_non_efficient_invariants.jl")
include("../data/NB_5_deg_6_invariants.jl")
# include("invariants_co_new.jl")

include("fastpolys.jl")
using FastPolys


function d_invariants_Q10_check(x)
   return ForwardDiff.gradient(invariants_Q10_check, x)
end

all_invariants(r) = vcat(invariants(r)...)
ad_invariants(r) = ForwardDiff.jacobian(all_invariants, r)
r = 1.0 + SVector(rand(10)...)
# dI1, dI2 = invariants_d(r)
dI1, dI2 = ad_invariants(r)
display(dI1)
display(dI2)
@test [dI1; dI2] ≈ ad_invariants(r)
print(".")


x = @SVector rand(10)
epsi = 1e-4;
h = epsi*(@SVector rand(10))

(Primary_slow, Sec_slow, Irr_sec_slow) = invariants_Q10_check(x)
(Primary_fast,Sec_fast) = invariants(x)
gradi = d_invariants_check(x)


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
