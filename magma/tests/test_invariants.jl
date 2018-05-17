
using StaticArrays, BenchmarkTools, Combinatorics, Base.Test
using ForwardDiff
using NBodyIPs

Nbody = 5;
Deg = 6;

NBlengths = Int(Nbody*(Nbody-1)/2);

include("../../src/fastpolys.jl")
using FastPolys

include("../data/NB_$Nbody"*"_deg_$Deg/NB_$Nbody"*"_deg_$Deg"*"_non_efficient_invariants.jl")
include("../data/NB_$Nbody"*"_deg_$Deg/NB_$Nbody"*"_deg_$Deg"*"_invariants.jl")



x = @SVector rand(NBlengths)

@btime invariants_ed_gen($x)

(Primary_slow, Sec_slow, Irr_sec_slow) = invariants_Q10_check(x)
(Primary_fast,Sec_fast) = invariants_gen(x)
(Primary_fast2,Sec_fast2,Prim_d,Sec_d) = invariants_ed_gen(x)

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
@btime invariants_ed_gen($x)
@btime invariants_Q10_check($x)
@btime invariants_gen($x)
@btime invariants_d_gen($x)


# function d_invariants_Q6_check(x)
#    return ForwardDiff.gradient(invariants_Q6_check, x)
# end
#
# all_invariants(r) = vcat(invariants_gen(r)...)
# ad_invariants(r) = ForwardDiff.jacobian(all_invariants, r)
# r = 1.0 + SVector(rand(6)...)
# # dI1, dI2 = invariants_d(r)
# dI1, dI2 = ad_invariants(r)
# display(dI1)
# display(dI2)
# @test [dI1; dI2] ≈ ad_invariants(r)
# print(".")
#
#
#
# epsi = 1e-4;
# h = epsi*(@SVector rand(10))
#
#
