using NBodyIPs, JuLIP, BenchmarkTools, StaticArrays
using Base.Test

using NBodyIPs: BondLengthDesc, transform, transform_d, fcut, fcut_d,
                invariants, invariants_d, invariants_ed
using NBodyIPs.Polys: NBPoly


nbasis = [0, 10, 10, 10, 10]
nbasis = [0, 10, 50, 100, 300]

println("Setting up the test systems ...")
r0 = rnn(:Cu)
TRANSFORM = "r -> exp( - 3 * ((r/$r0) - 1))"
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.05)
rcut3 = 3.1 * r0
D3 = BondLengthDesc(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
rcut4 = 2.1 * r0
D4 = BondLengthDesc(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
rcut5 = 1.5 * r0
D5 = BondLengthDesc(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )
DD = [nothing, D3, D3, D4, D5]

random_nbody(N, ntup) = (
   (N <= 3) ? NBPoly( [tuple( [rand(1:4, ((N*(N-1))÷2)); 0]... ) for n = 1:ntup],
                     (0.1+rand(ntup))/factorial(N), DD[N] )
           : NBPoly( [tuple( rand(0:N, ((N*(N-1))÷2+1))... ) for n = 1:ntup],
                   (0.1+rand(ntup))/factorial(N), DD[N] ) )

BB = [ nothing; [ [random_nbody(N, 1) for _=1:nbasis[N]] for N = 2:5 ] ]
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.05)

println("Performance Tests")
println("-----------------")
println("to get a rough estimate of the cost per atom take the `separate`")
println("timing and multiply by the factor provided")
for N = 2:5
   B = BB[N]
   println("[$N-body], length(B) == $(length(B)), factor = $(1/(length(at)*length(B)))")
   print(" E separate:" ); @btime ([energy(b, $at) for b in $B])
   print(" E combined:" ); @btime energy( $B, $at )
   print(" F separate:" ); @btime ([forces(b, $at) for b in $B])
   print(" F combined:" ); @btime forces( $B, $at )
   print(" S separate:" ); @btime ([virial(b, $at) for b in $B])
   print(" S combined:" ); @btime virial( $B, $at )
end

# B = BB[4]
#
# # @time [energy(b, at) for b in B]
# # @time [energy(b, at) for b in B]
# # @profile [energy(b, at) for b in B]
# # Base.Profile.print()
#
# @time [forces(b, at) for b in B];
# @time [forces(b, at) for b in B];
# @profile [forces(b, at) for b in B];
# Base.Profile.print()


# # THE CURRENT IMPLEMENTATION
# [2-body], length(B) == 10, factor = 0.003125
#  E separate:  26.563 ms (8693 allocations: 9.80 MiB)
#  E combined:  3.061 ms (804 allocations: 1002.56 KiB)
#  F separate:  44.766 ms (403043 allocations: 21.90 MiB)
#  F combined:  18.543 ms (395712 allocations: 13.08 MiB)
#  S separate:  30.726 ms (11363 allocations: 9.98 MiB)
#  S combined:  4.431 ms (5300 allocations: 1.22 MiB)
# [3-body], length(B) == 50, factor = 0.000625
#  E separate:  1.611 s (43403 allocations: 49.01 MiB)
#  E combined:  218.385 ms (804 allocations: 1002.94 KiB)
#  F separate:  2.663 s (2015203 allocations: 109.47 MiB)
#  F combined:  593.160 ms (1975312 allocations: 61.48 MiB)
#  S separate:  2.593 s (56803 allocations: 49.89 MiB)
#  S combined:  502.200 ms (23260 allocations: 2.19 MiB)
# [4-body], length(B) == 100, factor = 0.0003125
#  E separate:  4.788 s (86403 allocations: 34.84 MiB)
#  E combined:  312.957 ms (800 allocations: 356.28 KiB)
#  F separate:  7.037 s (1297203 allocations: 72.09 MiB)
#  F combined:  859.382 ms (1217008 allocations: 37.69 MiB)
#  S separate:  6.862 s (113203 allocations: 36.33 MiB)
#  S combined:  807.974 ms (45706 allocations: 2.50 MiB)
# [5-body], length(B) == 300, factor = 0.00010416666666666667
#  E separate:  2.188 s (259203 allocations: 78.31 MiB)
#  E combined:  28.564 ms (800 allocations: 268.41 KiB)
#  F separate:  2.774 s (1472403 allocations: 115.96 MiB)
#  F combined:  124.575 ms (1230208 allocations: 38.18 MiB)
#  S separate:  2.926 s (339604 allocations: 82.50 MiB)
#  S combined:  79.406 ms (135507 allocations: 6.44 MiB)
