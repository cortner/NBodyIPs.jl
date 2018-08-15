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
