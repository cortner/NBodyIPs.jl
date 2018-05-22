using NBodyIPs, StaticArrays, JuLIP, BenchmarkTools

using NBodyIPs.Polys: invariants, invariants_d, invariants_ed
using JuLIP.Potentials: evaluate, evaluate_d

println("----------------------------------------")
println("  [1] Single Invariants Evaluation:")
println("----------------------------------------")
for r in [ (@SVector rand(3)),
           (@SVector rand(6)),
           (@SVector rand(10)) ]
   println("dim = $(length(r))")
   print("     invariants: ")
   @btime invariants($r)
   print("   invariants_d: ")
   @btime invariants_d($r)
   print("  invariants_ed: ")
   @btime invariants_ed($r)
end




println("----------------------------------------")
println("   [2] `NBody` Evaluation")
println("----------------------------------------")

r0 = rnn(:Cu)
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
@show length(at)

TRANSFORM = let r0 = r0
   (@analytic r -> (r0/r)^3)
end
rcut3 = 3.1 * r0
D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
rcut4 = 2.1 * r0
D4 = Dictionary(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
rcut5 = 1.5 * r0
D5 = Dictionary(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )

n = 10
println("[1] Quick profiling for a N-body with $n basis functions")
r = 1.0 + rand(SVector{3, Float64})
V3 = NBody( [tuple([rand(0:4, 3);0]...) for n = 1:n], 1.0+rand(n), D3 )
print("     V3: "); @btime evaluate($V3, $r)
print("  @D V3: "); @btime evaluate_d($V3, $r)
# print("  nlist: "); @btime neighbourlist($at, $rcut3)
# print(" energy: "); @btime energy($V3, $at)
# print(" forces: "); @btime forces($V3, $at)

function rand_tuple_4b()
   t = rand(0:5, 7)
   t[rand(1:6, 5)] = 0
   return tuple(t...)
end

r = 1.0 + rand(SVector{6, Float64})
V4 = NBody( [rand_tuple_4b() for _ = 1:n], 1.0+rand(n), D4 )
print("     V4: "); @btime evaluate($V4, $r)
print("  @D V4: "); @btime evaluate_d($V4, $r)
# print("  nlist: "); @btime neighbourlist(at, rcut4)
# print(" energy: "); @btime energy(V4, at)
# print(" forces: "); @btime forces(V4, at)

function rand_tuple_5b()
   t = rand(0:6, 11)
   t[rand(1:10, 9)] = 0
   return tuple(t...)
end

r = 1.0 + rand(SVector{10, Float64})
V5 = NBody( [rand_tuple_5b() for _ = 1:n], 1.0+rand(n), D5 )
print("     V5: "); @btime evaluate($V5, $r)
print("  @D V5: "); @btime evaluate_d($V5, $r)
# print("  nlist: "); @btime neighbourlist(at, rcut4)
# print(" energy: "); @btime energy(V4, at)
# print(" forces: "); @btime forces(V4, at)
