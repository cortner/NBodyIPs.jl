using NBodyIPs
using JuLIP, Base.Test, StaticArrays, ForwardDiff, Combinatorics
using BenchmarkTools

using NBodyIPs.Polys: invariants, invariants_d, invariants_ed
using JuLIP.Potentials: evaluate, evaluate_d

all_invariants(r) = vcat(invariants(r)...)  # [I1; I2]
ad_invariants(r) = ForwardDiff.jacobian(all_invariants, r)

println("-------------------------------------------")
println("   Testing implementation of `invariants`")
println("-------------------------------------------")

println("[1] Quick profiling:")
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

# TODO: test correctness of the invariants implementation
#       against the MAGMA output

println("[2] Correctness of gradients")
for dim in [3, 6, 10]
   r = 1.0 + SVector(rand(dim)...)
   println("---------------")
   println("dim = $dim")
   println("---------------")

   I = all_invariants(r)
   dI1, dI2 = invariants_d(r)
   dI = [hcat(dI1...)'; hcat(dI2...)']
   dIh = zeros(size(dI))
   r0 = Vector(r)
   for p = 2:9
      h = .1^p
      dIh = zeros(size(dI))
      for j = 1:length(r)
         r0[j] += h
         Ih = all_invariants(SVector(r0...))
         dIh[:, j] = (Ih - I) / h
         r0[j] -= h
      end
      @printf(" %d | %.2e \n", p, vecnorm(dIh - dI, Inf))
   end
   println("---------------")
end
println()



println("[3] Symmetry")
for dim in [3, 6, 10]
   for n = 1:3
      r = 1.0 + SVector(rand(dim)...)
      I = all_invariants(r)
      for rπ in NBodyIPs.simplex_permutations(r)
         @test I ≈ all_invariants(SVector(rπ...))
      end
      print(".")
   end
end
println()



println("----------------------------------------")
println("   Testing Implementation of `NBody`")
println("----------------------------------------")

r0 = rnn(:Cu)
TRANSFORM = let r0 = r0
   (@analytic r -> (r0/r)^3)
end
rcut3 = 3.1 * r0
D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
rcut4 = 2.1 * r0
D4 = Dictionary(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )

n = 10
println("[1] Quick profiling for a N-body with $n basis functions")
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
@show length(at)
r = 1.0 + rand(SVector{3, Float64})
V3 = NBody( [tuple([rand(0:4, 3);0]...) for n = 1:n], 1.0+rand(n), D3 )
print("     V3: "); @btime evaluate($V3, $r)
print("  @D V3: "); @btime evaluate_d($V3, $r)
print("  nlist: "); @btime neighbourlist($at, $rcut3)
print(" energy: "); @btime energy($V3, $at)
print(" forces: "); @btime forces($V3, $at)

r = 1.0 + rand(SVector{6, Float64})
V4 = NBody( [tuple(rand(0:4, 7)...) for n = 1:n], 1.0+rand(n), D4 )
print("     V4: "); @btime evaluate($V4, $r)
print("  @D V4: "); @btime evaluate_d($V4, $r)
print("  nlist: "); @btime neighbourlist(at, rcut4)
print(" energy: "); @btime energy(V4, at)
print(" forces: "); @btime forces(V4, at)


println("[2] Gradient-test on simplices")
for n = [1, 3]
   V3 = NBody( [tuple([rand(0:3, 3); 0]...) for n = 1:n], 1.0 + rand(n), D3 )
   for _  = 1:10
      r = 1.0 + rand(SVector{3,Float64})
      @test (@D V3(r)) ≈ ForwardDiff.gradient(r_ -> V3(r_), r)
      print(".")
   end
end

for n = [1, 3]
   V4 = NBody( [tuple(rand(0:3, 7)...) for n = 1:n], 1.0 + rand(n), D4 )
   for _  = 1:10
      r = 1.0 + rand(SVector{6,Float64})
      @test evaluate_d(V4, r) ≈ ForwardDiff.gradient(r_ -> V4(r_), r)
      print(".")
   end
end
println()

println("[3] finite-difference test on configurations")
nb = 3
at1 = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
at2 = bulk(:Cu, cubic=true) * (1,1,2)
set_constraint!(at2, VariableCell(at2, free = []))
for at in [at1, at2]
   println("  3-body")
   V3 = NBody( [tuple([rand(0:5, 3); 0]...) for n = 1:nb], 1.0+rand(nb), D3 )
   @test JuLIP.Testing.fdtest(V3, at)

   println("  4-body")
   V4 = NBody( [tuple(rand(0:5, 7)...) for n = 1:nb], 1.0 + rand(nb), D4 )
   @test JuLIP.Testing.fdtest(V4, at)
end
