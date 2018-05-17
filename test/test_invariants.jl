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

# TODO: test correctness of the invariants implementation
#       against the MAGMA output

println("[1] Correctness of gradients")
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


println("[2] Symmetry")
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


println("[3] invariants_ed")
for dim in [3, 6, 10]
   for n = 1:3
      r = 1.0 + SVector(rand(dim)...)
      I1, I2 = invariants(r)
      dI1, dI2 = invariants_d(r)
      J1, J2, dJ1, dJ2 = invariants_ed(r)
      @test all(i ≈ j for (i,j) in zip(I1, J1))
      @test all(i ≈ j for (i,j) in zip(I2, J2))
      @test all(di ≈ dj for (di,dj) in zip(dI1, dJ1))
      @test all(di ≈ dj for (di,dj) in zip(dI2, dJ2))
      print(".")
   end
end
println()


# println("[4] `NBody` gradient-test on simplices")
# for n = [1, 3]
#    V3 = NBody( [tuple([rand(0:3, 3); 0]...) for n = 1:n], 1.0 + rand(n), D3 )
#    for _  = 1:10
#       r = 1.0 + rand(SVector{3,Float64})
#       @test (@D V3(r)) ≈ ForwardDiff.gradient(r_ -> V3(r_), r)
#       print(".")
#    end
# end
#
# for n = [1, 3]
#    V4 = NBody( [tuple(rand(0:3, 7)...) for n = 1:n], 1.0 + rand(n), D4 )
#    for _  = 1:10
#       r = 1.0 + rand(SVector{6,Float64})
#       @test evaluate_d(V4, r) ≈ ForwardDiff.gradient(r_ -> V4(r_), r)
#       print(".")
#    end
# end
# println()
#
# println("[5] `NBody` finite-difference test on configurations")
# nb = 3
# at1 = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
# at2 = bulk(:Cu, cubic=true) * (1,1,2)
# set_constraint!(at2, VariableCell(at2, free = []))
# for at in [at1, at2]
#    println("  3-body")
#    V3 = NBody( [tuple([rand(0:5, 3); 0]...) for n = 1:nb], 1.0+rand(nb), D3 )
#    @test JuLIP.Testing.fdtest(V3, at)
#
#    println("  4-body")
#    V4 = NBody( [tuple(rand(0:5, 7)...) for n = 1:nb], 1.0 + rand(nb), D4 )
#    @test JuLIP.Testing.fdtest(V4, at)
# end
