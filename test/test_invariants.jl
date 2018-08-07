using NBodyIPs
using JuLIP, Base.Test, StaticArrays, ForwardDiff, Combinatorics
using BenchmarkTools

using JuLIP.Potentials: evaluate, evaluate_d

const blinvariants = NBodyIPs.BLInvariants.invariants
const blinvariants_d = NBodyIPs.BLInvariants.invariants_d
const blinvariants_ed = NBodyIPs.BLInvariants.invariants_ed

include("aux_testing.jl")

all_blinvariants(r) = vcat(blinvariants(r)...)  # [I1; I2]
ad_blinvariants(r) = ForwardDiff.jacobian(all_blinvariants, r)

println("-------------------------------------------")
println("   Testing implementation of `blinvariants`")
println("-------------------------------------------")

# TODO: test correctness of the blinvariants implementation
#       against the MAGMA output

println("[1] Correctness of gradients")
for dim in [3, 6, 10]
   r = 1.0 + SVector(rand(dim)...)
   println("---------------")
   println("dim = $dim")
   println("---------------")

   I = all_blinvariants(r)
   dI1, dI2 = blinvariants_d(r)
   dI = [hcat(dI1...)'; hcat(dI2...)']
   dIh = zeros(size(dI))
   r0 = Vector(r)
   errs = []
   for p = 2:9
      h = .1^p
      dIh = zeros(size(dI))
      for j = 1:length(r)
         r0[j] += h
         Ih = all_blinvariants(SVector(r0...))
         dIh[:, j] = (Ih - I) / h
         r0[j] -= h
      end
      push!(errs, vecnorm(dIh - dI, Inf))
      @printf(" %d | %.2e \n", p, errs[end])
   end
   println("---------------")
   @test minimum(errs) <= 1e-3 * maximum(errs)
end
println()


println("[2] Symmetry")
for dim in [3, 6, 10]
   for n = 1:3
      r = 1.0 + SVector(rand(dim)...)
      I = all_blinvariants(r)
      for rπ in simplex_permutations(r)
         @test I ≈ all_blinvariants(SVector(rπ...))
      end
      print(".")
   end
end
println()


println("[3] blinvariants_ed")
for dim in [3, 6, 10]
   for n = 1:3
      r = 1.0 + SVector(rand(dim)...)
      I1, I2 = blinvariants(r)
      dI1, dI2 = blinvariants_d(r)
      J1, J2, dJ1, dJ2 = blinvariants_ed(r)
      @test all(i ≈ j for (i,j) in zip(I1, J1))
      @test all(i ≈ j for (i,j) in zip(I2, J2))
      @test all(di ≈ dj for (di,dj) in zip(dI1, dJ1))
      @test all(di ≈ dj for (di,dj) in zip(dI2, dJ2))
      print(".")
   end
end
println()
