using NBodyIPs
using JuLIP, Base.Test, StaticArrays, ForwardDiff, Combinatorics
using BenchmarkTools

using NBodyIPs.Polynomials: invariants, invariants_d
using JuLIP.Potentials: evaluate, evaluate_d

ad_invariants(r) = ForwardDiff.jacobian(invariants, r)
r = 0.5 + rand(SVector{3,Float64})

println("-------------------------------------------")
println("   Testing implementation of `invariants`")
println("-------------------------------------------")

println("[1] Quick profiling:")
for r in [(@SVector rand(3)), (@SVector rand(6))]
   println("dim = $(length(r))")
   print("     invariants: ")
   @btime invariants($r)
   print("   invariants_d: ")
   @btime invariants_d($r)
   print("  ad_inveriants: ")
   @btime ad_invariants($r)
end

println("[2] Correctness of gradients")
for n = 1:10
   r = 0.5 + rand(SVector{3,Float64})
   @test invariants_d(r) ≈ ad_invariants(r)
   print(".")
end
# for n = 1:10
#    r = 0.5 + rand(SVector{6,Float64})
#    @test hcat(grad_invariants(Inv, r)...)' ≈ ad_invariants(Inv, r)
# end
println()


println("[3] Symmetry")
for n = 1:10
   r = 0.5 + (@SVector rand(3))
   I = invariants(r)
   for rπ in NBodyIPs.simplex_permutations(r)
      @test I ≈ invariants(SVector{3}(rπ))
   end
   print(".")
end
for n = 1:10
   r = 1.0 + (@SVector rand(6))
   I = invariants(r)
   for rπ in NBodyIPs.simplex_permutations(r)
      Iπ = invariants(SVector{6}(rπ))
      @test I ≈ Iπ
   end
   print(".")
end
println()


println("----------------------------------------")
println("   Testing Implementation of `NBody{3}`")
println("----------------------------------------")

r0 = rnn(:Cu)
rcut3 = 3.1 * r0
D3 = Dictionary((:poly, -1), (:cos, 0.66*rcut3, rcut3))
rcut4 = 2.1 * r0
D4 = Dictionary( (:poly, -1), (:cos, 0.66*rcut4, rcut4))

n = 10
println("[1] Quick profiling for a 3-body with $n basis functions")
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
@show length(at)
r = 1.0 + rand(SVector{3, Float64})
V3 = NBody( [tuple(rand(0:4, 3)...) for n = 1:n], rand(n), D3 )
print("     V3: "); @btime evaluate($V3, $r)
print("  @D V3: "); @btime evaluate_d($V3, $r)
print("  nlist: "); @btime neighbourlist(at, rcut3)
print(" energy: "); @btime energy(V3, at)
print(" forces: "); @btime forces(V3, at)

r = 1.0 + rand(SVector{6, Float64})
V4 = NBody( [tuple(rand(0:4, 7)...) for n = 1:n], rand(n), D4 )
print("     V4: "); @btime evaluate($V4, $r)
print("  @D V4: "); @btime evaluate_d($V4, $r)
print("  nlist: "); @btime neighbourlist(at, rcut4)
print(" energy: "); @btime energy(V4, at)
print(" forces: "); @btime forces(V4, at)

quit()

println("[2] Gradient- test on triangles")
for n = [1, 3]
   V3 = NBody( [tuple(rand(0:3, 3)...) for n = 1:n], rand(n), D3 )
   for _  = 1:10
      r = 1.0 + rand(SVector{3,Float64})
      @test (@D V3(r)) ≈ ForwardDiff.gradient(r_ -> V3(r_), r)
   end
end

for n = [1, 3]
   V4 = NBody( [tuple(rand(0:3, 7)...) for n = 1:n], rand(n), D4 )
   for _  = 1:10
      r = 1.0 + rand(SVector{6,Float64})
      @test evaluate_d(V4, r) ≈ ForwardDiff.gradient(r_ -> V4(r_), r)
   end
end


println("[3] finite-difference test on configurations")
at = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
println("  3-body")
for n in [1, 3]
   V3 = NBody( [tuple(rand(0:3, 3)...) for n = 1:n], 0.01 * rand(n), D3 )
   @test JuLIP.Testing.fdtest(V3, at)
end
println("  4-body")
for n in [1, 3]
   V4 = NBody( [tuple(rand(0:3, 7)...) for n = 1:n], rand(n), D4 )
   @test JuLIP.Testing.fdtest(V4, at)
end

# r = 2.5 + (@SVector rand(6))
# V4(r)
# invariants(D4, r)
#


# using NBodyIPs
# using JuLIP, Base.Test, StaticArrays, ForwardDiff, Combinatorics
# using BenchmarkTools
# using NBodyIPs.Invariants
#
# using NBodyIPs.Invariants: invariants, grad_invariants
# using JuLIP.Potentials: evaluate, evaluate_d
#
# ad_invariants(inv, r) = ForwardDiff.jacobian(r_->invariants(inv, r_), r)
# Inv = InvInvariants()
#
# r = 1.0 + (@SVector rand(6))
# @btime invariants($Inv, $r)
# @btime ad_invariants($Inv, $r)
#
# const temp =  @MMatrix zeros(11, 6)
# const cfg = ForwardDiff.JacobianConfig(Invariants.invariants_inv, r,
#                ForwardDiff.Chunk{6}())
# dinv(r) = ForwardDiff.jacobian!(temp, Invariants.invariants_inv, r, cfg)
#
# @btime Invariants.invariants_inv($r)
# @btime dinv($r)
