using NBodyIPs
using JuLIP, Base.Test, StaticArrays, ForwardDiff, Combinatorics
using BenchmarkTools
using NBodyIPs.Invariants

using NBodyIPs.Invariants: invariants, grad_invariants
using JuLIP.Potentials: evaluate, evaluate_d

ad_invariants(inv, r) = ForwardDiff.jacobian(r_->invariants(inv, r_), r)
Inv = InvInvariants()
r = 0.5 + rand(SVector{3,Float64})

println("-------------------------------------------")
println("   Testing implementation of `invariants`")
println("-------------------------------------------")

println("[1] Quick profiling:")
print("     invariants: ")
@btime invariants($Inv, $r)
print("grad_inveriants: ")
@btime grad_invariants($Inv, $r)
print("  ad_inveriants: ")
@btime ad_invariants($Inv, $r)

println("[2] Correctness of gradients")
for n = 1:10
   r = 0.5 + rand(SVector{3,Float64})
   @test hcat(grad_invariants(Inv, r)...)' ≈ ad_invariants(Inv, r)
end

println("[3] Symmetry")
for n = 1:10
   r = 0.5 + rand(SVector{3,Float64})
   Q = invariants(Inv, r)
   for rπ in permutations(r)
      @test Q ≈ invariants(Inv, SVector{3}(rπ))
   end
end


println("----------------------------------------")
println("   Testing Implementation of `NBody{3}`")
println("----------------------------------------")

r0 = rnn(:Cu)
rcut = 3.1 * r0
D = Dictionary(InvInvariants, rcut)

n = 10
println("[1] Quick profiling for a 3-body with $n basis functions")
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
@show length(at)
r = 1.0 + rand(SVector{3, Float64})
V3 = NBody( [tuple(rand(0:4, 3)...) for n = 1:n], rand(n), D )
print("     V3: "); @btime evaluate($V3, $r)
print("  @D V3: "); @btime evaluate_d($V3, $r)
print("  nlist: "); @btime neighbourlist(at, rcut)
print(" energy: "); @btime energy(V3, at)
print(" forces: "); @btime forces(V3, at)

println("[2] finite-difference test on triangles")
for n = [1, 3]
   V3 = NBody( [tuple(rand(0:3, 3)...) for n = 1:n], rand(n), D )
   for _  = 1:10
      r = 1.0 + rand(SVector{3,Float64})
      @test (@D V3(r)) ≈ ForwardDiff.gradient(r_ -> V3(r_), r)
   end
end

println("[3] finite-difference test on configurations")
at = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
for n in [1, 3]
   V3 = NBody( [tuple(rand(0:3, 3)...) for n = 1:n], 0.01 * rand(n), D )
   @test JuLIP.Testing.fdtest(V3, at)
end




# using NBodyIPs
# using JuLIP, Base.Test, StaticArrays, ForwardDiff, Combinatorics
# using BenchmarkTools
# using NBodyIPs.Invariants
#
# using NBodyIPs.Invariants: invariants, grad_invariants
# using JuLIP.Potentials: evaluate, evaluate_d
#
# r0 = rnn(:Cu)
# rcut = 3.1 * r0
# D = Dictionary(InvInvariants, rcut)
# B = gen_basis(3, D, 12)
# at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
