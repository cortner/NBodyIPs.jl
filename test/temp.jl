

using NBodyIPs, JuLIP, NeighbourLists, StaticArrays, BenchmarkTools
using NBodyIPs.Invariants

const r0 = rnn(:Cu)
const rcut3 = 2.1 * r0
const D3 = Dictionary(InvInvariants, rcut3)
const rcut4 = 2.1 * r0
const D4 = Dictionary(InvInvariants, rcut4)
const at = rattle!(bulk(:Cu, cubic=true) * (2,2,3), 0.02)
@show length(at)

# warmup
V3 = NBody( [tuple(rand(0:4, 3)...)], rand(1), D3 )
energy(V3, at)
forces(V3, at)
V4 = NBody( [tuple(rand(0:4, 7)...)], rand(1), D4 )
energy(V4, at)
forces(V4, at)

# this allocates no memory and is the fastest possible 4-body term
# immaginable
f4(r) = prod(r)
df4(r) = (rand(1); ForwardDiff.gradient(f4, r))
energy_f4(at) = maptosites!( f4, zeros(length(at)),
                            nbodies(4, neighbourlist(at, rcut4)) ) |> sum_kbn
forces_f4(at) = maptosites_d!( df4, zeros(SVector{3,Float64}, length(at)),
                              nbodies(4, neighbourlist(at, rcut4)) )
energy_f4(at)
forces_f4(at)


for n in [1, 10]
   println("n = $n")
   println("   V3")
   V3 = NBody( [tuple(rand(0:4, 3)...) for n = 1:n], rand(n), D3 )
   print("     energy: ")
   @time energy(V3, at)
   print("     forces: ")
   @time forces(V3, at)
   println("   V4")
   V4 = NBody( [tuple(rand(0:4, 7)...) for n = 1:n], rand(n), D4 )
   print("     energy: ")
   @time energy(V4, at)
   print("     forces: ")
   @time forces(V4, at)
   println("   Naive 4-body term")
   print("     energy: ")
   @time energy_f4(at)
   print("     forces: ")
   @time forces_f4(at)
end



# V4 = NBody( [tuple(rand(0:4, 7)...)], rand(1), D4 )
# forces(V4, at)
