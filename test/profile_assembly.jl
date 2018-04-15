# TODO:
#  - rewrite for Invariants
#  - move test_invariants profiling to here

using JuLIP, NBodyIPs, PyCall, ProgressMeter, ASE
using NeighbourLists, BenchmarkTools, PyCall


function load_data(Nconfig = 411)   # (how can I load 411?)
   fname = "~/Dropbox/PIBmat/Ti_DFTB_Data/Ti_N54_T2000.xyz"
   @pyimport ase.io as ase_io
   data = Tuple{Atoms{Float64, Int}, Float64, JVecsF}[]
   @showprogress 0.1 "Loading Ti data ..." for n = 1:Nconfig
      atpy = try
         ase_io.read(fname, n-1)
      catch
         warn("Nconfig is apparently too large")
         break
      end
      E = atpy[:get_potential_energy]()
      F = atpy[:get_array]("force")' |> vecs
      at = Atoms(ASEAtoms(atpy))
      push!(data, (at, E, F))
   end
   return data
end

data = load_data(10)

RCUT = [4.1, 3.1, 2.1] * rnn(:Ti)
tibasis(ndict, RCUT) = get_basis(:inv2, [ndict+4, ndict+2, ndict], RCUT)

BB = tibasis(8, RCUT)
B = vcat(BB...)
c = rand(length(B))
IP = NBodyIP(B, c)

@time energy(IP, at)

D = dict(:inv2, 10, 2.1 * rnn(:Cu))
ex, f, df = parse([(3,3,3,3,3,3)], D...)
b0 = NBody(4, f, df, 2.1 * rnn(:Cu))

@btime b0(at);
@btime @D b0(at);

@time begin
   for b in vcat(BB...)
      b(at)
   end
end

at = data[1][1]
@show length(at)

for b in BB[1]   # 2-body
   @time (b(at), @D b(at))
end

for b in BB[2]   # 3-body
   @time (b(at), @D b(at))
end

@time begin for b in BB[3]   # 4-body
   b(at)
end
end


print("N-body energy and forces benchmark: ")
println("# Threads = ", Base.Threads.nthreads())

at = bulk(:Cu, cubic=true) * 5
println("bulk :Cu with nat = $(length(at))")
r0 = rnn(:Cu)
rcut = 2.1 * r0
X = positions(at)
C = cell(at)
f, f_d = gen_fnbody(rcut, r0)
nlist = PairList(X, rcut, C, (false, false, false), sorted = true)

for M in [2, 3, 4, 5]
   println("M = $M")
   print("  energy: ")
   @btime n_body($f, $M, $nlist)
   print("  forces: ")
   @btime grad_n_body($f_d, $M, $nlist)
end


lj = LennardJones(r0, 1.0) * C1Shift(cutoff)
@btime forces($lj, $at)




# -------- a profile code from a temp.jl

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
