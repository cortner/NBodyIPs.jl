
using JuLIP, NBodyIPs, PyCall, ProgressMeter, ASE

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

data = load_data(400)
train_data = data[1:300]
test_data = data[301:400]

# (3.5 x r0 is the cell dimension)
r0 = rnn(:Ti)

# some notes on the orders of magnitude
# E per atom ~ 6.0 eV, std ca 0.03
# mean force on each atom is ~ 1 eV / A
# but range of forces is between 1 and 12 eV / A

# construct a sequence of bases
# -----------------------------

println("generating basis functions")

BASES = []

B1 = [NBody(1.0)]

rcut2 = 9.2
rcut3 = 6.5

TRANSFORM = (@analytic r -> (2.9/r)^3)
CUTOFF2 = (:cos, 0.66*rcut2, rcut2)

D2 = Dictionary(TRANSFORM, CUTOFF2)
B2 = gen_basis(2, D2, 14)

for deg in [10,12,14,16]
   B = [B1; gen_basis(2, D2, deg)]
   push!(BASES, (B, D2, "2 / $(length(B)) / $rcut2"))
end

println("   first a few 3-body bases ...")
for (deg, rcut) in zip( [6, 8, 10],
                        [rcut3, rcut3, rcut3, rcut3] )
   CUTOFF3 = (:cos, 0.66*rcut, rcut)
   D = Dictionary(TRANSFORM, CUTOFF3)
   B3 = gen_basis(3, D, deg)
   B = [B1; B2; B3]
   push!(BASES, (B, D, "2+3 / $(length(B)) / $rcut2+$rcut"))
end

D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3))
B3 = gen_basis(3, D3, 8)
for (deg, rcut) in zip([4, 6, 8],
                       [4.5, 4.5, 4.5, 4.5])
   CUTOFF4 = (:cos, 0.66*rcut, rcut)
   D = Dictionary(TRANSFORM, CUTOFF4)
   B4 = gen_basis(4, D, deg)
   B = [B1; B2; B3; B4]
   push!(BASES, (B, D, "2+3+4 / $(length(B)) / $rcut2...$rcut") )
end

rmsE = Float64[]
rmsF = Float64[]
maeE = Float64[]
maeF = Float64[]

length_bases = [length(B[1]) for B in BASES]
@show length_bases

println("For each generated basis test the fit: " )
for (B, D, description) in BASES
   println(" - ", description)
   # standard least squares (see NBodyIPs/src/fitting.jl)
   # nforces = number of (randomly chosen) forces per configuration added
   #           to the LSQ problem
   ndata = min(length(train_data), 2 * length(B))
   c = regression(B, train_data[1:ndata], nforces = 50, stab = 0.0)
   @show norm(c, Inf)
   # construct an IP from the the basis and the weights
   IP = NBodyIP(B, c)
   # check error => the normalisation is w.r.t. natoms, not a genuine
   # relative error; we can discuss
   rE, rF = rms(IP, test_data)
   mE, mF = mae(IP, test_data)
   push!(rmsE, rE); push!(rmsF, rF)
   push!(maeE, mE); push!(maeF, mF)
   println("   E-rms, E-mae on testset = ", rE, ", ", mE)
   println("   F-rms, F-mae on testset = ", rF, ", ", mF)
end

using JLD
JLD.save("backup.jld", "desc", [B[3] for B in BASES],
         "rmsE", rmsE, "rmsF", rmsF, "maeE", maeE, "maeF", maeF)

using DataFrames
df = DataFrame(:desc => [B[3] for B in BASES])
df[Symbol("rms-E")] = rmsE
df[Symbol("rms-F")] = rmsF
df[Symbol("mae-E")] = maeE
df[Symbol("mae-F")] = maeF
println("Energy and Force Errors for Ti-DFTB Database: ")
println(df)






# using BenchmarkTools
#
# data = load_data(1)
# r0 = rnn(:Ti)
# at = data[1][1]
#
# D3 = Dictionary(InvInvariants,5.0)
# B3 = gen_basis(3, D3, 10)
# c3 = rand(length(B3))
# V3 = NBody(B3, c3, D3)
# @btime energy($V3, $at)
# @btime forces($V3, $at)
#
# D4 = Dictionary(InvInvariants, 5.0)
# B4 = gen_basis(4, D4, 6)
# c4 = rand(length(B4))
# V4 = NBody(B4, c4, D4)
# @btime energy($V4, $at)
# @btime forces($V4, $at)



# B, D, desc = BASES[end]
# at = data[1][1]
# B
# Es = [energy(b, at) for b in B]
# forces(B[10], at)
# IP = NBodyIP(B)
# IP
# energy(IP, at)
#
#
# B, D, desc = BASES[end]
# at = data[1][1]
# Es = [energy(b, at) for b in B2]
# display(Es)
#
# @show Es
