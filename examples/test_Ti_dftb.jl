
using JuLIP, NBodyIPs, PyCall, ProgressMeter, ASE
using NBodyIPs.Invariants

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

data = load_data(250)
train_data = data[1:220]
test_data = data[221:250]

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

rcut2 = 8.5
D2 = Dictionary(InvInvariants, rcut2)
B2 = gen_basis(2, D2, 12)

push!(BASES, ([B1; B2], D2, "2 / $(length(B2)) / $rcut2"))

println("   first a few 3-body bases ...")
for (deg, rcut) in zip( [6, 8, 10],
                        [4.5, 4.5, 4.5, 4.5] )
   D = Dictionary(InvInvariants, rcut)
   B3 = gen_basis(3, D, deg)
   B = [B1; B2; B3]
   push!(BASES, (B, D, "2+3 / $(length(B)) / 8.5+$(round(rcut,2))"))
end

# # add a few 4-body
# rcut3 = 2.3 * r0
# D3 = Dictionary(InvInvariants, rcut3)
# B3 = gen_basis(3, D3, 10)
# for (deg, rcut) in zip([4, 6, 8],
#                        [1.7, 1.7, 1.7] * r0)
#    D = Dictionary(InvInvariants, rcut)
#    B4 = gen_basis(4, D, deg)
#    B = [B3; B4]
#    push!(BASES, (B, D, "3+4 / $(length(B3))+$(length(B4)) / $(round(rcut3,2))+$(round(rcut,2))") )
# end

for (deg, rcut) in zip([4, 6, 8],
                       [4.5, 4.5, 4.5, 4.5])
   D = Dictionary(InvInvariants, rcut)
   B4 = gen_basis(4, D, deg)
   B = [B1; B2; B4]
   push!(BASES, (B, D, "2+4 / $(length(B2))+$(length(B4)) / $rcut2+$(round(rcut,2))") )
end
#

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
   ndata = min(length(train_data), length(B))
   c = regression(B, train_data[1:ndata], nforces = 20, stab = 0.0)
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

using DataFrames
df = DataFrame(:desc => [B[3] for B in BASES])
df[Symbol("rms-E")] = rmsE
df[Symbol("rms-F")] = rmsF
df[Symbol("mae-E")] = maeE
df[Symbol("mae-F")] = maeF
println("Energy and Force Errors for Ti-DFTB Database: ")
println(df)
