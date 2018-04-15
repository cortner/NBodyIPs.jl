
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
@show length(data)
train_data = data[1:220]
test_data = data[221:250]

# useful cutoffs to try
# (3.5 x r0 is the cell dimension)
r0 = rnn(:Ti)
rcuts = [2.1, 2.8, 3.5] * r0

# some notes on the orders of magnitude
# E per atom ~ 6.0 eV
# mean force on each atom is ~ 1 eV / A
# but range of forces is between 1 and 12 eV / A

# construct a sequence of bases
# -----------------------------

println("generating basis functions")

BASES = []

println("   first a few 3-body bases ...")
for (deg, rcut) in zip( [4, 6, 8], rcuts )
   D = Dictionary(InvInvariants, rcut)
   B = gen_basis(3, D, deg)
   push!(BASES, (B, D, "3 / $(length(B)) / $(round(rcut,2))"))
end

# # add a few 4-body
# for (deg, rcut) in zip( [4, 6, 8], rcuts )
#    D = Dictionary(InvInvariants, rcut)
#    B = gen_basis(3, D, deg)
#    push!(BASES, B, "3 / $(length(B)) / $(round(rcut,2))")
# end


errE = zeros(length(BASES))
errF = zeros(length(BASES))

println("For each generated basis test the fit: " )
for (B, D, description) in BASES
   println(" - ", description)
   # standard least squares (see NBodyIPs/src/fitting.jl)
   # nforces = number of (randomly chosen) forces per configuration added
   #           to the LSQ problem
   ndata = min(length(train_data), length(B))
   c = regression(B, train_data[1:ndata], nforces = 5)
   # construct an IP from the the basis and the weights
   # IP = NBodyIP(B, c, D)
   IP = NBody(B, c, D)
   # check error => the normalisation is w.r.t. natoms, not a genuine
   # relative error; we can discuss
   errE[in], errF[in] = rms(IP, test_data)
   println("   E-rms on testset = ", errE[in])
   println("   F-rms on testset = ", errF[in])
end

using DataFrames
df = DataFrame(:desc => [B[2] for B in BASES])
df[Symbol("E")] = errE
df[Symbol("F")] = errF
println("Energy and Force Errors for Ti-DFTB Database: ")
println(df)
