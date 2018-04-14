
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
train_data = data[1:100]
test_data = data[101:120]


# parameters for fitting
# -----------------------
rcut = 4.1 * rnn(:Ti)
DEG = 4:2:10
D = Dictionary(InvInvariants, rcut)

tibasis(deg) = gen_basis(3, D, deg)

errE = zeros(length(DEG))
errF = zeros(length(DEG))
nbasis = zeros(Int, length(DEG))

for (in, deg) = enumerate(DEG)
   # this generates the permutation symmetric polynomials then
   # wraps them into calculators that represent the actual basis functions
   # for total energies
   B = tibasis(deg)
   @show deg, length(B)
   nbasis[in] = length(B)
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
df = DataFrame(:nbasis => nbasis)
df[Symbol("E")] = errE
df[Symbol("F")] = errF
println("Energy and Force Errors: (float numbers are the cut-offs)")
println(df)
