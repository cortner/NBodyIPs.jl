
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

data = load_data(250)
@show length(data)
train_data = data[1:220]
test_data = data[221:250]


# parameters for fitting
# -----------------------
# [1] DICTTYPE: which dictionary to use, look at the top of `polynomials.jl`
# which symbols are implemented and which dictionary they generate, or
# how to generate new dictionaries
DICTTYPE = :inv2
# [2] RCUT : obviously the cut-off radius, I found for Si a good rule of
# thumb is to use twice the site-energy cutoff!
RCUT = [4.1, 3.1, 2.1] * rnn(:Ti)   # [2.1, 3.1, 4.1]
# [3] NICT : how many entries in the 1D basis (=dictionary), essentially
# the polynomial degree
NDICT = 4:2:8    #

tibasis(ndict) = vcat( get_basis(:inv2, [ndict+4, ndict+2, ndict], RCUT)... )

errE = zeros(length(NDICT))
errF = zeros(length(NDICT))
nbasis = zeros(Int, length(NDICT))

for (in, ndict) = enumerate(NDICT)
   # this generates the permutation symmetric polynomials then
   # wraps them into calculators that represent the actual basis functions
   # for total energies
   B = tibasis(ndict)
   @show ndict, length(B)
   nbasis[in] = length(B)
   # standard least squares (see NBodyIPs/src/fitting.jl)
   # nforces = number of (randomly chosen) forces per configuration added
   #           to the LSQ problem
   ndata = min(length(train_data), length(B))
   c = regression(B, train_data[1:ndata], nforces = 5)
   # construct an IP from the the basis and the weights
   IP = NBodyIP(B, c)
   # check error => the normalisation is w.r.t. natoms, not a genuine
   # relative error; we can discuss
   errE[in], errF[in] = rms(IP, test_data)
   println("   E-rms on testset = ", errE[in])
   println("   F-rms on testset = ", errF[in])
end

B = tibasis(4)
c = rand(length(B))
IP = NBodyIP(B, c)
rms(IP, test_data)

using DataFrames
df = DataFrame(:nbasis => nbasis)
df[Symbol("E")] = errE
df[Symbol("F")] = errF
println("Energy and Force Errors: (float numbers are the cut-offs)")
println(df)
