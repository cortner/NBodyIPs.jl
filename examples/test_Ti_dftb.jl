
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

data = load_data(150)
@show length(data)
train_data = data[1:100]
test_data = data[101:150]

basis(ndict::Integer, bord::Integer, rcut, sym=:inv2)  =
   get_basis(bord, dict(sym, ndict, rcut)..., rcut)

# parameters for fitting
# -----------------------
# [1] DICTTYPE: which dictionary to use, look at the top of `polynomials.jl`
# which symbols are implemented and which dictionary they generate, or
# how to generate new dictionaries
DICTTYPE = :inv2
# [2] RCUT : obviously the cut-off radius, I found for Si a good rule of
# thumb is to use twice the site-energy cutoff!
RCUT = [2.1, 3.1, 4.1] * rnn(:Ti)   # [2.1, 3.1, 4.1]
# [3] NICT : how many entries in the 1D basis (=dictionary), essentially
# the polynomial degree
NDICT = 4:2:12    # 4:2:12
# [4] BORD : just 3 for now, 4 is very slow, 5 is impossible. we need
# some optimisations first!
BORD = 3

errE = zeros(length(NDICT), length(RCUT))
errF = zeros(length(NDICT), length(RCUT))
nbasis = zeros(Int, length(NDICT))

for (in, ndict) = enumerate(NDICT), (ir, rcut) in enumerate(RCUT)
   # this generates the 3-dimensional permutation symmetric polynomials then
   # wraps them into calculators that represent the actual basis functions
   # for total energies
   B = basis(ndict, BORD, rcut, DICTTYPE)
   @show ndict, length(B), rcut
   nbasis[in] = length(B)
   # standard least squares (see NBodyIPs/src/fitting.jl)
   # nforces = number of (randomly chosen) forces per configuration added
   #           to the LSQ problem
   c = regression(B, train_data, nforces = 5)
   # construct an IP from the the basis and the weights
   IP = NBodyIP(B, c)
   # check error => the normalisation is w.r.t. natoms, not a genuine
   # relative error; we can discuss
   errE[in, ir], errF[in, ir] = rms(IP, test_data)
   println("   E-rms on testset = ", errE[in, ir])
   println("   F-rms on testset = ", errF[in, ir])
end


using DataFrames
df = DataFrame(:nbasis => nbasis)
for (ir, rcut) in enumerate(RCUT)
   df[Symbol("E($(round(rcut,2)))")] = errE[:, ir]
   df[Symbol("F($(round(rcut,2)))")] = errF[:, ir]
end
println("Energy and Force Errors: (float numbers are the cut-offs)")
println(df)
