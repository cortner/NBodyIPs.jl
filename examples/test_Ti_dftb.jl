
using JuLIP, ManyBodyIPs, PyCall, ProgressMeter, ASE

function load_data()
   fname = "/Users/ortner/Dropbox/PIBmat/Ti_DFTB_Data/Ti_N54_T2000.xyz"
   Nconfig = 411  # (how can I load this?)
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

data = load_data()
@show length(data)
train_data = data[1:350]
test_data = data[351:end]

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
RCUT = [2.1, 3.1, 4.1] * rnn(:Ti)
# [3] NICT : how many entries in the 1D basis (=dictionary), essentially
# the polynomial degree
NDICT = 4:2:10
# [4] BORD : just 3 for now, 4 is very slow, 5 is impossible. we need
# some optimisations first!
BORD = 3

errors = zeros(length(NDICT), length(RCUT))
nbasis = zeros(Int, length(NDICT))

for (in, ndict) = enumerate(NDICT), (ir, rcut) in enumerate(RCUT)
   # this generates the 3-dimensional permutation symmetric polynomials then
   # wraps them into calculators that represent the actual basis functions
   # for total energies
   B = basis(ndict, BORD, rcut, DICTTYPE)
   @show ndict, length(B)
   nbasis[in] = length(B)
   # standard least squares (see ManyBodyIPs/src/fitting.jl)
   c = regression(B, train_data)
   # check error => the normalisation 54 / 300 = natoms / typical total energy
   # that is, the errors stored are roughly the relative error per atom
   # if this is not what you are after, then we can discuss.
   errors[in, ir] = ManyBodyIPs.rms(c, B, test_data) * sqrt(54 / 300)
   println("rms on testset = ", errors[in, ir])
end


using DataFrames
df = DataFrame(:nbasis => nbasis)
for (ir, rcut) in enumerate(RCUT)
   df[Symbol("$(round(rcut,2))")] = errors[:, ir]
end
println("Errors: (float numbers are the cut-offs)")
println(df)
