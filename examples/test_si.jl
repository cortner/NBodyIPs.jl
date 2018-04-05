
using JuLIP, NBodyIPs, JSON, NeighbourLists

si_data = "/Users/ortner/Dropbox/PIBmat/si_data.json"

function load_json(fname)
   data_any = JSON.parsefile(fname)
   data = [ ( Atoms(:Si, convert(Vector{JVecF}, d[1])),
              d[2]::Float64,
              convert(Vector{JVecF}, d[3]) )    for d in data_any ]
   return data
end

basis(ndict::Integer, bord::Integer)  =
   get_basis(bord, dict(:inv2, ndict, rcut)..., rcut)

data = load_json(si_data)
length(data)

r0 = rnn(:Si)
rcut = 3.2 * r0

# test_params =
#       [ (8, 3), (10, 3), (12, 3), (7, 4), (8, 4), (6, 5) ]
#
# for (ndict, bord) in test_params
#    @show (ndict, bord)
#    @show length(basis(ndict, bord))
# end



test_params =
      [ (8, 3), (10, 3), (12, 3), (7, 4), (8, 4) ]

for (ndict, bord) in test_params
   B = basis(ndict, bord)
   @show (ndict, bord, length(B))
   c = NBodyIPs.regression(B, data[1:400])
end

at = data[1][1]::Atoms
for b in basis(8, 4)
   @time b(at)
end

# ENV["JULIA_NUM_THREADS"] = 1
# NeighbourLists.set_maxthreads!(1)
# NeighbourLists.MAX_THREADS
