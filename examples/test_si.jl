
using JuLIP, ManyBodyIPs, JSON

si_data = "/Users/ortner/Dropbox/PIBmat/si_data.json"

function load_json(fname)
   data_any = JSON.parsefile(fname)
   data = [ ( Atoms(:Si, convert(Vector{JVecF}, d[1])),
              d[2]::Float64,
              convert(Vector{JVecF}, d[3]) )    for d in data_any ]
   return data
end

data = load_json(si_data)
length(data)

r0 = rnn(:Si)
rcut = 3.2 * r0
D = dict(:inv2, 5, rcut)
B = get_basis(3, D..., rcut)
