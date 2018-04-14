
using JuLIP, NBodyIPs, JSON, NeighbourLists
using NBodyIPs.Invariants

si_data = "/Users/ortner/Dropbox/PIBmat/si_data.json"

function load_json(fname)
   data_any = JSON.parsefile(fname)
   data = [ ( Atoms(:Si, convert(Vector{JVecF}, d[1])),
              d[2]::Float64,
              convert(Vector{JVecF}, d[3]) )    for d in data_any ]
   return data
end

data = load_json(si_data)
train_data = data[1:100]
test_data =  data[101:120]

r0 = rnn(:Si)
rcutN = 4.3*r0   # (ca 2 x SW cutoff)
DEG = 4:2:10
errE = zeros(length(DEG))
errF = zeros(length(DEG))
nbasis = zeros(Int, length(DEG))

D = Dictionary(InvInvariants, rcutN)
basis(deg) = gen_basis(3, D, deg)

for (in, deg) in enumerate(DEG)
   B = basis(deg)
   nbasis[in] = length(B)
   @show (deg, length(B))
   c = regression(B, train_data, nforces = 5)
   IP = NBody(B, c, D)
   errE[in], errF[in] = rms(IP, test_data)
   println("   E-rms on testset = ", errE[in])
   println("   F-rms on testset = ", errF[in])
end

using DataFrames
df = DataFrame(:nbasis => nbasis, :errE => errE, :errF => errF)
println(df)
