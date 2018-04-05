
using JuLIP, NBodyIPs, NeighbourLists


function gen_data(N, rnd=0.1)
   sw = StillingerWeber()
   r0 = rnn(:Si)
   rcut = cutoff(sw)
   data = Tuple{typeof(bulk(:Si)), Float64}[]
   for n = 1:N
      at = bulk(:Si, cubic=true) * 2
      rattle!(at, rnd * r0)
      push!(data, (at, energy(sw, at)))
   end
   return data
end

train_data = gen_data(1_000, 0.1)
test_data =  gen_data(100, 0.1)

sw = StillingerWeber()
r0 = rnn(:Si)
rcut = cutoff(sw)
rcutN = 2 * rcut

basis(ndict::Integer)  =
   get_basis(3, dict(:inv2, ndict, rcutN)..., rcutN)

NDICT = 4:2:12
err = zeros(length(NDICT))
nbasis = zeros(Int, length(NDICT))

for (in, ndict) in enumerate(NDICT)
   B = basis(ndict)
   nbasis[in] = length(B)
   @show (ndict, length(B))
   c = NBodyIPs.regression(B, train_data[1:5*length(B)])
   err[in] = NBodyIPs.rms(c, B, test_data)
   println("rms on testset = ", err[in])
end

using DataFrames
df = DataFrame(:nbasis => nbasis, :err_inv2 => err)
println(df)
