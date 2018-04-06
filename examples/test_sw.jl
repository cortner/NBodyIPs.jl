
using JuLIP, NBodyIPs


function gen_data(N, rnd=0.1)
   sw = StillingerWeber()
   r0 = rnn(:Si)
   rcut = cutoff(sw)
   data = Tuple{Atoms{Float64, Int}, Float64, JVecsF}[]
   for n = 1:N
      at = bulk(:Si, cubic=true) * 2
      rattle!(at, rnd * r0)
      push!(data, (at, energy(sw, at), forces(sw, at)) )
   end
   return data
end

train_data = gen_data(200, 0.1)
test_data =  gen_data(100, 0.1)

sw = StillingerWeber()
r0 = rnn(:Si)
rcut = cutoff(sw)
rcutN = 2 * rcut

basis(ndict::Integer)  =
   get_basis(3, dict(:poly, ndict, rcutN)..., rcutN)

NDICT = 4:2:12
err = zeros(length(NDICT))
nbasis = zeros(Int, length(NDICT))

for (in, ndict) in enumerate(NDICT)
   B = basis(ndict)
   nbasis[in] = length(B)
   @show (ndict, length(B))
   c = NBodyIPs.regression(B, train_data, nforces = 5)
   IP = NBodyIP([NBodies(3, c, [b.f for b in B], [b.d_f for b in B], rcutN)])
   err[in] = NBodyIPs.rms(c, B, test_data)
   println("rms on testset = ", err[in])
end

using DataFrames
df = DataFrame(:nbasis => nbasis, :err_inv2 => err)
println(df)
