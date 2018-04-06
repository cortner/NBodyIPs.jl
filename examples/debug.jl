
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

train_data = gen_data(50, 0.1)
test_data =  gen_data(50, 0.1)

sw = StillingerWeber()
r0 = rnn(:Si)
rcut = cutoff(sw)
rcutN = 2 * rcut

basis(ndict::Integer)  =
   get_basis(3, dict(:inv2, ndict, rcutN)..., rcutN)

B = basis(4)
c = NBodyIPs.regression(B, train_data, nforces = 0)
IP = NBodyIP(B, c)

at = bulk(:Si, cubic=true) * 3
energy(IP, at)

errE, errF = rms(IP, test_data)
println("E-rms on testset = ", errE[in])
println("F-rms on testset = ", errF[in])
