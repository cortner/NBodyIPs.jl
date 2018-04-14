
using JuLIP, NBodyIPs
using NBodyIPs.Invariants


function gen_data(N, rnd=0.1)
   sw = StillingerWeber()
   r0 = rnn(:Si)
   rcut = cutoff(sw)
   data = Tuple{Atoms{Float64, Int}, Float64, JVecsF}[]
   for n = 1:N
      at = bulk(:Si, cubic=true) * 2
      rattle!(at, rnd * r0)
      push!(data, (at, energy(sw, at), forces(sw, at)))
   end
   return data
end

train_data = gen_data(100, 0.1)
test_data =  gen_data(20, 0.1)

sw = StillingerWeber()
r0 = rnn(:Si)
rcut = cutoff(sw)
rcutN = 2 * rcut

DEG = 4:2:12
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
