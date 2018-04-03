
using JuLIP, ManyBodyIPs, NeighbourLists

# generate data
sw = StillingerWeber()
r0 = rnn(:Si)
rcut = cutoff(sw)
data = Tuple{typeof(bulk(:Si)), Float64}[]
for n = 1:1_000
   at = bulk(:Si, cubic=true) * 2
   rattle!(at, 0.2 * r0)
   push!(data, (at, energy(sw, at)))
end

ManyBodyIPs.dict(n) =
   ["r^$n * (exp(-4*(r/$r0-1)) - $(exp(-4*(rcut/r0-1))) + $(4/r0 * exp(-4*(rcut/r0-1))) * (r-$rcut)"
      for n = 0:(n-1)], "r"

# basis(ndict::Integer)  =
#    get_basis(3, dict(ndict)..., rcut)

basis(ndict::Integer)  =
   get_basis(3, dict(:inv1, ndict, rcut)..., rcut)

for ndict = 6:14
   B = basis(ndict)
   @show (ndict, length(B))
   c = ManyBodyIPs.regression(B, data[1:3*length(B)])
   println("rms on testset = ", ManyBodyIPs.rms(c, B, data[801:900]))
end


# at = data[1][1]::Atoms
# for b in basis(12)
#    @time b(at)
# end
