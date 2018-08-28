
using NBodyIPs, JuLIP, Base.Test
using NBodyIPs.Data: Dat
using NBodyIPs: LsqSys

# generate random data
function generate_data(species, L, rmax, N, calc)
   data = Dat[]
   for n = 1:N
      at = bulk(species; cubic=true, pbc=true) * L
      rattle!(at, rand() * rmax)
      E = energy(calc, at)
      F = forces(calc, at)
      push!(data, Dat(at, E, F, nothing, 1.0, "rand"))
   end
   return data
end

r0 = rnn(:Si)
calc = StillingerWeber()
data = generate_data(:Si, 2, 0.2*r0, 30, calc)


TRANSFORM = "(@analytic r -> exp( - 2 * (r/$r0 - 1) ))"

rcut2 = 2 * cutoff(calc)
CUTOFF2 = "(:cos, $(rcut2-1), $rcut2)"
D2 = Dictionary(TRANSFORM, CUTOFF2)

rcut3 = 2 * cutoff(calc)
CUTOFF3 = "(:cos, $(rcut3-1), $rcut3)"
D3 = Dictionary(TRANSFORM, CUTOFF3)

deg3 = 10
B = poly_basis(3, D3, deg3)
@show length(B)
LsqSys(data[1:2], B)
@time LsqSys(data, B)
