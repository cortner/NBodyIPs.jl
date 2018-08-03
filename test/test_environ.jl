using Base.Test

using JuLIP, NBodyIPs
using NBodyIPs.BLPolys: BLDictionary, BLNBody
const Env = NBodyIPs.EnvBLs

println("-------------------")
println("Tests for EnvBLs")
println("-------------------")
println("Setting up the test systems ...")
r0 = rnn(:Cu)
rcut3 = 2.1 * r0
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
D3 = BLDictionary("exp( - 3 * ((r/$r0) - 1))",
                 "(:cos, $(0.66*rcut3), $rcut3)" )

random_3body(ntup=1) = BLNBody( [tuple( [rand(1:4, 3); 0]... ) for n = 1:ntup],
                                 (0.1+rand(ntup))/factorial(3), D3 )

Vn = ("exp(- 3 * ((r/$r0)-1))", 1.8*r0)


V3env_0 = Env.EnvBL(0, random_3body(), Vn...)
V3env_1 = Env.EnvBL(1, random_3body(), Vn...)
V3env_2 = Env.EnvBL(2, random_3body(), Vn...)

println("Finite-difference test on configurations")
println("------------------------------------------------")
set_constraint!(at, VariableCell(at))
for Venv in [V3env_0, V3env_1, V3env_2]
   (@test JuLIP.Testing.fdtest(Venv, at)) |> println
end



println("Testing combined assembly of basis")
B = Env.envbl_basis(3, D3, 5, Vn, 2)

# nb: the , false is to set typewarn=false, since eltype(B) is not a leaftype
print("energy: ")
E1 = energy(B, at, false)
E2 = [ energy(b, at) for b in B ]
println(@test E1 ≈ E2)

print("forces: ")
F1 = forces(B, at, false)
F2 = [ forces(b, at) for b in B ]
println(@test F1 ≈ F2)

print("virial: ")
V1 = virial(B, at, false)
V2 = [ virial(b, at) for b in B ]
println(@test V1 ≈ V2)
