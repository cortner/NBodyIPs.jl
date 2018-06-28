using Base.Test

using NBodyIPs
const Env = NBodyIPs.EnvBLs

println("Setting up the test systems ...")
r0 = rnn(:Cu)
TRANSFORM = let r0 = r0
   # (@analytic r -> (r0/r)^3)
   (@analytic r -> exp( - 3 * ((r/r0) - 1)))
end
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
rcut3 = 2.1 * r0
D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )

random_3body(ntup=1) = NBody( [tuple( [rand(1:4, 3); 0]... ) for n = 1:ntup],
                              (0.1+rand(ntup))/factorial(3), D3 )

Vn = let r0=r0, rcut=1.5*r0
   (@analytic r -> exp(- 3 * ((r/r0)-1))) * C1Shift(rcut)
end

V3env_0 = Env.EnvBL([0,], [1.0+rand()], random_3body(), Vn)
V3env_1 = Env.EnvBL([1,], [1.0+rand()], random_3body(), Vn)
V3env_2 = Env.EnvBL([2,], [1.0+rand()], random_3body(), Vn)

println("Finite-difference test on configurations")
println("------------------------------------------------")
set_constraint!(at, VariableCell(at))
for Venv in [V3env_0, V3env_1, V3env_2]
   (@test JuLIP.Testing.fdtest(Venv, at)) |> println
end


B = Env.envbl_basis(3, D3, 5, Vn, 2)

E1 = energy(B, at)
E2 = [ energy(b, at) for b in B ]
@test E1 ≈ E2

F1 = forces(B, at)
F2 = [ forces(b, at) for b in B ]
@test F1 ≈ F2

V1 = virial(B, at)
V2 = [ virial(b, at) for b in B ]
@test V1 ≈ V2
