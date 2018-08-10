
using NBodyIPs, JuLIP, BenchmarkTools, StaticArrays, ForwardDiff
using Base.Test

using NBodyIPs: BondLengthDesc, BondAngleDesc,
                transform, transform_d, fcut, fcut_d,
                invariants, invariants_d, invariants_ed
using NBodyIPs.Polys: NBPoly

Desc = BondAngleDesc

println("Setting up the test systems ...")
r0 = rnn(:Cu)
TRANSFORM = "r -> exp( - 3 * ((r/$r0) - 1))"
rcut4 = 2.1 * r0
desc = Desc(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )

const Ir = @SVector [1,2,3]
const Iθ = @SVector [4,5,6]
all_invs(rθ) = vcat(invariants(desc, (rθ[Ir], rθ[Iθ]))...)
jac_all_invs(rθ) = hcat(vcat(invariants_ed(desc, (rθ[Ir], rθ[Iθ]))[3:4]...)...)'
# ad_all_invs(rθ) = ForwardDiff.jacobian(all_invs, rθ)

function fdjac(rθ; h=1e-5)
   v = Vector(rθ)
   f = all_invs(rθ)
   J = zeros(length(f), length(v))
   for i = 1:length(v)
      v[i] += h
      fh = all_invs(SVector(v...))
      J[:, i] = (fh-f) / h
      v[i] -= h
   end
   return J
end

rθ = (@SVector rand(6)) .+ r0
all_invs(rθ)
jac_all_invs(rθ)
# ad_all_invs(rθ)

for h in [1e-3, 1e-4, 1e-5, 1e-6]
   @show h
   err = (jac_all_invs(rθ) - fdjac(rθ; h=h))
   @show vecnorm(err, Inf)
end
