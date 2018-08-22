using NBodyIPs, JuLIP, BenchmarkTools, StaticArrays
using Base.Test

using NBodyIPs: BondLengthDesc, BondAngleDesc,
                transform, transform_d, fcut, fcut_d,
                invariants, invariants_d, invariants_ed
using NBodyIPs.Polys: NBPoly

profile = false

if profile
   nbasis = [0, 10, 50, 100, 300]
else
   nbasis = [0, 10, 10, 10, 10]
end

const Ir = @SVector [1,2,3]
const Iθ = @SVector [4,5,6]
all_invs(desc::BondLengthDesc, r) = vcat(invariants(desc, r)...)
jac_all_invs(desc::BondLengthDesc, r) = hcat(vcat(invariants_ed(desc, r)[3:4]...)...)'
all_invs(desc::BondAngleDesc, rθ) = vcat(invariants(desc, (rθ[Ir], rθ[Iθ]))...)
jac_all_invs(desc::BondAngleDesc, rθ) = hcat(vcat(invariants_ed(desc, (rθ[Ir], rθ[Iθ]))[3:4]...)...)'

function fdjac(desc, rθ; h=1e-5)
   v = Vector(rθ)
   f = all_invs(desc, rθ)
   J = zeros(length(f), length(v))
   for i = 1:length(v)
      v[i] += h
      fh = all_invs(desc, SVector(v...))
      J[:, i] = (fh-f) / h
      v[i] -= h
   end
   return J
end


for Desc in [BondLengthDesc, BondAngleDesc]
# Desc = BondAngleDesc

println("------------------------------------------------")
println(" Testing Polys based on $(Desc) descriptors")
println("------------------------------------------------")

println("Setting up the test systems ...")
r0 = rnn(:Cu)
TRANSFORM = "r -> exp( - 3 * ((r/$r0) - 1))"
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.05)
rcut3 = 3.1 * r0
D3 = Desc(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
rcut4 = 2.1 * r0
D4 = Desc(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
rcut5 = 1.5 * r0
D5 = Desc(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )
DD = [nothing, D3, D3, D4, D5]

random_nbody(N, ntup) = (
   (N <= 3) ? NBPoly( [tuple( [rand(1:4, ((N*(N-1))÷2)); 0]... ) for n = 1:ntup],
                     (0.1+rand(ntup))/factorial(N), DD[N] )
           : NBPoly( [tuple( rand(0:N, ((N*(N-1))÷2+1))... ) for n = 1:ntup],
                   (0.1+rand(ntup))/factorial(N), DD[N] ) )

rand_rθ(::BondLengthDesc, N) = SVector( (1.0 + rand((N*(N-1))÷2))... )
rand_rθ(::BondAngleDesc, N) = SVector((1.0+rand(N-1))...),
                               SVector( (-0.5+rand(((N-1)*(N-2))÷2))... )


println("testing the transform")
println("---------------------")
D = D3
rr = linspace(r0, r0+1, 10)
tdh = (transform.(D, rr + 1e-5) - transform.(D, rr - 1e-5)) / (2e-5)
(@test norm(tdh - transform_d.(D, rr), Inf) < 1e-8) |> println

println("testing the cutoff")
println("------------------")
D = D3
rr = linspace(rcut3-0.9, rcut3+0.1, 100)
dfh = (fcut.(D.cutoff, rr + 1e-5) - fcut.(D.cutoff, rr - 1e-5)) / (2e-5)
(@test norm(dfh - fcut_d.(D.cutoff, rr), Inf) < 1e-8) |> println

println(" testing a transformed 4B invariant")
println("------------------------------------")
r = SVector( (r0 + rand(6))... )
DI = jac_all_invs(D4, r)
errs = []
for p = 2:11
   h = 0.1^p
   DIh = fdjac(D4, r; h = h)
   push!(errs, vecnorm(DIh - DI, Inf))
   @printf(" %.2e | %.2e \n", h, errs[end])
end
(@test minimum(errs) <= 1e-3 * maximum(errs)) |> println

if D isa BondLengthDesc
   println("`NBPoly` gradient-test on simplices")
   println("----------------------------------")
   for N in 2:5, ntup = [1,3]
      VN = random_nbody(N, ntup)
      println("[$N-body, ntup=$ntup, t[1] = $(VN.t[1])]")
      # r = SVector( (r0 + rand((N*(N-1))÷2))... )
      r = rand_rθ(D, N)
      dvN = @D VN(r)
      vN = VN(r)
      rv = Vector(r)
      errs = []
      println("      h   |    err")
      for p = 2:10
         h = 0.1^p
         dvh = zeros(length(dvN))
         for i = 1:length(rv)
            rv[i] += h
            dvh[i] = (VN( SVector(rv...) ) - vN) / h
            rv[i] -=h
         end
         push!(errs, vecnorm(dvN - dvh, Inf))
         @printf(" %.2e | %.2e \n", h, errs[end])
      end
      (@test minimum(errs) <= 1e-3 * maximum(errs)) |> println
   end
end

##
println("`NBPoly` finite-difference test on configurations")
println("------------------------------------------------")
at1 = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
at2 = bulk(:Cu)
set_constraint!(at2, VariableCell(at2, free = []))
for at in [at1, at2]
   for N = 3:4, nbas = [1,3]
      println("[$N-body, nbasis = $nbas]")
      VN = random_nbody(N, 1)
      (@test JuLIP.Testing.fdtest(VN, at)) |> println
   end
end

##
println("Testing Collective Assembly across a Basis Set")
println("----------------------------------------------")

BB = [ nothing; [ [random_nbody(N, 1) for _=1:nbasis[N]] for N = 2:5 ] ]
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.05)

for N = 2:5
   B = BB[N]
   println("[$N-body]")
   E1 = [energy(b, at) for b in B]
   E2 = energy( B, at )
   print("E: "); (@test E1 ≈ E2) |> println
   F1 = [forces(b, at) for b in B]
   F2 = forces(B, at)
   print("F: "); (@test F1 ≈ F2) |> println
   S1 = [virial(b, at) for b in B]
   S2 = virial(B, at)
   print("S: "); (@test S1 ≈ S2) |> println;
end

if profile
   println("Performance Tests")
   println("-----------------")
   println("to get a rough estimate of the cost per atom take the `separate`")
   println("timing and multiply by the factor provided")
   for N = 2:5
      B = BB[N]
      println("[$N-body], length(B) == $(length(B)), factor = $(1/(length(at)*length(B)))")
      print(" E separate:" ); @btime ([energy(b, $at) for b in $B])
      print(" E combined:" ); @btime energy( $B, $at )
      print(" F separate:" ); @btime ([forces(b, $at) for b in $B])
      print(" F combined:" ); @btime forces( $B, $at )
      print(" S separate:" ); @btime ([virial(b, $at) for b in $B])
      print(" S combined:" ); @btime virial( $B, $at )
   end
end


end
