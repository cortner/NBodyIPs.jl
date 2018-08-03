using NBodyIPs, JuLIP, BenchmarkTools, StaticArrays
using Base.Test


using NBodyIPs.BLPolys: BLNBody, BLDictionary, transform, transform_d, fcut, fcut_d
using NBodyIPs.BLPolys.BLInvariants: invariants, invariants_d, invariants_ed

profile = false

if profile
   nbasis = [0, 10, 50, 100, 300]
else
   nbasis = [0, 10, 10, 10, 10]
end

println("Setting up the test systems ...")
r0 = rnn(:Cu)
TRANSFORM = "r -> exp( - 3 * ((r/$r0) - 1))"
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
rcut3 = 3.1 * r0
D3 = BLDictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
rcut4 = 2.1 * r0
D4 = BLDictionary(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
rcut5 = 1.5 * r0
D5 = BLDictionary(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )
DD = [nothing, D3, D3, D4, D5]

random_nbody(N, ntup) = (
   (N <= 3) ? BLNBody( [tuple( [rand(1:4, ((N*(N-1))÷2)); 0]... ) for n = 1:ntup],
                     (0.1+rand(ntup))/factorial(N), DD[N] )
           : BLNBody( [tuple( rand(0:N, ((N*(N-1))÷2+1))... ) for n = 1:ntup],
                   (0.1+rand(ntup))/factorial(N), DD[N] ) )

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
dfh = (fcut.(D, rr + 1e-5) - fcut.(D, rr - 1e-5)) / (2e-5)
(@test norm(dfh - fcut_d.(D, rr), Inf) < 1e-8) |> println


println("testing a transformed invariant")
println("-------------------------------")
r = SVector( (r0 + rand(3))... )
rv = Vector(r)
I1, _, dI1, _ = invariants_ed(D, r)
I1 = I1[2]
dI1 = dI1[2]
errs = []
for p = 2:11
   h = 0.1^p
   dI1h = zeros(3)
   for n = 1:3
      rv[n] += h
      dI1h[n] = (invariants(D, SVector(rv...))[1][2] - I1) / h
      rv[n] -= h
   end
   push!(errs, norm(dI1h - dI1, Inf))
   @printf(" %.2e | %.2e \n", h, errs[end])
end
(@test minimum(errs) <= 1e-3 * maximum(errs)) |> println


println("`BLNBody` gradient-test on simplices")
println("----------------------------------")
for N in 2:5, ntup = [1,3]
   VN = random_nbody(N, ntup)
   println("[$N-body, ntup=$ntup, t[1] = $(VN.t[1])]")
   r = SVector( (r0 + rand((N*(N-1))÷2))... )
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


println("`BLNBody` finite-difference test on configurations")
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


println("Testing Collective Assembly across a Basis Set")
println("----------------------------------------------")

BB = [ nothing; [ [random_nbody(N, 1) for _=1:nbasis[N]] for N = 2:5 ] ]

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