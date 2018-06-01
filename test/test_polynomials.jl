using NBodyIPs, JuLIP, BenchmarkTools, StaticArrays
using Base.Test

profile = false

if profile
   nbasis = [0, 10, 50, 100, 300]
else
   nbasis = [0, 10, 10, 10, 10]
end

println("Setting up the test systems ...")
r0 = rnn(:Cu)
TRANSFORM = let r0 = r0
   # (@analytic r -> (r0/r)^3)
   (@analytic r -> exp( - 3 * ((r/r0) - 1)))
end
at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
rcut3 = 3.1 * r0
D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
rcut4 = 2.1 * r0
D4 = Dictionary(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
rcut5 = 1.5 * r0
D5 = Dictionary(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )
DD = [nothing, D3, D3, D4, D5]

random_nbody(N, ntup) = (
   (N <= 3) ? NBody( [tuple( [rand(0:4, ((N*(N-1))÷2)); 0]... ) for n = 1:ntup],
                     (0.1+rand(ntup))/factorial(N), DD[N] )
           : NBody( [tuple( rand(0:N, ((N*(N-1))÷2+1))... ) for n = 1:ntup],
                   (0.1+rand(ntup))/factorial(N), DD[N] ) )


println("`NBody` gradient-test on simplices")
println("----------------------------------")
for N = 2:5, ntup = [1, 3]
   VN = random_nbody(N, ntup)
   r = SVector( (r0 + rand((N*(N-1))÷2))... )
   dvN = @D VN(r)
   vN = VN(r)
   rv = Vector(r)
   errs = []
   println("      h   |    err   [N=$N, ntup=$ntup]")
   for p = 2:10
      h = 0.1^p
      dvh = zeros(length(dvN))
      for i = 1:length(r0)
         rv[i] += h
         dvh[i] = (VN( SVector(rv...) ) - vN) / h
         rv[i] -=h
      end
      push!(errs, vecnorm(dvN - dvh, Inf))
      @printf(" %.2e | %.2e \n", h, errs[end])
   end
   # @test minimum(errs) <= 1e-3 * maximum(errs)
end



println("`NBody` finite-difference test on configurations")
println("------------------------------------------------")
nb = 3
at1 = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
at2 = bulk(:Cu, cubic=true) * (1,1,2)
set_constraint!(at2, VariableCell(at2, free = []))
for at in [at1, at2]
   for N = 2:5 #nbas = [1,3]
      println("  $N-body")
      VN = random_nbody(N, 1)
      # @test
      JuLIP.Testing.fdtest(VN, at)
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

profile = true
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
