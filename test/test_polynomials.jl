using NBodyIPs, JuLIP
using BenchmarkTools
using Base.Test

profile = true
if profile
   nbasis3 = 50
   nbasis4 = 100
   nbasis5 = 300
else
   nbasis3 = 10
   nbasis4 = 10
   nbasis5 = 10
end

println("Testing Collective Assembly across a Basis Set")
r0 = rnn(:Cu)
TRANSFORM = let r0 = r0
   # (@analytic r -> (r0/r)^3)
   (@analytic r -> exp( - 3 * ((r/r0) - 1)))
end
at = rattle!(bulk(:Cu, cubic=true) * 3, 0.02)

println("3-body")
rcut3 = 3.1 * r0
D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
B3 = [ NBody( [tuple([rand(0:4, 3);0]...)], [1.0+rand()], D3 )
      for n = 1:nbasis3 ]
E1 = [energy(b, at) for b in B3]
E2 = energy( B3, at )
(@test E1 ≈ E2) |> println
F1 = [forces(b, at) for b in B3]
F2 = forces(B3, at)
(@test F1 ≈ F2) |> println
S1 = [stress(b, at) for b in B3]
S2 = stress(B3, at)
(@test S1 ≈ S2) |> println;


println("4-body")
rcut4 = 2.1 * r0
D4 = Dictionary(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
B4 = [ NBody( [tuple(rand(0:3, 7)...)], [1.0+rand()], D4 )
      for n = 1:nbasis4 ]
E1 = [energy(b, at) for b in B4]
E2 = energy( B4, at )
(@test E1 ≈ E2) |> println;
F1 = [forces(b, at) for b in B4]
F2 = forces(B4, at)
(@test F1 ≈ F2) |> println;
S1 = [stress(b, at) for b in B4]
S2 = stress(B4, at)
(@test S1 ≈ S2) |> println;

println("5-body")
rcut5 = 1.5 * r0
D5 = Dictionary(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )
B5 = [ NBody( [tuple(rand(0:5, 11)...)], [1.0+rand()], D5 )
      for n = 1:nbasis5 ]
E1 = [energy(b, at) for b in B5]
E2 = energy( B5, at )
(@test E1 ≈ E2) |> println;
F1 = [forces(b, at) for b in B5]
F2 = forces(B5, at)
(@test F1 ≈ F2) |> println;
S1 = [stress(b, at) for b in B5]
S2 = stress(B5, at)
(@test S1 ≈ S2) |> println;

if profile
   println("Performance")
   print(" E 3-body old:" ); @btime ([energy(b, at) for b in B3])
   print(" E 3-body new:" ); @btime energy( B3, at )
   print(" F 3-body old:" ); @btime ([forces(b, at) for b in B3])
   print(" F 3-body new:" ); @btime forces( B3, at )
   print(" S 3-body old:" ); @btime ([stress(b, at) for b in B3])
   print(" S 3-body new:" ); @btime stress( B3, at )
   print(" E 4-body old:" ); @btime ([energy(b, at) for b in B4])
   print(" E 4-body new:" ); @btime energy( $B4, $at )
   print(" F 4-body old:" ); @btime ([forces(b, at) for b in B4])
   print(" F 4-body new:" ); @btime forces( $B4, $at )
   print(" S 4-body old:" ); @btime ([stress(b, at) for b in B4])
   print(" S 4-body new:" ); @btime stress( $B4, $at )
   print(" E 5-body old:" ); @btime ([energy(b, at) for b in B5])
   print(" E 5-body new:" ); @btime energy( $B5, $at )
   print(" F 5-body old:" ); @btime ([forces(b, at) for b in B5])
   print(" F 5-body new:" ); @btime forces( $B5, $at )
   print(" S 5-body old:" ); @btime ([stress(b, at) for b in B5])
   print(" S 5-body new:" ); @btime stress( $B5, $at )
end




# println("`NBody` gradient-test on simplices")
# for n = [1, 3]
#    V3 = NBody( [tuple([rand(0:3, 3); 0]...) for n = 1:n], 1.0 + rand(n), D3 )
#    for _  = 1:10
#       r = 1.0 + rand(SVector{3,Float64})
#       @test (@D V3(r)) ≈ ForwardDiff.gradient(r_ -> V3(r_), r)
#       print(".")
#    end
# end

# for n = [1, 3]
#    V4 = NBody( [tuple(rand(0:3, 7)...) for n = 1:n], 1.0 + rand(n), D4 )
#    for _  = 1:10
#       r = 1.0 + rand(SVector{6,Float64})
#       @test evaluate_d(V4, r) ≈ ForwardDiff.gradient(r_ -> V4(r_), r)
#       print(".")
#    end
# end
# println()
#
# println("[5] `NBody` finite-difference test on configurations")
# nb = 3
# at1 = rattle!(bulk(:Cu, cubic=true) * (1,2,2), 0.02)
# at2 = bulk(:Cu, cubic=true) * (1,1,2)
# set_constraint!(at2, VariableCell(at2, free = []))
# for at in [at1, at2]
#    println("  3-body")
#    V3 = NBody( [tuple([rand(0:5, 3); 0]...) for n = 1:nb], 1.0+rand(nb), D3 )
#    @test JuLIP.Testing.fdtest(V3, at)
#
#    println("  4-body")
#    V4 = NBody( [tuple(rand(0:5, 7)...) for n = 1:nb], 1.0 + rand(nb), D4 )
#    @test JuLIP.Testing.fdtest(V4, at)
# end




# at = rattle!(bulk(:Cu, cubic=true) * 2, 0.02)
#
# @btime ([stress(b, $at) for b in $B3])
# @btime stress( $B3, $at )
#
# Profile.clear()
# @profile stress( B3, at );
#
# Profile.print()
#
# stress(B3, at)
