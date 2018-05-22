using NBodyIPs, JuLIP
using BenchmarkTools
using Base.Test

println("Testing NBBasis")
r0 = rnn(:Cu)
TRANSFORM = let r0 = r0
   # (@analytic r -> (r0/r)^3)
   (@analytic r -> exp( - 3 * ((r/r0) - 1)))
end
at = rattle!(bulk(:Cu, cubic=true) * 3, 0.02)

println("3-body")
rcut3 = 3.1 * r0
D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
nbasis = 50
B3 = [ NBody( [tuple([rand(0:4, 3);0]...)], [1.0+rand()], D3 )
      for n = 1:nbasis ]
E1 = [energy(b, at) for b in B3]
E2 = energy( B3, at )
@test E1 ≈ E2

println("4-body")
rcut4 = 2.1 * r0
D4 = Dictionary(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )
nbasis4 = 100
B4 = [ NBody( [tuple(rand(0:3, 7)...)], [1.0+rand()], D4 )
      for n = 1:nbasis4 ]
E1 = [energy(b, at) for b in B4]
E2 = energy( B4, at )
@test E1 ≈ E2

println("5-body")
rcut5 = 1.5 * r0
D5 = Dictionary(TRANSFORM, (:cos, 0.66*rcut5, rcut5) )
nbasis5 = 300
B5 = [ NBody( [tuple(rand(0:5, 11)...)], [1.0+rand()], D5 )
      for n = 1:nbasis5 ]
E1 = [energy(b, at) for b in B5]
E2 = energy( B5, at )
@test E1 ≈ E2

# println("Performance")
# print(" 3-body old:" ); @btime ([energy(b, at) for b in B3])
# print(" 3-body new:" ); @btime energy( B3, at )
# print(" 4-body old:" ); @btime ([energy(b, at) for b in B4])
# print(" 4-body new:" ); @btime energy( $B4, $at )
# print(" 5-body old:" ); @btime ([energy(b, at) for b in B5])
# print(" 5-body new:" ); @btime energy( $B5, $at )


# println("[4] `NBody` gradient-test on simplices")
# for n = [1, 3]
#    V3 = NBody( [tuple([rand(0:3, 3); 0]...) for n = 1:n], 1.0 + rand(n), D3 )
#    for _  = 1:10
#       r = 1.0 + rand(SVector{3,Float64})
#       @test (@D V3(r)) ≈ ForwardDiff.gradient(r_ -> V3(r_), r)
#       print(".")
#    end
# end
#
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
