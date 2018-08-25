using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP

using JuLIP.Potentials: evaluate, evaluate_d
using NBodyIPs: BondLengthDesc, BondAngleDesc, invariants, invariants_d, descriptor,
               bo2angles, bo2edges
using NBodyIPs.Polys


rcut = 7.0
TRANSFORM = "r -> (2.9/r)^3"
CUTOFF3 = (:cos, 0.66*rcut, rcut)
Dbl = BondLengthDesc(TRANSFORM, CUTOFF3)
Dba = BondAngleDesc(TRANSFORM, CUTOFF3)

randri(N, D::BondLengthDesc) = SVector((rand(bo2edges(N)) + 3.0)...)
randri(N, D::BondAngleDesc) = SVector((rand(N-1) + 3.0)...),
                               SVector((3*rand(bo2angles(N))-1.5)...)

println("================")
println(" StNBPoly Tests ")
println("================")

for D in [Dbl, Dba]
   println("---------------------------------------------------")
   println("Testing StaticPolynmials with D = ")
   println("           $(typeof(D))")
   println("---------------------------------------------------")

   # 3-body potential
   B3 = nbpolys(3, D, 10)
   c = rand(length(B3))
   V3 = NBPoly(B3, c, D)
   V3sp = StNBPoly(V3)

   # 4-body potential
   B4 = nbpolys(4, D, 8)
   c = rand(length(B4))
   V4 = NBPoly(B4, c, D)
   V4sp = StNBPoly(V4)

   ri3 = randri(3, D)
   ri4 = randri(4, D)

   println("3-body test:")
   @show errV3 = V3(ri3) - V3sp(ri3)
   println(@test abs(errV3) < 1e-10)
   @show errdV3 = norm((@D V3(ri3)) - (@D V3sp(ri3)))
   println(@test abs(errdV3) < 1e-10)
   print(" evaluate old: "); @btime evaluate($V3, $ri3)
   print(" evaluate new: "); @btime evaluate($V3sp, $ri3)
   print(" gradient old: "); @btime evaluate_d($V3, $ri3)
   print(" gradient new: "); @btime evaluate_d($V3sp, $ri3)

   println("4-body test:")
   @show errV4 = V4(ri4) - V4sp(ri4)
   println(@test abs(errV4) < 1e-8)
   @show errdV4 = norm((@D V4(ri4)) - (@D V4sp(ri4)))
   println(@test abs(errdV4) < 1e-8)
   print(" evaluate old: "); @btime evaluate($V4, $ri4)
   print(" evaluate new: "); @btime evaluate($V4sp, $ri4)
   print(" gradient old: "); @btime evaluate_d($V4, $ri4)
   print(" gradient new: "); @btime evaluate_d($V4sp, $ri4)
end
