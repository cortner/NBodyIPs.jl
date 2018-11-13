
using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP, Base.Test

# using JuLIP.Potentials: evaluate, evaluate_d
# using NBodyIPs: BondLengthDesc, BondAngleDesc, invariants, invariants_d, descriptor,
#                bo2angles, bo2edges
using NBodyIPs.Polys

rcut = 7.0
TRANSFORM = "r -> (2.9/r)^3"
DEGREES = [14, 10, 8]   # [16, 12, 10]
RCUT = [8.0, 6.5, 4.8]

# test correctness first

at = bulk(:W, cubic=true) * 10


for DT in [BondLengthDesc, ClusterBLDesc]
   println("-----------------------------------------------------------")
   println(" Testing StaticPolynomials with D = $(DT)")
   println("-----------------------------------------------------------")

   DD = [ DT(TRANSFORM, (:cos, (0.66*rcut), rcut))  for rcut in RCUT ]

   # 2-body potential
   B2 = nbpolys(2, DD[1], DEGREES[1])
   c = rand(length(B2))
   V2 = NBPoly(B2, c, DD[1])
   V2sp = StNBPoly(V2)

   # 3-body potential
   B3 = nbpolys(3, DD[2], DEGREES[2])
   c = rand(length(B3))
   V3 = NBPoly(B3, c, DD[2])
   V3sp = StNBPoly(V3)

   # 4-body potential
   B4 = nbpolys(4, DD[3], DEGREES[3])
   c = rand(length(B4))
   V4 = NBPoly(B4, c, DD[3])
   V4sp = StNBPoly(V4)

   IP = NBodyIP([V2, V3, V4])
   IPf = NBodyIP([V2sp, V3sp, V4sp])

   println("Neighbourlist")
   @time neighbourlist(at, cutoff(V4))
   @time neighbourlist(at, cutoff(V4))
   println("Dynamic Polynomials")
   @time energy(IP, at)
   @time energy(IP, at)
   println("Static Polynomials")
   @time energy(IPf, at)
   @time energy(IPf, at)
end
