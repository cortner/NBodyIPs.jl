
using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP, Test, Profile

# using JuLIP.Potentials: evaluate, evaluate_d
# using NBodyIPs: BondLengthDesc, BondAngleDesc, invariants, invariants_d, descriptor,
#                bo2angles, bo2edges
using NBodyIPs.Polys

TRANSFORM = PolyTransform(3, 2.9) # "r -> (2.9/r)^3"
DEGREES = [18, 14, 12, 10]
RCUT = [7.0, 5.80, 4.5, 4.1]

function random_ip(DT)
   DD = [ DT(TRANSFORM, CosCut(rcut-1.0, rcut))  for rcut in RCUT ]

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

   # 5-body potential
   B5 = nbpolys(5, DD[4], DEGREES[4])
   c = rand(length(B5))
   V5 = NBPoly(B5, c, DD[4])
   V5sp = StNBPoly(V5)

   IP = NBodyIP([V2, V3, V4, V5])
   IPf = NBodyIP([V2sp, V3sp, V4sp, V5sp])

   return IP, IPf
end

_, IPf = random_ip(BondLengthDesc)
at = rattle!(bulk(:W, cubic=true) * 3, 0.01)
energy(IPf, at)
forces(IPf, at)
at = rattle!(bulk(:W, cubic=true) * 10, 0.01)
@info("Energy")
@time energy(IPf, at)
@time energy(IPf, at)
@info("Forces")
@time forces(IPf, at)
@time forces(IPf, at)

# @profile forces(IPf, at)
# Profile.print()

# F64WRAP
#   1.172545 seconds (21.05 M allocations: 685.697 MiB, 10.16% gc time)
#   1.182909 seconds (21.05 M allocations: 685.697 MiB, 9.84% gc time)
#   1.722956 seconds (33.95 M allocations: 729.967 MiB, 7.34% gc time)
#   1.799147 seconds (33.95 M allocations: 729.967 MiB, 10.56% gc time)

# ORIGINAL
# 81c72ac038fc913362498daabd3c10b7c454aead
#   0.876188 seconds (6.16 M allocations: 458.036 MiB, 7.11% gc time)
#   1.105670 seconds (4.16 M allocations: 275.011 MiB, 3.40% gc time)
#
# REWRITE
#  0.905021 seconds (6.16 M allocations: 458.402 MiB, 7.75% gc time)
#  1.163815 seconds (4.16 M allocations: 275.378 MiB, 3.88% gc time)


# test correctness first
IP_bl, _ = random_ip(BondLengthDesc)
IP_cl, _ = random_ip(ClusterBLDesc)
for (n, (Vbl, Vcl)) in enumerate(zip(IP_bl.components, IP_cl.components))
   Vcl.c[:] = Vbl.c[:] * (n+1)
end

at = rattle!(bulk(:W, cubic=true) * 10, 0.01)
# set_pbc!(at, false)

V2_bl = IP_bl.components[1]
V2_cl = IP_cl.components[1]
@test energy(V2_bl, at) ≈ energy(V2_cl, at)

V3_bl = IP_bl.components[2]
V3_cl = IP_cl.components[2]
@test energy(V3_bl, at) ≈ energy(V3_cl, at)

V4_bl = IP_bl.components[3]
V4_cl = IP_cl.components[3]
@test energy(V4_bl, at) ≈ energy(V4_cl, at)



for DT in [BondLengthDesc, ClusterBLDesc]
   println("-----------------------------------------------------------")
   println(" Testing StaticPolynomials with D = $(DT)")
   println("-----------------------------------------------------------")

   DD = [ DT(TRANSFORM, CosCut(rcut-1, rcut))  for rcut in RCUT ]

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
   @time neighbourlist(at, cutoff(V2_bl))
   @time neighbourlist(at, cutoff(V2_bl))
   println("Dynamic Polynomials")
   @time energy(IP, at)
   @time energy(IP, at)
   @time forces(IP, at)
   @time forces(IP, at)
   println("Static Polynomials")
   @time energy(IPf, at)
   @time energy(IPf, at)
   @time forces(IPf, at)
   @time forces(IPf, at)
end
