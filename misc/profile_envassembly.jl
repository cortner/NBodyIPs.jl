info("loading libraries...")
using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP, Base.Test

# using JuLIP.Potentials: evaluate, evaluate_d
# using NBodyIPs: BondLengthDesc, BondAngleDesc, invariants, invariants_d, descriptor,
#                bo2angles, bo2edges
using NBodyIPs.Polys
using NBodyIPs.EnvIPs
using NBodyIPs.EnvIPs: EnvPoly

TRANSFORM = "r -> (2.9/r)^3"
DEGREES = [14, 12, 10]
ENVDEGS = [3, 2, 1]
RCUT = [7.0, 5.80, 4.5]


function random_ip(DT)
   DD = [ DT(TRANSFORM, (:cos, (0.66*rcut), rcut))  for rcut in RCUT ]

   rn = 1.5 * rnn(:W)
   Vn = ("1.0/(1.0 + exp(1.0 / ($rn - r + 1e-2)))", rn)

   # 2-body potential
   BB = [ envpolys(n+1, DD[n], DEGREES[n], Vn, ENVDEGS[n]) for n = 1:3 ]
   B = vcat(BB...)
   c = rand(length(B))
   IP = NBodyIP(B, c)

   return IP, fast(IP)
end

info("Generate random EnvIP...")
IP, IPf = random_ip(BondLengthDesc)

info("Generate faster EnvIP => EnvPoly ...")
V2 = IPf.components[1:4]
V3 = IPf.components[5:7]
V4 = IPf.components[8:9]

IPff = NBodyIP( [ EnvPoly([v.Vr for v in V2], V2[1].Vn, V2[1].str_Vn),
                  EnvPoly([v.Vr for v in V3], V3[1].Vn, V3[1].str_Vn),
                  EnvPoly([v.Vr for v in V4], V4[1].Vn, V4[1].str_Vn) ] )

info("Test Correctness...")
at = rattle!(bulk(:W, cubic=true) * 3, 0.02)
@show energy(IP, at) - energy(IPf, at)
@show energy(IPff, at) - energy(IPf, at)

info("Test evaluation time...")
at = rattle!(bulk(:W, cubic=true) * 8, 0.02)
print("IP  : ")
energy(IP, at); @time energy(IP, at);
print("IPf : ")
energy(IPf, at); @time energy(IPf, at);
print("IPff: ")
energy(IPff, at); @time energy(IPff, at);

info("For comparison some simpler potentials:")
IP0 = NBodyIP( [V2[1], V3[1], V4[1]] )
IP1 = NBodyIP( [V2[1].Vr, V3[1].Vr, V4[1].Vr] )
print("EnvIP deg = 0 : ")
energy(IP0, at); @time energy(IP0, at);
print(" 4B IP        : ")
energy(IP1, at); @time energy(IP1, at);








##
quit()
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
