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

info("Test Error in Forces")
@show maximum(norm.(forces(IP, at) - forces(IPf, at)))
@show maximum(norm.(forces(IPff, at) - forces(IPf, at)))

info("Test evaluation time for forces...")
at = rattle!(bulk(:W, cubic=true) * 8, 0.02)
print("IP  : ")
forces(IP, at); @time forces(IP, at);
print("IPf  : ")
forces(IPf, at); @time forces(IPf, at);
print("IPff  : ")
forces(IPff, at); @time forces(IPff, at);
print("IP0  : ")
forces(IP0, at); @time forces(IP0, at);
print("IP1  : ")
forces(IP1, at); @time forces(IP1, at);

at = rattle!(bulk(:W, cubic=true) * 7, 0.02)
Profile.clear(); @profile forces(IPff, at);
Profile.print()

at = rattle!(bulk(:W, cubic=true) * 3, 0.02)
V4 = IPff.components[end]
@code_warntype forces(V4, at);
