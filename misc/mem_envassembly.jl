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

at = rattle!(bulk(:W, cubic=true) * 7, 0.02)
forces(IPff, at);
Profile.clear(); @profile forces(IPff, at);
Profile.print()
