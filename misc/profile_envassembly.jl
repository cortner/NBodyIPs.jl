@info("loading libraries...")
using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP, Test, LinearAlgebra

using NBodyIPs.Polys
using NBodyIPs.EnvIPs
using NBodyIPs.EnvIPs: EnvPoly

# TRANSFORM = PolyTransform(3, 2.9)
# DEGREES = [14, 12, 10]
# ENVDEGS = [3, 2, 1]
# r0 = rnn(:W)
# RCUT = [3.1*r0, 2.5*r0, 2.1*r0]
# RCUT = [2.6*r0, 2.2*r0, 1.8*r0]

## fast tests
TRANSFORM = PolyTransform(3, 2.9)
DEGREES = [14, 12, 10]
ENVDEGS = [3, 2, 1]
r0 = rnn(:W)
RCUT = [2.6*r0, 2.2*r0, 1.8*r0]


function random_ip(DT)
   DD = [ DT(TRANSFORM, CosCut(rcut-1.0, rcut))  for rcut in RCUT ]
   rn = 1.45 * rnn(:W)
   Vn = ("1.0/(1.0 + exp(1.0 / ($rn - r + 1e-2)))", rn)
   BB = [ envpolys(n+1, DD[n], DEGREES[n], Vn, ENVDEGS[n]) for n = 1:3 ]
   B = vcat(BB...)
   c = rand(length(B))
   IP = NBodyIP(B, c)
   return IP
end

##
@info("Generate random EnvIP...")
IP = random_ip(BondLengthDesc)
IPf = fast(IP)
IPcl = random_ip(ClusterBLDesc)
IPclf = fast(IPcl)

@info("Generate faster EnvIP => EnvPoly ...")
V2 = IPf.components[1:4]
V3 = IPf.components[5:7]
V4 = IPf.components[8:9]
IPff = NBodyIP( [ EnvPoly([v.Vr for v in V2], V2[1].Vn, V2[1].str_Vn),
                  EnvPoly([v.Vr for v in V3], V3[1].Vn, V3[1].str_Vn),
                  EnvPoly([v.Vr for v in V4], V4[1].Vn, V4[1].str_Vn) ] )

V2cl = IPclf.components[1:4]
V3cl = IPclf.components[5:7]
V4cl = IPclf.components[8:9]
IPclff = NBodyIP( [ EnvPoly([v.Vr for v in V2cl], V2cl[1].Vn, V2cl[1].str_Vn),
                    EnvPoly([v.Vr for v in V3cl], V3cl[1].Vn, V3cl[1].str_Vn),
                    EnvPoly([v.Vr for v in V4cl], V4cl[1].Vn, V4cl[1].str_Vn) ] )


@info("Test Correctness...")
at = rattle!(bulk(:W, cubic=true) * 3, 0.02)
@show energy(IP, at)
@show energy(IP, at) - energy(IPf, at)
@show energy(IPff, at) - energy(IPf, at)
@show energy(IPcl, at)
@show energy(IPcl, at) - energy(IPclf, at)
@show energy(IPclff, at) - energy(IPclf, at)

@info("Test evaluation time...")
at = rattle!(bulk(:W, cubic=true) * 8, 0.02)
GC.enable(false)
print("nlist : ")
neighbourlist(at, 3.0); @time begin
      neighbourlist(at, RCUT[1]);
      neighbourlist(at, RCUT[2]);
      neighbourlist(at, RCUT[3])
end 
GC.gc()
print("IP    : ")
energy(IP, at); @time energy(IP, at); GC.gc()
print("IPf   : ")
energy(IPf, at); @time energy(IPf, at);GC.gc()
print("IPff  : ")
energy(IPff, at); @time energy(IPff, at);GC.gc()
print("IPcl  : ")
energy(IPcl, at); @time energy(IPcl, at);GC.gc()
print("IPclf : ")
energy(IPclf, at); @time energy(IPclf, at);GC.gc()
print("IPclff: ")
energy(IPclff, at); @time energy(IPclff, at);GC.gc()
GC.enable(true)

exit()

@info("For comparison some simpler potentials:")
IP0 = NBodyIP( [V2[1], V3[1], V4[1]] )
IP1 = NBodyIP( [V2[1].Vr, V3[1].Vr, V4[1].Vr] )
print("IP0: EnvIP deg = 0 : ")
energy(IP0, at); @time energy(IP0, at);
print("IP1:  4B IP        : ")
energy(IP1, at); @time energy(IP1, at);

@info("Test Error in Forces")
@show maximum(norm.(forces(IP, at) - forces(IPf, at)))
@show maximum(norm.(forces(IPff, at) - forces(IPf, at)))

@info("Test evaluation time for forces...")
at = rattle!(bulk(:W, cubic=true) * 8, 0.02)
GC.enable(false)
print("IP  : ")
forces(IP, at); @time forces(IP, at);
GC.gc()
print("IPf  : ")
forces(IPf, at); @time forces(IPf, at);
GC.gc()
print("IPff  : ")
forces(IPff, at); @time forces(IPff, at);
GC.gc()
print("IP0  : ")
forces(IP0, at); @time forces(IP0, at);
GC.gc()
print("IP1  : ")
forces(IP1, at); @time forces(IP1, at);
GC.enable(true)

# BACKUP PROFILING CODE
# at = rattle!(bulk(:W, cubic=true) * 7, 0.02)
# Profile.clear(); @profile forces(IPff, at);
# Profile.print()
#
# at = rattle!(bulk(:W, cubic=true) * 3, 0.02)
# V4 = IPff.components[end]
# @code_warntype forces(V4, at);
