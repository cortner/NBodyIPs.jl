using NBodyIPs, JuLIP, BenchmarkTools, StaticArrays
using Base.Test

using NBodyIPs: ClusterBLDesc, transform, transform_d, fcut, fcut_d,
                invariants, invariants_d, invariants_ed
using NBodyIPs.Polys: NBPoly


info("Setting up the test systems ...")
r0 = rnn(:W)
# TRANSFORM = "r -> exp( - 3 * ((r/$r0) - 1))"
TRANSFORM = "r -> 1/r"
rcut2 = 3.1 * r0
D2 = ClusterBLDesc(TRANSFORM, (:cos, 0.66*rcut2, rcut2) )
rcut3 = 2.3 * r0
D3 = ClusterBLDesc(TRANSFORM, (:cos, 0.66*rcut3, rcut3) )
rcut4 = 1.9 * r0
D4 = ClusterBLDesc(TRANSFORM, (:cos, 0.66*rcut4, rcut4) )

DD = [nothing, D2, D3, D4]


# TEST 1:
B2 = nbpolys(2, D2, 23)
c2 = rand(length(B2))
B3 = nbpolys(3, D3, 12)
c3 = rand(length(B3))
B4 = nbpolys(4, D4, 8)
c4 = rand(length(B4))

B4l = nbpolys(4, D4, 14)
c4l = rand(length(B4l))


IP2 = NBodyIP(B2, c2)
IP3 = NBodyIP([B2;B3], [c2;c3])
IP4 = NBodyIP([B2;B3;B4], [c2;c3;c4])
IP4l = NBodyIP([B2;B3;B4l], [c2;c3;c4l])

IP2f = fast(IP2)
IP3f = fast(IP3)
IP4f = fast(IP4)
IP4lf = fast(IP4l)

at = bulk(:W, cubic=true)*4
@btime forces($IP2f, $at)
@btime forces($IP3f, $at)
@btime forces($IP4f, $at)
@btime forces($IP4lf, $at)

# println((@elapsed forces(IP2f, at)) / length(at) * 1_000)
# println((@elapsed forces(IP3f, at)) / length(at) * 1_000)
# println((@elapsed forces(IP4f, at)) / length(at) * 1_000)
# println((@elapsed forces(IP4lf, at)) / length(at) * 1_000)


println(9.7490 / length(at))
println(20.673 / length(at))
println(28.727 / length(at))
println(37.033 / length(at))
