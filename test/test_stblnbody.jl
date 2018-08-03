using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP
using NBodyIPs.BLPolys
using JuLIP.Potentials: evaluate, evaluate_d

rcut = 7.0
TRANSFORM = "r -> (2.9/r)^3"
CUTOFF3 = (:cos, 0.66*rcut, rcut)
D = BLDictionary(TRANSFORM, CUTOFF3)

# 3-body potential
B3 = bl_basis(3, D, 10)
c = rand(length(B3))
V3 = BLNBody(B3, c, D)
V3sp = StBLNBody(V3)

# 4-body potential
B4 = bl_basis(4, D, 8)
c = rand(length(B4))
V4 = BLNBody(B4, c, D)
V4sp = StBLNBody(V4)

r3 = (@SVector rand(3)) + 3.0
r4 = (@SVector rand(6)) + 3.0

println("3-body test:")
@show errV3 = V3(r3) - V3sp(r3)
println(@test abs(errV3) < 1e-10)
@show errdV3 = norm((@D V3(r3)) - (@D V3sp(r3)))
println(@test abs(errdV3) < 1e-10)
print(" evaluate old: "); @btime evaluate($V3, $r3)
print(" evaluate new: "); @btime evaluate($V3sp, $r3)
print(" gradient old: "); @btime evaluate_d($V3, $r3)
print(" gradient new: "); @btime evaluate_d($V3sp, $r3)

println("4-body test:")
@show errV4 = V4(r4) - V4sp(r4)
println(@test abs(errV4) < 1e-8)
@show errdV4 = norm((@D V4(r4)) - (@D V4sp(r4)))
println(@test abs(errdV4) < 1e-8)
print(" evaluate old: "); @btime evaluate($V4, $r4)
print(" evaluate new: "); @btime evaluate($V4sp, $r4)
print(" gradient old: "); @btime evaluate_d($V4, $r4)
print(" gradient new: "); @btime evaluate_d($V4sp, $r4)
