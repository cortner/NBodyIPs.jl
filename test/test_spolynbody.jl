using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP
using NBodyIPs.Polys.SPolyNBody
using JuLIP.Potentials: evaluate, evaluate_d


rcut = 7.0
TRANSFORM = (@analytic r -> (2.9/r)^3)
CUTOFF3 = (:cos, 0.66*rcut, rcut)
D = Dictionary(TRANSFORM, CUTOFF3)

# 3-body potential
B3 = gen_basis(3, D, 10)
c = rand(length(B3))
V3 = NBody(B3, c, D)
V3sp = SPolyNBody(V3)

# 4-body potential
B4 = gen_basis(4, D, 8)
c = rand(length(B4))
V4 = NBody(B4, c, D)
V4sp = SPolyNBody(V4)

r3 = (@SVector rand(3)) + 3.0
r4 = (@SVector rand(6)) + 3.0


println("3-body test:")
@show V3(r3) - V3sp(r3)
@show norm((@D V3(r3)) - (@D V3sp(r3)))
print(" evaluate old: "); @btime evaluate($V3, $r3)
print(" evaluate new: "); @btime evaluate($V3sp, $r3)
print(" gradient old: "); @btime evaluate_d($V3, $r3)
print(" gradient new: "); @btime evaluate_d($V3sp, $r3)

println("4-body test:")
@show V4(r4) - V4sp(r4)
@show norm((@D V4(r4)) - (@D V4sp(r4)))
print(" evaluate old: "); @btime evaluate($V4, $r4)
print(" evaluate new: "); @btime evaluate($V4sp, $r4)
print(" gradient old: "); @btime evaluate_d($V4, $r4)
print(" gradient new: "); @btime evaluate_d($V4sp, $r4)
