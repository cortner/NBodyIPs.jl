using NBodyIPs, StaticArrays, BenchmarkTools
using NBodyIPs.Polys.SPolyNBody
using JuLIP.Potentials: evaluate, evaluate_d


rcut = 7.0
TRANSFORM = (@analytic r -> (2.9/r)^3)
CUTOFF3 = (:cos, 0.66*rcut, rcut)
D = Dictionary(TRANSFORM, CUTOFF3)
B3 = gen_basis(3, D, 6)
c = rand(length(B3))

V3 = NBody(B3, c, D)
V3sp = SPolyNBody(V3)

r = (@SVector rand(3)) + 3.0
@show V3(r), V3sp(r)

print(" evaluate old: "); @btime evaluate($V3, $r)
print(" evaluate new: "); @btime evaluate($V3sp, $r)
print(" gradient old: "); @btime evaluate_d($V3, $r)
print(" gradient new: "); @btime evaluate_d($V3sp, $r)
