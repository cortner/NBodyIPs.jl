
using NBodyIPs, JuLIP, LinearAlgebra
using NBodyIPs: combinebasis, BASIS

e1 = rand()
e2 = rand()
b1 = OneBody(e1)
b2 = OneBody(e2)
at = bulk(:Cu, cubic=true) * 3
println(@test(energy(b1, at) ≈ length(at) * b1.E0))
println(@test(maximum(norm.(forces(b1, at))) == 0.0))
println(@test(OneBody(2*e1+3*e2) == combinebasis([b1, b2], [2,3])))
println(@test(hash(BASIS(), b1) == hash(BASIS(), b2)))
