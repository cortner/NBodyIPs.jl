
using NBodyIPs
using NBodyIPs.Polynomials: ispure

v3 = Val(3)
@show ispure(v3, (1,0,0,0))
@show ispure(v3, (1,1,0,0))
@show ispure(v3, (0,1,0,0))
@show ispure(v3, (0,0,1,0))


v4 = Val(4)
@show ispure(v4, (1,0,0,0,0))
@show ispure(v4, (1,1,0,0,0))
@show ispure(v4, (0,1,0,0,0))
@show ispure(v4, (0,0,1,0,1))


NBodyIPs.gen_tuples(3, 5)
