
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


NBodyIPs.gen_tuples(3, 10, purify=false)
NBodyIPs.gen_tuples(3, 10, purify=true)

NBodyIPs.gen_tuples(4, 6, purify=false)
NBodyIPs.gen_tuples(4, 6, purify=true)

NBodyIPs.gen_tuples(4, 8, purify=false)
NBodyIPs.gen_tuples(4, 8, purify=true)

NBodyIPs.gen_tuples(4, 10, purify=false)
NBodyIPs.gen_tuples(4, 10, purify=true)
