
using NBodyIPs
import Base: hash
using NBodyIPs: OneBody
using NBodyIPs.Polys: NBPoly
##
println("generate some basis functions")
rcuts = [9.2, 6.2, 4.5]   # 9.2
TRANSFORM = "r -> (2.9/r)^3"
CUTOFF = ["(:cos, $(0.66*rcut), $(rcut))" for rcut in rcuts]
B1 = [OneBody(1.0)]
B2 = blpolys(2, TRANSFORM, CUTOFF[1], 12)
B3 = blpolys(3, TRANSFORM, CUTOFF[2], 10)
B4 = blpolys(4, TRANSFORM, CUTOFF[3], 6)
B = [B1; B2; B3; B4]

hash(typeof(B1))
b = B2[1]

hash((NBPoly, B2[1].D.transform.id, B2[1].D.cutoff))
hash((NBPoly, B2[2].D.transform.id, B2[2].D.cutoff))

function hash(::Val{:basis}, V::NBPoly)
   return hash((NBPoly, hash(V.D)))
end


using JuLIP
JuLIP.Potentials.SplineCutoff
