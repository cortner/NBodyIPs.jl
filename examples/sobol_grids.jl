
using NBodyIPs, StaticArrays, LinearAlgebra
Sob = NBodyIPs.Sobol

r0 = 1.0
r1 = 2.0
x0 = r0 * (@SVector ones(3))
x1 = r1 * (@SVector ones(3))

Xbl = Sob.filtered_sobol(x0, x1, Sob.bl_is_simplex,
                      npoints=1_000, nfailed=100_000)

const J3 = SVector(2, 3)
cart2r(R) = NBodyIPs.edge_lengths(R, J3)
filt(r) = (r0 <= minimum(r)) && (maximum(r) <= r1)
Xcart = Sob.filtered_cart_sobol(r1, 3, cart2r, filt,
                                npoints = 1_000, nfailed = 100_000
                       )
