
# TODO:
# * energy and forces for ASEAtoms ????
#   -> maybe define these for AbstractAtoms, then remove the
#      conversion from ASE???
# * poly_regularise
# * Remove NBPoly{1} altogether and replace it with an OneBody?
#   that is available for all sub-modules?

__precompile__()

"""
# `NBodyIPs.jl`

Package for specifying interatomic potentials based on the
N-Body expansion (ANOVA, HDMR, ...).

See `NBodyIPFitting` for the associated fitting and testing framework.
"""
module NBodyIPs

using Reexport

# generic types and function prototypes
include("common.jl")

# the machinery for evaluating the invariants as fast as possible
include("fastpolys.jl")

# bond-length invariants
include("blinvariants.jl")

# bond-angle invariants
include("bainvariants.jl")

include("descriptors.jl")

# Polynomials of an invariant coordinate system
include("polys.jl")
import NBodyIPs.Polys: blpolys # , bapolys
export blpolys # , bapolys


include("environ.jl")
import NBodyIPs.EnvIPs: envblpolys # , envbapolys
export envblpolys # , envbapolys

end # module
