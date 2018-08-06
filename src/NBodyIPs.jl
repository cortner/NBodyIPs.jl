
# TODO:
# * energy and forces for ASEAtoms ????
#   -> maybe define these for AbstractAtoms, then remove the
#      conversion from ASE???
# * poly_regularise
# * Remove NBPoly{1} altogether and replace it with an OneBody?
#   that is available for all sub-modules?
# * reformulate the basis assembly for EnvIP to allow arbitrary inner functions?

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

# # bond-angle invariants
# include("bainvariants.jl")

include("descriptors.jl")

# Polynomials of an invariant coordinate system
include("polys.jl")
import NBodyIPs.Polys.bl_basis
export bl_basis


# include("environ.jl")
# import NBodyIPs.EnvIPs: envbl_basis
# export envbl_basis

end # module
