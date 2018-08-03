
# TODO:
# * energy and forces for ASEAtoms ????
#   -> maybe define these for AbstractAtoms, then remove the
#      conversion from ASE???
# * poly_regularise
# * environ
# * Remove BLNBody{1} altogether and replace it with an OneBody?
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

# Bond-length polynomials
include("blpolys.jl")
import NBodyIPs.BLPolys.bl_basis
export bl_basis

# include("bapolys.jl")

# include("environ.jl")
# import NBodyIPs.EnvBLs: envbl_basis
# export envbl_basis

end # module
