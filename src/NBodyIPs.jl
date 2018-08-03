
# TODO:
#  - energy and forces for ASEAtoms ????
#   -> maybe define these for AbstractAtoms, then remove the
#      conversion from ASE???

__precompile__()

"""
# `NBodyIPs.jl`

Package for specifying interatomic potentials based on the
N-Body expansion (ANOVA, HDMR, ...).

See `NBodyIPFitting` for the associated fitting and testing framework.
"""
module NBodyIPs

# generic types and function prototypes
include("common.jl")

# the machinery for evaluating the invariants as fast as possible
include("fastpolys.jl")

# Bond-length polynomials
include("blpolys.jl")

# include("bapolys.jl")

# some generically useful code that
# could be used across different n-body basis function implementations
# TODO: move some codes from Invariants submodule to here
#       or maybe other parts of the package
# include("misc.jl")


# include("environ.jl")
# import NBodyIPs.EnvBLs: envbl_basis
# export envbl_basis

end # module
