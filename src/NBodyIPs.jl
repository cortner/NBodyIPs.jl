
# TODO:
#  - energy and forces for ASEAtoms ????

"""
# `NBodyIPs.jl`

Package for specifying interatomic potentials based on the
N-Body expansion (ANOVA, HDMR, ...).

See `NBodyIPFitting` for the associated fitting and testing framework.
"""
module NBodyIPs

using Reexport

@reexport using StaticArrays
@reexport using JuLIP

# two auxiliary functions to make for easier assembly of the code
# TODO: move these somewhere else, or better get rid of them
push_str!(ex::Vector{Expr}, s::String) = push!(ex, parse(s))
append_str!(ex::Vector{Expr}, s::Vector{String}) = append!(ex, parse.(s))


include("fastpolys.jl")

include("invariants.jl")

include("common.jl")


# some generically useful code that
# could be used across different n-body basis function implementations
# TODO: move some codes from Invariants submodule to here
#       or maybe other parts of the package
include("misc.jl")


# describe basis functions in terms of symmetry invariants
include("polynomials.jl")
@reexport using NBodyIPs.Polys

include("environ.jl")
import NBodyIPs.EnvBLs: envbl_basis
export envbl_basis

end # module
