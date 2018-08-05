

module BAPolys


"""
`module Polys`

The exported symbols are
* `BLDictionary`: collects all the information about a basis
* `NBody`: an N-body function wrapped into a JuLIP calculator
* `poly_basis` : generate a basis of N-body functions

## Usage

## Notation

* `N` : N-body
* `M` : number of edges, i.e., M = N (N-1)/2
* `K` : length of the tuples defining polynomials, K = M+1
"""
module Polys

import StaticPolynomials
using StaticArrays

using NBodyIPs: NBodyFunction,
                SpaceTransform,
                Cutoff

import Base:              length,
                          Dict,
                          ==
import JuLIP:             cutoff
import JuLIP.Potentials:  @pot
import NBodyIPs:          bodyorder,
                          fast,
                          degree,
                          combinebasis

const Tup{M} = NTuple{M, Int}
const VecTup{M} = Vector{NTuple{M, Int}}

export NBPoly,
       StNBPoly,
       BLDictionary,
       bl_basis


# import the raw invariant polynomials
include("bainvariants.jl")
using NBodyIPs.BAPolys.BAInvariants: invariants,
                                     invariants_d,
                                     invariants_ed



end
