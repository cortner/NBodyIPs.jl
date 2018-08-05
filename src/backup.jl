# ------- distribution functions on the invariants --------------

function rdf(at::Atoms, rcut, transform=identity)
   if transform != identity
      return idf(2, at, rcut, transform)[1][1]
   end
   nlist = neighbourlist(at, rcut)
   return nlist.r
end

"""
`idf(N::Integer, at, rcut, transform)`

invariants distribution function => accumulates all the invariants values
arising during an N-body assembly over at.
"""
idf(N::Integer, at, rcut, transform) =
      _idf(Val(N), Val(bo2edges(N)), at, rcut, transform)

function _idf(valN::Val{N}, valM::Val{M}, at::Atoms{T}, rcut::T,
              transform) where {N, M, T}
   # compute invariants vectors to learn how many there are
   x = rand(SVector{M,T})
   I1, I2 = invariants(x)

   I1acc = [ T[] for n = 1:length(I1) ]
   I2acc = [ T[] for n = 1:length(I2) ]
   Iacc = (I1acc, I2acc)

   for (i, j, r, R) in sites(at, rcut)
      eval_site_nbody!(valN, R, rcut,
                  (Iacc, s, _1, _2, _3) -> idf_accumulator(Iacc, transform.(s)),
                  Iacc, nothing)
   end

   return I1acc, I2acc
end

function idf_accumulator(Iacc, x)
   I1, I2 = invariants(x)
   for n = 1:length(I1)
      push!(Iacc[1][n], I1[n])
   end
   for n = 1:length(I2)
      push!(Iacc[2][n], I2[n])
   end
   return Iacc
end



# testing the custom JLD serialization

using JuLIP, NBodyIPs, DataFrames, Plots, JLD
using NBodyIPs.Data: Dat

include(homedir() * "/Dropbox/PIBmat/Ti_DFTB_Data/Ti.jl")
data = Ti.read_Ti(; exclude = ["wire", "surface", "omega", "hcp"])

# reference energy
B1 = [ NBody(Ti.get_E0()) ]

# generate long-range 2B basis
r0 = rnn(:Ti)
TLONG = "@analytic r -> exp( - 2.5 * (r/$r0 - 1))"
TSHORT = "@analytic r -> ($r0/r)^8"
CUT2 = "(:cos, 5.5, 7.5)"
B2 = [ gen_basis(2, Dictionary(TLONG, CUT2), 9);
       gen_basis(2, Dictionary(TSHORT, CUT2), 6) ]

# 3B BASIS
TRANS3 = "@analytic r -> exp( - 3.0 * (r/$r0 - 1))"
CUT3 = "(:cos, 5.0, 6.5)"
B3 = gen_basis(3, Dictionary(TRANS3, CUT3), 8)

B = [B1; B2; B3]
@show length(B)

lsq = kron(data[1:10], B)
JLD.save("temp.jld", "d1", lsq)
d1 = JLD.load("temp.jld", "d1")


JLD.@load "temp.jld" lsq



"""
`struct BLDictionary` : specifies all details about the basis functions

## Constructor

`BLDictionary(ftrans, fcut)`, where `ftrans`, `fcut` specify in one of several
ways how the dictionary is defined. If
the radius is not provided then `BLDictionary` will try to infer it from the
`fcut` argument. E.g.,
```julia
D = BLDictionary("r -> 1/r", (:cos, rcut1, rcut2))
```

Known symbols for the cutoff are
```
[:cos, :sw, :spline, :square, :cos2s]
```

## Developer Doc: Methods associated with a `D::BLDictionary`:

* `invariants`, `invariants_d`, `invariants_ed`: compute invariants and jacobian
   in transformed coordinates defined by `D.transform`
* `evaluate`, `evaluate_d`: evaluate the (univariate) basis
   function associated with this dictionary; at the moment only
   standard polynomials are admitted
* `fcut, fcut_d`: evulate the cut-off function and its derivative / gradient
   when interpreted as a product of cut-off functions
"""
struct BLDictionary{TT, TC}
   transform::TT
   cutoff::TC
end

function ==(D1::BLDictionary, D2::BLDictionary)
   return (D1.transform == D2.transform) && (D1.cutoff == D2.cutoff)
end

Dict(D::BLDictionary) = Dict("__id__" => "BLDictionary",
                             "transform" => Dict(D.transform),
                             "cutoff" => Dict(D.cutoff))

BLDictionary(D::Dict) = BLDictionary( SpaceTransform(D["transform"]),
                                      Cutoff(D["cutoff"]) )

Base.convert(::Val{:BLDictionary}, D::Dict) = BLDictionary(D)


BLDictionary(transform::String, cutoff::Union{String, Tuple}) =
         BLDictionary(SpaceTransform(transform), Cutoff(cutoff))


@inline transform(D::BLDictionary, r::Number) = D.transform.f(r)
@inline transform_d(D::BLDictionary, r::Number) = D.transform.f_d(r)
@inline fcut(D::BLDictionary, r::Number) = D.cutoff.f(r)
@inline fcut_d(D::BLDictionary, r::Number) = D.cutoff.f_d(r)
@inline cutoff(D::BLDictionary) = D.cutoff.rcut
