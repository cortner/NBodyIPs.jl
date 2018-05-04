"""
# `module IO`

Writing and reading of `NBodyIPs`.

## Usage

If `V <: NBodyIP` then
```
NBodyIPs.IO.write(fname, V)
```
will write `V` into a file, while
```
V = NBodyIPs.IO.read(fname)
```
will load a previously saved `NBodyIP`. One can then simply call
`energy(V), forces(V)`, etc, as with any `JuLIP.AbstractCalculator`.

## Restrictions

`NBodyIP` simply stores a list of N-body terms, where each of these
N-body terms must implement `serialize` and `deserialize`. At the moment,
only `NBodyIPs.Polynomials.NBody` is available for this. Currently,
and `NBody` can only be serialised if its dictionary constructed using
strings as parameters, e.g.,
```
Dictionary("@analytic r -> (\$(r0)/r)^3", "(:cos, \$(0.66*rcut), \$rcut)")
```
since this allows the `Dictionary` to store the string from which it is
constructed, which can be easily (de-)serialized.
"""
module IO

import Base: read, write

using JLD
using NBodyIPs: NBodyIP
using NBodyIPs.Polynomials: NBody, Dictionary

function write(fname::AbstractString, IP::NBodyIP)
   @assert fname[end-3:end] == ".jld"
   IPs = serialize.(IP.orders)
   JLD.save(fname, "IP", IPs)
end

function read(fname::AbstractString)
   IPs = JLD.load(fname, "IP")
   return NBodyIP(deserialize.(NBody, IPs))
end

end
