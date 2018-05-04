
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
