using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

for args in [ (2, 12), (3,8), (4, 6), (4, 10), (5,6), (5,7), (5,8) ]
   println("args = $args")
   @btime NBodyIPs.Polys.gen_tuples($args[1], $args[2])
end
