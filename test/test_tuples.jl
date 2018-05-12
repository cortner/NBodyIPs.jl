using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

for args in [ (2, 12), (3,8), (4, 6), (4, 10) ]
   Aold = NBodyIPs.Polys.gen_tuples(args...)
   Anew = NBodyIPs.Polys.gen_tuples_new2(args...)
   println("args = $args")
   @show sort(Aold) == sort(Anew)
   @test sort(Aold) == sort(Anew)
   @btime NBodyIPs.Polys.gen_tuples($args[1], $args[2])
   @btime NBodyIPs.Polys.gen_tuples_new2($args[1], $args[2])
end

println(" Runtime for 5-body:")
for args in [ (5,6), (5,7), (5,8) ]
   println("args = $args")
   @btime NBodyIPs.Polys.gen_tuples_new2($args[1], $args[2])
end
