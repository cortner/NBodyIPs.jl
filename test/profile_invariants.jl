include("../src/fastpolys.jl")
using FastPolys
using StaticArrays, BenchmarkTools, StaticPolynomials


IS15_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,2,2,2,2,3,3,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,5,5,6,7,6,7,5,5,5,5,6,6,6,7,5,5,7,6,8,9,8,8,9,8,5,5,6,7,6,7,5,5,5,5,6,6,6,7,5,5,7,6,8,9,8,8,9,8,)
IS15_2 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,2,2,3,4,3,4,2,2,3,4,3,4,3,4,3,4,4,4,3,4,3,4,4,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,10,10,9,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,10,10,9,9,10,10,)
IS15_3 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,7,6,7,6,5,5,9,8,10,10,8,9,9,8,10,10,8,9,5,5,6,7,6,7,)
const IS15 = Val((IS15_1, IS15_2, IS15_3))

function mypoly(x)
   x2 = x.*x
   return FastPolys.fpoly((x2,x2,x), IS15)
end

function mypoly_d(x)
   x2 = x.*x
   return FastPolys.fpoly_d((x2,x2,x),(2*x,2*x,(@SVector ones(10))), IS15)
end


C = ones(length(IS15_1))
E = zeros(Int, length(IS15_1), 10)
for (j, (i1,i2,i3)) in enumerate(zip(IS15_1, IS15_2, IS15_3))
   E[j, i1] = 2
   E[j, i2] = 2
   E[j, i3] = 1
end

SP = StaticPolynomials.Polynomial(C, E')

x = @SVector rand(10)

println("fpoly with overhead")
@btime mypoly($x)
@btime mypoly_d($x)
println("fpoly without overhead")
@btime FastPolys.fpoly($((x.*x,x.*x,x)), $IS15)
@btime FastPolys.fpoly_d($((x.*x,x.*x,x)), $((2*x,2*x,(@SVector ones(10)))), $IS15)
println("StaticPolynomials")
@btime StaticPolynomials.evaluate($SP, $x)
@btime StaticPolynomials.gradient($SP, $x)
