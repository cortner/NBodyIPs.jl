#
# This set of tests indicates that maybe we should go back
# to using StaticPolynomials everywhere. Their performance has
# improved very much
#


include("../src/fastpolys.jl")
using FastPolys
using StaticArrays, BenchmarkTools, StaticPolynomials

IS17_1 = (1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10)
IS17_2 = (2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 8, 9, 1, 2, 4, 6, 8, 10, 1, 2, 3, 7, 9, 10, 1, 2, 6, 7, 8, 9, 1, 3, 5, 7, 8, 10, 1, 4, 5, 6, 9, 10, 2,3, 5, 6, 9, 10, 2, 4, 5, 7, 8, 10, 3, 4, 6, 7, 8, 9)
const IS17 = Val((IS17_1, IS17_2))

IS15_1 = (1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,2,2,2,2,3,3,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,5,5,6,7,6,7,5,5,5,5,6,6,6,7,5,5,7,6,8,9,8,8,9,8,5,5,6,7,6,7,5,5,5,5,6,6,6,7,5,5,7,6,8,9,8,8,9,8,)
IS15_2 = (5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,2,2,3,4,3,4,2,2,3,4,3,4,3,4,3,4,4,4,3,4,3,4,4,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,10,10,9,9,10,10,8,9,8,9,10,10,6,7,6,7,7,7,8,9,8,9,10,10,10,10,9,9,10,10,)
IS15_3 = (3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,7,6,7,6,5,5,9,8,10,10,8,9,9,8,10,10,8,9,5,5,6,7,6,7,)
const IS15 = Val((IS15_1, IS15_2, IS15_3))

P6_1 = (1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,1,1,1,2,2,3,3,2,2,1,1,1,4,3,2,4,3,4,3,2,2,1,1,1,)
P6_2 = (2,2,2,2,3,3,2,2,3,3,3,3,5,5,5,5,6,6,5,5,6,7,6,7,5,5,6,6,5,5,7,6,7,6,5,5,4,3,2,4,3,4,4,3,4,2,3,4,5,5,5,5,5,6,7,6,7,5,6,7,)
P6_3 = (3,4,3,4,4,4,3,4,4,4,4,4,6,7,6,7,7,7,8,9,8,8,9,8,7,6,7,7,6,7,8,8,8,9,8,9,5,5,6,5,5,6,7,6,7,8,8,9,6,6,6,8,8,8,9,8,9,8,8,9,)
P6_4 = (10,10,9,8,9,8,7,6,5,7,6,5,10,10,9,8,9,8,10,10,9,9,10,10,8,9,8,9,10,10,9,9,10,10,10,10,6,7,7,8,9,8,9,10,10,9,10,10,7,7,7,9,9,10,10,10,10,9,10,10,)
const P6 = Val((P6_1, P6_2, P6_3, P6_4))


function exps17()
   E = zeros(Int, length(IS17_1), 10)
   for (j, (i1,i2)) in enumerate(zip(IS17_1, IS17_2))
      E[j, i1] = 4
      E[j, i2] = 2
   end
   return E
end


function exps15()
   E = zeros(Int, length(IS15_1), 10)
   for (j, (i1,i2,i3)) in enumerate(zip(IS15_1, IS15_2, IS15_3))
      E[j, i1] = 2
      E[j, i2] = 2
      E[j, i3] = 1
   end
   return E
end

function exps6()
   E = zeros(Int, length(P6_1), 10)
   for (j, (i1,i2,i3,i4)) in enumerate(zip(P6_1, P6_2, P6_3, P6_4))
      E[j, i1] = 1
      E[j, i2] = 1
      E[j, i3] = 1
      E[j, i4] = 1
   end
   return E
end


x = @SVector rand(10)
x2 = x .* x
x3 = x2 .* x
x4 = x2 .* x2
dx4 = 3 * x3
o = @SVector ones(10)


println("Test 1: 2nd-order")
C17 = rand(length(IS17_1))
E17 = exps17()
SP17 = StaticPolynomials.Polynomial(C17, E17')

println("fpoly without overhead")
@btime FastPolys.fpoly($((x4,x2)), $IS17)
@btime FastPolys.fpoly_d($((x2,x4)), $((2*x,dx4)), $IS17)
println("StaticPolynomials")
@btime StaticPolynomials.evaluate($SP17, $x)
@btime StaticPolynomials.gradient($SP17, $x)

println("Test 2: 3-rd order")

C15 = ones(length(IS15_1))
E15 = exps15()
SP15 = StaticPolynomials.Polynomial(C15, E15')

println("fpoly without overhead")
@btime FastPolys.fpoly($((x.*x,x.*x,x)), $IS15)
@btime FastPolys.fpoly_d($((x.*x,x.*x,x)), $((2*x,2*x,o)), $IS15)
println("StaticPolynomials")
@btime StaticPolynomials.evaluate($SP15, $x)
@btime StaticPolynomials.gradient($SP15, $x)


println("Test 3: 4-th order")

C6 = ones(length(P6_1))
E6 = exps6()
SP6 = StaticPolynomials.Polynomial(C6, E6')

println("fpoly without overhead")
@btime FastPolys.fpoly($((x,x,x,x)), $P6)
@btime FastPolys.fpoly_d($((x,x,x,x)), $((o,o,o,o)), $P6)
println("StaticPolynomials")
@btime StaticPolynomials.evaluate($SP6, $x)
@btime StaticPolynomials.gradient($SP6, $x)
