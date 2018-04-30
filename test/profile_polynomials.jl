
using DynamicPolynomials: @polyvar
using StaticPolynomials, StaticArrays, BenchmarkTools

f3(r::SVector{3}) = SVector(
         r[1] + r[2] + r[3],
         r[1]*r[2] + r[2]*r[3] + r[1]*r[3],
         r[1]*r[2]*r[3] )

Df3(r::SVector{3}) =
   @SMatrix [ 1.0 1.0 1.0;
             r[2]*r[3] r[1]*r[3] r[1]*r[2];
             r[2]*r[3] r[1]*r[3] r[1]*r[2] ]


f4(r::SVector{4}) = SVector(
         r[1] + r[2] + r[3] + r[4],
         r[1]*r[2] + r[1]*r[3] + r[1]*r[4] + r[2]*r[3] + r[2]*r[4] + r[3]*r[4],
         r[1]*r[2]*r[3] + r[1]*r[2]*r[4] + r[1]*r[3]*r[4] + r[2]*r[3]*r[4],
         r[1]*r[2]*r[3]*r[4] )

Df4(r::SVector{4}) =
   @SMatrix [ 1.0 1.0 1.0 1.0;
         r[2]+r[3]+r[4]  r[1]+r[3]+r[4]  r[1]+r[2]+r[4]  r[1]+r[2]+r[3];
         r[2]*r[3]+r[2]*r[4]+r[3]*r[4]  r[1]*r[3]+r[1]*r[4]+r[3]*r[4]  r[1]*r[2]+r[1]*r[4]+r[2]*r[4]  r[1]*r[2]+r[1]*r[3]+r[2]*r[3];
         r[2]*r[3]*r[4]  r[1]*r[3]*r[4]  r[1]*r[2]*r[4]  r[1]*r[2]*r[3] ]

@polyvar r1 r2 r3 r4

P3 = system([   r1 + r2 + r3,
                  r1*r2 + r2*r3 + r1*r3,
                  r1*r2*r3 ])


P4 = system([ r1 + r2 + r3 + r4,
              r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4,
              r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4,
              r1*r2*r3*r4   ])


r_3 = @SVector rand(3)
r_4 = @SVector rand(4)

cfg3 = ForwardDiff.JacobianConfig(f3, r_3)
cfg4 = ForwardDiff.JacobianConfig(f4, r_4)

@btime f3($r_3)    # 2.634 ns (0 allocations: 0 bytes)
@btime StaticPolynomials.evaluate($P3, $r_3)   # 17.045 ns (0 allocations: 0 bytes)
@btime f4($r_4)
@btime StaticPolynomials.evaluate($P4, $r_4)   # 17.045 ns (0 allocations: 0 bytes)

@btime Df3($r_3)  #  3.701 ns
@btime StaticPolynomials.jacobian($P3, $r_3) # 26.666 ns
@btime ForwardDiff.jacobian($f3, $r_3, $cfg3)   # 80.462 ns
@btime Df4($r_4)   #  12.988 ns
@btime StaticPolynomials.jacobian($P4, $r_4)  #  51.952 ns
@btime ForwardDiff.jacobian($f4, $r_4, $cfg4)   #   199.193 ns


@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10

q10 = Polynomial(
   48*(x2^3 + x3^3 + x4^3 + x5^3 + x6^3 + x7^3 + x8^3 + x9^3 + x10^3) +
   72*(x1^2*(x2 + x3 + x4 + x5 + x7 + x8) +
       x2^2*(x1 + x3 + x4 + x6 + x7 + x9) +
       x3^2*(x1 + x2 + x5 + x6 + x8 + x9) +
       x4^2*(x1 + x2 + x5 + x6 + x7 + x10) +
       x5^2*(x1 + x3 + x4 + x6 + x8 + x10) +
       x6^2*(x2 + x3 + x4 + x5 + x9 + x10) +
       x7^2*(x1 + x2 + x4 + x8 + x9 + x10) +
       x8^2*(x1 + x3 + x5 + x7 + x9 + x10) +
       x9^2*(x2 + x3 + x6 + x7 + x8 + x10) +
       x10^2*(x4 + x5 + x6 + x7 + x8 + x9)) +
   144*(x1*x2*(x4 + x7) +
        x1*x3*(x5 + x8) +
        x2*x3*(x6+x9) +
        (x1+x3)*x5*x8 +
        (x1+x2)*x4*x7 +
        (x2+x3)*x6*x9 +
        x4*x10*(x5+x6) +
        x5*x6*(x4+x10) +
        x7*x9*(x8+x10) +
        (x7+x9)*x8*x10)
   )

p10 = Polynomial(
   48*x1^3 + 72*x1^2*x2 + 72*x1^2*x3 + 72*x1^2*x4 + 72*x1^2*x5 + 72*x1^2*x7 +
   72*x1^2*x8 + 72*x1*x2^2 + 144*x1*x2*x4 + 144*x1*x2*x7 + 72*x1*x3^2 +
   144*x1*x3*x5 + 144*x1*x3*x8 + 72*x1*x4^2 + 144*x1*x4*x7 + 72*x1*x5^2 +
   144*x1*x5*x8 + 72*x1*x7^2 + 72*x1*x8^2 + 48*x2^3 + 72*x2^2*x3 +
   72*x2^2*x4 + 72*x2^2*x6 + 72*x2^2*x7 + 72*x2^2*x9 + 72*x2*x3^2 +
   144*x2*x3*x6 + 144*x2*x3*x9 + 72*x2*x4^2 + 144*x2*x4*x7 + 72*x2*x6^2 +
   144*x2*x6*x9 + 72*x2*x7^2 + 72*x2*x9^2 + 48*x3^3 + 72*x3^2*x5 +
   72*x3^2*x6 + 72*x3^2*x8 + 72*x3^2*x9 + 72*x3*x5^2 + 144*x3*x5*x8 +
   72*x3*x6^2 + 144*x3*x6*x9 + 72*x3*x8^2 + 72*x3*x9^2 + 48*x4^3 +
   72*x4^2*x5 + 72*x4^2*x6 + 72*x4^2*x7 + 72*x4^2*x10 + 72*x4*x5^2 +
   144*x4*x5*x6 + 144*x4*x5*x10 + 72*x4*x6^2 + 144*x4*x6*x10 + 72*x4*x7^2 +
   72*x4*x10^2 + 48*x5^3 + 72*x5^2*x6 + 72*x5^2*x8 + 72*x5^2*x10 +
   72*x5*x6^2 + 144*x5*x6*x10 + 72*x5*x8^2 + 72*x5*x10^2 + 48*x6^3 +
   72*x6^2*x9 + 72*x6^2*x10 + 72*x6*x9^2 + 72*x6*x10^2 + 48*x7^3 +
   72*x7^2*x8 + 72*x7^2*x9 + 72*x7^2*x10 + 72*x7*x8^2 + 144*x7*x8*x9 +
   144*x7*x8*x10 + 72*x7*x9^2 + 144*x7*x9*x10 + 72*x7*x10^2 + 48*x8^3 +
   72*x8^2*x9 + 72*x8^2*x10 + 72*x8*x9^2 + 144*x8*x9*x10 + 72*x8*x10^2 +
   48*x9^3 + 72*x9^2*x10 + 72*x9*x10^2 + 48*x10^3 )

function f10(x)
   f  = 48*x[1]^3 + 72*x[1]^2*x[2] + 72*x[1]^2*x[3] + 72*x[1]^2*x[4] + 72*x[1]^2*x[5] + 72*x[1]^2*x[7]
   f += 72*x[1]^2*x[8] + 72*x[1]*x[2]^2 + 144*x[1]*x[2]*x[4] + 144*x[1]*x[2]*x[7] + 72*x[1]*x[3]^2
   f += 144*x[1]*x[3]*x[5] + 144*x[1]*x[3]*x[8] + 72*x[1]*x[4]^2 + 144*x[1]*x[4]*x[7] + 72*x[1]*x[5]^2
   f += 144*x[1]*x[5]*x[8] + 72*x[1]*x[7]^2 + 72*x[1]*x[8]^2 + 48*x[2]^3 + 72*x[2]^2*x[3]
   f += 72*x[2]^2*x[4] + 72*x[2]^2*x[6] + 72*x[2]^2*x[7] + 72*x[2]^2*x[9] + 72*x[2]*x[3]^2
   f += 144*x[2]*x[3]*x[6] + 144*x[2]*x[3]*x[9] + 72*x[2]*x[4]^2 + 144*x[2]*x[4]*x[7] + 72*x[2]*x[6]^2
   f += 144*x[2]*x[6]*x[9] + 72*x[2]*x[7]^2 + 72*x[2]*x[9]^2 + 48*x[3]^3 + 72*x[3]^2*x[5]
   f += 72*x[3]^2*x[6] + 72*x[3]^2*x[8] + 72*x[3]^2*x[9] + 72*x[3]*x[5]^2 + 144*x[3]*x[5]*x[8]
   f += 72*x[3]*x[6]^2 + 144*x[3]*x[6]*x[9] + 72*x[3]*x[8]^2 + 72*x[3]*x[9]^2 + 48*x[4]^3
   f += 72*x[4]^2*x[5] + 72*x[4]^2*x[6] + 72*x[4]^2*x[7] + 72*x[4]^2*x[10] + 72*x[4]*x[5]^2
   f += 144*x[4]*x[5]*x[6] + 144*x[4]*x[5]*x[10] + 72*x[4]*x[6]^2 + 144*x[4]*x[6]*x[10] + 72*x[4]*x[7]^2
   f += 72*x[4]*x[10]^2 + 48*x[5]^3 + 72*x[5]^2*x[6] + 72*x[5]^2*x[8] + 72*x[5]^2*x[10]
   f += 72*x[5]*x[6]^2 + 144*x[5]*x[6]*x[10] + 72*x[5]*x[8]^2 + 72*x[5]*x[10]^2 + 48*x[6]^3
   f += 72*x[6]^2*x[9] + 72*x[6]^2*x[10] + 72*x[6]*x[9]^2 + 72*x[6]*x[10]^2 + 48*x[7]^3
   f += 72*x[7]^2*x[8] + 72*x[7]^2*x[9] + 72*x[7]^2*x[10] + 72*x[7]*x[8]^2 + 144*x[7]*x[8]*x[9]
   f += 144*x[7]*x[8]*x[10] + 72*x[7]*x[9]^2 + 144*x[7]*x[9]*x[10] + 72*x[7]*x[10]^2 + 48*x[8]^3
   f += 72*x[8]^2*x[9] + 72*x[8]^2*x[10] + 72*x[8]*x[9]^2 + 144*x[8]*x[9]*x[10] + 72*x[8]*x[10]^2
   f += 48*x[9]^3 + 72*x[9]^2*x[10] + 72*x[9]*x[10]^2 + 48*x[10]^3
   return f
end

function g10(x)
   g = sum(x)
   x2 = x .* x
   x3 = x2 .* x
   f = -24 * sum(x3) + 72 * g * sum(x2)
   f -= 72*(
       x2[1]*(x[6]+x[9]+x[10]) +
       x2[2]*(x[5]+x[8]+x[10]) +
       x2[3]*(x[4]+x[7]+x[10]) +
       x2[4]*(x[3]+x[8]+x[9]) +
       x2[5]*(x[2]+x[7]+x[9]))
   f += 72 *(
       x2[6]* (x[1]+x[7]+x[8]) +
       x2[7]* (x[3]+x[5]+x[6]) +
       x2[8]* (x[1]+x[4]+x[6]) +
       x2[9]* (x[2]+x[4]+x[5]) +
       x2[10]*(x[1]+x[2]+x[3]) )
   f += 144*(
         x[1]*(x[2]*(x[4] + x[7]) + x[3]*(x[5] + x[8])) +
         x[2]*x[3]*(x[6]+x[9]) +
         x[4]*(x[10]*(x[5]+x[6]) + (x[1]+x[2])*x[7]) +
         x[5]*((x[1]+x[3])*x[8] + x[6]*(x[4]+x[10])) +
         x[9]*((x[2]+x[3])*x[6] + x[7]*(x[8]+x[10])) +
         (x[7]+x[9])*x[8]*x[10] )
   return f
end


x = @SVector rand(10)

@btime f10($x)
@btime g10($x)
@btime StaticPolynomials.evaluate($p10, $x)
@btime StaticPolynomials.evaluate($q10, $x)


@btime f10($x)   # 34.261 ns (0 allocations: 0 bytes)
@btime StaticPolynomials.evaluate($p10, $x)   # 47.131 ns (0 allocations: 0 bytes)
@btime ForwardDiff.gradient($f10, $x)   # 2.017 μs (73 allocations: 6.84 KiB)
@btime StaticPolynomials.gradient($p10, $x)   # 102.675 ns (0 allocations: 0 bytes)



using StaticArrays, BenchmarkTools
import StaticPolynomials
using DynamicPolynomials: @polyvar

@polyvar Q1 Q2 Q3 Q4 Q5 Q6

const INV6Q = StaticPolynomials.system(
   [  Q1,
      Q2^2 + Q3^2 + Q4^2,
      Q2 * Q3 * Q4,
      Q3^2 * Q4^2 + Q2^2 * Q4^2 + Q2^2 * Q3^2,
      Q5^2 + Q6^2,
      Q6^3 - 3*Q5^2 * Q6,
      1.0,
      2*Q6*Q2^2 - 2*Q6*Q3^2 - 2*Q6*Q4^2 + √3*Q5*Q3^2 - √3*Q5*Q4^2,
      (Q6^2 - Q5^2) * (2*Q2^2 - Q3^2 - Q4^2) - (2*√3)*Q5*Q6*Q3^2 - (2*√3)*Q5*Q6*Q4^2,
      2*Q6*Q3^2*Q4^2 - Q6*Q2^2*Q4^2 - Q6*Q2^2*Q3^2 + √3*Q2*(Q2^2 * Q4^2 - Q2^2 * Q3^2),
      (Q6^2 - Q5^2)*(2*Q3^2*Q4^2 - Q2^2*Q4^2 -Q2^2*Q3^2) - 2*√3 * Q5 * Q6 * (Q2^2*Q4^2 - Q2^2*Q3^2),
      (Q3^2 - Q4^2) * (Q4^2 - Q2^2) * (Q2^2 - Q3^2) * Q5 * (3*Q6^2 - Q5^2)
   ])

Q = @SVector rand(6)
@btime StaticPolynomials.evaluate($INV6Q, $Q)
@btime StaticPolynomials.jacobian($INV6Q, $Q)


const _2 = 2.0^(-0.5)
const _3 = 3.0^(-0.5)
const _6 = 6.0^(-0.5)
const _12 = 12.0^(-0.5)

const R2Q = @SMatrix [ _6     _6     _6    _6     _6     _6
                        _2      0      0   -_2      0      0
                         0     _2      0     0    -_2      0
                         0      0     _2     0      0    -_2
                         0    0.5   -0.5     0    0.5   -0.5
                        _3   -_12   -_12    _3   -_12   -_12 ]

r2Q(r) =
   SVector( 6.0^(-0.5) * sum(r),
            2.0^(-0.5) * (r[1]-r[4]),
            2.0^(-0.5) * (r[2]-r[5]),
            2.0^(-0.5) * (r[3]-r[6]),
            0.5 * (r[2]-r[3]+r[5]-r[6]),
            3.0^(-0.5) * (r[1]+r[4]) - 12.0^(-0.5)*(r[2]+r[3]+r[5]+r[6]))

f(r) = StaticPolynomials.evaluate(INV6Q, r2Q(r))
df(r) = StaticPolynomials.jacobian(INV6Q, r2Q(r)) * R2Q
r = @SVector rand(6)
@btime f($r)
@btime df($r)

using StaticArrays, BenchmarkTools
using DynamicPolynomials: @polyvar
import XGrad, ForwardDiff, StaticPolynomials
r = @SVector rand(6)
ra = rand(6)
@polyvar x1 x2 x3 x4 x5 x6
P = StaticPolynomials.Polynomial(x1^3*(x2+x3+x4+x5) + x1*(x2^3+x3^3+x4^3+x5^3) + x2^3*(x3+x4+x6) + x2*(x3^3+x4^3+x6^3) + (x3^3+x4^3)*(x5+x6) + (x3+x4)*(x5^3+x6^3) +  x5^3*x6 + x5*x6^3)
p(x) = (x[1]^3*(x[2]+x[3]+x[4]+x[5]) + x[1]*(x[2]^3+x[3]^3+x[4]^3+x[5]^3) + x[2]^3*(x[3]+x[4]+x[6])) + (x[2]*(x[3]^3+x[4]^3+x[6]^3) + (x[3]^3+x[4]^3)*(x[5]+x[6]) + (x[3]+x[4])*(x[5]^3+x[6]^3) +  x[5]^3*x[6] + x[5]*x[6]^3)
∇p(x) = XGrad.xdiff(p; x=ra)
@btime StaticPolynomials.evaluate($P, $r)
@btime StaticPolynomials.gradient($P, $r)
@btime p($r)
@btime ForwardDiff.gradient($p, $r)
∇p(r)



f(x) = ((x[1]*x[2])*x[3] + (x[1]*x[4])*x[5]) + ((x[2]*x[4])*x[6] + (x[3]*x[5])*x[6])
r = @SVector rand(6)
ra = rand(6)
XGrad.xdiff(f, x=ra)
