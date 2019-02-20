
import JuLIP
using FunctionWrappers,  BenchmarkTools, StaticPolynomials, StaticArrays

import NBodyIPs
using DynamicPolynomials: @polyvar
using JuLIP.Potentials: coscut, F64fun

##
@info "[1] Profile FunctionWrappers"
println("    1000 evaluations of coscut")
f(r) = coscut(r, 1.123, 3.456)
fwrap = F64fun(f)
r = 1.0 .+ 3*rand(1000)
test1(f, r) = (s=0.0; for n = 1:length(r); s += f(r[n]); end; s)
print(" Plain Function: "); @btime test1($f, $r)
print("FunctionWrapper: "); @btime test1($fwrap, $r)


@info "[2] Profile StaticPolynomials"

##
@info "[2.1] From the StaticPolynomials README"

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

@polyvar x[1:10]
p10 = Polynomial(f10(x))
y = @SVector rand(10)
test2(f, x) = (s=0.0; for n=1:1000; s+= f(x); end; s)
print("     Manual Implementation: "); @btime test2($f10, $y)
print("         StaticPolynomials: "); @btime test2($p10, $y)
print("StaticPolynomials.gradient: "); @btime StaticPolynomials.gradient($p10, $y)


##
@info "[2.2] A BL5 Invariant"
f10I(x) = (x[1]^2*x[2]*x[5]*x[3] + x[1]*x[2]^2*x[5]*x[3] + x[1]^2*x[2]*x[5]*x[6] + x[1]*x[2]*x[5]^2*x[6] + x[1]^2*x[2]*x[3]*x[6] + x[1]^2*x[5]*x[3]*x[6] + x[1]*x[2]*x[3]^2*x[6] + x[1]*x[5]*x[3]*x[6]^2 + x[1]*x[2]^2*x[5]*x[8] + x[1]*x[2]*x[5]^2*x[8] + x[1]*x[2]^2*x[3]*x[8] + x[2]^2*x[5]*x[3]*x[8] + x[1]*x[2]*x[3]^2*x[8] + x[1]*x[5]^2*x[6]*x[8] + x[2]*x[5]^2*x[6]*x[8] + x[1]*x[3]^2*x[6]*x[8] + x[2]*x[3]^2*x[6]*x[8] + x[1]*x[5]*x[6]^2*x[8] + x[1]*x[3]*x[6]^2*x[8] + x[5]*x[3]*x[6]^2*x[8] + x[2]*x[5]*x[3]*x[8]^2 + x[2]*x[5]*x[6]*x[8]^2 + x[2]*x[3]*x[6]*x[8]^2 + x[5]*x[3]*x[6]*x[8]^2 + x[1]^2*x[2]*x[5]*x[4] + x[1]*x[2]^2*x[5]*x[4] + x[1]^2*x[3]*x[6]*x[4] + x[1]*x[3]^2*x[6]*x[4] + x[2]^2*x[3]*x[8]*x[4] + x[2]*x[3]^2*x[8]*x[4] + x[1]^2*x[2]*x[5]*x[7] + x[1]*x[2]*x[5]^2*x[7] + x[1]^2*x[3]*x[6]*x[7] + x[1]*x[3]*x[6]^2*x[7] + x[5]^2*x[6]*x[8]*x[7] + x[5]*x[6]^2*x[8]*x[7] + x[1]^2*x[2]*x[4]*x[7] + x[1]^2*x[5]*x[4]*x[7] + x[1]^2*x[3]*x[4]*x[7] + x[1]^2*x[6]*x[4]*x[7] + x[1]*x[2]*x[4]^2*x[7] + x[1]*x[3]*x[4]^2*x[7] + x[1]*x[5]*x[4]*x[7]^2 + x[1]*x[6]*x[4]*x[7]^2 + x[1]*x[2]^2*x[5]*x[9] + x[1]*x[2]*x[5]^2*x[9] + x[2]^2*x[3]*x[8]*x[9] + x[5]^2*x[6]*x[8]*x[9] + x[2]*x[3]*x[8]^2*x[9] + x[5]*x[6]*x[8]^2*x[9] + x[1]*x[2]^2*x[4]*x[9] + x[2]^2*x[5]*x[4]*x[9] + x[2]^2*x[3]*x[4]*x[9] + x[2]^2*x[8]*x[4]*x[9] + x[1]*x[2]*x[4]^2*x[9] + x[2]*x[3]*x[4]^2*x[9] + x[1]*x[5]^2*x[7]*x[9] + x[2]*x[5]^2*x[7]*x[9] + x[5]^2*x[6]*x[7]*x[9] + x[5]^2*x[8]*x[7]*x[9] + x[1]*x[4]^2*x[7]*x[9] + x[2]*x[4]^2*x[7]*x[9] + x[1]*x[5]*x[7]^2*x[9] + x[5]*x[6]*x[7]^2*x[9] + x[1]*x[4]*x[7]^2*x[9] + x[5]*x[4]*x[7]^2*x[9] + x[2]*x[5]*x[4]*x[9]^2 + x[2]*x[8]*x[4]*x[9]^2 + x[2]*x[5]*x[7]*x[9]^2 + x[5]*x[8]*x[7]*x[9]^2 + x[2]*x[4]*x[7]*x[9]^2 + x[5]*x[4]*x[7]*x[9]^2 + x[1]*x[3]^2*x[6]*x[10] + x[1]*x[3]*x[6]^2*x[10] + x[2]*x[3]^2*x[8]*x[10] + x[5]*x[6]^2*x[8]*x[10] + x[2]*x[3]*x[8]^2*x[10] + x[5]*x[6]*x[8]^2*x[10] + x[1]*x[3]^2*x[4]*x[10] + x[2]*x[3]^2*x[4]*x[10] + x[3]^2*x[6]*x[4]*x[10] + x[3]^2*x[8]*x[4]*x[10] + x[1]*x[3]*x[4]^2*x[10] + x[2]*x[3]*x[4]^2*x[10] + x[1]*x[6]^2*x[7]*x[10] + x[5]*x[6]^2*x[7]*x[10] + x[3]*x[6]^2*x[7]*x[10] + x[6]^2*x[8]*x[7]*x[10] + x[1]*x[4]^2*x[7]*x[10] + x[3]*x[4]^2*x[7]*x[10] + x[1]*x[6]*x[7]^2*x[10] + x[5]*x[6]*x[7]^2*x[10] + x[1]*x[4]*x[7]^2*x[10] + x[6]*x[4]*x[7]^2*x[10] + x[2]*x[8]^2*x[9]*x[10] + x[5]*x[8]^2*x[9]*x[10] + x[3]*x[8]^2*x[9]*x[10] + x[6]*x[8]^2*x[9]*x[10] + x[2]*x[4]^2*x[9]*x[10] + x[3]*x[4]^2*x[9]*x[10] + x[5]*x[7]^2*x[9]*x[10] + x[6]*x[7]^2*x[9]*x[10] + x[2]*x[8]*x[9]^2*x[10] + x[5]*x[8]*x[9]^2*x[10] + x[2]*x[4]*x[9]^2*x[10] + x[8]*x[4]*x[9]^2*x[10] + x[5]*x[7]*x[9]^2*x[10] + x[8]*x[7]*x[9]^2*x[10] + x[3]*x[6]*x[4]*x[10]^2 + x[3]*x[8]*x[4]*x[10]^2 + x[3]*x[6]*x[7]*x[10]^2 + x[6]*x[8]*x[7]*x[10]^2 + x[3]*x[4]*x[7]*x[10]^2 + x[6]*x[4]*x[7]*x[10]^2 + x[3]*x[8]*x[9]*x[10]^2 + x[6]*x[8]*x[9]*x[10]^2 + x[3]*x[4]*x[9]*x[10]^2 + x[8]*x[4]*x[9]*x[10]^2 + x[6]*x[7]*x[9]*x[10]^2 + x[8]*x[7]*x[9]*x[10]^2)
p10I = Polynomial(f10I(x))
# print("     Manual Implementation: "); @btime test2($f10I, $y)
print("         StaticPolynomials: "); @btime test2($p10I, $y)
print("StaticPolynomials.gradient: "); @btime StaticPolynomials.gradient($p10I, $y)


##
@info "[2.3] A generic 4B total-degree 10 polynomial"
function sp_exps(A)
   exps = zeros(Int, (length(A), 12))
   for (i, a) in enumerate(A)
      exps[i,1:6] .= a[1:6]
      exps[i,7+a[end]] = 1
   end
   return exps
end
desc = NBodyIPs.BondLengthDesc("r->r", (:cos, 1.0, 2.0))
A = NBodyIPs.PolyBasis.gen_tuples(desc, 4, 10)
C = rand(Float64, length(A))
p4B = StaticPolynomials.Polynomial(C, collect(sp_exps(A)'))
x23 = @SVector rand(12)
print("         StaticPolynomials: "); @btime $p4B($x23)
print("StaticPolynomials.gradient: "); @btime StaticPolynomials.gradient($p4B, $x23)
