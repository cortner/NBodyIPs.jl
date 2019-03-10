#
# To plot all cutoff functions, uncomment the code at the end of th e
# file and then execute the entire file.
#

using Test, NBodyIPs, BenchmarkTools
using NBodyIPs.Cutoffs
using NBodyIPs: fcut, fcut_d
import JuLIP

r0 = 0.712
rn = 1.234
rc = 2.890

cutoffs = [ CosCut(rc-1.0, rc),
            CosCut2s(r0, r0+0.4, rc-0.9, rc),
            PolyCut(2, rc),
            PolyCutSym(2, rc),
            PolyCut2sA(r0, rn, rc) ]
cnames = ["CosCut", "CosCut2s", "PolyCut", "PolyCutSym", "PolyCut2sA"]

##
@info("Testing Correctness of Derivative Implementations")
xt = [range(r0-0.5, rc+0.5, length=50); r0; rc]
for (C, name) in zip(cutoffs, cnames)
   println("Testing Derivatives of $name")
   f = r -> fcut(C, r)
   df = r -> fcut_d(C, r)
   println(@test JuLIP.Testing.fdtest_R2R(f, df, xt))
   Cd = Dict(C)
   C1 = JuLIP.decode_dict(Cd)
   Cd1 = Dict(Cd)
   println((@test C1 == C))
   println((@test Cd == Cd1))
end



##
@info("Testing Performance of fcut and fcut_d")
@info("One test includes 1_000 fcut/fcut_d calls")
function runN(f, x, dx, N)
   s = 0.0
   for n = 1:N
      x += dx
      s += f(x)
   end
   return s
end

N = 1_000
x0 = r0
dx = (rc+0.3-r0)/N

for (C, name) in zip(cutoffs, cnames)
   println(name)
   print("   fcut  : ")
   @btime runN($(r->fcut(C,r)), $x0, $dx, $N)
   print("   fcut_d: ")
   @btime runN($(r->fcut_d(C,r)), $x0, $dx, $N)
end


##
# using Plots
# xx = range(r0-0.5, rc+0.5, length=200)
# P = plot(; ylims = [-0.2, 1.2])
# for (C, name) in zip(cutoffs, cnames)
#    plot!(xx, fcut.(Ref(C), xx), label=name)
# end
# P = scatter!([r0, rc], [0, 0], label = "r0, rcut")
# P = vline!([rn], c=:black, ls=:dash, label = "rn")
# display(P)
