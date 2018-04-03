using ManyBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools
using FunctionWrappers: FunctionWrapper
const FWrap{N, T} = FunctionWrapper{Float64, Tuple{SVector{N,T}}}
const F64Wrap = FunctionWrapper{Float64, Tuple{Float64}}

ex, f, df = ManyBodyIPs.psym_monomial([2,1,2], dict(:inv2, 5, 3.0))
b = NBody(3, f, df, 3.0)
R = 1.0 + @SVector rand(3)
@btime f($R)
@btime df($R)
@btime b.f($R)
@btime b.f_d($R)

@code_warntype f(R)

using Calculus
simplify(ex)

ex, f, df = ManyBodyIPs.psym_monomial([2,3,1,2,1], dict(:inv1, 8, 3.0))
simplify(ex)
R = 1.0 + @SVector rand(5)
@btime f(R)
fwrap = FWrap{5, Float64}(f)
@btime fwrap($R)

g = r -> 3.1*r^5 - 1.2 * r^4 + r^3 - 3.2 * r^1 + 1.123 * r + 1.234
gwrap = F64Wrap(g)
r = rand()
@btime g($r)
@btime gwrap($r)

function Ncalls(f::F, N, x) where F
   s = 0.0
   for n = 1:N
      s += f(x)
   end
   return s
end

Ncalls(g, 10, r);
Ncalls(gwrap, 10, r);

@time Ncalls(g, 1_000_000, r)
@time Ncalls(gwrap, 1_000_000, r)

@time Ncalls(f, 100_000, R)
@time Ncalls(fwrap, 100_000, R)
