
using FunctionWrappers, JuLIP, BenchmarkTools

using JuLIP.Potentials: coscut, F64fun

@info "[1] Profile FunctionWrappers"
println("    1000 evaluations of coscut")
f = let rc1 = 1.123, rc2 = 3.456
   r -> coscut(r, rc1, rc2)
end
fwrap = F64fun(f)
r = 1.0 .+ 3*rand(1000)
test1(f, r) = (s=0.0; for n = 1:length(r); s += f(r[n]); end; s)
print(" Plain Function: "); @btime test1($f, $r)
print("FunctionWrapper: "); @btime test1($fwrap, $r)
