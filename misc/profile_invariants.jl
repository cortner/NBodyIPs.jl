
using  BenchmarkTools, StaticArrays
import NBodyIPs

test(f, x) = (s = 0.0; for n = 1:10; s += f(x)[1][1]; end; s)

@info "Runtimes of single invariants computation"
@info "Devide runtimes by 10"
for N in 2:4
   x = @SVector rand((N*(N-1))÷2)
   print("$N-body BL-invariants: ");
   @btime test(NBodyIPs.BLInvariants.invariants_ed, $x)
   print("$N-body BA-invariants: ");
   @btime test(NBodyIPs.BAInvariants.invariants_ed, $x)
end


x = @SVector rand(10)
@btime NBodyIPs.BLInvariants.invariants($x)
@btime NBodyIPs.BLInvariants.invariants_d($x)
@btime NBodyIPs.BLInvariants.invariants_ed($x)
