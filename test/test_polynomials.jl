using ManyBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

# ndict = 5
# for n = 2:5 # dimension / body-order
#    for _ = 1:2 # number of random tests
#       alpha = rand(0:4, n)
#       d, sym = dict(:poly, ndict)
#       ex, f, df = ManyBodyIPs.psym_monomial(alpha, d, sym)
#       R = @SVector rand(n)
#       @show f(R)
#       @show df(R)
#       # df(R)
#    end
# end

ex, f, df = ManyBodyIPs.psym_monomial([1,3,2], dict(:inv1, 6, 3.0);
                                       simplify=false);
R = 1.0 + @SVector rand(3)
@btime f($R)
@btime df($R)



df1 = ManyBodyIPs.gen_psym_grad(3, ex, :x)
@btime df1($R)
