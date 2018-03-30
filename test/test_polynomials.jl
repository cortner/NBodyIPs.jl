using ManyBodyIPs
using JuLIP, Base.Test, StaticArrays, Combinatorics
using BenchmarkTools

println("testing symmetry for different body-orders and dictionaries")
ndict = 6
for n = 2:5 # dimension / body-order
   for _ = 1:2 # number of random tests
      alpha = rand(0:(ndict-1), n)
      ex, f, df = ManyBodyIPs.psym_monomial(alpha, dict(:inv2, ndict, 3.0))
      R = 1.0 + @SVector rand(n)
      f_ = f(R)
      issym = true
      for S in Combinatorics.permutations(R)
         if abs(f(S) - f(R)) > 1e-10
            issym = false
            break
         end
      end
      if issym
         print("✓")
      else
         print("x")
      end
      @test issym
   end
end

for dict_sym in [:poly, :poly1, :poly2, :inv1, :inv2]
   println("testing `psym_monomial` with $dict_sym dictionary")
   α = [2,1,3,1]
   dim = length(α)
   ex, f, df = ManyBodyIPs.psym_monomial(α, dict(:inv2, 8, 3.0));
   R = 1.0 + @SVector rand(4)
   # @btime f($R)
   # @btime df($R)
   println("checking symmetry")
   f_ = f(R)
   issym = true
   for S in Combinatorics.permutations(R)
      if abs(f(S) - f(R)) > 1e-10
         issym = false
         break
      end
   end
   @test issym

   err = []
   f_ = f(R)
   df_ = df(R)
   println("finite-difference test")
   for p = 2:12
      h = 0.1^p
      dfh = zeros(length(df_))
      Rh = Vector(R)
      for i = 1:length(R)
         Rh[i] += h
         dfh[i] = (f(Rh) - f_) / h
         Rh[i] -= h
      end
      push!(err, norm(df_ - dfh, Inf))
      @printf(" %.2e : %.3e\n", h, err[end])
   end
   @test minimum(err) < 1e-3 * maximum(err)
end
