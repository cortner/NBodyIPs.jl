

using JSON, NBodyIPs, BenchmarkTools, StaticArrays

function eval_mpoly(P::Vector{SVector{6, Int}}, x::SVector{6, T}) where {T}
   out = 0.0
   for p in P
      out += x[1]^p[1] * x[2]^p[2] * x[3]^p[3] * x[4]^p[4] * x[5]^p[5] * x[6]^p[6]
   end
   return out
end

function eval_mpoly(P::Vector{SVector{10, Int}}, x::SVector{10, T}) where {T}
   out = 0.0
   for p in P
      out += (x[1]^p[1] * x[2]^p[2] * x[3]^p[3] * x[4]^p[4] * x[5]^p[5] *
              x[6]^p[6] * x[7]^p[7] * x[8]^p[8] * x[9]^p[9] * x[10]^p[10] )
   end
   return out
end

p = JSON.parsefile(@__DIR__() * "/poly6.json")
p = Vector{Vector{Int}}(p)
ps = Vector{SVector{6, Int}}(p)
D = dict(:poly, 8)
ex, f, _ = NBodyIPs.gen_fun(p, D...)

x = @SVector rand(6)

@btime f($x)
@btime eval_mpoly($ps, $x)
