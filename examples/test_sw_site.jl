
using JuLIP, NBodyIPs, NeighbourLists, ForwardDiff, StaticArrays

Base.isvalid(r) = ( (r[1] + r[2] > r[3]) &&
                    (r[2] + r[3] > r[1]) &&
                    (r[1] + r[3] > r[2]) )

bondangle(S1, S2) = (dot(S1, S2) + 1.0/3.0)^2

function eval3body(calc::StillingerWeber, r)
   @assert isvalid(r)
   R1 = SVector(0.0, 0.0)
   R2 = SVector(r[1], 0.0)
   a = (r[1]^2 + r[2]^2 - r[3]^2) / (2*r[1])
   @assert a < r[2]
   R3 = SVector(a, sqrt(r[2]^2 - a^2))
   @assert (norm(R2) ≈ r[1]) && (norm(R3) ≈ r[2]) && (norm(R3-R2) ≈ r[3])
   S1 = R2/norm(R2)
   S2 = R3/norm(R3)
   S3 = (R3-R2)/norm(R3-R2)
   V3 = [ calc.V3(r1) for r1 in r ]
   return (   V3[1] * V3[2] * bondangle(S1, S2)
            + V3[1] * V3[3] * bondangle(-S1, S3)
            + V3[2] * V3[3] * bondangle(-S2, -S3) )
end

eval3body_d(calc::StillingerWeber, r) =
      ForwardDiff.gradient(r->eval3body(calc, r), r)

function gen_data(N, r0, r1)
   sw = StillingerWeber()
   data = Tuple{SVector{3, Float64}, Float64}[]
   for n = 1:N
      r = r0 + (r1-r0) * (@SVector rand(3))
      while !isvalid(r)
         r = r0 + (r1-r0) * (@SVector rand(3))
      end
      push!(data, (r, eval3body(sw, r)))
   end
   return data
end

r0 = rnn(:Si)
sw = StillingerWeber()
rcut = cutoff(sw)
rcutN = 2 * rcut

train_data = gen_data(2_000, 0.8*r0, 1.1*rcutN)
test_data =  gen_data(100, 0.9*r0, rcutN)

Es = Float64[d[2] for d in train_data]
Eavg = mean([extrema(Es)...])

NDICT = [3, 5, 7, 9, 11]
DSYM = [:poly1, :poly2, :inv1, :inv2]

err_rms = zeros(length(DSYM), length(NDICT))
nbasis = zeros(Int, length(NDICT))
for (is, dsym) in enumerate(DSYM), (in, ndict) in enumerate(NDICT)
   B = get_basis(3, dict(dsym, ndict, rcutN)..., rcutN)
   nbasis[in] = length(B)
   @show (dsym, ndict, length(B))
   BB = [b.f for b in B]
   c = NBodyIPs.regression(BB, train_data; verbose=false)
   err_rms[is, in] = NBodyIPs.rms(c, BB, test_data) * sqrt(3) / Eavg
   println(" => rms(testset) = ", err_rms[is, in])
end


using DataFrames
df = DataFrame()
df[:nbasis] = nbasis
for n = 1:length(DSYM)
   df[DSYM[n]] = err_rms[n, :]
end
print(df)

# at = data[1][1]::Atoms
# for b in basis(12)
#    @time b(at)
# end
