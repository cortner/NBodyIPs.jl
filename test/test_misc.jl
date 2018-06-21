using NBodyIPs, JuLIP, StaticArrays

R = rand(JVecF, 10)
J = @SVector [3,5,7]
s = Array(norm.(R[J]))
for i = 1:2, j = i+1:3
   push!(s, norm(R[J[i]]-R[J[j]]))
end
s2, _ = NBodyIPs.simplex_edges(R, J)
s2 == s
