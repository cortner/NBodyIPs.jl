
using StaticArrays

simplex(r::SVector{1,T}) = @SVector T[r[1], 0, 0]

function simplex(r::SVector{3,T})
   R1 = @SVector T[r[1], 0, 0]
   R2 = 
