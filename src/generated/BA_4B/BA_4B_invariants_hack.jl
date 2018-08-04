module BA_4B

using NBodyIPs.FastPolys
using StaticArrays
using BenchmarkTools: @btime

import NBodyIPs.tdegrees

const G_BA_4B = [
[ 1, 2, 3, 4, 5, 6 ]
,[ 2, 3, 1, 6, 4, 5 ]
,[ 3, 1, 2, 5, 6, 4 ]
,]
simplex_permutations(x::SVector{6}) = [x[G_BA_4B[i]] for i=1:3]
# Primary invariants for BA_4B
 # : definitions at the beginning of the file
const P1_1 = (1,3,2,)

const P2_1 = (4,5,6,)

const P3_1 = (1,3,2,)

const P4_1 = (1,3,2,)
const P4_2 = (4,5,6,)

const P5_1 = (4,5,6,)

const P6_1 = (1,3,2,4,6,5)
const P6_2 = (1,3,2,5,4,6)

# Irreducible secondaries for group BA_4B
 # : definitions at the beginning of the file
const IS1_1 = (1,3,2,)
const IS1_2 = (5,6,4,)

const IS2_1 = (4,5,6,)

const IS3_1 = (1,3,2,)
const IS3_2 = (3,2,1,)

const IS4_1 = (1,1,2,)
const IS4_2 = (2,3,3,)
const IS4_3 = (4,5,6,)

const IS5_1 = (4,5,6,)
const IS5_2 = (1,3,2,)

const IS6_1 = (1,3,2,)
const IS6_2 = (5,6,4,)

const IS7_1 = (1,3,2,)
const IS7_2 = (4,5,4,)
const IS7_3 = (5,6,6,)

const IS8_1 = (4,5,6,)
const IS8_2 = (6,4,5,)


# Primary invariants for BA_4B
 # : definitions of the types at the beginning of the file
const P1 = Val((P1_1,))
const P2 = Val((P2_1,))
const P3 = Val((P3_1,))
const P4 = Val((P4_1,P4_2,))
const P5 = Val((P5_1,))
const P6 = Val((P6_1,P6_2))
# Irreducible secondaries for group BA_4B
 # : definitions of the types at the beginning of the file
const IS1 = Val((IS1_1,IS1_2,))
const IS2 = Val((IS2_1,))
const IS3 = Val((IS3_1,IS3_2,))
const IS4 = Val((IS4_1,IS4_2,IS4_3,))
const IS5 = Val((IS5_1,IS5_2,))
const IS6 = Val((IS6_1,IS6_2,))
const IS7 = Val((IS7_1,IS7_2,IS7_3,))
const IS8 = Val((IS8_1,IS8_2,))


function invariants(x1::SVector{6, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_4B
 # : what goes in the function for the evaluation
P1 = fpoly((x1,) , BA_4B.P1)
P2 = fpoly((x1,) , BA_4B.P2)
P3 = fpoly((x2,) , BA_4B.P3)
P4 = fpoly((x1,x1,) , BA_4B.P4)
P5 = fpoly((x3,) , BA_4B.P5)
P6 = fpoly((x2,x1,) , BA_4B.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_4B
 # : what goes in the function for the evaluation
IS1 = fpoly((x1,x1,) , BA_4B.IS1)
IS2 = fpoly((x2,) , BA_4B.IS2)
IS3 = fpoly((x2,x1,) , BA_4B.IS3)
IS4 = fpoly((x1,x1,x1,) , BA_4B.IS4)
IS5 = fpoly((x2,x1,) , BA_4B.IS5)
IS6 = fpoly((x2,x1,) , BA_4B.IS6)
IS7 = fpoly((x1,x1,x1,) , BA_4B.IS7)
IS8 = fpoly((x2,x1,) , BA_4B.IS8)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1
SEC3  = IS2
SEC4  = IS3
SEC5  = IS4
SEC6  = IS5
SEC7  = IS6
SEC8  = IS7
SEC9  = IS8
SEC10  = IS1*IS2
SEC11  = IS2^2
SEC12  = IS4^2


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,SEC7,SEC8,SEC9,SEC10,SEC11,SEC12,])
 end



function invariants_d(x1::SVector{6, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1

   dx1 = @SVector ones(6)
   dx2 = 2 * x1
   dx3 = 3 * x2
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_4B
 # : what goes in the function for the derivatives
dP1 = fpoly_d((x1,),(dx1,) , BA_4B.P1)
dP2 = fpoly_d((x1,),(dx1,) , BA_4B.P2)
dP3 = fpoly_d((x2,),(dx2,) , BA_4B.P3)
dP4 = fpoly_d((x1,x1,),(dx1,dx1,) , BA_4B.P4)
dP5 = fpoly_d((x3,),(dx3,) , BA_4B.P5)
dP6 = fpoly_d((x2,x1,),(dx2,dx1,) , BA_4B.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_4B
 # : what goes in the function for the evaluation
IS1 = fpoly((x1,x1,) , BA_4B.IS1)
IS2 = fpoly((x2,) , BA_4B.IS2)
IS3 = fpoly((x2,x1,) , BA_4B.IS3)
IS4 = fpoly((x1,x1,x1,) , BA_4B.IS4)
IS5 = fpoly((x2,x1,) , BA_4B.IS5)
IS6 = fpoly((x2,x1,) , BA_4B.IS6)
IS7 = fpoly((x1,x1,x1,) , BA_4B.IS7)
IS8 = fpoly((x2,x1,) , BA_4B.IS8)


# Irreducible secondaries for group BA_4B
 # : what goes in the function for the derivatives
dIS1 = fpoly_d((x1,x1,),(dx1,dx1,) , BA_4B.IS1)
dIS2 = fpoly_d((x2,),(dx2,) , BA_4B.IS2)
dIS3 = fpoly_d((x2,x1,),(dx2,dx1,) , BA_4B.IS3)
dIS4 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , BA_4B.IS4)
dIS5 = fpoly_d((x2,x1,),(dx2,dx1,) , BA_4B.IS5)
dIS6 = fpoly_d((x2,x1,),(dx2,dx1,) , BA_4B.IS6)
dIS7 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , BA_4B.IS7)
dIS8 = fpoly_d((x2,x1,),(dx2,dx1,) , BA_4B.IS8)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


dSEC1   = @SVector zeros(6)
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  = dIS4
dSEC6  = dIS5
dSEC7  = dIS6
dSEC8  = dIS7
dSEC9  = dIS8
dSEC10  =  + dIS1*IS2 + dIS2*IS1
dSEC11  =  + dIS2*2IS2
dSEC12  =  + dIS4*2IS4


return (dP1,dP2,dP3,dP4,dP5,dP6,), (dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,dSEC7,dSEC8,dSEC9,dSEC10,dSEC11,dSEC12,)
 end



function invariants_ed(x1::SVector{6, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1

   dx1 = @SVector ones(6)
   dx2 = 2 * x1
   dx3 = 3 * x2
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_4B
 # : what goes in the function for the evaluation and derivatives
P1, dP1 = fpoly_ed((x1,),(dx1,) , BA_4B.P1)
P2, dP2 = fpoly_ed((x1,),(dx1,) , BA_4B.P2)
P3, dP3 = fpoly_ed((x2,),(dx2,) , BA_4B.P3)
P4, dP4 = fpoly_ed((x1,x1,),(dx1,dx1,) , BA_4B.P4)
P5, dP5 = fpoly_ed((x3,),(dx3,) , BA_4B.P5)
P6, dP6 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_4B.P6)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_4B
 # : what goes in the function for the evaluation and derivatives
IS1, dIS1 = fpoly_ed((x1,x1,),(dx1,dx1,) , BA_4B.IS1)
IS2, dIS2 = fpoly_ed((x2,),(dx2,) , BA_4B.IS2)
IS3, dIS3 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_4B.IS3)
IS4, dIS4 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , BA_4B.IS4)
IS5, dIS5 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_4B.IS5)
IS6, dIS6 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_4B.IS6)
IS7, dIS7 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , BA_4B.IS7)
IS8, dIS8 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_4B.IS8)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


SEC1  = 1
SEC2  = IS1
SEC3  = IS2
SEC4  = IS3
SEC5  = IS4
SEC6  = IS5
SEC7  = IS6
SEC8  = IS7
SEC9  = IS8
SEC10  = IS1*IS2
SEC11  = IS2^2
SEC12  = IS4^2


dSEC1   = @SVector zeros(6)
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  = dIS4
dSEC6  = dIS5
dSEC7  = dIS6
dSEC8  = dIS7
dSEC9  = dIS8
dSEC10  =  + dIS1*IS2 + dIS2*IS1
dSEC11  =  + dIS2*2IS2
dSEC12  =  + dIS4*2IS4


return (@SVector [P1,P2,P3,P4,P5,P6,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,SEC7,SEC8,SEC9,SEC10,SEC11,SEC12,]), (@SVector [dP1,dP2,dP3,dP4,dP5,dP6,]), (@SVector [dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,dSEC7,dSEC8,dSEC9,dSEC10,dSEC11,dSEC12,])
 end

end
