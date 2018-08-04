module BA_5B

using NBodyIPs.FastPolys
using StaticArrays
using BenchmarkTools: @btime

import NBodyIPs.tdegrees

const G_BA_5B = [
[ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
,[ 3, 1, 4, 2, 6, 10, 8, 7, 5, 9 ]
,[ 2, 4, 1, 3, 9, 5, 8, 7, 10, 6 ]
,[ 4, 1, 2, 3, 7, 9, 10, 5, 6, 8 ]
,[ 1, 3, 4, 2, 6, 7, 5, 10, 8, 9 ]
,[ 3, 4, 2, 1, 10, 8, 6, 9, 7, 5 ]
,[ 2, 1, 3, 4, 5, 8, 9, 6, 7, 10 ]
,[ 4, 2, 3, 1, 9, 10, 7, 8, 5, 6 ]
,[ 1, 4, 2, 3, 7, 5, 6, 9, 10, 8 ]
,[ 3, 2, 1, 4, 8, 6, 10, 5, 9, 7 ]
,[ 2, 3, 4, 1, 8, 9, 5, 10, 6, 7 ]
,[ 4, 3, 1, 2, 10, 7, 9, 6, 8, 5 ]
,[ 1, 3, 2, 4, 6, 5, 7, 8, 10, 9 ]
,[ 3, 4, 1, 2, 10, 6, 8, 7, 9, 5 ]
,[ 2, 1, 4, 3, 5, 9, 8, 7, 6, 10 ]
,[ 4, 2, 1, 3, 9, 7, 10, 5, 8, 6 ]
,[ 1, 4, 3, 2, 7, 6, 5, 10, 9, 8 ]
,[ 3, 2, 4, 1, 8, 10, 6, 9, 5, 7 ]
,[ 2, 3, 1, 4, 8, 5, 9, 6, 10, 7 ]
,[ 4, 3, 2, 1, 10, 9, 7, 8, 6, 5 ]
,[ 1, 2, 4, 3, 5, 7, 6, 9, 8, 10 ]
,[ 3, 1, 2, 4, 6, 8, 10, 5, 7, 9 ]
,[ 2, 4, 3, 1, 9, 8, 5, 10, 7, 6 ]
,[ 4, 1, 3, 2, 7, 10, 9, 6, 5, 8 ]
,]
simplex_permutations(x::SVector{10}) = [x[G_BA_5B[i]] for i=1:24]
# Primary invariants for BA_5B
 # : definitions at the beginning of the file
const P1_1 = (1,2,3,4,)

const P2_1 = (5,9,6,8,7,10,)

const P3_1 = (1,2,3,4,)

const P4_1 = (5,9,6,8,7,10,)

const P5_1 = (1,2,3,2,1,4,2,4,1,3,4,3,)
const P5_2 = (5,9,6,8,7,10,5,9,6,8,7,10,)

const P6_1 = (1,2,3,4,)

const P7_1 = (5,5,6,7,)
const P7_2 = (6,8,8,9,)
const P7_3 = (7,9,10,10,)

const P8_1 = (5,9,6,8,7,10,)

const P9_1 = (5,9,6,8,7,10,)

const P10_1 = (1,2,3,4,5,5,5,5,6,7,8,9,6,6,6 ,7,8,10,7,7 ,9,10,8,8 ,9,10,9 ,10)
const P10_2 = (1,2,3,4,6,7,8,9,5,5,5,5,7,8,10,6,6,6 ,9,10,7,7 ,9,10,8,8 ,10,9)

# Irreducible secondaries for group BA_5B
 # : definitions at the beginning of the file
const IS1_1 = (5,5,6,8,5,7,5,9,6,6,7,8,)
const IS1_2 = (6,9,10,9,7,10,8,10,7,8,9,10,)

const IS2_1 = (1,2,1,2,1,3,)
const IS2_2 = (2,4,3,3,4,4,)
const IS2_3 = (10,6,9,7,8,5,)

const IS3_1 = (1,2,1,2,1,3,)
const IS3_2 = (2,4,3,3,4,4,)
const IS3_3 = (5,9,6,8,7,10,)

const IS4_1 = (8,7,7,10,9,6,6,8,10,5,5,9,)
const IS4_2 = (1,2,3,2,1,4,2,4,1,3,4,3,)

const IS5_1 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,)
const IS5_2 = (5,7,6,8,7,6,5,8,6,5,5,9,6,7,5,8,7,5,5,8,5,6,6,9,)
const IS5_3 = (8,9,7,10,9,10,6,9,10,8,7,10,8,10,7,9,10,9,6,10,9,8,7,10,)

const IS6_1 = (5,9,6,8,7,10,5,9,6,8,7,10,6,10,5,8,7,9,6,10,5,8,7,9,)
const IS6_2 = (6,5,10,9,5,7,8,10,7,6,9,8,5,6,9,10,6,7,8,9,7,5,10,8,)

const IS7_1 = (1,2,3,2,1,4,2,4,1,3,4,3,)
const IS7_2 = (2,1,1,3,2,1,1,2,3,1,1,2,)
const IS7_3 = (3,4,4,4,4,3,3,3,4,2,2,4,)

const IS8_1 = (1,2,3,2,1,4,2,4,1,3,4,3,)
const IS8_2 = (2,4,1,3,4,3,1,2,3,2,1,4,)
const IS8_3 = (10,6,9,7,8,5,10,6,9,7,8,5,)

const IS9_1 = (1,2,1,2,1,3,)
const IS9_2 = (2,4,3,3,4,4,)
const IS9_3 = (5,6,6,7,7,5,)
const IS9_4 = (10,9,9,8,8,10,)

const IS10_1 = (1,2,1,2,1,3,1,2,1,2,1,3,1,3,1,2,1,2,1,3,1,2,1,2,)
const IS10_2 = (2,4,3,3,4,4,2,4,3,3,4,4,3,4,2,3,4,4,3,4,2,3,4,4,)
const IS10_3 = (6,5,9,7,5,5,8,6,7,6,8,5,5,5,9,7,6,6,8,5,7,5,8,6,)
const IS10_4 = (10,6,10,9,8,7,10,10,9,7,9,8,9,6,10,10,8,7,9,9,10,7,10,8,)

const IS11_1 = (5,9,6,8,7,10,)
const IS11_2 = (1,2,1,2,1,3,)
const IS11_3 = (2,4,3,3,4,4,)

const IS12_1 = (10,6,9,7,8,5,10,6,9,7,8,5,)
const IS12_2 = (1,2,3,2,1,4,2,4,1,3,4,3,)
const IS12_3 = (5,9,6,8,7,10,5,9,6,8,7,10,)

const IS13_1 = (8,7,7,10,9,6,6,8,10,5,5,9,8,7,7,9,10,5,5,8,9,6,6,10,)
const IS13_2 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,)
const IS13_3 = (5,9,6,8,7,10,5,9,6,8,7,10,6,10,5,8,7,9,6,10,5,8,7,9,)

const IS14_1 = (8,7,7,10,9,6,6,8,10,5,5,9,8,7,7,9,10,5,5,8,9,6,6,10,)
const IS14_2 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,)
const IS14_3 = (9,10,5,6,10,8,7,5,8,9,6,7,10,9,6,5,9,8,7,6,8,10,5,7,)

const IS15_1 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,)
const IS15_2 = (2,1,1,3,2,1,1,2,3,1,1,2,2,1,1,2,3,1,1,2,2,1,1,3,)
const IS15_3 = (3,4,4,4,4,3,3,3,4,2,2,4,3,4,4,4,4,2,2,3,4,3,3,4,)
const IS15_4 = (5,9,6,8,7,10,5,9,6,8,7,10,6,10,5,8,7,9,6,10,5,8,7,9,)

const IS16_1 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,)
const IS16_2 = (2,4,1,3,4,3,1,2,3,2,1,4,3,4,1,2,4,2,1,3,2,3,1,4,)
const IS16_3 = (8,6,7,7,8,5,6,6,9,5,5,5,8,5,7,7,8,5,5,5,9,6,6,6,)
const IS16_4 = (10,7,9,10,9,6,10,8,10,7,8,9,9,7,10,9,10,6,9,8,10,7,8,10,)

const IS17_1 = (5,9,6,8,7,10,)
const IS17_2 = (1,2,1,2,1,3,)
const IS17_3 = (2,4,3,3,4,4,)
const IS17_4 = (10,6,9,7,8,5,)

const IS18_1 = (1,2,1,2,1,3,1,2,1,2,1,3,1,3,1,2,1,2,1,3,1,2,1,2,)
const IS18_2 = (2,4,3,3,4,4,2,4,3,3,4,4,3,4,2,3,4,4,3,4,2,3,4,4,)
const IS18_3 = (5,5,6,7,5,5,5,6,6,6,7,5,5,5,5,7,6,6,6,5,5,5,7,6,)
const IS18_4 = (6,6,9,8,7,7,8,9,7,7,8,8,6,6,9,8,7,7,8,9,7,7,8,8,)
const IS18_5 = (10,9,10,9,8,10,10,10,9,8,9,10,9,10,10,10,8,9,9,10,10,8,10,9,)

const IS19_1 = (9,10,5,6,10,8,7,5,8,9,6,7,10,9,6,5,9,8,7,6,8,10,5,7,)
const IS19_2 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,)
const IS19_3 = (5,5,6,8,5,7,5,9,6,6,7,8,5,6,5,8,6,7,6,9,5,5,7,8,)
const IS19_4 = (6,9,10,9,7,10,8,10,7,8,9,10,6,10,9,10,7,9,8,10,7,8,10,9,)

const IS20_1 = (10,6,9,7,8,5,10,6,9,7,8,5,9,5,10,7,8,6,9,5,10,7,8,6,)
const IS20_2 = (1,2,3,2,1,4,2,4,1,3,4,3,1,3,2,3,1,4,3,4,1,2,4,2,)
const IS20_3 = (5,7,6,8,7,6,5,8,6,5,5,9,6,7,5,8,7,5,5,8,5,6,6,9,)
const IS20_4 = (8,9,7,10,9,10,6,9,10,8,7,10,8,10,7,9,10,9,6,10,9,8,7,10,)

const IS21_1 = (5,6,6,7,7,5,5,6,7,5,7,6,)
const IS21_2 = (10,9,9,8,8,10,10,9,8,10,8,9,)
const IS21_3 = (6,5,10,9,5,7,8,7,6,9,10,8,)


# Primary invariants for BA_5B
 # : definitions of the types at the beginning of the file
const P1 = Val((P1_1,))
const P2 = Val((P2_1,))
const P3 = Val((P3_1,))
const P4 = Val((P4_1,))
const P5 = Val((P5_1,P5_2,))
const P6 = Val((P6_1,))
const P7 = Val((P7_1,P7_2,P7_3,))
const P8 = Val((P8_1,))
const P9 = Val((P9_1,))
const P10 = Val((P10_1,P10_2))

# Irreducible secondaries for group BA_5B
 # : definitions of the types at the beginning of the file
const IS1 = Val((IS1_1,IS1_2,))
const IS2 = Val((IS2_1,IS2_2,IS2_3,))
const IS3 = Val((IS3_1,IS3_2,IS3_3,))
const IS4 = Val((IS4_1,IS4_2,))
const IS5 = Val((IS5_1,IS5_2,IS5_3,))
const IS6 = Val((IS6_1,IS6_2,))
const IS7 = Val((IS7_1,IS7_2,IS7_3,))
const IS8 = Val((IS8_1,IS8_2,IS8_3,))
const IS9 = Val((IS9_1,IS9_2,IS9_3,IS9_4,))
const IS10 = Val((IS10_1,IS10_2,IS10_3,IS10_4,))
const IS11 = Val((IS11_1,IS11_2,IS11_3,))
const IS12 = Val((IS12_1,IS12_2,IS12_3,))
const IS13 = Val((IS13_1,IS13_2,IS13_3,))
const IS14 = Val((IS14_1,IS14_2,IS14_3,))
const IS15 = Val((IS15_1,IS15_2,IS15_3,IS15_4,))
const IS16 = Val((IS16_1,IS16_2,IS16_3,IS16_4,))
const IS17 = Val((IS17_1,IS17_2,IS17_3,IS17_4,))
const IS18 = Val((IS18_1,IS18_2,IS18_3,IS18_4,IS18_5,))
const IS19 = Val((IS19_1,IS19_2,IS19_3,IS19_4,))
const IS20 = Val((IS20_1,IS20_2,IS20_3,IS20_4,))
const IS21 = Val((IS21_1,IS21_2,IS21_3,))


function invariants(x1::SVector{10, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1
   x4 = x3.*x1
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_5B
 # : what goes in the function for the evaluation
P1 = fpoly((x1,) , BA_5B.P1)
P2 = fpoly((x1,) , BA_5B.P2)
P3 = fpoly((x2,) , BA_5B.P3)
P4 = fpoly((x2,) , BA_5B.P4)
P5 = fpoly((x1,x1,) , BA_5B.P5)
P6 = fpoly((x3,) , BA_5B.P6)
P7 = fpoly((x1,x1,x1,) , BA_5B.P7)
P8 = fpoly((x3,) , BA_5B.P8)
P9 = fpoly((x4,) , BA_5B.P9)
P10 = fpoly((x3,x1) , BA_5B.P10)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_5B
 # : what goes in the function for the evaluation
IS1 = fpoly((x1,x1,) , BA_5B.IS1)
IS2 = fpoly((x1,x1,x1,) , BA_5B.IS2)
IS3 = fpoly((x1,x1,x1,) , BA_5B.IS3)
IS4 = fpoly((x2,x1,) , BA_5B.IS4)
IS5 = fpoly((x1,x1,x1,) , BA_5B.IS5)
IS6 = fpoly((x2,x1,) , BA_5B.IS6)
IS7 = fpoly((x2,x1,x1,) , BA_5B.IS7)
IS8 = fpoly((x2,x1,x1,) , BA_5B.IS8)
IS9 = fpoly((x1,x1,x1,x1,) , BA_5B.IS9)
IS10 = fpoly((x1,x1,x1,x1,) , BA_5B.IS10)
IS11 = fpoly((x2,x1,x1,) , BA_5B.IS11)
IS12 = fpoly((x2,x1,x1,) , BA_5B.IS12)
IS13 = fpoly((x2,x1,x1,) , BA_5B.IS13)
IS14 = fpoly((x2,x1,x1,) , BA_5B.IS14)
IS15 = fpoly((x2,x1,x1,x1,) , BA_5B.IS15)
IS16 = fpoly((x2,x1,x1,x1,) , BA_5B.IS16)
IS17 = fpoly((x2,x1,x1,x1,) , BA_5B.IS17)
IS18 = fpoly((x1,x1,x1,x1,x1,) , BA_5B.IS18)
IS19 = fpoly((x2,x1,x1,x1,) , BA_5B.IS19)
IS20 = fpoly((x2,x1,x1,x1,) , BA_5B.IS20)
IS21 = fpoly((x2,x2,x1,) , BA_5B.IS21)



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
SEC8  = IS1^2
SEC9  = IS7
SEC10  = IS8
SEC11  = IS9
SEC12  = IS10
SEC13  = IS11
SEC14  = IS12
SEC15  = IS13
SEC16  = IS14
SEC17  = IS1*IS5
SEC18  = IS1*IS4
SEC19  = IS1*IS3
SEC20  = IS1*IS6
SEC21  = IS1*IS2
SEC22  = IS15
SEC23  = IS16
SEC24  = IS17
SEC25  = IS18
SEC26  = IS19
SEC27  = IS20
SEC28  = IS21
SEC29  = IS4^2
SEC30  = IS3*IS6
SEC31  = IS5*IS6
SEC32  = IS1*IS7
SEC33  = IS4*IS6
SEC34  = IS2*IS6
SEC35  = IS1*IS9
SEC36  = IS6^2
SEC37  = IS1*IS8
SEC38  = IS1*IS11
SEC39  = IS1^3
SEC40  = IS3^2
SEC41  = IS3*IS5
SEC42  = IS1*IS12
SEC43  = IS1*IS14
SEC44  = IS3*IS4
SEC45  = IS2*IS3
SEC46  = IS1*IS13
SEC47  = IS5*IS7
SEC48  = IS6*IS13
SEC49  = IS4*IS13
SEC50  = IS5*IS13
SEC51  = IS5*IS10
SEC52  = IS3*IS11
SEC53  = IS3*IS8
SEC54  = IS3*IS9
SEC55  = IS1^2*IS3
SEC56  = IS6*IS11
SEC57  = IS6*IS8
SEC58  = IS2*IS11
SEC59  = IS1*IS18
SEC60  = IS4*IS11
SEC61  = IS1*IS16
SEC62  = IS6*IS9
SEC63  = IS3*IS12
SEC64  = IS1*IS21
SEC65  = IS6*IS12
SEC66  = IS6*IS14
SEC67  = IS1*IS20
SEC68  = IS12*IS13
SEC69  = IS10*IS14
SEC70  = IS6*IS20
SEC71  = IS10^2
SEC72  = IS12*IS14
SEC73  = IS4*IS17
SEC74  = IS3*IS21
SEC75  = IS2*IS17
SEC76  = IS7*IS10
SEC77  = IS9*IS14
SEC78  = IS1*IS3*IS5
SEC79  = IS7*IS12
SEC80  = IS11*IS13
SEC81  = IS1^2*IS13
SEC82  = IS6*IS17
SEC83  = IS7*IS9
SEC84  = IS5*IS15
SEC85  = IS1*IS3*IS4
SEC86  = IS9*IS10
SEC87  = IS5*IS21
SEC88  = IS11^2
SEC89  = IS9^2
SEC90  = IS4*IS21
SEC91  = IS6*IS21
SEC92  = IS8*IS17
SEC93  = IS14*IS18
SEC94  = IS1*IS4*IS7
SEC95  = IS1*IS3*IS11
SEC96  = IS1*IS6*IS14
SEC97  = IS1^3*IS3
SEC98  = IS7*IS19
SEC99  = IS3*IS6^2
SEC100  = IS7*IS18
SEC101  = IS10*IS21
SEC102  = IS1*IS3*IS8
SEC103  = IS11*IS15
SEC104  = IS1*IS6*IS7
SEC105  = IS1*IS5*IS10
SEC106  = IS1^2*IS15
SEC107  = IS12*IS21
SEC108  = IS1*IS5*IS12
SEC109  = IS13*IS20
SEC110  = IS8*IS15
SEC111  = IS9*IS21
SEC112  = IS4^2*IS5
SEC113  = IS1*IS4*IS10
SEC114  = IS10*IS16
SEC115  = IS1^2*IS21
SEC116  = IS1*IS6*IS12
SEC117  = IS1*IS5*IS11
SEC118  = IS1*IS9*IS13
SEC119  = IS1*IS10^2
SEC120  = IS16*IS18
SEC121  = IS4*IS6*IS10
SEC122  = IS4*IS6*IS11
SEC123  = IS1^3*IS8
SEC124  = IS1*IS6*IS20
SEC125  = IS15*IS20
SEC126  = IS1*IS3*IS21
SEC127  = IS5*IS6*IS12
SEC128  = IS1*IS7*IS9
SEC129  = IS17*IS18
SEC130  = IS2*IS4*IS11
SEC131  = IS1*IS2*IS20
SEC132  = IS4^2*IS7
SEC133  = IS2*IS3*IS17
SEC134  = IS1^2*IS6*IS12
SEC135  = IS5*IS8*IS13
SEC136  = IS4^2*IS21
SEC137  = IS1*IS9*IS19
SEC138  = IS3*IS4*IS17
SEC139  = IS1^3*IS17
SEC140  = IS1^3*IS21
SEC141  = IS5*IS11*IS21
SEC142  = IS1*IS4*IS6*IS7
SEC143  = IS2*IS11*IS18
SEC144  = IS5*IS13*IS20


return (@SVector [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,SEC7,SEC8,SEC9,SEC10,SEC11,SEC12,SEC13,SEC14,SEC15,SEC16,SEC17,SEC18,SEC19,SEC20,SEC21,SEC22,SEC23,SEC24,SEC25,SEC26,SEC27,SEC28,SEC29,SEC30,SEC31,SEC32,SEC33,SEC34,SEC35,SEC36,SEC37,SEC38,SEC39,SEC40,SEC41,SEC42,SEC43,SEC44,SEC45,SEC46,SEC47,SEC48,SEC49,SEC50,SEC51,SEC52,SEC53,SEC54,SEC55,SEC56,SEC57,SEC58,SEC59,SEC60,SEC61,SEC62,SEC63,SEC64,SEC65,SEC66,SEC67,SEC68,SEC69,SEC70,SEC71,SEC72,SEC73,SEC74,SEC75,SEC76,SEC77,SEC78,SEC79,SEC80,SEC81,SEC82,SEC83,SEC84,SEC85,SEC86,SEC87,SEC88,SEC89,SEC90,SEC91,SEC92,SEC93,SEC94,SEC95,SEC96,SEC97,SEC98,SEC99,SEC100,SEC101,SEC102,SEC103,SEC104,SEC105,SEC106,SEC107,SEC108,SEC109,SEC110,SEC111,SEC112,SEC113,SEC114,SEC115,SEC116,SEC117,SEC118,SEC119,SEC120,SEC121,SEC122,SEC123,SEC124,SEC125,SEC126,SEC127,SEC128,SEC129,SEC130,SEC131,SEC132,SEC133,SEC134,SEC135,SEC136,SEC137,SEC138,SEC139,SEC140,SEC141,SEC142,SEC143,SEC144,])
 end



function invariants_d(x1::SVector{10, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1
   x4 = x3.*x1

   dx1 = @SVector ones(10)
   dx2 = 2 * x1
   dx3 = 3 * x2
   dx4 = 4 * x3
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_5B
 # : what goes in the function for the derivatives
dP1 = fpoly_d((x1,),(dx1,) , BA_5B.P1)
dP2 = fpoly_d((x1,),(dx1,) , BA_5B.P2)
dP3 = fpoly_d((x2,),(dx2,) , BA_5B.P3)
dP4 = fpoly_d((x2,),(dx2,) , BA_5B.P4)
dP5 = fpoly_d((x1,x1,),(dx1,dx1,) , BA_5B.P5)
dP6 = fpoly_d((x3,),(dx3,) , BA_5B.P6)
dP7 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.P7)
dP8 = fpoly_d((x3,),(dx3,) , BA_5B.P8)
dP9 = fpoly_d((x4,),(dx4,) , BA_5B.P9)
dP10 = fpoly_d((x3,x1,),(dx3,dx1,) , BA_5B.P10)



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_5B
 # : what goes in the function for the evaluation
IS1 = fpoly((x1,x1,) , BA_5B.IS1)
IS2 = fpoly((x1,x1,x1,) , BA_5B.IS2)
IS3 = fpoly((x1,x1,x1,) , BA_5B.IS3)
IS4 = fpoly((x2,x1,) , BA_5B.IS4)
IS5 = fpoly((x1,x1,x1,) , BA_5B.IS5)
IS6 = fpoly((x2,x1,) , BA_5B.IS6)
IS7 = fpoly((x2,x1,x1,) , BA_5B.IS7)
IS8 = fpoly((x2,x1,x1,) , BA_5B.IS8)
IS9 = fpoly((x1,x1,x1,x1,) , BA_5B.IS9)
IS10 = fpoly((x1,x1,x1,x1,) , BA_5B.IS10)
IS11 = fpoly((x2,x1,x1,) , BA_5B.IS11)
IS12 = fpoly((x2,x1,x1,) , BA_5B.IS12)
IS13 = fpoly((x2,x1,x1,) , BA_5B.IS13)
IS14 = fpoly((x2,x1,x1,) , BA_5B.IS14)
IS15 = fpoly((x2,x1,x1,x1,) , BA_5B.IS15)
IS16 = fpoly((x2,x1,x1,x1,) , BA_5B.IS16)
IS17 = fpoly((x2,x1,x1,x1,) , BA_5B.IS17)
IS18 = fpoly((x1,x1,x1,x1,x1,) , BA_5B.IS18)
IS19 = fpoly((x2,x1,x1,x1,) , BA_5B.IS19)
IS20 = fpoly((x2,x1,x1,x1,) , BA_5B.IS20)
IS21 = fpoly((x2,x2,x1,) , BA_5B.IS21)


# Irreducible secondaries for group BA_5B
 # : what goes in the function for the derivatives
dIS1 = fpoly_d((x1,x1,),(dx1,dx1,) , BA_5B.IS1)
dIS2 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.IS2)
dIS3 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.IS3)
dIS4 = fpoly_d((x2,x1,),(dx2,dx1,) , BA_5B.IS4)
dIS5 = fpoly_d((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.IS5)
dIS6 = fpoly_d((x2,x1,),(dx2,dx1,) , BA_5B.IS6)
dIS7 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS7)
dIS8 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS8)
dIS9 = fpoly_d((x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,) , BA_5B.IS9)
dIS10 = fpoly_d((x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,) , BA_5B.IS10)
dIS11 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS11)
dIS12 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS12)
dIS13 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS13)
dIS14 = fpoly_d((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS14)
dIS15 = fpoly_d((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS15)
dIS16 = fpoly_d((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS16)
dIS17 = fpoly_d((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS17)
dIS18 = fpoly_d((x1,x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,dx1,) , BA_5B.IS18)
dIS19 = fpoly_d((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS19)
dIS20 = fpoly_d((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS20)
dIS21 = fpoly_d((x2,x2,x1,),(dx2,dx2,dx1,) , BA_5B.IS21)



   #------------------------------------------------
   # All secondaries
   #------------------------------------------------


dSEC1   = @SVector zeros(10)
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  = dIS4
dSEC6  = dIS5
dSEC7  = dIS6
dSEC8  =  + dIS1*2IS1
dSEC9  = dIS7
dSEC10  = dIS8
dSEC11  = dIS9
dSEC12  = dIS10
dSEC13  = dIS11
dSEC14  = dIS12
dSEC15  = dIS13
dSEC16  = dIS14
dSEC17  =  + dIS1*IS5 + dIS5*IS1
dSEC18  =  + dIS1*IS4 + dIS4*IS1
dSEC19  =  + dIS1*IS3 + dIS3*IS1
dSEC20  =  + dIS1*IS6 + dIS6*IS1
dSEC21  =  + dIS1*IS2 + dIS2*IS1
dSEC22  = dIS15
dSEC23  = dIS16
dSEC24  = dIS17
dSEC25  = dIS18
dSEC26  = dIS19
dSEC27  = dIS20
dSEC28  = dIS21
dSEC29  =  + dIS4*2IS4
dSEC30  =  + dIS3*IS6 + dIS6*IS3
dSEC31  =  + dIS5*IS6 + dIS6*IS5
dSEC32  =  + dIS1*IS7 + dIS7*IS1
dSEC33  =  + dIS4*IS6 + dIS6*IS4
dSEC34  =  + dIS2*IS6 + dIS6*IS2
dSEC35  =  + dIS1*IS9 + dIS9*IS1
dSEC36  =  + dIS6*2IS6
dSEC37  =  + dIS1*IS8 + dIS8*IS1
dSEC38  =  + dIS1*IS11 + dIS11*IS1
dSEC39  =  + dIS1*3 * IS1 ^ 2
dSEC40  =  + dIS3*2IS3
dSEC41  =  + dIS3*IS5 + dIS5*IS3
dSEC42  =  + dIS1*IS12 + dIS12*IS1
dSEC43  =  + dIS1*IS14 + dIS14*IS1
dSEC44  =  + dIS3*IS4 + dIS4*IS3
dSEC45  =  + dIS2*IS3 + dIS3*IS2
dSEC46  =  + dIS1*IS13 + dIS13*IS1
dSEC47  =  + dIS5*IS7 + dIS7*IS5
dSEC48  =  + dIS6*IS13 + dIS13*IS6
dSEC49  =  + dIS4*IS13 + dIS13*IS4
dSEC50  =  + dIS5*IS13 + dIS13*IS5
dSEC51  =  + dIS5*IS10 + dIS10*IS5
dSEC52  =  + dIS3*IS11 + dIS11*IS3
dSEC53  =  + dIS3*IS8 + dIS8*IS3
dSEC54  =  + dIS3*IS9 + dIS9*IS3
dSEC55  =  + dIS1*(2IS1) * IS3 + dIS3*IS1 ^ 2
dSEC56  =  + dIS6*IS11 + dIS11*IS6
dSEC57  =  + dIS6*IS8 + dIS8*IS6
dSEC58  =  + dIS2*IS11 + dIS11*IS2
dSEC59  =  + dIS1*IS18 + dIS18*IS1
dSEC60  =  + dIS4*IS11 + dIS11*IS4
dSEC61  =  + dIS1*IS16 + dIS16*IS1
dSEC62  =  + dIS6*IS9 + dIS9*IS6
dSEC63  =  + dIS3*IS12 + dIS12*IS3
dSEC64  =  + dIS1*IS21 + dIS21*IS1
dSEC65  =  + dIS6*IS12 + dIS12*IS6
dSEC66  =  + dIS6*IS14 + dIS14*IS6
dSEC67  =  + dIS1*IS20 + dIS20*IS1
dSEC68  =  + dIS12*IS13 + dIS13*IS12
dSEC69  =  + dIS10*IS14 + dIS14*IS10
dSEC70  =  + dIS6*IS20 + dIS20*IS6
dSEC71  =  + dIS10*2IS10
dSEC72  =  + dIS12*IS14 + dIS14*IS12
dSEC73  =  + dIS4*IS17 + dIS17*IS4
dSEC74  =  + dIS3*IS21 + dIS21*IS3
dSEC75  =  + dIS2*IS17 + dIS17*IS2
dSEC76  =  + dIS7*IS10 + dIS10*IS7
dSEC77  =  + dIS9*IS14 + dIS14*IS9
dSEC78  =  + dIS1*IS3 * IS5 + dIS3*IS1 * IS5 + dIS5*IS1 * IS3
dSEC79  =  + dIS7*IS12 + dIS12*IS7
dSEC80  =  + dIS11*IS13 + dIS13*IS11
dSEC81  =  + dIS1*(2IS1) * IS13 + dIS13*IS1 ^ 2
dSEC82  =  + dIS6*IS17 + dIS17*IS6
dSEC83  =  + dIS7*IS9 + dIS9*IS7
dSEC84  =  + dIS5*IS15 + dIS15*IS5
dSEC85  =  + dIS1*IS3 * IS4 + dIS3*IS1 * IS4 + dIS4*IS1 * IS3
dSEC86  =  + dIS9*IS10 + dIS10*IS9
dSEC87  =  + dIS5*IS21 + dIS21*IS5
dSEC88  =  + dIS11*2IS11
dSEC89  =  + dIS9*2IS9
dSEC90  =  + dIS4*IS21 + dIS21*IS4
dSEC91  =  + dIS6*IS21 + dIS21*IS6
dSEC92  =  + dIS8*IS17 + dIS17*IS8
dSEC93  =  + dIS14*IS18 + dIS18*IS14
dSEC94  =  + dIS1*IS4 * IS7 + dIS4*IS1 * IS7 + dIS7*IS1 * IS4
dSEC95  =  + dIS1*IS3 * IS11 + dIS3*IS1 * IS11 + dIS11*IS1 * IS3
dSEC96  =  + dIS1*IS6 * IS14 + dIS6*IS1 * IS14 + dIS14*IS1 * IS6
dSEC97  =  + dIS1*(3 * IS1 ^ 2) * IS3 + dIS3*IS1 ^ 3
dSEC98  =  + dIS7*IS19 + dIS19*IS7
dSEC99  =  + dIS3*IS6 ^ 2 + dIS6*IS3 * (2IS6)
dSEC100  =  + dIS7*IS18 + dIS18*IS7
dSEC101  =  + dIS10*IS21 + dIS21*IS10
dSEC102  =  + dIS1*IS3 * IS8 + dIS3*IS1 * IS8 + dIS8*IS1 * IS3
dSEC103  =  + dIS11*IS15 + dIS15*IS11
dSEC104  =  + dIS1*IS6 * IS7 + dIS6*IS1 * IS7 + dIS7*IS1 * IS6
dSEC105  =  + dIS1*IS5 * IS10 + dIS5*IS1 * IS10 + dIS10*IS1 * IS5
dSEC106  =  + dIS1*(2IS1) * IS15 + dIS15*IS1 ^ 2
dSEC107  =  + dIS12*IS21 + dIS21*IS12
dSEC108  =  + dIS1*IS5 * IS12 + dIS5*IS1 * IS12 + dIS12*IS1 * IS5
dSEC109  =  + dIS13*IS20 + dIS20*IS13
dSEC110  =  + dIS8*IS15 + dIS15*IS8
dSEC111  =  + dIS9*IS21 + dIS21*IS9
dSEC112  =  + dIS4*(2IS4) * IS5 + dIS5*IS4 ^ 2
dSEC113  =  + dIS1*IS4 * IS10 + dIS4*IS1 * IS10 + dIS10*IS1 * IS4
dSEC114  =  + dIS10*IS16 + dIS16*IS10
dSEC115  =  + dIS1*(2IS1) * IS21 + dIS21*IS1 ^ 2
dSEC116  =  + dIS1*IS6 * IS12 + dIS6*IS1 * IS12 + dIS12*IS1 * IS6
dSEC117  =  + dIS1*IS5 * IS11 + dIS5*IS1 * IS11 + dIS11*IS1 * IS5
dSEC118  =  + dIS1*IS9 * IS13 + dIS9*IS1 * IS13 + dIS13*IS1 * IS9
dSEC119  =  + dIS1*IS10 ^ 2 + dIS10*IS1 * (2IS10)
dSEC120  =  + dIS16*IS18 + dIS18*IS16
dSEC121  =  + dIS4*IS6 * IS10 + dIS6*IS4 * IS10 + dIS10*IS4 * IS6
dSEC122  =  + dIS4*IS6 * IS11 + dIS6*IS4 * IS11 + dIS11*IS4 * IS6
dSEC123  =  + dIS1*(3 * IS1 ^ 2) * IS8 + dIS8*IS1 ^ 3
dSEC124  =  + dIS1*IS6 * IS20 + dIS6*IS1 * IS20 + dIS20*IS1 * IS6
dSEC125  =  + dIS15*IS20 + dIS20*IS15
dSEC126  =  + dIS1*IS3 * IS21 + dIS3*IS1 * IS21 + dIS21*IS1 * IS3
dSEC127  =  + dIS5*IS6 * IS12 + dIS6*IS5 * IS12 + dIS12*IS5 * IS6
dSEC128  =  + dIS1*IS7 * IS9 + dIS7*IS1 * IS9 + dIS9*IS1 * IS7
dSEC129  =  + dIS17*IS18 + dIS18*IS17
dSEC130  =  + dIS2*IS4 * IS11 + dIS4*IS2 * IS11 + dIS11*IS2 * IS4
dSEC131  =  + dIS1*IS2 * IS20 + dIS2*IS1 * IS20 + dIS20*IS1 * IS2
dSEC132  =  + dIS4*(2IS4) * IS7 + dIS7*IS4 ^ 2
dSEC133  =  + dIS2*IS3 * IS17 + dIS3*IS2 * IS17 + dIS17*IS2 * IS3
dSEC134  =  + dIS1*(2IS1) * IS6 * IS12 + dIS6*IS1 ^ 2 * IS12 + dIS12*IS1 ^ 2 * IS6
dSEC135  =  + dIS5*IS8 * IS13 + dIS8*IS5 * IS13 + dIS13*IS5 * IS8
dSEC136  =  + dIS4*(2IS4) * IS21 + dIS21*IS4 ^ 2
dSEC137  =  + dIS1*IS9 * IS19 + dIS9*IS1 * IS19 + dIS19*IS1 * IS9
dSEC138  =  + dIS3*IS4 * IS17 + dIS4*IS3 * IS17 + dIS17*IS3 * IS4
dSEC139  =  + dIS1*(3 * IS1 ^ 2) * IS17 + dIS17*IS1 ^ 3
dSEC140  =  + dIS1*(3 * IS1 ^ 2) * IS21 + dIS21*IS1 ^ 3
dSEC141  =  + dIS5*IS11 * IS21 + dIS11*IS5 * IS21 + dIS21*IS5 * IS11
dSEC142  =  + dIS1*IS4 * IS6 * IS7 + dIS4*IS1 * IS6 * IS7 + dIS6*IS1 * IS4 * IS7 + dIS7*IS1 * IS4 * IS6
dSEC143  =  + dIS2*IS11 * IS18 + dIS11*IS2 * IS18 + dIS18*IS2 * IS11
dSEC144  =  + dIS5*IS13 * IS20 + dIS13*IS5 * IS20 + dIS20*IS5 * IS13


return (dP1,dP2,dP3,dP4,dP5,dP6,dP7,dP8,dP9,dP10,), (dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,dSEC7,dSEC8,dSEC9,dSEC10,dSEC11,dSEC12,dSEC13,dSEC14,dSEC15,dSEC16,dSEC17,dSEC18,dSEC19,dSEC20,dSEC21,dSEC22,dSEC23,dSEC24,dSEC25,dSEC26,dSEC27,dSEC28,dSEC29,dSEC30,dSEC31,dSEC32,dSEC33,dSEC34,dSEC35,dSEC36,dSEC37,dSEC38,dSEC39,dSEC40,dSEC41,dSEC42,dSEC43,dSEC44,dSEC45,dSEC46,dSEC47,dSEC48,dSEC49,dSEC50,dSEC51,dSEC52,dSEC53,dSEC54,dSEC55,dSEC56,dSEC57,dSEC58,dSEC59,dSEC60,dSEC61,dSEC62,dSEC63,dSEC64,dSEC65,dSEC66,dSEC67,dSEC68,dSEC69,dSEC70,dSEC71,dSEC72,dSEC73,dSEC74,dSEC75,dSEC76,dSEC77,dSEC78,dSEC79,dSEC80,dSEC81,dSEC82,dSEC83,dSEC84,dSEC85,dSEC86,dSEC87,dSEC88,dSEC89,dSEC90,dSEC91,dSEC92,dSEC93,dSEC94,dSEC95,dSEC96,dSEC97,dSEC98,dSEC99,dSEC100,dSEC101,dSEC102,dSEC103,dSEC104,dSEC105,dSEC106,dSEC107,dSEC108,dSEC109,dSEC110,dSEC111,dSEC112,dSEC113,dSEC114,dSEC115,dSEC116,dSEC117,dSEC118,dSEC119,dSEC120,dSEC121,dSEC122,dSEC123,dSEC124,dSEC125,dSEC126,dSEC127,dSEC128,dSEC129,dSEC130,dSEC131,dSEC132,dSEC133,dSEC134,dSEC135,dSEC136,dSEC137,dSEC138,dSEC139,dSEC140,dSEC141,dSEC142,dSEC143,dSEC144,)
 end



function invariants_ed(x1::SVector{10, T}) where {T}
   x2 = x1.*x1
   x3 = x2.*x1
   x4 = x3.*x1

   dx1 = @SVector ones(10)
   dx2 = 2 * x1
   dx3 = 3 * x2
   dx4 = 4 * x3
   #------------------------------------------------
   # Primaries
   #------------------------------------------------

# Primary invariants for BA_5B
 # : what goes in the function for the evaluation and derivatives
P1, dP1 = fpoly_ed((x1,),(dx1,) , BA_5B.P1)
P2, dP2 = fpoly_ed((x1,),(dx1,) , BA_5B.P2)
P3, dP3 = fpoly_ed((x2,),(dx2,) , BA_5B.P3)
P4, dP4 = fpoly_ed((x2,),(dx2,) , BA_5B.P4)
P5, dP5 = fpoly_ed((x1,x1,),(dx1,dx1,) , BA_5B.P5)
P6, dP6 = fpoly_ed((x3,),(dx3,) , BA_5B.P6)
P7, dP7 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.P7)
P8, dP8 = fpoly_ed((x3,),(dx3,) , BA_5B.P8)
P9, dP9 = fpoly_ed((x4,),(dx4,) , BA_5B.P9)
P10, dP10 = fpoly_ed((x3,x1),(dx3,dx1) , BA_5B.P10) 



   #------------------------------------------------
   # Irreducible secondaries
   #------------------------------------------------


# Irreducible secondaries for group BA_5B
 # : what goes in the function for the evaluation and derivatives
IS1, dIS1 = fpoly_ed((x1,x1,),(dx1,dx1,) , BA_5B.IS1)
IS2, dIS2 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.IS2)
IS3, dIS3 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.IS3)
IS4, dIS4 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_5B.IS4)
IS5, dIS5 = fpoly_ed((x1,x1,x1,),(dx1,dx1,dx1,) , BA_5B.IS5)
IS6, dIS6 = fpoly_ed((x2,x1,),(dx2,dx1,) , BA_5B.IS6)
IS7, dIS7 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS7)
IS8, dIS8 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS8)
IS9, dIS9 = fpoly_ed((x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,) , BA_5B.IS9)
IS10, dIS10 = fpoly_ed((x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,) , BA_5B.IS10)
IS11, dIS11 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS11)
IS12, dIS12 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS12)
IS13, dIS13 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS13)
IS14, dIS14 = fpoly_ed((x2,x1,x1,),(dx2,dx1,dx1,) , BA_5B.IS14)
IS15, dIS15 = fpoly_ed((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS15)
IS16, dIS16 = fpoly_ed((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS16)
IS17, dIS17 = fpoly_ed((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS17)
IS18, dIS18 = fpoly_ed((x1,x1,x1,x1,x1,),(dx1,dx1,dx1,dx1,dx1,) , BA_5B.IS18)
IS19, dIS19 = fpoly_ed((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS19)
IS20, dIS20 = fpoly_ed((x2,x1,x1,x1,),(dx2,dx1,dx1,dx1,) , BA_5B.IS20)
IS21, dIS21 = fpoly_ed((x2,x2,x1,),(dx2,dx2,dx1,) , BA_5B.IS21)



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
SEC8  = IS1^2
SEC9  = IS7
SEC10  = IS8
SEC11  = IS9
SEC12  = IS10
SEC13  = IS11
SEC14  = IS12
SEC15  = IS13
SEC16  = IS14
SEC17  = IS1*IS5
SEC18  = IS1*IS4
SEC19  = IS1*IS3
SEC20  = IS1*IS6
SEC21  = IS1*IS2
SEC22  = IS15
SEC23  = IS16
SEC24  = IS17
SEC25  = IS18
SEC26  = IS19
SEC27  = IS20
SEC28  = IS21
SEC29  = IS4^2
SEC30  = IS3*IS6
SEC31  = IS5*IS6
SEC32  = IS1*IS7
SEC33  = IS4*IS6
SEC34  = IS2*IS6
SEC35  = IS1*IS9
SEC36  = IS6^2
SEC37  = IS1*IS8
SEC38  = IS1*IS11
SEC39  = IS1^3
SEC40  = IS3^2
SEC41  = IS3*IS5
SEC42  = IS1*IS12
SEC43  = IS1*IS14
SEC44  = IS3*IS4
SEC45  = IS2*IS3
SEC46  = IS1*IS13
SEC47  = IS5*IS7
SEC48  = IS6*IS13
SEC49  = IS4*IS13
SEC50  = IS5*IS13
SEC51  = IS5*IS10
SEC52  = IS3*IS11
SEC53  = IS3*IS8
SEC54  = IS3*IS9
SEC55  = IS1^2*IS3
SEC56  = IS6*IS11
SEC57  = IS6*IS8
SEC58  = IS2*IS11
SEC59  = IS1*IS18
SEC60  = IS4*IS11
SEC61  = IS1*IS16
SEC62  = IS6*IS9
SEC63  = IS3*IS12
SEC64  = IS1*IS21
SEC65  = IS6*IS12
SEC66  = IS6*IS14
SEC67  = IS1*IS20
SEC68  = IS12*IS13
SEC69  = IS10*IS14
SEC70  = IS6*IS20
SEC71  = IS10^2
SEC72  = IS12*IS14
SEC73  = IS4*IS17
SEC74  = IS3*IS21
SEC75  = IS2*IS17
SEC76  = IS7*IS10
SEC77  = IS9*IS14
SEC78  = IS1*IS3*IS5
SEC79  = IS7*IS12
SEC80  = IS11*IS13
SEC81  = IS1^2*IS13
SEC82  = IS6*IS17
SEC83  = IS7*IS9
SEC84  = IS5*IS15
SEC85  = IS1*IS3*IS4
SEC86  = IS9*IS10
SEC87  = IS5*IS21
SEC88  = IS11^2
SEC89  = IS9^2
SEC90  = IS4*IS21
SEC91  = IS6*IS21
SEC92  = IS8*IS17
SEC93  = IS14*IS18
SEC94  = IS1*IS4*IS7
SEC95  = IS1*IS3*IS11
SEC96  = IS1*IS6*IS14
SEC97  = IS1^3*IS3
SEC98  = IS7*IS19
SEC99  = IS3*IS6^2
SEC100  = IS7*IS18
SEC101  = IS10*IS21
SEC102  = IS1*IS3*IS8
SEC103  = IS11*IS15
SEC104  = IS1*IS6*IS7
SEC105  = IS1*IS5*IS10
SEC106  = IS1^2*IS15
SEC107  = IS12*IS21
SEC108  = IS1*IS5*IS12
SEC109  = IS13*IS20
SEC110  = IS8*IS15
SEC111  = IS9*IS21
SEC112  = IS4^2*IS5
SEC113  = IS1*IS4*IS10
SEC114  = IS10*IS16
SEC115  = IS1^2*IS21
SEC116  = IS1*IS6*IS12
SEC117  = IS1*IS5*IS11
SEC118  = IS1*IS9*IS13
SEC119  = IS1*IS10^2
SEC120  = IS16*IS18
SEC121  = IS4*IS6*IS10
SEC122  = IS4*IS6*IS11
SEC123  = IS1^3*IS8
SEC124  = IS1*IS6*IS20
SEC125  = IS15*IS20
SEC126  = IS1*IS3*IS21
SEC127  = IS5*IS6*IS12
SEC128  = IS1*IS7*IS9
SEC129  = IS17*IS18
SEC130  = IS2*IS4*IS11
SEC131  = IS1*IS2*IS20
SEC132  = IS4^2*IS7
SEC133  = IS2*IS3*IS17
SEC134  = IS1^2*IS6*IS12
SEC135  = IS5*IS8*IS13
SEC136  = IS4^2*IS21
SEC137  = IS1*IS9*IS19
SEC138  = IS3*IS4*IS17
SEC139  = IS1^3*IS17
SEC140  = IS1^3*IS21
SEC141  = IS5*IS11*IS21
SEC142  = IS1*IS4*IS6*IS7
SEC143  = IS2*IS11*IS18
SEC144  = IS5*IS13*IS20


dSEC1   = @SVector zeros(10)
dSEC2  = dIS1
dSEC3  = dIS2
dSEC4  = dIS3
dSEC5  = dIS4
dSEC6  = dIS5
dSEC7  = dIS6
dSEC8  =  + dIS1*2IS1
dSEC9  = dIS7
dSEC10  = dIS8
dSEC11  = dIS9
dSEC12  = dIS10
dSEC13  = dIS11
dSEC14  = dIS12
dSEC15  = dIS13
dSEC16  = dIS14
dSEC17  =  + dIS1*IS5 + dIS5*IS1
dSEC18  =  + dIS1*IS4 + dIS4*IS1
dSEC19  =  + dIS1*IS3 + dIS3*IS1
dSEC20  =  + dIS1*IS6 + dIS6*IS1
dSEC21  =  + dIS1*IS2 + dIS2*IS1
dSEC22  = dIS15
dSEC23  = dIS16
dSEC24  = dIS17
dSEC25  = dIS18
dSEC26  = dIS19
dSEC27  = dIS20
dSEC28  = dIS21
dSEC29  =  + dIS4*2IS4
dSEC30  =  + dIS3*IS6 + dIS6*IS3
dSEC31  =  + dIS5*IS6 + dIS6*IS5
dSEC32  =  + dIS1*IS7 + dIS7*IS1
dSEC33  =  + dIS4*IS6 + dIS6*IS4
dSEC34  =  + dIS2*IS6 + dIS6*IS2
dSEC35  =  + dIS1*IS9 + dIS9*IS1
dSEC36  =  + dIS6*2IS6
dSEC37  =  + dIS1*IS8 + dIS8*IS1
dSEC38  =  + dIS1*IS11 + dIS11*IS1
dSEC39  =  + dIS1*3 * IS1 ^ 2
dSEC40  =  + dIS3*2IS3
dSEC41  =  + dIS3*IS5 + dIS5*IS3
dSEC42  =  + dIS1*IS12 + dIS12*IS1
dSEC43  =  + dIS1*IS14 + dIS14*IS1
dSEC44  =  + dIS3*IS4 + dIS4*IS3
dSEC45  =  + dIS2*IS3 + dIS3*IS2
dSEC46  =  + dIS1*IS13 + dIS13*IS1
dSEC47  =  + dIS5*IS7 + dIS7*IS5
dSEC48  =  + dIS6*IS13 + dIS13*IS6
dSEC49  =  + dIS4*IS13 + dIS13*IS4
dSEC50  =  + dIS5*IS13 + dIS13*IS5
dSEC51  =  + dIS5*IS10 + dIS10*IS5
dSEC52  =  + dIS3*IS11 + dIS11*IS3
dSEC53  =  + dIS3*IS8 + dIS8*IS3
dSEC54  =  + dIS3*IS9 + dIS9*IS3
dSEC55  =  + dIS1*(2IS1) * IS3 + dIS3*IS1 ^ 2
dSEC56  =  + dIS6*IS11 + dIS11*IS6
dSEC57  =  + dIS6*IS8 + dIS8*IS6
dSEC58  =  + dIS2*IS11 + dIS11*IS2
dSEC59  =  + dIS1*IS18 + dIS18*IS1
dSEC60  =  + dIS4*IS11 + dIS11*IS4
dSEC61  =  + dIS1*IS16 + dIS16*IS1
dSEC62  =  + dIS6*IS9 + dIS9*IS6
dSEC63  =  + dIS3*IS12 + dIS12*IS3
dSEC64  =  + dIS1*IS21 + dIS21*IS1
dSEC65  =  + dIS6*IS12 + dIS12*IS6
dSEC66  =  + dIS6*IS14 + dIS14*IS6
dSEC67  =  + dIS1*IS20 + dIS20*IS1
dSEC68  =  + dIS12*IS13 + dIS13*IS12
dSEC69  =  + dIS10*IS14 + dIS14*IS10
dSEC70  =  + dIS6*IS20 + dIS20*IS6
dSEC71  =  + dIS10*2IS10
dSEC72  =  + dIS12*IS14 + dIS14*IS12
dSEC73  =  + dIS4*IS17 + dIS17*IS4
dSEC74  =  + dIS3*IS21 + dIS21*IS3
dSEC75  =  + dIS2*IS17 + dIS17*IS2
dSEC76  =  + dIS7*IS10 + dIS10*IS7
dSEC77  =  + dIS9*IS14 + dIS14*IS9
dSEC78  =  + dIS1*IS3 * IS5 + dIS3*IS1 * IS5 + dIS5*IS1 * IS3
dSEC79  =  + dIS7*IS12 + dIS12*IS7
dSEC80  =  + dIS11*IS13 + dIS13*IS11
dSEC81  =  + dIS1*(2IS1) * IS13 + dIS13*IS1 ^ 2
dSEC82  =  + dIS6*IS17 + dIS17*IS6
dSEC83  =  + dIS7*IS9 + dIS9*IS7
dSEC84  =  + dIS5*IS15 + dIS15*IS5
dSEC85  =  + dIS1*IS3 * IS4 + dIS3*IS1 * IS4 + dIS4*IS1 * IS3
dSEC86  =  + dIS9*IS10 + dIS10*IS9
dSEC87  =  + dIS5*IS21 + dIS21*IS5
dSEC88  =  + dIS11*2IS11
dSEC89  =  + dIS9*2IS9
dSEC90  =  + dIS4*IS21 + dIS21*IS4
dSEC91  =  + dIS6*IS21 + dIS21*IS6
dSEC92  =  + dIS8*IS17 + dIS17*IS8
dSEC93  =  + dIS14*IS18 + dIS18*IS14
dSEC94  =  + dIS1*IS4 * IS7 + dIS4*IS1 * IS7 + dIS7*IS1 * IS4
dSEC95  =  + dIS1*IS3 * IS11 + dIS3*IS1 * IS11 + dIS11*IS1 * IS3
dSEC96  =  + dIS1*IS6 * IS14 + dIS6*IS1 * IS14 + dIS14*IS1 * IS6
dSEC97  =  + dIS1*(3 * IS1 ^ 2) * IS3 + dIS3*IS1 ^ 3
dSEC98  =  + dIS7*IS19 + dIS19*IS7
dSEC99  =  + dIS3*IS6 ^ 2 + dIS6*IS3 * (2IS6)
dSEC100  =  + dIS7*IS18 + dIS18*IS7
dSEC101  =  + dIS10*IS21 + dIS21*IS10
dSEC102  =  + dIS1*IS3 * IS8 + dIS3*IS1 * IS8 + dIS8*IS1 * IS3
dSEC103  =  + dIS11*IS15 + dIS15*IS11
dSEC104  =  + dIS1*IS6 * IS7 + dIS6*IS1 * IS7 + dIS7*IS1 * IS6
dSEC105  =  + dIS1*IS5 * IS10 + dIS5*IS1 * IS10 + dIS10*IS1 * IS5
dSEC106  =  + dIS1*(2IS1) * IS15 + dIS15*IS1 ^ 2
dSEC107  =  + dIS12*IS21 + dIS21*IS12
dSEC108  =  + dIS1*IS5 * IS12 + dIS5*IS1 * IS12 + dIS12*IS1 * IS5
dSEC109  =  + dIS13*IS20 + dIS20*IS13
dSEC110  =  + dIS8*IS15 + dIS15*IS8
dSEC111  =  + dIS9*IS21 + dIS21*IS9
dSEC112  =  + dIS4*(2IS4) * IS5 + dIS5*IS4 ^ 2
dSEC113  =  + dIS1*IS4 * IS10 + dIS4*IS1 * IS10 + dIS10*IS1 * IS4
dSEC114  =  + dIS10*IS16 + dIS16*IS10
dSEC115  =  + dIS1*(2IS1) * IS21 + dIS21*IS1 ^ 2
dSEC116  =  + dIS1*IS6 * IS12 + dIS6*IS1 * IS12 + dIS12*IS1 * IS6
dSEC117  =  + dIS1*IS5 * IS11 + dIS5*IS1 * IS11 + dIS11*IS1 * IS5
dSEC118  =  + dIS1*IS9 * IS13 + dIS9*IS1 * IS13 + dIS13*IS1 * IS9
dSEC119  =  + dIS1*IS10 ^ 2 + dIS10*IS1 * (2IS10)
dSEC120  =  + dIS16*IS18 + dIS18*IS16
dSEC121  =  + dIS4*IS6 * IS10 + dIS6*IS4 * IS10 + dIS10*IS4 * IS6
dSEC122  =  + dIS4*IS6 * IS11 + dIS6*IS4 * IS11 + dIS11*IS4 * IS6
dSEC123  =  + dIS1*(3 * IS1 ^ 2) * IS8 + dIS8*IS1 ^ 3
dSEC124  =  + dIS1*IS6 * IS20 + dIS6*IS1 * IS20 + dIS20*IS1 * IS6
dSEC125  =  + dIS15*IS20 + dIS20*IS15
dSEC126  =  + dIS1*IS3 * IS21 + dIS3*IS1 * IS21 + dIS21*IS1 * IS3
dSEC127  =  + dIS5*IS6 * IS12 + dIS6*IS5 * IS12 + dIS12*IS5 * IS6
dSEC128  =  + dIS1*IS7 * IS9 + dIS7*IS1 * IS9 + dIS9*IS1 * IS7
dSEC129  =  + dIS17*IS18 + dIS18*IS17
dSEC130  =  + dIS2*IS4 * IS11 + dIS4*IS2 * IS11 + dIS11*IS2 * IS4
dSEC131  =  + dIS1*IS2 * IS20 + dIS2*IS1 * IS20 + dIS20*IS1 * IS2
dSEC132  =  + dIS4*(2IS4) * IS7 + dIS7*IS4 ^ 2
dSEC133  =  + dIS2*IS3 * IS17 + dIS3*IS2 * IS17 + dIS17*IS2 * IS3
dSEC134  =  + dIS1*(2IS1) * IS6 * IS12 + dIS6*IS1 ^ 2 * IS12 + dIS12*IS1 ^ 2 * IS6
dSEC135  =  + dIS5*IS8 * IS13 + dIS8*IS5 * IS13 + dIS13*IS5 * IS8
dSEC136  =  + dIS4*(2IS4) * IS21 + dIS21*IS4 ^ 2
dSEC137  =  + dIS1*IS9 * IS19 + dIS9*IS1 * IS19 + dIS19*IS1 * IS9
dSEC138  =  + dIS3*IS4 * IS17 + dIS4*IS3 * IS17 + dIS17*IS3 * IS4
dSEC139  =  + dIS1*(3 * IS1 ^ 2) * IS17 + dIS17*IS1 ^ 3
dSEC140  =  + dIS1*(3 * IS1 ^ 2) * IS21 + dIS21*IS1 ^ 3
dSEC141  =  + dIS5*IS11 * IS21 + dIS11*IS5 * IS21 + dIS21*IS5 * IS11
dSEC142  =  + dIS1*IS4 * IS6 * IS7 + dIS4*IS1 * IS6 * IS7 + dIS6*IS1 * IS4 * IS7 + dIS7*IS1 * IS4 * IS6
dSEC143  =  + dIS2*IS11 * IS18 + dIS11*IS2 * IS18 + dIS18*IS2 * IS11
dSEC144  =  + dIS5*IS13 * IS20 + dIS13*IS5 * IS20 + dIS20*IS5 * IS13


return (@SVector [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,]), (@SVector [SEC1,SEC2,SEC3,SEC4,SEC5,SEC6,SEC7,SEC8,SEC9,SEC10,SEC11,SEC12,SEC13,SEC14,SEC15,SEC16,SEC17,SEC18,SEC19,SEC20,SEC21,SEC22,SEC23,SEC24,SEC25,SEC26,SEC27,SEC28,SEC29,SEC30,SEC31,SEC32,SEC33,SEC34,SEC35,SEC36,SEC37,SEC38,SEC39,SEC40,SEC41,SEC42,SEC43,SEC44,SEC45,SEC46,SEC47,SEC48,SEC49,SEC50,SEC51,SEC52,SEC53,SEC54,SEC55,SEC56,SEC57,SEC58,SEC59,SEC60,SEC61,SEC62,SEC63,SEC64,SEC65,SEC66,SEC67,SEC68,SEC69,SEC70,SEC71,SEC72,SEC73,SEC74,SEC75,SEC76,SEC77,SEC78,SEC79,SEC80,SEC81,SEC82,SEC83,SEC84,SEC85,SEC86,SEC87,SEC88,SEC89,SEC90,SEC91,SEC92,SEC93,SEC94,SEC95,SEC96,SEC97,SEC98,SEC99,SEC100,SEC101,SEC102,SEC103,SEC104,SEC105,SEC106,SEC107,SEC108,SEC109,SEC110,SEC111,SEC112,SEC113,SEC114,SEC115,SEC116,SEC117,SEC118,SEC119,SEC120,SEC121,SEC122,SEC123,SEC124,SEC125,SEC126,SEC127,SEC128,SEC129,SEC130,SEC131,SEC132,SEC133,SEC134,SEC135,SEC136,SEC137,SEC138,SEC139,SEC140,SEC141,SEC142,SEC143,SEC144,]), (@SVector [dP1,dP2,dP3,dP4,dP5,dP6,dP7,dP8,dP9,dP10,]), (@SVector [dSEC1,dSEC2,dSEC3,dSEC4,dSEC5,dSEC6,dSEC7,dSEC8,dSEC9,dSEC10,dSEC11,dSEC12,dSEC13,dSEC14,dSEC15,dSEC16,dSEC17,dSEC18,dSEC19,dSEC20,dSEC21,dSEC22,dSEC23,dSEC24,dSEC25,dSEC26,dSEC27,dSEC28,dSEC29,dSEC30,dSEC31,dSEC32,dSEC33,dSEC34,dSEC35,dSEC36,dSEC37,dSEC38,dSEC39,dSEC40,dSEC41,dSEC42,dSEC43,dSEC44,dSEC45,dSEC46,dSEC47,dSEC48,dSEC49,dSEC50,dSEC51,dSEC52,dSEC53,dSEC54,dSEC55,dSEC56,dSEC57,dSEC58,dSEC59,dSEC60,dSEC61,dSEC62,dSEC63,dSEC64,dSEC65,dSEC66,dSEC67,dSEC68,dSEC69,dSEC70,dSEC71,dSEC72,dSEC73,dSEC74,dSEC75,dSEC76,dSEC77,dSEC78,dSEC79,dSEC80,dSEC81,dSEC82,dSEC83,dSEC84,dSEC85,dSEC86,dSEC87,dSEC88,dSEC89,dSEC90,dSEC91,dSEC92,dSEC93,dSEC94,dSEC95,dSEC96,dSEC97,dSEC98,dSEC99,dSEC100,dSEC101,dSEC102,dSEC103,dSEC104,dSEC105,dSEC106,dSEC107,dSEC108,dSEC109,dSEC110,dSEC111,dSEC112,dSEC113,dSEC114,dSEC115,dSEC116,dSEC117,dSEC118,dSEC119,dSEC120,dSEC121,dSEC122,dSEC123,dSEC124,dSEC125,dSEC126,dSEC127,dSEC128,dSEC129,dSEC130,dSEC131,dSEC132,dSEC133,dSEC134,dSEC135,dSEC136,dSEC137,dSEC138,dSEC139,dSEC140,dSEC141,dSEC142,dSEC143,dSEC144,])
 end

end
