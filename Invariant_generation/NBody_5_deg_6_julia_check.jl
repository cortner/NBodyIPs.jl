function invariants_Q10_check(x)
pv=zeros(27,1);

v=zeros(31,1);
 v[1] = 1
 pv[1] = x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^2*x[5] +
x[2]^2*x[5] + x[1]*x[5]^2 + x[2]*x[5]^2 +
x[1]^2*x[3] + x[2]^2*x[3] + x[1]*x[3]^2 +
x[2]*x[3]^2 + x[1]^2*x[6] + x[5]^2*x[6] +
x[3]^2*x[6] + x[1]*x[6]^2 + x[5]*x[6]^2 +
x[3]*x[6]^2 + x[2]^2*x[8] + x[5]^2*x[8] +
x[3]^2*x[8] + x[6]^2*x[8] + x[2]*x[8]^2 +
x[5]*x[8]^2 + x[3]*x[8]^2 + x[6]*x[8]^2 +
x[1]^2*x[4] + x[2]^2*x[4] + x[3]^2*x[4] +
x[1]*x[4]^2 + x[2]*x[4]^2 + x[3]*x[4]^2 +
x[1]^2*x[7] + x[5]^2*x[7] + x[6]^2*x[7] +
x[4]^2*x[7] + x[1]*x[7]^2 + x[5]*x[7]^2 +
x[6]*x[7]^2 + x[4]*x[7]^2 + x[2]^2*x[9] +
x[5]^2*x[9] + x[8]^2*x[9] + x[4]^2*x[9] +
x[7]^2*x[9] + x[2]*x[9]^2 + x[5]*x[9]^2 +
x[8]*x[9]^2 + x[4]*x[9]^2 + x[7]*x[9]^2 +
x[3]^2*x[10] + x[6]^2*x[10] + x[8]^2*x[10] +
x[4]^2*x[10] + x[7]^2*x[10] + x[9]^2*x[10] +
x[3]*x[10]^2 + x[6]*x[10]^2 + x[8]*x[10]^2 +
x[4]*x[10]^2 + x[7]*x[10]^2 + x[9]*x[10]^2
 pv[2] = x[1]*x[2]*x[5] + x[1]*x[3]*x[6] +
x[2]*x[3]*x[8] + x[5]*x[6]*x[8] +
x[1]*x[4]*x[7] + x[2]*x[4]*x[9] +
x[5]*x[7]*x[9] + x[3]*x[4]*x[10] +
x[6]*x[7]*x[10] + x[8]*x[9]*x[10]
 v[2:3] = pv[1:2]
 pv[3] = x[1]^3*x[2] + x[1]*x[2]^3 + x[1]^3*x[5] +
x[2]^3*x[5] + x[1]*x[5]^3 + x[2]*x[5]^3 +
x[1]^3*x[3] + x[2]^3*x[3] + x[1]*x[3]^3 +
x[2]*x[3]^3 + x[1]^3*x[6] + x[5]^3*x[6] +
x[3]^3*x[6] + x[1]*x[6]^3 + x[5]*x[6]^3 +
x[3]*x[6]^3 + x[2]^3*x[8] + x[5]^3*x[8] +
x[3]^3*x[8] + x[6]^3*x[8] + x[2]*x[8]^3 +
x[5]*x[8]^3 + x[3]*x[8]^3 + x[6]*x[8]^3 +
x[1]^3*x[4] + x[2]^3*x[4] + x[3]^3*x[4] +
x[1]*x[4]^3 + x[2]*x[4]^3 + x[3]*x[4]^3 +
x[1]^3*x[7] + x[5]^3*x[7] + x[6]^3*x[7] +
x[4]^3*x[7] + x[1]*x[7]^3 + x[5]*x[7]^3 +
x[6]*x[7]^3 + x[4]*x[7]^3 + x[2]^3*x[9] +
x[5]^3*x[9] + x[8]^3*x[9] + x[4]^3*x[9] +
x[7]^3*x[9] + x[2]*x[9]^3 + x[5]*x[9]^3 +
x[8]*x[9]^3 + x[4]*x[9]^3 + x[7]*x[9]^3 +
x[3]^3*x[10] + x[6]^3*x[10] + x[8]^3*x[10] +
x[4]^3*x[10] + x[7]^3*x[10] + x[9]^3*x[10] +
x[3]*x[10]^3 + x[6]*x[10]^3 + x[8]*x[10]^3 +
x[4]*x[10]^3 + x[7]*x[10]^3 + x[9]*x[10]^3
 pv[4] = x[1]^2*x[2]^2 + x[1]^2*x[5]^2 + x[2]^2*x[5]^2 +
x[1]^2*x[3]^2 + x[2]^2*x[3]^2 + x[1]^2*x[6]^2 +
x[5]^2*x[6]^2 + x[3]^2*x[6]^2 + x[2]^2*x[8]^2 +
x[5]^2*x[8]^2 + x[3]^2*x[8]^2 + x[6]^2*x[8]^2 +
x[1]^2*x[4]^2 + x[2]^2*x[4]^2 + x[3]^2*x[4]^2 +
x[1]^2*x[7]^2 + x[5]^2*x[7]^2 + x[6]^2*x[7]^2 +
x[4]^2*x[7]^2 + x[2]^2*x[9]^2 + x[5]^2*x[9]^2 +
x[8]^2*x[9]^2 + x[4]^2*x[9]^2 + x[7]^2*x[9]^2 +
x[3]^2*x[10]^2 + x[6]^2*x[10]^2 + x[8]^2*x[10]^2 +
x[4]^2*x[10]^2 + x[7]^2*x[10]^2 + x[9]^2*x[10]^2
 pv[5] = x[1]^2*x[2]*x[5] + x[1]*x[2]^2*x[5] +
x[1]*x[2]*x[5]^2 + x[1]^2*x[3]*x[6] +
x[1]*x[3]^2*x[6] + x[1]*x[3]*x[6]^2 +
x[2]^2*x[3]*x[8] + x[2]*x[3]^2*x[8] +
x[5]^2*x[6]*x[8] + x[5]*x[6]^2*x[8] +
x[2]*x[3]*x[8]^2 + x[5]*x[6]*x[8]^2 +
x[1]^2*x[4]*x[7] + x[1]*x[4]^2*x[7] +
x[1]*x[4]*x[7]^2 + x[2]^2*x[4]*x[9] +
x[2]*x[4]^2*x[9] + x[5]^2*x[7]*x[9] +
x[5]*x[7]^2*x[9] + x[2]*x[4]*x[9]^2 +
x[5]*x[7]*x[9]^2 + x[3]^2*x[4]*x[10] +
x[3]*x[4]^2*x[10] + x[6]^2*x[7]*x[10] +
x[6]*x[7]^2*x[10] + x[8]^2*x[9]*x[10] +
x[8]*x[9]^2*x[10] + x[3]*x[4]*x[10]^2 +
x[6]*x[7]*x[10]^2 + x[8]*x[9]*x[10]^2
 pv[6] = x[1]^2*x[2]*x[3] + x[1]*x[2]^2*x[3] +
x[1]*x[2]*x[3]^2 + x[1]^2*x[5]*x[6] +
x[1]*x[5]^2*x[6] + x[1]*x[5]*x[6]^2 +
x[2]^2*x[5]*x[8] + x[2]*x[5]^2*x[8] +
x[3]^2*x[6]*x[8] + x[3]*x[6]^2*x[8] +
x[2]*x[5]*x[8]^2 + x[3]*x[6]*x[8]^2 +
x[1]^2*x[2]*x[4] + x[1]*x[2]^2*x[4] +
x[1]^2*x[3]*x[4] + x[2]^2*x[3]*x[4] +
x[1]*x[3]^2*x[4] + x[2]*x[3]^2*x[4] +
x[1]*x[2]*x[4]^2 + x[1]*x[3]*x[4]^2 +
x[2]*x[3]*x[4]^2 + x[1]^2*x[5]*x[7] +
x[1]*x[5]^2*x[7] + x[1]^2*x[6]*x[7] +
x[5]^2*x[6]*x[7] + x[1]*x[6]^2*x[7] +
x[5]*x[6]^2*x[7] + x[1]*x[5]*x[7]^2 +
x[1]*x[6]*x[7]^2 + x[5]*x[6]*x[7]^2 +
x[2]^2*x[5]*x[9] + x[2]*x[5]^2*x[9] +
x[2]^2*x[8]*x[9] + x[5]^2*x[8]*x[9] +
x[2]*x[8]^2*x[9] + x[5]*x[8]^2*x[9] +
x[4]^2*x[7]*x[9] + x[4]*x[7]^2*x[9] +
x[2]*x[5]*x[9]^2 + x[2]*x[8]*x[9]^2 +
x[5]*x[8]*x[9]^2 + x[4]*x[7]*x[9]^2 +
x[3]^2*x[6]*x[10] + x[3]*x[6]^2*x[10] +
x[3]^2*x[8]*x[10] + x[6]^2*x[8]*x[10] +
x[3]*x[8]^2*x[10] + x[6]*x[8]^2*x[10] +
x[4]^2*x[7]*x[10] + x[4]*x[7]^2*x[10] +
x[4]^2*x[9]*x[10] + x[7]^2*x[9]*x[10] +
x[4]*x[9]^2*x[10] + x[7]*x[9]^2*x[10] +
x[3]*x[6]*x[10]^2 + x[3]*x[8]*x[10]^2 +
x[6]*x[8]*x[10]^2 + x[4]*x[7]*x[10]^2 +
x[4]*x[9]*x[10]^2 + x[7]*x[9]*x[10]^2
 pv[7] = x[1]^2*x[5]*x[3] + x[2]^2*x[5]*x[3] +
x[1]^2*x[2]*x[6] + x[2]*x[5]^2*x[6] +
x[2]*x[3]^2*x[6] + x[5]*x[3]*x[6]^2 +
x[1]*x[2]^2*x[8] + x[1]*x[5]^2*x[8] +
x[1]*x[3]^2*x[8] + x[1]*x[6]^2*x[8] +
x[5]*x[3]*x[8]^2 + x[2]*x[6]*x[8]^2 +
x[1]^2*x[5]*x[4] + x[2]^2*x[5]*x[4] +
x[1]^2*x[6]*x[4] + x[3]^2*x[6]*x[4] +
x[2]^2*x[8]*x[4] + x[3]^2*x[8]*x[4] +
x[1]^2*x[2]*x[7] + x[2]*x[5]^2*x[7] +
x[1]^2*x[3]*x[7] + x[3]*x[6]^2*x[7] +
x[5]^2*x[8]*x[7] + x[6]^2*x[8]*x[7] +
x[2]*x[4]^2*x[7] + x[3]*x[4]^2*x[7] +
x[5]*x[4]*x[7]^2 + x[6]*x[4]*x[7]^2 +
x[1]*x[2]^2*x[9] + x[1]*x[5]^2*x[9] +
x[2]^2*x[3]*x[9] + x[5]^2*x[6]*x[9] +
x[3]*x[8]^2*x[9] + x[6]*x[8]^2*x[9] +
x[1]*x[4]^2*x[9] + x[3]*x[4]^2*x[9] +
x[1]*x[7]^2*x[9] + x[6]*x[7]^2*x[9] +
x[5]*x[4]*x[9]^2 + x[8]*x[4]*x[9]^2 +
x[2]*x[7]*x[9]^2 + x[8]*x[7]*x[9]^2 +
x[1]*x[3]^2*x[10] + x[2]*x[3]^2*x[10] +
x[1]*x[6]^2*x[10] + x[5]*x[6]^2*x[10] +
x[2]*x[8]^2*x[10] + x[5]*x[8]^2*x[10] +
x[1]*x[4]^2*x[10] + x[2]*x[4]^2*x[10] +
x[1]*x[7]^2*x[10] + x[5]*x[7]^2*x[10] +
x[2]*x[9]^2*x[10] + x[5]*x[9]^2*x[10] +
x[6]*x[4]*x[10]^2 + x[8]*x[4]*x[10]^2 +
x[3]*x[7]*x[10]^2 + x[8]*x[7]*x[10]^2 +
x[3]*x[9]*x[10]^2 + x[6]*x[9]*x[10]^2
 v[4:8] = pv[3:7]
 pv[8] = x[1]^4*x[2] + x[1]*x[2]^4 + x[1]^4*x[5] +
x[2]^4*x[5] + x[1]*x[5]^4 + x[2]*x[5]^4 +
x[1]^4*x[3] + x[2]^4*x[3] + x[1]*x[3]^4 +
x[2]*x[3]^4 + x[1]^4*x[6] + x[5]^4*x[6] +
x[3]^4*x[6] + x[1]*x[6]^4 + x[5]*x[6]^4 +
x[3]*x[6]^4 + x[2]^4*x[8] + x[5]^4*x[8] +
x[3]^4*x[8] + x[6]^4*x[8] + x[2]*x[8]^4 +
x[5]*x[8]^4 + x[3]*x[8]^4 + x[6]*x[8]^4 +
x[1]^4*x[4] + x[2]^4*x[4] + x[3]^4*x[4] +
x[1]*x[4]^4 + x[2]*x[4]^4 + x[3]*x[4]^4 +
x[1]^4*x[7] + x[5]^4*x[7] + x[6]^4*x[7] +
x[4]^4*x[7] + x[1]*x[7]^4 + x[5]*x[7]^4 +
x[6]*x[7]^4 + x[4]*x[7]^4 + x[2]^4*x[9] +
x[5]^4*x[9] + x[8]^4*x[9] + x[4]^4*x[9] +
x[7]^4*x[9] + x[2]*x[9]^4 + x[5]*x[9]^4 +
x[8]*x[9]^4 + x[4]*x[9]^4 + x[7]*x[9]^4 +
x[3]^4*x[10] + x[6]^4*x[10] + x[8]^4*x[10] +
x[4]^4*x[10] + x[7]^4*x[10] + x[9]^4*x[10] +
x[3]*x[10]^4 + x[6]*x[10]^4 + x[8]*x[10]^4 +
x[4]*x[10]^4 + x[7]*x[10]^4 + x[9]*x[10]^4
 pv[9] = x[1]^3*x[2]^2 + x[1]^2*x[2]^3 + x[1]^3*x[5]^2 +
x[2]^3*x[5]^2 + x[1]^2*x[5]^3 + x[2]^2*x[5]^3 +
x[1]^3*x[3]^2 + x[2]^3*x[3]^2 + x[1]^2*x[3]^3 +
x[2]^2*x[3]^3 + x[1]^3*x[6]^2 + x[5]^3*x[6]^2 +
x[3]^3*x[6]^2 + x[1]^2*x[6]^3 + x[5]^2*x[6]^3 +
x[3]^2*x[6]^3 + x[2]^3*x[8]^2 + x[5]^3*x[8]^2 +
x[3]^3*x[8]^2 + x[6]^3*x[8]^2 + x[2]^2*x[8]^3 +
x[5]^2*x[8]^3 + x[3]^2*x[8]^3 + x[6]^2*x[8]^3 +
x[1]^3*x[4]^2 + x[2]^3*x[4]^2 + x[3]^3*x[4]^2 +
x[1]^2*x[4]^3 + x[2]^2*x[4]^3 + x[3]^2*x[4]^3 +
x[1]^3*x[7]^2 + x[5]^3*x[7]^2 + x[6]^3*x[7]^2 +
x[4]^3*x[7]^2 + x[1]^2*x[7]^3 + x[5]^2*x[7]^3 +
x[6]^2*x[7]^3 + x[4]^2*x[7]^3 + x[2]^3*x[9]^2 +
x[5]^3*x[9]^2 + x[8]^3*x[9]^2 + x[4]^3*x[9]^2 +
x[7]^3*x[9]^2 + x[2]^2*x[9]^3 + x[5]^2*x[9]^3 +
x[8]^2*x[9]^3 + x[4]^2*x[9]^3 + x[7]^2*x[9]^3 +
x[3]^3*x[10]^2 + x[6]^3*x[10]^2 + x[8]^3*x[10]^2 +
x[4]^3*x[10]^2 + x[7]^3*x[10]^2 + x[9]^3*x[10]^2 +
x[3]^2*x[10]^3 + x[6]^2*x[10]^3 + x[8]^2*x[10]^3 +
x[4]^2*x[10]^3 + x[7]^2*x[10]^3 + x[9]^2*x[10]^3
 pv[10] = x[1]^3*x[2]*x[5] + x[1]*x[2]^3*x[5] +
x[1]*x[2]*x[5]^3 + x[1]^3*x[3]*x[6] +
x[1]*x[3]^3*x[6] + x[1]*x[3]*x[6]^3 +
x[2]^3*x[3]*x[8] + x[2]*x[3]^3*x[8] +
x[5]^3*x[6]*x[8] + x[5]*x[6]^3*x[8] +
x[2]*x[3]*x[8]^3 + x[5]*x[6]*x[8]^3 +
x[1]^3*x[4]*x[7] + x[1]*x[4]^3*x[7] +
x[1]*x[4]*x[7]^3 + x[2]^3*x[4]*x[9] +
x[2]*x[4]^3*x[9] + x[5]^3*x[7]*x[9] +
x[5]*x[7]^3*x[9] + x[2]*x[4]*x[9]^3 +
x[5]*x[7]*x[9]^3 + x[3]^3*x[4]*x[10] +
x[3]*x[4]^3*x[10] + x[6]^3*x[7]*x[10] +
x[6]*x[7]^3*x[10] + x[8]^3*x[9]*x[10] +
x[8]*x[9]^3*x[10] + x[3]*x[4]*x[10]^3 +
x[6]*x[7]*x[10]^3 + x[8]*x[9]*x[10]^3
 pv[11] = x[1]^2*x[2]^2*x[5] + x[1]^2*x[2]*x[5]^2 +
x[1]*x[2]^2*x[5]^2 + x[1]^2*x[3]^2*x[6] +
x[1]^2*x[3]*x[6]^2 + x[1]*x[3]^2*x[6]^2 +
x[2]^2*x[3]^2*x[8] + x[5]^2*x[6]^2*x[8] +
x[2]^2*x[3]*x[8]^2 + x[2]*x[3]^2*x[8]^2 +
x[5]^2*x[6]*x[8]^2 + x[5]*x[6]^2*x[8]^2 +
x[1]^2*x[4]^2*x[7] + x[1]^2*x[4]*x[7]^2 +
x[1]*x[4]^2*x[7]^2 + x[2]^2*x[4]^2*x[9] +
x[5]^2*x[7]^2*x[9] + x[2]^2*x[4]*x[9]^2 +
x[2]*x[4]^2*x[9]^2 + x[5]^2*x[7]*x[9]^2 +
x[5]*x[7]^2*x[9]^2 + x[3]^2*x[4]^2*x[10] +
x[6]^2*x[7]^2*x[10] + x[8]^2*x[9]^2*x[10] +
x[3]^2*x[4]*x[10]^2 + x[3]*x[4]^2*x[10]^2 +
x[6]^2*x[7]*x[10]^2 + x[6]*x[7]^2*x[10]^2 +
x[8]^2*x[9]*x[10]^2 + x[8]*x[9]^2*x[10]^2
 pv[12] = x[1]^3*x[2]*x[3] + x[1]*x[2]^3*x[3] +
x[1]*x[2]*x[3]^3 + x[1]^3*x[5]*x[6] +
x[1]*x[5]^3*x[6] + x[1]*x[5]*x[6]^3 +
x[2]^3*x[5]*x[8] + x[2]*x[5]^3*x[8] +
x[3]^3*x[6]*x[8] + x[3]*x[6]^3*x[8] +
x[2]*x[5]*x[8]^3 + x[3]*x[6]*x[8]^3 +
x[1]^3*x[2]*x[4] + x[1]*x[2]^3*x[4] +
x[1]^3*x[3]*x[4] + x[2]^3*x[3]*x[4] +
x[1]*x[3]^3*x[4] + x[2]*x[3]^3*x[4] +
x[1]*x[2]*x[4]^3 + x[1]*x[3]*x[4]^3 +
x[2]*x[3]*x[4]^3 + x[1]^3*x[5]*x[7] +
x[1]*x[5]^3*x[7] + x[1]^3*x[6]*x[7] +
x[5]^3*x[6]*x[7] + x[1]*x[6]^3*x[7] +
x[5]*x[6]^3*x[7] + x[1]*x[5]*x[7]^3 +
x[1]*x[6]*x[7]^3 + x[5]*x[6]*x[7]^3 +
x[2]^3*x[5]*x[9] + x[2]*x[5]^3*x[9] +
x[2]^3*x[8]*x[9] + x[5]^3*x[8]*x[9] +
x[2]*x[8]^3*x[9] + x[5]*x[8]^3*x[9] +
x[4]^3*x[7]*x[9] + x[4]*x[7]^3*x[9] +
x[2]*x[5]*x[9]^3 + x[2]*x[8]*x[9]^3 +
x[5]*x[8]*x[9]^3 + x[4]*x[7]*x[9]^3 +
x[3]^3*x[6]*x[10] + x[3]*x[6]^3*x[10] +
x[3]^3*x[8]*x[10] + x[6]^3*x[8]*x[10] +
x[3]*x[8]^3*x[10] + x[6]*x[8]^3*x[10] +
x[4]^3*x[7]*x[10] + x[4]*x[7]^3*x[10] +
x[4]^3*x[9]*x[10] + x[7]^3*x[9]*x[10] +
x[4]*x[9]^3*x[10] + x[7]*x[9]^3*x[10] +
x[3]*x[6]*x[10]^3 + x[3]*x[8]*x[10]^3 +
x[6]*x[8]*x[10]^3 + x[4]*x[7]*x[10]^3 +
x[4]*x[9]*x[10]^3 + x[7]*x[9]*x[10]^3
 pv[13] = x[1]^2*x[2]^2*x[3] + x[1]^2*x[2]*x[3]^2 +
x[1]*x[2]^2*x[3]^2 + x[1]^2*x[5]^2*x[6] +
x[1]^2*x[5]*x[6]^2 + x[1]*x[5]^2*x[6]^2 +
x[2]^2*x[5]^2*x[8] + x[3]^2*x[6]^2*x[8] +
x[2]^2*x[5]*x[8]^2 + x[2]*x[5]^2*x[8]^2 +
x[3]^2*x[6]*x[8]^2 + x[3]*x[6]^2*x[8]^2 +
x[1]^2*x[2]^2*x[4] + x[1]^2*x[3]^2*x[4] +
x[2]^2*x[3]^2*x[4] + x[1]^2*x[2]*x[4]^2 +
x[1]*x[2]^2*x[4]^2 + x[1]^2*x[3]*x[4]^2 +
x[2]^2*x[3]*x[4]^2 + x[1]*x[3]^2*x[4]^2 +
x[2]*x[3]^2*x[4]^2 + x[1]^2*x[5]^2*x[7] +
x[1]^2*x[6]^2*x[7] + x[5]^2*x[6]^2*x[7] +
x[1]^2*x[5]*x[7]^2 + x[1]*x[5]^2*x[7]^2 +
x[1]^2*x[6]*x[7]^2 + x[5]^2*x[6]*x[7]^2 +
x[1]*x[6]^2*x[7]^2 + x[5]*x[6]^2*x[7]^2 +
x[2]^2*x[5]^2*x[9] + x[2]^2*x[8]^2*x[9] +
x[5]^2*x[8]^2*x[9] + x[4]^2*x[7]^2*x[9] +
x[2]^2*x[5]*x[9]^2 + x[2]*x[5]^2*x[9]^2 +
x[2]^2*x[8]*x[9]^2 + x[5]^2*x[8]*x[9]^2 +
x[2]*x[8]^2*x[9]^2 + x[5]*x[8]^2*x[9]^2 +
x[4]^2*x[7]*x[9]^2 + x[4]*x[7]^2*x[9]^2 +
x[3]^2*x[6]^2*x[10] + x[3]^2*x[8]^2*x[10] +
x[6]^2*x[8]^2*x[10] + x[4]^2*x[7]^2*x[10] +
x[4]^2*x[9]^2*x[10] + x[7]^2*x[9]^2*x[10] +
x[3]^2*x[6]*x[10]^2 + x[3]*x[6]^2*x[10]^2 +
x[3]^2*x[8]*x[10]^2 + x[6]^2*x[8]*x[10]^2 +
x[3]*x[8]^2*x[10]^2 + x[6]*x[8]^2*x[10]^2 +
x[4]^2*x[7]*x[10]^2 + x[4]*x[7]^2*x[10]^2 +
x[4]^2*x[9]*x[10]^2 + x[7]^2*x[9]*x[10]^2 +
x[4]*x[9]^2*x[10]^2 + x[7]*x[9]^2*x[10]^2
 pv[14] = x[1]^3*x[5]*x[3] + x[2]^3*x[5]*x[3] +
x[1]^3*x[2]*x[6] + x[2]*x[5]^3*x[6] +
x[2]*x[3]^3*x[6] + x[5]*x[3]*x[6]^3 +
x[1]*x[2]^3*x[8] + x[1]*x[5]^3*x[8] +
x[1]*x[3]^3*x[8] + x[1]*x[6]^3*x[8] +
x[5]*x[3]*x[8]^3 + x[2]*x[6]*x[8]^3 +
x[1]^3*x[5]*x[4] + x[2]^3*x[5]*x[4] +
x[1]^3*x[6]*x[4] + x[3]^3*x[6]*x[4] +
x[2]^3*x[8]*x[4] + x[3]^3*x[8]*x[4] +
x[1]^3*x[2]*x[7] + x[2]*x[5]^3*x[7] +
x[1]^3*x[3]*x[7] + x[3]*x[6]^3*x[7] +
x[5]^3*x[8]*x[7] + x[6]^3*x[8]*x[7] +
x[2]*x[4]^3*x[7] + x[3]*x[4]^3*x[7] +
x[5]*x[4]*x[7]^3 + x[6]*x[4]*x[7]^3 +
x[1]*x[2]^3*x[9] + x[1]*x[5]^3*x[9] +
x[2]^3*x[3]*x[9] + x[5]^3*x[6]*x[9] +
x[3]*x[8]^3*x[9] + x[6]*x[8]^3*x[9] +
x[1]*x[4]^3*x[9] + x[3]*x[4]^3*x[9] +
x[1]*x[7]^3*x[9] + x[6]*x[7]^3*x[9] +
x[5]*x[4]*x[9]^3 + x[8]*x[4]*x[9]^3 +
x[2]*x[7]*x[9]^3 + x[8]*x[7]*x[9]^3 +
x[1]*x[3]^3*x[10] + x[2]*x[3]^3*x[10] +
x[1]*x[6]^3*x[10] + x[5]*x[6]^3*x[10] +
x[2]*x[8]^3*x[10] + x[5]*x[8]^3*x[10] +
x[1]*x[4]^3*x[10] + x[2]*x[4]^3*x[10] +
x[1]*x[7]^3*x[10] + x[5]*x[7]^3*x[10] +
x[2]*x[9]^3*x[10] + x[5]*x[9]^3*x[10] +
x[6]*x[4]*x[10]^3 + x[8]*x[4]*x[10]^3 +
x[3]*x[7]*x[10]^3 + x[8]*x[7]*x[10]^3 +
x[3]*x[9]*x[10]^3 + x[6]*x[9]*x[10]^3
 pv[15] = x[1]^2*x[5]^2*x[3] + x[2]^2*x[5]^2*x[3] +
x[1]^2*x[5]*x[3]^2 + x[2]^2*x[5]*x[3]^2 +
x[1]^2*x[2]^2*x[6] + x[2]^2*x[5]^2*x[6] +
x[2]^2*x[3]^2*x[6] + x[1]^2*x[2]*x[6]^2 +
x[2]*x[5]^2*x[6]^2 + x[5]^2*x[3]*x[6]^2 +
x[2]*x[3]^2*x[6]^2 + x[5]*x[3]^2*x[6]^2 +
x[1]^2*x[2]^2*x[8] + x[1]^2*x[5]^2*x[8] +
x[1]^2*x[3]^2*x[8] + x[1]^2*x[6]^2*x[8] +
x[1]*x[2]^2*x[8]^2 + x[1]*x[5]^2*x[8]^2 +
x[5]^2*x[3]*x[8]^2 + x[1]*x[3]^2*x[8]^2 +
x[5]*x[3]^2*x[8]^2 + x[2]^2*x[6]*x[8]^2 +
x[1]*x[6]^2*x[8]^2 + x[2]*x[6]^2*x[8]^2 +
x[1]^2*x[5]^2*x[4] + x[2]^2*x[5]^2*x[4] +
x[1]^2*x[6]^2*x[4] + x[3]^2*x[6]^2*x[4] +
x[2]^2*x[8]^2*x[4] + x[3]^2*x[8]^2*x[4] +
x[1]^2*x[5]*x[4]^2 + x[2]^2*x[5]*x[4]^2 +
x[1]^2*x[6]*x[4]^2 + x[3]^2*x[6]*x[4]^2 +
x[2]^2*x[8]*x[4]^2 + x[3]^2*x[8]*x[4]^2 +
x[1]^2*x[2]^2*x[7] + x[2]^2*x[5]^2*x[7] +
x[1]^2*x[3]^2*x[7] + x[3]^2*x[6]^2*x[7] +
x[5]^2*x[8]^2*x[7] + x[6]^2*x[8]^2*x[7] +
x[2]^2*x[4]^2*x[7] + x[3]^2*x[4]^2*x[7] +
x[1]^2*x[2]*x[7]^2 + x[2]*x[5]^2*x[7]^2 +
x[1]^2*x[3]*x[7]^2 + x[3]*x[6]^2*x[7]^2 +
x[5]^2*x[8]*x[7]^2 + x[6]^2*x[8]*x[7]^2 +
x[5]^2*x[4]*x[7]^2 + x[6]^2*x[4]*x[7]^2 +
x[2]*x[4]^2*x[7]^2 + x[5]*x[4]^2*x[7]^2 +
x[3]*x[4]^2*x[7]^2 + x[6]*x[4]^2*x[7]^2 +
x[1]^2*x[2]^2*x[9] + x[1]^2*x[5]^2*x[9] +
x[2]^2*x[3]^2*x[9] + x[5]^2*x[6]^2*x[9] +
x[3]^2*x[8]^2*x[9] + x[6]^2*x[8]^2*x[9] +
x[1]^2*x[4]^2*x[9] + x[3]^2*x[4]^2*x[9] +
x[1]^2*x[7]^2*x[9] + x[6]^2*x[7]^2*x[9] +
x[1]*x[2]^2*x[9]^2 + x[1]*x[5]^2*x[9]^2 +
x[2]^2*x[3]*x[9]^2 + x[5]^2*x[6]*x[9]^2 +
x[3]*x[8]^2*x[9]^2 + x[6]*x[8]^2*x[9]^2 +
x[5]^2*x[4]*x[9]^2 + x[8]^2*x[4]*x[9]^2 +
x[1]*x[4]^2*x[9]^2 + x[5]*x[4]^2*x[9]^2 +
x[3]*x[4]^2*x[9]^2 + x[8]*x[4]^2*x[9]^2 +
x[2]^2*x[7]*x[9]^2 + x[8]^2*x[7]*x[9]^2 +
x[1]*x[7]^2*x[9]^2 + x[2]*x[7]^2*x[9]^2 +
x[6]*x[7]^2*x[9]^2 + x[8]*x[7]^2*x[9]^2 +
x[1]^2*x[3]^2*x[10] + x[2]^2*x[3]^2*x[10] +
x[1]^2*x[6]^2*x[10] + x[5]^2*x[6]^2*x[10] +
x[2]^2*x[8]^2*x[10] + x[5]^2*x[8]^2*x[10] +
x[1]^2*x[4]^2*x[10] + x[2]^2*x[4]^2*x[10] +
x[1]^2*x[7]^2*x[10] + x[5]^2*x[7]^2*x[10] +
x[2]^2*x[9]^2*x[10] + x[5]^2*x[9]^2*x[10] +
x[1]*x[3]^2*x[10]^2 + x[2]*x[3]^2*x[10]^2 +
x[1]*x[6]^2*x[10]^2 + x[5]*x[6]^2*x[10]^2 +
x[2]*x[8]^2*x[10]^2 + x[5]*x[8]^2*x[10]^2 +
x[6]^2*x[4]*x[10]^2 + x[8]^2*x[4]*x[10]^2 +
x[1]*x[4]^2*x[10]^2 + x[2]*x[4]^2*x[10]^2 +
x[6]*x[4]^2*x[10]^2 + x[8]*x[4]^2*x[10]^2 +
x[3]^2*x[7]*x[10]^2 + x[8]^2*x[7]*x[10]^2 +
x[1]*x[7]^2*x[10]^2 + x[5]*x[7]^2*x[10]^2 +
x[3]*x[7]^2*x[10]^2 + x[8]*x[7]^2*x[10]^2 +
x[3]^2*x[9]*x[10]^2 + x[6]^2*x[9]*x[10]^2 +
x[2]*x[9]^2*x[10]^2 + x[5]*x[9]^2*x[10]^2 +
x[3]*x[9]^2*x[10]^2 + x[6]*x[9]^2*x[10]^2
 v[9:16] = pv[8:15]
 v[17] = pv[1]*pv[1]
 v[18] = pv[1]*pv[2]
 v[19] = pv[2]*pv[2]
 pv[16] = x[1]^5*x[2] + x[1]*x[2]^5 + x[1]^5*x[5] +
x[2]^5*x[5] + x[1]*x[5]^5 + x[2]*x[5]^5 +
x[1]^5*x[3] + x[2]^5*x[3] + x[1]*x[3]^5 +
x[2]*x[3]^5 + x[1]^5*x[6] + x[5]^5*x[6] +
x[3]^5*x[6] + x[1]*x[6]^5 + x[5]*x[6]^5 +
x[3]*x[6]^5 + x[2]^5*x[8] + x[5]^5*x[8] +
x[3]^5*x[8] + x[6]^5*x[8] + x[2]*x[8]^5 +
x[5]*x[8]^5 + x[3]*x[8]^5 + x[6]*x[8]^5 +
x[1]^5*x[4] + x[2]^5*x[4] + x[3]^5*x[4] +
x[1]*x[4]^5 + x[2]*x[4]^5 + x[3]*x[4]^5 +
x[1]^5*x[7] + x[5]^5*x[7] + x[6]^5*x[7] +
x[4]^5*x[7] + x[1]*x[7]^5 + x[5]*x[7]^5 +
x[6]*x[7]^5 + x[4]*x[7]^5 + x[2]^5*x[9] +
x[5]^5*x[9] + x[8]^5*x[9] + x[4]^5*x[9] +
x[7]^5*x[9] + x[2]*x[9]^5 + x[5]*x[9]^5 +
x[8]*x[9]^5 + x[4]*x[9]^5 + x[7]*x[9]^5 +
x[3]^5*x[10] + x[6]^5*x[10] + x[8]^5*x[10] +
x[4]^5*x[10] + x[7]^5*x[10] + x[9]^5*x[10] +
x[3]*x[10]^5 + x[6]*x[10]^5 + x[8]*x[10]^5 +
x[4]*x[10]^5 + x[7]*x[10]^5 + x[9]*x[10]^5
 pv[17] = x[1]^4*x[2]^2 + x[1]^2*x[2]^4 + x[1]^4*x[5]^2 +
x[2]^4*x[5]^2 + x[1]^2*x[5]^4 + x[2]^2*x[5]^4 +
x[1]^4*x[3]^2 + x[2]^4*x[3]^2 + x[1]^2*x[3]^4 +
x[2]^2*x[3]^4 + x[1]^4*x[6]^2 + x[5]^4*x[6]^2 +
x[3]^4*x[6]^2 + x[1]^2*x[6]^4 + x[5]^2*x[6]^4 +
x[3]^2*x[6]^4 + x[2]^4*x[8]^2 + x[5]^4*x[8]^2 +
x[3]^4*x[8]^2 + x[6]^4*x[8]^2 + x[2]^2*x[8]^4 +
x[5]^2*x[8]^4 + x[3]^2*x[8]^4 + x[6]^2*x[8]^4 +
x[1]^4*x[4]^2 + x[2]^4*x[4]^2 + x[3]^4*x[4]^2 +
x[1]^2*x[4]^4 + x[2]^2*x[4]^4 + x[3]^2*x[4]^4 +
x[1]^4*x[7]^2 + x[5]^4*x[7]^2 + x[6]^4*x[7]^2 +
x[4]^4*x[7]^2 + x[1]^2*x[7]^4 + x[5]^2*x[7]^4 +
x[6]^2*x[7]^4 + x[4]^2*x[7]^4 + x[2]^4*x[9]^2 +
x[5]^4*x[9]^2 + x[8]^4*x[9]^2 + x[4]^4*x[9]^2 +
x[7]^4*x[9]^2 + x[2]^2*x[9]^4 + x[5]^2*x[9]^4 +
x[8]^2*x[9]^4 + x[4]^2*x[9]^4 + x[7]^2*x[9]^4 +
x[3]^4*x[10]^2 + x[6]^4*x[10]^2 + x[8]^4*x[10]^2 +
x[4]^4*x[10]^2 + x[7]^4*x[10]^2 + x[9]^4*x[10]^2 +
x[3]^2*x[10]^4 + x[6]^2*x[10]^4 + x[8]^2*x[10]^4 +
x[4]^2*x[10]^4 + x[7]^2*x[10]^4 + x[9]^2*x[10]^4
 pv[18] = x[1]^3*x[2]^3 + x[1]^3*x[5]^3 + x[2]^3*x[5]^3+
x[1]^3*x[3]^3 + x[2]^3*x[3]^3 + x[1]^3*x[6]^3 +
x[5]^3*x[6]^3 + x[3]^3*x[6]^3 + x[2]^3*x[8]^3 +
x[5]^3*x[8]^3 + x[3]^3*x[8]^3 + x[6]^3*x[8]^3 +
x[1]^3*x[4]^3 + x[2]^3*x[4]^3 + x[3]^3*x[4]^3 +
x[1]^3*x[7]^3 + x[5]^3*x[7]^3 + x[6]^3*x[7]^3 +
x[4]^3*x[7]^3 + x[2]^3*x[9]^3 + x[5]^3*x[9]^3 +
x[8]^3*x[9]^3 + x[4]^3*x[9]^3 + x[7]^3*x[9]^3 +
x[3]^3*x[10]^3 + x[6]^3*x[10]^3 + x[8]^3*x[10]^3 +
x[4]^3*x[10]^3 + x[7]^3*x[10]^3 + x[9]^3*x[10]^3
 pv[19] = x[1]^4*x[2]*x[5] + x[1]*x[2]^4*x[5] +
x[1]*x[2]*x[5]^4 + x[1]^4*x[3]*x[6] +
x[1]*x[3]^4*x[6] + x[1]*x[3]*x[6]^4 +
x[2]^4*x[3]*x[8] + x[2]*x[3]^4*x[8] +
x[5]^4*x[6]*x[8] + x[5]*x[6]^4*x[8] +
x[2]*x[3]*x[8]^4 + x[5]*x[6]*x[8]^4 +
x[1]^4*x[4]*x[7] + x[1]*x[4]^4*x[7] +
x[1]*x[4]*x[7]^4 + x[2]^4*x[4]*x[9] +
x[2]*x[4]^4*x[9] + x[5]^4*x[7]*x[9] +
x[5]*x[7]^4*x[9] + x[2]*x[4]*x[9]^4 +
x[5]*x[7]*x[9]^4 + x[3]^4*x[4]*x[10] +
x[3]*x[4]^4*x[10] + x[6]^4*x[7]*x[10] +
x[6]*x[7]^4*x[10] + x[8]^4*x[9]*x[10] +
x[8]*x[9]^4*x[10] + x[3]*x[4]*x[10]^4 +
x[6]*x[7]*x[10]^4 + x[8]*x[9]*x[10]^4
 pv[20] = x[1]^3*x[2]^2*x[5] + x[1]^2*x[2]^3*x[5] +
x[1]^3*x[2]*x[5]^2 + x[1]*x[2]^3*x[5]^2 +
x[1]^2*x[2]*x[5]^3 + x[1]*x[2]^2*x[5]^3 +
x[1]^3*x[3]^2*x[6] + x[1]^2*x[3]^3*x[6] +
x[1]^3*x[3]*x[6]^2 + x[1]*x[3]^3*x[6]^2 +
x[1]^2*x[3]*x[6]^3 + x[1]*x[3]^2*x[6]^3 +
x[2]^3*x[3]^2*x[8] + x[2]^2*x[3]^3*x[8] +
x[5]^3*x[6]^2*x[8] + x[5]^2*x[6]^3*x[8] +
x[2]^3*x[3]*x[8]^2 + x[2]*x[3]^3*x[8]^2 +
x[5]^3*x[6]*x[8]^2 + x[5]*x[6]^3*x[8]^2 +
x[2]^2*x[3]*x[8]^3 + x[2]*x[3]^2*x[8]^3 +
x[5]^2*x[6]*x[8]^3 + x[5]*x[6]^2*x[8]^3 +
x[1]^3*x[4]^2*x[7] + x[1]^2*x[4]^3*x[7] +
x[1]^3*x[4]*x[7]^2 + x[1]*x[4]^3*x[7]^2 +
x[1]^2*x[4]*x[7]^3 + x[1]*x[4]^2*x[7]^3 +
x[2]^3*x[4]^2*x[9] + x[2]^2*x[4]^3*x[9] +
x[5]^3*x[7]^2*x[9] + x[5]^2*x[7]^3*x[9] +
x[2]^3*x[4]*x[9]^2 + x[2]*x[4]^3*x[9]^2 +
x[5]^3*x[7]*x[9]^2 + x[5]*x[7]^3*x[9]^2 +
x[2]^2*x[4]*x[9]^3 + x[2]*x[4]^2*x[9]^3 +
x[5]^2*x[7]*x[9]^3 + x[5]*x[7]^2*x[9]^3 +
x[3]^3*x[4]^2*x[10] + x[3]^2*x[4]^3*x[10] +
x[6]^3*x[7]^2*x[10] + x[6]^2*x[7]^3*x[10] +
x[8]^3*x[9]^2*x[10] + x[8]^2*x[9]^3*x[10] +
x[3]^3*x[4]*x[10]^2 + x[3]*x[4]^3*x[10]^2 +
x[6]^3*x[7]*x[10]^2 + x[6]*x[7]^3*x[10]^2 +
x[8]^3*x[9]*x[10]^2 + x[8]*x[9]^3*x[10]^2 +
x[3]^2*x[4]*x[10]^3 + x[3]*x[4]^2*x[10]^3 +
x[6]^2*x[7]*x[10]^3 + x[6]*x[7]^2*x[10]^3 +
x[8]^2*x[9]*x[10]^3 + x[8]*x[9]^2*x[10]^3
 pv[21] = x[1]^2*x[2]^2*x[5]^2 + x[1]^2*x[3]^2*x[6]^2 +
x[2]^2*x[3]^2*x[8]^2 + x[5]^2*x[6]^2*x[8]^2 +
x[1]^2*x[4]^2*x[7]^2 + x[2]^2*x[4]^2*x[9]^2 +
x[5]^2*x[7]^2*x[9]^2 + x[3]^2*x[4]^2*x[10]^2 +
x[6]^2*x[7]^2*x[10]^2 + x[8]^2*x[9]^2*x[10]^2
 pv[22] = x[1]^4*x[2]*x[3] + x[1]*x[2]^4*x[3] +
x[1]*x[2]*x[3]^4 + x[1]^4*x[5]*x[6] +
x[1]*x[5]^4*x[6] + x[1]*x[5]*x[6]^4 +
x[2]^4*x[5]*x[8] + x[2]*x[5]^4*x[8] +
x[3]^4*x[6]*x[8] + x[3]*x[6]^4*x[8] +
x[2]*x[5]*x[8]^4 + x[3]*x[6]*x[8]^4 +
x[1]^4*x[2]*x[4] + x[1]*x[2]^4*x[4] +
x[1]^4*x[3]*x[4] + x[2]^4*x[3]*x[4] +
x[1]*x[3]^4*x[4] + x[2]*x[3]^4*x[4] +
x[1]*x[2]*x[4]^4 + x[1]*x[3]*x[4]^4 +
x[2]*x[3]*x[4]^4 + x[1]^4*x[5]*x[7] +
x[1]*x[5]^4*x[7] + x[1]^4*x[6]*x[7] +
x[5]^4*x[6]*x[7] + x[1]*x[6]^4*x[7] +
x[5]*x[6]^4*x[7] + x[1]*x[5]*x[7]^4 +
x[1]*x[6]*x[7]^4 + x[5]*x[6]*x[7]^4 +
x[2]^4*x[5]*x[9] + x[2]*x[5]^4*x[9] +
x[2]^4*x[8]*x[9] + x[5]^4*x[8]*x[9] +
x[2]*x[8]^4*x[9] + x[5]*x[8]^4*x[9] +
x[4]^4*x[7]*x[9] + x[4]*x[7]^4*x[9] +
x[2]*x[5]*x[9]^4 + x[2]*x[8]*x[9]^4 +
x[5]*x[8]*x[9]^4 + x[4]*x[7]*x[9]^4 +
x[3]^4*x[6]*x[10] + x[3]*x[6]^4*x[10] +
x[3]^4*x[8]*x[10] + x[6]^4*x[8]*x[10] +
x[3]*x[8]^4*x[10] + x[6]*x[8]^4*x[10] +
x[4]^4*x[7]*x[10] + x[4]*x[7]^4*x[10] +
x[4]^4*x[9]*x[10] + x[7]^4*x[9]*x[10] +
x[4]*x[9]^4*x[10] + x[7]*x[9]^4*x[10] +
x[3]*x[6]*x[10]^4 + x[3]*x[8]*x[10]^4 +
x[6]*x[8]*x[10]^4 + x[4]*x[7]*x[10]^4 +
x[4]*x[9]*x[10]^4 + x[7]*x[9]*x[10]^4
 pv[23] = x[1]^3*x[2]^2*x[3] + x[1]^2*x[2]^3*x[3] +
x[1]^3*x[2]*x[3]^2 + x[1]*x[2]^3*x[3]^2 +
x[1]^2*x[2]*x[3]^3 + x[1]*x[2]^2*x[3]^3 +
x[1]^3*x[5]^2*x[6] + x[1]^2*x[5]^3*x[6] +
x[1]^3*x[5]*x[6]^2 + x[1]*x[5]^3*x[6]^2 +
x[1]^2*x[5]*x[6]^3 + x[1]*x[5]^2*x[6]^3 +
x[2]^3*x[5]^2*x[8] + x[2]^2*x[5]^3*x[8] +
x[3]^3*x[6]^2*x[8] + x[3]^2*x[6]^3*x[8] +
x[2]^3*x[5]*x[8]^2 + x[2]*x[5]^3*x[8]^2 +
x[3]^3*x[6]*x[8]^2 + x[3]*x[6]^3*x[8]^2 +
x[2]^2*x[5]*x[8]^3 + x[2]*x[5]^2*x[8]^3 +
x[3]^2*x[6]*x[8]^3 + x[3]*x[6]^2*x[8]^3 +
x[1]^3*x[2]^2*x[4] + x[1]^2*x[2]^3*x[4] +
x[1]^3*x[3]^2*x[4] + x[2]^3*x[3]^2*x[4] +
x[1]^2*x[3]^3*x[4] + x[2]^2*x[3]^3*x[4] +
x[1]^3*x[2]*x[4]^2 + x[1]*x[2]^3*x[4]^2 +
x[1]^3*x[3]*x[4]^2 + x[2]^3*x[3]*x[4]^2 +
x[1]*x[3]^3*x[4]^2 + x[2]*x[3]^3*x[4]^2 +
x[1]^2*x[2]*x[4]^3 + x[1]*x[2]^2*x[4]^3 +
x[1]^2*x[3]*x[4]^3 + x[2]^2*x[3]*x[4]^3 +
x[1]*x[3]^2*x[4]^3 + x[2]*x[3]^2*x[4]^3 +
x[1]^3*x[5]^2*x[7] + x[1]^2*x[5]^3*x[7] +
x[1]^3*x[6]^2*x[7] + x[5]^3*x[6]^2*x[7] +
x[1]^2*x[6]^3*x[7] + x[5]^2*x[6]^3*x[7] +
x[1]^3*x[5]*x[7]^2 + x[1]*x[5]^3*x[7]^2 +
x[1]^3*x[6]*x[7]^2 + x[5]^3*x[6]*x[7]^2 +
x[1]*x[6]^3*x[7]^2 + x[5]*x[6]^3*x[7]^2 +
x[1]^2*x[5]*x[7]^3 + x[1]*x[5]^2*x[7]^3 +
x[1]^2*x[6]*x[7]^3 + x[5]^2*x[6]*x[7]^3 +
x[1]*x[6]^2*x[7]^3 + x[5]*x[6]^2*x[7]^3 +
x[2]^3*x[5]^2*x[9] + x[2]^2*x[5]^3*x[9] +
x[2]^3*x[8]^2*x[9] + x[5]^3*x[8]^2*x[9] +
x[2]^2*x[8]^3*x[9] + x[5]^2*x[8]^3*x[9] +
x[4]^3*x[7]^2*x[9] + x[4]^2*x[7]^3*x[9] +
x[2]^3*x[5]*x[9]^2 + x[2]*x[5]^3*x[9]^2 +
x[2]^3*x[8]*x[9]^2 + x[5]^3*x[8]*x[9]^2 +
x[2]*x[8]^3*x[9]^2 + x[5]*x[8]^3*x[9]^2 +
x[4]^3*x[7]*x[9]^2 + x[4]*x[7]^3*x[9]^2 +
x[2]^2*x[5]*x[9]^3 + x[2]*x[5]^2*x[9]^3 +
x[2]^2*x[8]*x[9]^3 + x[5]^2*x[8]*x[9]^3 +
x[2]*x[8]^2*x[9]^3 + x[5]*x[8]^2*x[9]^3 +
x[4]^2*x[7]*x[9]^3 + x[4]*x[7]^2*x[9]^3 +
x[3]^3*x[6]^2*x[10] + x[3]^2*x[6]^3*x[10] +
x[3]^3*x[8]^2*x[10] + x[6]^3*x[8]^2*x[10] +
x[3]^2*x[8]^3*x[10] + x[6]^2*x[8]^3*x[10] +
x[4]^3*x[7]^2*x[10] + x[4]^2*x[7]^3*x[10] +
x[4]^3*x[9]^2*x[10] + x[7]^3*x[9]^2*x[10] +
x[4]^2*x[9]^3*x[10] + x[7]^2*x[9]^3*x[10] +
x[3]^3*x[6]*x[10]^2 + x[3]*x[6]^3*x[10]^2 +
x[3]^3*x[8]*x[10]^2 + x[6]^3*x[8]*x[10]^2 +
x[3]*x[8]^3*x[10]^2 + x[6]*x[8]^3*x[10]^2 +
x[4]^3*x[7]*x[10]^2 + x[4]*x[7]^3*x[10]^2 +
x[4]^3*x[9]*x[10]^2 + x[7]^3*x[9]*x[10]^2 +
x[4]*x[9]^3*x[10]^2 + x[7]*x[9]^3*x[10]^2 +
x[3]^2*x[6]*x[10]^3 + x[3]*x[6]^2*x[10]^3 +
x[3]^2*x[8]*x[10]^3 + x[6]^2*x[8]*x[10]^3 +
x[3]*x[8]^2*x[10]^3 + x[6]*x[8]^2*x[10]^3 +
x[4]^2*x[7]*x[10]^3 + x[4]*x[7]^2*x[10]^3 +
x[4]^2*x[9]*x[10]^3 + x[7]^2*x[9]*x[10]^3 +
x[4]*x[9]^2*x[10]^3 + x[7]*x[9]^2*x[10]^3
 pv[24] = x[1]^4*x[5]*x[3] + x[2]^4*x[5]*x[3] +
x[1]^4*x[2]*x[6] + x[2]*x[5]^4*x[6] +
x[2]*x[3]^4*x[6] + x[5]*x[3]*x[6]^4 +
x[1]*x[2]^4*x[8] + x[1]*x[5]^4*x[8] +
x[1]*x[3]^4*x[8] + x[1]*x[6]^4*x[8] +
x[5]*x[3]*x[8]^4 + x[2]*x[6]*x[8]^4 +
x[1]^4*x[5]*x[4] + x[2]^4*x[5]*x[4] +
x[1]^4*x[6]*x[4] + x[3]^4*x[6]*x[4] +
x[2]^4*x[8]*x[4] + x[3]^4*x[8]*x[4] +
x[1]^4*x[2]*x[7] + x[2]*x[5]^4*x[7] +
x[1]^4*x[3]*x[7] + x[3]*x[6]^4*x[7] +
x[5]^4*x[8]*x[7] + x[6]^4*x[8]*x[7] +
x[2]*x[4]^4*x[7] + x[3]*x[4]^4*x[7] +
x[5]*x[4]*x[7]^4 + x[6]*x[4]*x[7]^4 +
x[1]*x[2]^4*x[9] + x[1]*x[5]^4*x[9] +
x[2]^4*x[3]*x[9] + x[5]^4*x[6]*x[9] +
x[3]*x[8]^4*x[9] + x[6]*x[8]^4*x[9] +
x[1]*x[4]^4*x[9] + x[3]*x[4]^4*x[9] +
x[1]*x[7]^4*x[9] + x[6]*x[7]^4*x[9] +
x[5]*x[4]*x[9]^4 + x[8]*x[4]*x[9]^4 +
x[2]*x[7]*x[9]^4 + x[8]*x[7]*x[9]^4 +
x[1]*x[3]^4*x[10] + x[2]*x[3]^4*x[10] +
x[1]*x[6]^4*x[10] + x[5]*x[6]^4*x[10] +
x[2]*x[8]^4*x[10] + x[5]*x[8]^4*x[10] +
x[1]*x[4]^4*x[10] + x[2]*x[4]^4*x[10] +
x[1]*x[7]^4*x[10] + x[5]*x[7]^4*x[10] +
x[2]*x[9]^4*x[10] + x[5]*x[9]^4*x[10] +
x[6]*x[4]*x[10]^4 + x[8]*x[4]*x[10]^4 +
x[3]*x[7]*x[10]^4 + x[8]*x[7]*x[10]^4 +
x[3]*x[9]*x[10]^4 + x[6]*x[9]*x[10]^4
 pv[25] = x[1]^3*x[2]*x[5]*x[3] +
x[1]*x[2]^3*x[5]*x[3] + x[1]^3*x[2]*x[5]*x[6] +
x[1]*x[2]*x[5]^3*x[6] + x[1]^3*x[2]*x[3]*x[6] +
x[1]^3*x[5]*x[3]*x[6] + x[1]*x[2]*x[3]^3*x[6] +
x[1]*x[5]*x[3]*x[6]^3 + x[1]*x[2]^3*x[5]*x[8] +
x[1]*x[2]*x[5]^3*x[8] + x[1]*x[2]^3*x[3]*x[8] +
x[2]^3*x[5]*x[3]*x[8] + x[1]*x[2]*x[3]^3*x[8] +
x[1]*x[5]^3*x[6]*x[8] + x[2]*x[5]^3*x[6]*x[8] +
x[1]*x[3]^3*x[6]*x[8] + x[2]*x[3]^3*x[6]*x[8] +
x[1]*x[5]*x[6]^3*x[8] + x[1]*x[3]*x[6]^3*x[8] +
x[5]*x[3]*x[6]^3*x[8] + x[2]*x[5]*x[3]*x[8]^3 +
x[2]*x[5]*x[6]*x[8]^3 + x[2]*x[3]*x[6]*x[8]^3 +
x[5]*x[3]*x[6]*x[8]^3 + x[1]^3*x[2]*x[5]*x[4] +
x[1]*x[2]^3*x[5]*x[4] + x[1]^3*x[3]*x[6]*x[4] +
x[1]*x[3]^3*x[6]*x[4] + x[2]^3*x[3]*x[8]*x[4] +
x[2]*x[3]^3*x[8]*x[4] + x[1]^3*x[2]*x[5]*x[7] +
x[1]*x[2]*x[5]^3*x[7] + x[1]^3*x[3]*x[6]*x[7] +
x[1]*x[3]*x[6]^3*x[7] + x[5]^3*x[6]*x[8]*x[7] +
x[5]*x[6]^3*x[8]*x[7] + x[1]^3*x[2]*x[4]*x[7] +
x[1]^3*x[5]*x[4]*x[7] + x[1]^3*x[3]*x[4]*x[7] +
x[1]^3*x[6]*x[4]*x[7] + x[1]*x[2]*x[4]^3*x[7] +
x[1]*x[3]*x[4]^3*x[7] + x[1]*x[5]*x[4]*x[7]^3 +
x[1]*x[6]*x[4]*x[7]^3 + x[1]*x[2]^3*x[5]*x[9] +
x[1]*x[2]*x[5]^3*x[9] + x[2]^3*x[3]*x[8]*x[9] +
x[5]^3*x[6]*x[8]*x[9] + x[2]*x[3]*x[8]^3*x[9] +
x[5]*x[6]*x[8]^3*x[9] + x[1]*x[2]^3*x[4]*x[9] +
x[2]^3*x[5]*x[4]*x[9] + x[2]^3*x[3]*x[4]*x[9] +
x[2]^3*x[8]*x[4]*x[9] + x[1]*x[2]*x[4]^3*x[9] +
x[2]*x[3]*x[4]^3*x[9] + x[1]*x[5]^3*x[7]*x[9] +
x[2]*x[5]^3*x[7]*x[9] + x[5]^3*x[6]*x[7]*x[9] +
x[5]^3*x[8]*x[7]*x[9] + x[1]*x[4]^3*x[7]*x[9] +
x[2]*x[4]^3*x[7]*x[9] + x[1]*x[5]*x[7]^3*x[9] +
x[5]*x[6]*x[7]^3*x[9] + x[1]*x[4]*x[7]^3*x[9] +
x[5]*x[4]*x[7]^3*x[9] + x[2]*x[5]*x[4]*x[9]^3 +
x[2]*x[8]*x[4]*x[9]^3 + x[2]*x[5]*x[7]*x[9]^3 +
x[5]*x[8]*x[7]*x[9]^3 + x[2]*x[4]*x[7]*x[9]^3 +
x[5]*x[4]*x[7]*x[9]^3 + x[1]*x[3]^3*x[6]*x[10] +
x[1]*x[3]*x[6]^3*x[10] + x[2]*x[3]^3*x[8]*x[10] +
x[5]*x[6]^3*x[8]*x[10] + x[2]*x[3]*x[8]^3*x[10] +
x[5]*x[6]*x[8]^3*x[10] + x[1]*x[3]^3*x[4]*x[10] +
x[2]*x[3]^3*x[4]*x[10] + x[3]^3*x[6]*x[4]*x[10] +
x[3]^3*x[8]*x[4]*x[10] + x[1]*x[3]*x[4]^3*x[10] +
x[2]*x[3]*x[4]^3*x[10] + x[1]*x[6]^3*x[7]*x[10] +
x[5]*x[6]^3*x[7]*x[10] + x[3]*x[6]^3*x[7]*x[10] +
x[6]^3*x[8]*x[7]*x[10] + x[1]*x[4]^3*x[7]*x[10] +
x[3]*x[4]^3*x[7]*x[10] + x[1]*x[6]*x[7]^3*x[10] +
x[5]*x[6]*x[7]^3*x[10] + x[1]*x[4]*x[7]^3*x[10] +
x[6]*x[4]*x[7]^3*x[10] + x[2]*x[8]^3*x[9]*x[10] +
x[5]*x[8]^3*x[9]*x[10] + x[3]*x[8]^3*x[9]*x[10] +
x[6]*x[8]^3*x[9]*x[10] + x[2]*x[4]^3*x[9]*x[10] +
x[3]*x[4]^3*x[9]*x[10] + x[5]*x[7]^3*x[9]*x[10] +
x[6]*x[7]^3*x[9]*x[10] + x[2]*x[8]*x[9]^3*x[10] +
x[5]*x[8]*x[9]^3*x[10] + x[2]*x[4]*x[9]^3*x[10] +
x[8]*x[4]*x[9]^3*x[10] + x[5]*x[7]*x[9]^3*x[10] +
x[8]*x[7]*x[9]^3*x[10] + x[3]*x[6]*x[4]*x[10]^3 +
x[3]*x[8]*x[4]*x[10]^3 + x[3]*x[6]*x[7]*x[10]^3 +
x[6]*x[8]*x[7]*x[10]^3 + x[3]*x[4]*x[7]*x[10]^3 +
x[6]*x[4]*x[7]*x[10]^3 + x[3]*x[8]*x[9]*x[10]^3 +
x[6]*x[8]*x[9]*x[10]^3 + x[3]*x[4]*x[9]*x[10]^3 +
x[8]*x[4]*x[9]*x[10]^3 + x[6]*x[7]*x[9]*x[10]^3 +
x[8]*x[7]*x[9]*x[10]^3
 pv[26] = x[1]^2*x[2]^2*x[5]*x[3] +
x[1]^2*x[2]*x[5]^2*x[6] +
x[1]^2*x[2]*x[3]^2*x[6] +
x[1]^2*x[5]*x[3]*x[6]^2 +
x[1]*x[2]^2*x[5]^2*x[8] +
x[1]*x[2]^2*x[3]^2*x[8] +
x[1]*x[5]^2*x[6]^2*x[8] +
x[1]*x[3]^2*x[6]^2*x[8] +
x[2]^2*x[5]*x[3]*x[8]^2 +
x[2]*x[5]^2*x[6]*x[8]^2 +
x[2]*x[3]^2*x[6]*x[8]^2 +
x[5]*x[3]*x[6]^2*x[8]^2 +
x[1]^2*x[2]^2*x[5]*x[4] +
x[1]^2*x[3]^2*x[6]*x[4] +
x[2]^2*x[3]^2*x[8]*x[4] +
x[1]^2*x[2]*x[5]^2*x[7] +
x[1]^2*x[3]*x[6]^2*x[7] +
x[5]^2*x[6]^2*x[8]*x[7] +
x[1]^2*x[2]*x[4]^2*x[7] +
x[1]^2*x[3]*x[4]^2*x[7] +
x[1]^2*x[5]*x[4]*x[7]^2 +
x[1]^2*x[6]*x[4]*x[7]^2 +
x[1]*x[2]^2*x[5]^2*x[9] +
x[2]^2*x[3]*x[8]^2*x[9] +
x[5]^2*x[6]*x[8]^2*x[9] +
x[1]*x[2]^2*x[4]^2*x[9] +
x[2]^2*x[3]*x[4]^2*x[9] +
x[1]*x[5]^2*x[7]^2*x[9] +
x[5]^2*x[6]*x[7]^2*x[9] +
x[1]*x[4]^2*x[7]^2*x[9] +
x[2]^2*x[5]*x[4]*x[9]^2 +
x[2]^2*x[8]*x[4]*x[9]^2 +
x[2]*x[5]^2*x[7]*x[9]^2 +
x[5]^2*x[8]*x[7]*x[9]^2 +
x[2]*x[4]^2*x[7]*x[9]^2 +
x[5]*x[4]*x[7]^2*x[9]^2 +
x[1]*x[3]^2*x[6]^2*x[10] +
x[2]*x[3]^2*x[8]^2*x[10] +
x[5]*x[6]^2*x[8]^2*x[10] +
x[1]*x[3]^2*x[4]^2*x[10] +
x[2]*x[3]^2*x[4]^2*x[10] +
x[1]*x[6]^2*x[7]^2*x[10] +
x[5]*x[6]^2*x[7]^2*x[10] +
x[1]*x[4]^2*x[7]^2*x[10] +
x[2]*x[8]^2*x[9]^2*x[10] +
x[5]*x[8]^2*x[9]^2*x[10] +
x[2]*x[4]^2*x[9]^2*x[10] +
x[5]*x[7]^2*x[9]^2*x[10] +
x[3]^2*x[6]*x[4]*x[10]^2 +
x[3]^2*x[8]*x[4]*x[10]^2 +
x[3]*x[6]^2*x[7]*x[10]^2 +
x[6]^2*x[8]*x[7]*x[10]^2 +
x[3]*x[4]^2*x[7]*x[10]^2 +
x[6]*x[4]*x[7]^2*x[10]^2 +
x[3]*x[8]^2*x[9]*x[10]^2 +
x[6]*x[8]^2*x[9]*x[10]^2 +
x[3]*x[4]^2*x[9]*x[10]^2 +
x[6]*x[7]^2*x[9]*x[10]^2 +
x[8]*x[4]*x[9]^2*x[10]^2 +
x[8]*x[7]*x[9]^2*x[10]^2
 pv[27] = x[1]^3*x[5]^2*x[3] + x[2]^3*x[5]^2*x[3] +
x[1]^3*x[5]*x[3]^2 + x[2]^3*x[5]*x[3]^2 +
x[1]^3*x[2]^2*x[6] + x[2]^2*x[5]^3*x[6] +
x[2]^2*x[3]^3*x[6] + x[1]^3*x[2]*x[6]^2 +
x[2]*x[5]^3*x[6]^2 + x[2]*x[3]^3*x[6]^2 +
x[5]^2*x[3]*x[6]^3 + x[5]*x[3]^2*x[6]^3 +
x[1]^2*x[2]^3*x[8] + x[1]^2*x[5]^3*x[8] +
x[1]^2*x[3]^3*x[8] + x[1]^2*x[6]^3*x[8] +
x[1]*x[2]^3*x[8]^2 + x[1]*x[5]^3*x[8]^2 +
x[1]*x[3]^3*x[8]^2 + x[1]*x[6]^3*x[8]^2 +
x[5]^2*x[3]*x[8]^3 + x[5]*x[3]^2*x[8]^3 +
x[2]^2*x[6]*x[8]^3 + x[2]*x[6]^2*x[8]^3 +
x[1]^3*x[5]^2*x[4] + x[2]^3*x[5]^2*x[4] +
x[1]^3*x[6]^2*x[4] + x[3]^3*x[6]^2*x[4] +
x[2]^3*x[8]^2*x[4] + x[3]^3*x[8]^2*x[4] +
x[1]^3*x[5]*x[4]^2 + x[2]^3*x[5]*x[4]^2 +
x[1]^3*x[6]*x[4]^2 + x[3]^3*x[6]*x[4]^2 +
x[2]^3*x[8]*x[4]^2 + x[3]^3*x[8]*x[4]^2 +
x[1]^3*x[2]^2*x[7] + x[2]^2*x[5]^3*x[7] +
x[1]^3*x[3]^2*x[7] + x[3]^2*x[6]^3*x[7] +
x[5]^3*x[8]^2*x[7] + x[6]^3*x[8]^2*x[7] +
x[2]^2*x[4]^3*x[7] + x[3]^2*x[4]^3*x[7] +
x[1]^3*x[2]*x[7]^2 + x[2]*x[5]^3*x[7]^2 +
x[1]^3*x[3]*x[7]^2 + x[3]*x[6]^3*x[7]^2 +
x[5]^3*x[8]*x[7]^2 + x[6]^3*x[8]*x[7]^2 +
x[2]*x[4]^3*x[7]^2 + x[3]*x[4]^3*x[7]^2 +
x[5]^2*x[4]*x[7]^3 + x[6]^2*x[4]*x[7]^3 +
x[5]*x[4]^2*x[7]^3 + x[6]*x[4]^2*x[7]^3 +
x[1]^2*x[2]^3*x[9] + x[1]^2*x[5]^3*x[9] +
x[2]^3*x[3]^2*x[9] + x[5]^3*x[6]^2*x[9] +
x[3]^2*x[8]^3*x[9] + x[6]^2*x[8]^3*x[9] +
x[1]^2*x[4]^3*x[9] + x[3]^2*x[4]^3*x[9] +
x[1]^2*x[7]^3*x[9] + x[6]^2*x[7]^3*x[9] +
x[1]*x[2]^3*x[9]^2 + x[1]*x[5]^3*x[9]^2 +
x[2]^3*x[3]*x[9]^2 + x[5]^3*x[6]*x[9]^2 +
x[3]*x[8]^3*x[9]^2 + x[6]*x[8]^3*x[9]^2 +
x[1]*x[4]^3*x[9]^2 + x[3]*x[4]^3*x[9]^2 +
x[1]*x[7]^3*x[9]^2 + x[6]*x[7]^3*x[9]^2 +
x[5]^2*x[4]*x[9]^3 + x[8]^2*x[4]*x[9]^3 +
x[5]*x[4]^2*x[9]^3 + x[8]*x[4]^2*x[9]^3 +
x[2]^2*x[7]*x[9]^3 + x[8]^2*x[7]*x[9]^3 +
x[2]*x[7]^2*x[9]^3 + x[8]*x[7]^2*x[9]^3 +
x[1]^2*x[3]^3*x[10] + x[2]^2*x[3]^3*x[10] +
x[1]^2*x[6]^3*x[10] + x[5]^2*x[6]^3*x[10] +
x[2]^2*x[8]^3*x[10] + x[5]^2*x[8]^3*x[10] +
x[1]^2*x[4]^3*x[10] + x[2]^2*x[4]^3*x[10] +
x[1]^2*x[7]^3*x[10] + x[5]^2*x[7]^3*x[10] +
x[2]^2*x[9]^3*x[10] + x[5]^2*x[9]^3*x[10] +
x[1]*x[3]^3*x[10]^2 + x[2]*x[3]^3*x[10]^2 +
x[1]*x[6]^3*x[10]^2 + x[5]*x[6]^3*x[10]^2 +
x[2]*x[8]^3*x[10]^2 + x[5]*x[8]^3*x[10]^2 +
x[1]*x[4]^3*x[10]^2 + x[2]*x[4]^3*x[10]^2 +
x[1]*x[7]^3*x[10]^2 + x[5]*x[7]^3*x[10]^2 +
x[2]*x[9]^3*x[10]^2 + x[5]*x[9]^3*x[10]^2 +
x[6]^2*x[4]*x[10]^3 + x[8]^2*x[4]*x[10]^3 +
x[6]*x[4]^2*x[10]^3 + x[8]*x[4]^2*x[10]^3 +
x[3]^2*x[7]*x[10]^3 + x[8]^2*x[7]*x[10]^3 +
x[3]*x[7]^2*x[10]^3 + x[8]*x[7]^2*x[10]^3 +
x[3]^2*x[9]*x[10]^3 + x[6]^2*x[9]*x[10]^3 +
x[3]*x[9]^2*x[10]^3 + x[6]*x[9]^2*x[10]^3
 v[20:31] = pv[16:27]
 Primary_invariants= [ x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[10], x[1]*x[2] + x[1]*x[3] + x[1]*x[4] + x[1]*x[5] + x[1]*x[7] + x[1]*x[8] + x[2]*x[3] + x[2]*x[4] + x[2]*x[6] + x[2]*x[7] + x[2]*x[9] + x[3]*x[5] + x[3]*x[6] + x[3]*x[8] + x[3]*x[9] + x[4]*x[5] + x[4]*x[6] + x[4]*x[7] + x[4]*x[10] + x[5]*x[6] + x[5]*x[8] + x[5]*x[10] + x[6]*x[9] + x[6]*x[10] + x[7]*x[8] + x[7]*x[9] + x[7]*x[10] + x[8]*x[9] + x[8]*x[10] + x[9]*x[10], x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[5]^2 + x[6]^2 + x[7]^2 + x[8]^2 + x[9]^2 + x[10]^2, x[1]*x[2]*x[4] + x[1]*x[2]*x[7] + x[1]*x[3]*x[5] + x[1]*x[3]*x[8] + x[1]*x[4]*x[7] + x[1]*x[5]*x[8] + x[2]*x[3]*x[6] + x[2]*x[3]*x[9] + x[2]*x[4]*x[7] + x[2]*x[6]*x[9] + x[3]*x[5]*x[8] + x[3]*x[6]*x[9] + x[4]*x[5]*x[6] + x[4]*x[5]*x[10] + x[4]*x[6]*x[10] + x[5]*x[6]*x[10] + x[7]*x[8]*x[9] + x[7]*x[8]*x[10] + x[7]*x[9]*x[10] + x[8]*x[9]*x[10], x[1]^3 + x[2]^3 + x[3]^3 + x[4]^3 + x[5]^3 + x[6]^3 + x[7]^3 + x[8]^3 + x[9]^3 + x[10]^3, x[1]*x[2]*x[3]*x[10] + x[1]*x[4]*x[5]*x[9] + x[1]*x[6]*x[7]*x[8] + x[1]*x[6]*x[9]*x[10] + x[2]*x[4]*x[6]*x[8] + x[2]*x[5]*x[7]*x[9] + x[2]*x[5]*x[8]*x[10] + x[3]*x[4]*x[7]*x[10] + x[3]*x[4]*x[8]*x[9] + x[3]*x[5]*x[6]*x[7], x[1]^4 + x[2]^4 + x[3]^4 + x[4]^4 + x[5]^4 + x[6]^4 + x[7]^4 + x[8]^4 + x[9]^4 + x[10]^4, 2*x[1]*x[2]*x[3]*x[4]*x[8] + 2*x[1]*x[2]*x[3]*x[4]*x[9] + 2*x[1]*x[2]*x[3]*x[5]*x[7] + 2*x[1]*x[2]*x[3]*x[5]*x[9] + 2*x[1]*x[2]*x[3]*x[6]*x[7] + 2*x[1]*x[2]*x[3]*x[6]*x[8] + 2*x[1]*x[2]*x[4]*x[5]*x[8] + 2*x[1]*x[2]*x[4]*x[5]*x[10] + 2*x[1]*x[2]*x[4]*x[6]*x[9] + 2*x[1]*x[2]*x[4]*x[6]*x[10] + 2*x[1]*x[2]*x[5]*x[7]*x[8] + 2*x[1]*x[2]*x[6]*x[7]*x[9] + 2*x[1]*x[2]*x[7]*x[8]*x[10] + 2*x[1]*x[2]*x[7]*x[9]*x[10] + 2*x[1]*x[3]*x[4]*x[5]*x[7] + 2*x[1]*x[3]*x[4]*x[5]*x[10] + 2*x[1]*x[3]*x[4]*x[7]*x[8] + 2*x[1]*x[3]*x[5]*x[6]*x[9] + 2*x[1]*x[3]*x[5]*x[6]*x[10] + 2*x[1]*x[3]*x[6]*x[8]*x[9] + 2*x[1]*x[3]*x[7]*x[8]*x[10] + 2*x[1]*x[3]*x[8]*x[9]*x[10] + 2*x[1]*x[4]*x[5]*x[6]*x[7] + 2*x[1]*x[4]*x[5]*x[6]*x[8] + 2*x[1]*x[4]*x[6]*x[7]*x[10] + 2*x[1]*x[4]*x[7]*x[8]*x[9] + 2*x[1]*x[4]*x[7]*x[9]*x[10] + 2*x[1]*x[5]*x[6]*x[8]*x[10] + 2*x[1]*x[5]*x[7]*x[8]*x[9] + 2*x[1]*x[5]*x[8]*x[9]*x[10] + 2*x[2]*x[3]*x[4]*x[6]*x[7] + 2*x[2]*x[3]*x[4]*x[6]*x[10] + 2*x[2]*x[3]*x[4]*x[7]*x[9] + 2*x[2]*x[3]*x[5]*x[6]*x[8] + 2*x[2]*x[3]*x[5]*x[6]*x[10] + 2*x[2]*x[3]*x[5]*x[8]*x[9] + 2*x[2]*x[3]*x[7]*x[9]*x[10] + 2*x[2]*x[3]*x[8]*x[9]*x[10] + 2*x[2]*x[4]*x[5]*x[6]*x[7] + 2*x[2]*x[4]*x[5]*x[6]*x[9] + 2*x[2]*x[4]*x[5]*x[7]*x[10] + 2*x[2]*x[4]*x[7]*x[8]*x[9] + 2*x[2]*x[4]*x[7]*x[8]*x[10] + 2*x[2]*x[5]*x[6]*x[9]*x[10] + 2*x[2]*x[6]*x[7]*x[8]*x[9] + 2*x[2]*x[6]*x[8]*x[9]*x[10] + 2*x[3]*x[4]*x[5]*x[6]*x[8] + 2*x[3]*x[4]*x[5]*x[6]*x[9] + 2*x[3]*x[4]*x[5]*x[8]*x[10] + 2*x[3]*x[4]*x[6]*x[9]*x[10] + 2*x[3]*x[5]*x[7]*x[8]*x[9] + 2*x[3]*x[5]*x[7]*x[8]*x[10] + 2*x[3]*x[6]*x[7]*x[8]*x[9] + 2*x[3]*x[6]*x[7]*x[9]*x[10] + 2*x[4]*x[5]*x[7]*x[9]*x[10] + 2*x[4]*x[5]*x[8]*x[9]*x[10] + 2*x[4]*x[6]*x[7]*x[8]*x[10] + 2*x[4]*x[6]*x[8]*x[9]*x[10] + 2*x[5]*x[6]*x[7]*x[8]*x[10] + 2*x[5]*x[6]*x[7]*x[9]*x[10], x[1]^5 + x[2]^5 + x[3]^5 + x[4]^5 + x[5]^5 + x[6]^5 + x[7]^5 + x[8]^5 + x[9]^5 + x[10]^5, x[1]^6 + x[2]^6 + x[3]^6 + x[4]^6 + x[5]^6 + x[6]^6 + x[7]^6 + x[8]^6 + x[9]^6 + x[10]^6 ]
return Primary_invariants, v, pv

end
x = rand(10)
display(invariants_Q10_check(x))
