
//Define the field
K := RationalField();
//Define the symmetry group (line that can vary)
// G := PermutationGroup<6| (2,3,1)(5,6,4)>; //4B
GROUP_DEF;
// G := PermutationGroup<10 | (2,3)(5,6)(10,9), (2,3,4)(5,6,7)(8,10,9), (1,3,4,2)(5,6,10,9)(7,8)>; //5B
//Define the invariant ring
R := InvariantRing(G, K);
R;

R0 := PrimaryInvariants(R);
R0;
[TotalDegree(R0[i]): i in [1..#R0]];

// R1 := SecondaryInvariants(R);
// R1;
// [TotalDegree(R1[i]): i in [1..#R1]];

R2 := IrreducibleSecondaryInvariants(R);
R2;
[TotalDegree(R2[i]): i in [1..#R2]];

// F := FundamentalInvariants(R);
// F;
// [TotalDegree(F[i]): i in [1..#F]];

RR2 := IrreducibleSecondaryInvariants(R);
RR2;

A, Q := Algebra(R);
A;
Q;
