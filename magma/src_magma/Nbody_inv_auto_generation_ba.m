
//Group symmetry - body order ([4] for 4-body, [5] for 5-body, etc.)
// nbody := NBODY;

//Degree (up to)
// deg := DEGREE;

//Outputfile name
// outputfile := "Nbody_output.txt";
//Logfile name
logfile := "logNbody_output.txt";

//Attach package from Braams
// Attach("pack_opt_primaries.m");

//Prescribe outputfile and logfile
// SetOutputFile(outputfile: Overwrite := true);
SetLogFile(logfile: Overwrite := true);



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

printf "nb_primary_invariants_begin\n";
#R0;
printf "nb_primary_invariants_end\n";

printf "primary_invariants_begin\n";
for i:=1 to #R0 do
  printf "prim[";
  printf IntegerToString(i);
  printf "]=";
  R0[i];
  printf "\n";
end for;
printf "primary_invariants_end\n";

printf "prim_inv_tot_degree_begin\n";
// printf "degree_prim=";
[TotalDegree(R0[i]): i in [1..#R0]];
printf "prim_inv_tot_degree_end\n";

printf "prim_inv_mon_begin\n";
for i:=1 to #R0 do
  printf "prim[";
  printf IntegerToString(i);
  printf "]=";
  LeadingMonomial(R0[i]);
  printf "\n";
end for;
printf "prim_inv_mon_end\n";


R2 := IrreducibleSecondaryInvariants(R);

printf "nb_irr_sec_invariants_begin\n";
#R2;
printf "nb_irr_sec_invariants_end\n";

printf "irreducible_sec_invariants_begin\n";
for i:=1 to #R2 do
  printf "pv[";
  printf IntegerToString(i);
  printf "]=";
  R2[i];
  printf "\n";
end for;
printf "irreducible_sec_invariants_end\n";

printf "irr_sec_degrees_begin\n";
[TotalDegree(R2[i]): i in [1..#R2]];
printf "irr_sec_degrees_end\n";


printf "irrsec_inv_mon_begin\n";
for i:=1 to #R2 do
  printf "pv[";
  printf IntegerToString(i);
  printf "]=";
  LeadingMonomial(R2[i]);
  printf "\n";
end for;
printf "irrsec_inv_mon_end\n";

// F := FundamentalInvariants(R);
// F;
// [TotalDegree(F[i]): i in [1..#F]];

// RR2 := IrreducibleSecondaryInvariants(R);
// RR2;

A, Q := Algebra(R);
A;
printf "inv_relations_begin\n";
for i:=1 to #Q do
  printf "v[";
  printf IntegerToString(i);
  printf "] = ";
  Q[i];
end for;
printf "inv_relations_end\n";


R1 := SecondaryInvariants(R);

printf "nb_secondaries_begin\n";
#Q;
printf "nb_secondaries_end\n";

printf "degrees_secondaries_begin\n";
[TotalDegree(R1[i]): i in [1..#R1]];
printf "degrees_secondaries_end\n";





// // printf " Nbody group where N="*IntegerToString(nbody)*"\n";
// //  " Maximal polynomial degree="*IntegerToString(deg)*"\n";
//
// // Write(outputfile, "N-body group where N=" cat IntegerToString(nbody));
// // Write(outputfile, "Maximal polynomial degree:" cat IntegerToString(deg));
//
// //Define the field
// K := RationalField();
// //Define the symmetry group
// G := MolSymGen([nbody]);
// //Define the invariant ring
// R := InvariantRing(G, K);
//
// //Primary and secondary invariants
// PSI :=MolInvRngGen([nbody],deg);
//
// // S<t>:=PowerSeriesRing(Rationals(),deg+1);
// // Write(outputfile, "Group degree and order:");
// // Write(outputfile, Degree(Group(PSI))) ;
// // Write(outputfile, #Group(PSI)) ;
// // dnpr:=[0:i in [0..Precision(S)-1]];
// // for f in PSI`PrimaryInvariants do
// //  i:=TotalDegree(f);
// //  if i lt Precision(S) then
// //   dnpr[i+1]:=dnpr[i+1]+1;
// //  end if;
// // end for;
// // Write(outputfile, "Degrees of Primaries:");
// // Write(outputfile, [TotalDegree(f): f in PSI`PrimaryInvariants]);
// // Write(outputfile, "Number of primaries at each degree, and sums:");
// // Write(outputfile, [dnpr[1+i]:i in [0..Precision(S)-1]]);
// // Write(outputfile, [&+[dnpr[1+j]:j in [0..i]]:i in [0..Precision(S)-1]]);
// // dnb:=S!PSI`HilbertSeries;
// // dnpb:=1/&*[1-t^TotalDegree(f):f in PSI`PrimaryInvariants];
// // Write(outputfile, "Dimensions of the Primaries Ring:");
// // Write(outputfile, [Coefficient(dnpb,i):i in [0..Precision(S)-1]]);
// // Write(outputfile, [&+[Coefficient(dnpb,j):j in [0..i]]:i in [0..Precision(S)-1]]);
// // dnsc:=dnb*&*[1-t^TotalDegree(f):f in PSI`PrimaryInvariants];
// // Write(outputfile, "Expected Numbers of Secondaries, from Degree 0:");
// // Write(outputfile, [Coefficient(dnsc,i):i in [0..Precision(S)-1]]);
// // Write(outputfile, [&+[Coefficient(dnsc,j):j in [0..i]]:i in [0..Precision(S)-1]]);
// // Write(outputfile, "Molien/Hilbert Series from Degree 0, and sums:");
// // Write(outputfile, [Coefficient(dnb,i):i in [0..Precision(S)-1]]);
// // Write(outputfile, [&+[Coefficient(dnb,j):j in [0..i]]:i in [0..Precision(S)-1]]);
//
//
// // Write(outputfile, "Primary invariants");
// // Write(outputfile, PSI`PrimaryInvariants);
//
// printf " Names_of_variables=";
// Name := MolSymCGNames(nbody);
// Name;
// printf "\n end_names_of_variables \n";
//
//
// // printf " Primary_invariants=\n";
// // PSI`PrimaryInvariants;
//
// // Write(outputfile, "Secondary invariants");
// // Write(outputfile, PSI`SecondaryInvariants);
// //
// // Write(outputfile, "Degrees of the secondary invariants");
// // Write(outputfile, [TotalDegree(PSI`SecondaryInvariants[i]): i in [1..#PSI`SecondaryInvariants]]);
//
//
// // //Fundamental invariants
// // F := FundamentalInvariants(R);
// // Write(outputfile, "Fundamental invariants");
// // Fd := [];
// // count := 0;
// // for f in F do
// //  i:=TotalDegree(f);
// //  if i le deg then
// //  count := count +1;
// //   Fd[count]:=f;
// //  end if;
// // end for;
// // Fd;
// // #Fd;
//
// // Write(outputfile, Fd);
// // Write(outputfile, "Degrees of the fundamental invariants");
// // Write(outputfile, [TotalDegree(Fd[i]): i in [1..#Fd]]);
//
// // //Molien series
// // M<t> := MolienSeries(G);
// // Write(outputfile, "Molien series");
// // Write(outputfile, M);
// Sec := PSI`SecondaryInvariants;
// IrrSec := PSI`IrreducibleSecondaries;
//
// printf " Nb_secondary_invariants=";
// #Sec;
// printf "\n";
//
// printf " Nb_irr_sec_invariants=";
// #IrrSec;
//
//
//
// // P := PolynomialRing(R);
// // AssignNames(~PSI, Name);
// // R0 := MolInvPrimsRngCG(nbody);
// // AssignNames(~R0, Name);
// // R0;
// // P;
// //
// // R0 := PrimaryInvariants(R);
// // R0;
//
// // PSI:Names := MolSymCGNames(nbody);
// // printf " Primary invariants=\n";
// // PSI`PrimaryInvariants;

exit;
