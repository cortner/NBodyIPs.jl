#!/bin/bash

# Script generating files containing primary and secondary invariants for a given group

# -------------------------------------------
# Parameters
# -------------------------------------------
GROUP_DEF="PermutationGroup<10 | (2,3)(5,6)(10,9), (2,3,4)(5,6,7)(8,10,9), (1,3,4,2)(5,6,10,9)(7,8)>"
GROUP_NAME="BA_5B"

PREFSEC="SEC" #prefix for the secondary invariants
PREFIRRSEC="IS" #prefix for the irreducible secondary invariants

MAGMA_RUN=1

#------------------------------------------------------------
#Move files in data folder and delete useless files
#------------------------------------------------------------
# create directory
mkdir -p ../data/$GROUP_NAME

#Define the file names used later and print out their names
filename_log="../data/$GROUP_NAME/log.txt"
fn_jl_check="../data/$GROUP_NAME/$GROUP_NAME""_non_efficient_invariants.jl"
fn_jl_irr_inv="../data/$GROUP_NAME/$GROUP_NAME""_irr_invariants.jl"
fn_jl_prim_inv="../data/$GROUP_NAME/$GROUP_NAME""_prim_invariants.jl"
fn_jl_sec_rel_inv="../data/$GROUP_NAME/$GROUP_NAME""_relations_invariants.jl"
fn_jl_group_elts="../data/$GROUP_NAME/$GROUP_NAME""_group_elements.jl"

ECHO Output files:

ECHO $filename_log
ECHO $fn_jl_check
ECHO $fn_jl_irr_inv
ECHO $fn_jl_prim_inv
ECHO $fn_jl_sec_rel_inv
ECHO $fn_jl_group_elts

# -------------------------------------------
# Run of the Magma computation
# -------------------------------------------
if [ $MAGMA_RUN -eq 1 ]
then
#put the paramters into the input file
cp Nbody_inv_auto_generation_ba.m Nbody_run.m;

ECHO $GROUP_DEF

sed -i -e "s/GROUP_DEF/$GROUP_DEF/g" Nbody_run.m;

#connect to galois and copy the input files
# scp pack_opt_primaries.m dusson@galois.warwick.ac.uk: ;
scp Nbody_run.m dusson@galois.warwick.ac.uk: ;

#run the magma computation
ssh dusson@galois.warwick.ac.uk << EOF
magma Nbody_run.m
EOF

#remove now useless file
rm Nbody_run.m;

#copy the output on the local machine
scp dusson@galois.warwick.ac.uk:logNbody_output.txt .;

# change the name of the output file
mv logNbody_output.txt $filename_log

#Remove useless file
rm NBody_run.m-e
fi


#------------------------------------------------------------
#extract number of primaries, secondaries and irreducible secondaries
#------------------------------------------------------------

#get number of primaries
NBprimaries=$(gsed '0,/^nb_primary_invariants_begin/d;/^nb_primary_invariants_end/,$d' $filename_log)
ECHO "Nb of primaries="$NBprimaries

#get number of irreducible secondaries
NBirrsec=$(gsed '0,/^nb_irr_sec_invariants_begin/d;/^nb_irr_sec_invariants_end/,$d' $filename_log)
ECHO "Nb of irreducible secondaries="$NBirrsec

#get number of secondaries
NBsec=$(gsed '0,/^nb_secondaries_begin/d;/^nb_secondaries_end/,$d' $filename_log)
ECHO "Nb of secondaries="$NBsec

#------------------------------------------------------------
#extract degrees for primaries, secondaries and irreducible secondaries
#------------------------------------------------------------

#get degrees of primaries
Degprimaries=$(gsed '0,/^prim_inv_tot_degree_begin/d;/^prim_inv_tot_degree_end/,$d' $filename_log)
ECHO "Degrees of primaries="$Degprimaries

#get degrees of irreducible secondaries
Degirrsec=$(gsed '0,/^irr_sec_degrees_begin/d;/^irr_sec_degrees_end/,$d' $filename_log)
ECHO "Degrees of irreducible secondaries="$Degirrsec

#get degrees of secondaries
Degsec=$(gsed '0,/^degrees_secondaries_begin/d;/^degrees_secondaries_end/,$d' $filename_log)
ECHO "Degrees of secondaries="$Degsec


#------------------------------------------------------------
#write the julia file for the non-efficient invariants function
#------------------------------------------------------------
# Generate julia file with function computing primary and secondary invariants (not efficient but hopefully correct)
cp $filename_log $fn_jl_check
# Pick lines with primaries, irreducible secondaries and secondaries
gsed -i '0,/^primary_invariants_begin$/d' $fn_jl_check
gsed -i '/^primary_invariants_end$/,/^irreducible_sec_invariants_begin$/d' $fn_jl_check
gsed -i '/^irreducible_sec_invariants_end$/,/^inv_relations_begin$/d' $fn_jl_check
gsed -i '/^inv_relations_end$/,$d' $fn_jl_check


echo "v=zeros($NBsec"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "pv=zeros($NBirrsec"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "prim=zeros($NBprimaries"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "function invariants_$GROUP_NAME""_check(x)" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check


echo "return prim, v, pv"  >> $fn_jl_check
echo "" >> $fn_jl_check
echo "end" >> $fn_jl_check

# replace variables for the primaries
# xi -> x[i]
ECHO "replacing variables for the primaries"
for a in `seq $(($NBprimaries)) -1 1`; do
		OLD="x$a" ;
		NEW="x\[$a\]" ;
		sed -i '' "s/$OLD/$NEW/g" $fn_jl_check
done

# replace variables for the secondaries
# hi -> pv[i]
ECHO "replacing variables for the secondaries relations"
for a in `seq $(($NBirrsec)) -1 1`; do
		OLD="h$a" ;
		NEW="pv\[$a\]" ;
		sed -i '' "s/$OLD/$NEW/g" $fn_jl_check
done

# Moving lines starting with +
# ---------------------------------------------------------
gsed -i '$!N;s/\n\s*+/ +/;P;D' $fn_jl_check

#------------------------------------------------------------
#Generate a file with only monomials of primaries
#------------------------------------------------------------
echo "generating file with monomials of primaries"

cp $filename_log $fn_jl_prim_inv

# gsed -i '0,/^prim_inv_mon_begin/d;/^prim_inv_mon_end/,$d' $fn_jl_prim_inv
#
# ECHO "replacing variables for the primaries"
# for a in `seq $(($NBprimaries)) -1 1`; do
# 		OLD="x$a" ;
# 		NEW="x\[$a\]" ;
# 		sed -i '' "s/$OLD/$NEW/g" $fn_jl_prim_inv
# done

gsed -i '0,/^prim_inv_mon_list_begin/d;/^prim_inv_mon_list_end/,$d' $fn_jl_prim_inv


#------------------------------------------------------------
#Generate a file with only monomials of irreducible secondaries
#------------------------------------------------------------
echo "generating file with monomials of irreducible secondaries"

cp $filename_log $fn_jl_irr_inv
# gsed -i '0,/^irrsec_inv_mon_begin/d;/^irrsec_inv_mon_end/,$d' $fn_jl_irr_inv

gsed -i '0,/^irrsec_mon_list_begin/d;/^irrsec_mon_list_end/,$d' $fn_jl_irr_inv

#------------------------------------------------------------
#Generate a file with relations between irreducible and secondary invariants
#------------------------------------------------------------
echo "generating file with relations between secondaries and irreducible secondaries"

cp $filename_log $fn_jl_sec_rel_inv
gsed -i '0,/^inv_relations_begin/d;/^inv_relations_end/,$d' $fn_jl_sec_rel_inv
sed -i '' "s/h/$PREFIRRSEC/g" $fn_jl_sec_rel_inv
sed -i '' "s/v\[/$PREFSEC/g" $fn_jl_sec_rel_inv
sed -i '' "s/\]/ /g" $fn_jl_sec_rel_inv

#------------------------------------------------------------
#Generate a file with group elements
#------------------------------------------------------------
#get number of elements in the group
NB_elts=$(gsed '0,/^nb_group_elements_begin/d;/^nb_group_elements_end/,$d' $filename_log)
ECHO "Number of elements="$NB_elts


echo "generating file with group elements"

cp $filename_log $fn_jl_group_elts
gsed -i '0,/^group_def_begin/d;/^group_def_end/,$d' $fn_jl_group_elts
gsed -i 's/Gr_elts/G_'$GROUP_NAME'/g' $fn_jl_group_elts

echo "simplex_permutations(x::SVector{$NBprimaries}) = [x[G_$GROUP_NAME[i]] for i=1:$NB_elts]" >>  $fn_jl_group_elts
