using Combinatorics, StaticArrays, NBodyIPs

include("invariants_generator.jl")
include("misc.jl")

# Parameters
NBody = 5;
Deg = 6;
# --------------
NBlengths = Int(NBody*(NBody-1)/2);

# -------------------------------------------
#
# Generate irreducible secondaries
#
# -------------------------------------------
# filename = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text.jl";
filename1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text1.jl";
filename2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text2.jl";
filename3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text3.jl";
filenamedata = "magma/data/NB_$NBody""_deg_$Deg""_irr_invariants.jl";
preword = "# Irreducible secondaries for NBody=$NBody"*"and deg=$Deg \n"
pref = "IS"

generate_invariants(filenamedata,filename1,filename2,filename3,NBlengths,Deg,preword,pref)


# -------------------------------------------
#
# Generate primary invariants
#
# -------------------------------------------
filename1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text1.jl";
filename2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text2.jl";
filename3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text3.jl";
filenamedata = "magma/data/NB_$NBody""_deg_$Deg""_prim_invariants.jl";
preword = "# Primary invariants for NBody=$NBody"*"and deg=$Deg \n"
pref = "P"

generate_invariants(filenamedata,filename1,filename2,filename3,NBlengths,Deg,preword,pref)


# -------------------------------------------
#
# Generate function with all invariants
#
# -------------------------------------------
prefix_file = "using StaticArrays \n using BenchmarkTools: @btime \n\n include(\"fastpolys.jl\") \n using FastPolys \n"
