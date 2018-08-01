using Combinatorics, StaticArrays

include("invariants_generator_ba.jl")
include("misc.jl")
include("inv_monomials.jl")

# Parameters
#TODO: check that prefix are the same as the ones in generate_invariants.sh
GROUP_NAME="BA_5B"
prefsec = "SEC" #prefix for the secondaries
prefirrsec = "IS" #prefix for the irreducible secondaries
prefprim = "P" #prefix for the primaries
# --------------

# -------------------------------------------
#
# Generate irreducible secondaries
#
# -------------------------------------------
filenameirrsec1 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/"*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text1.jl";
filenameirrsec2 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/"*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text2.jl";
filenameirrsec3 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/"*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text3.jl";
filenameirrsec4 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/"*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text4.jl";
filenameirrsec5 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/"*GROUP_NAME*"/"*GROUP_NAME*"_irr_sec_text5.jl";

# -------------------------------------------
#
# Define data files
#
# -------------------------------------------
filenameirrsecdata = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/"*GROUP_NAME*"/"*GROUP_NAME*"_irr_invariants.jl";
preword = "# Irreducible secondaries for group "*GROUP_NAME*"\n"

NB_irrsec = countlines(filenameirrsecdata)

max_exp_irrsec = generate_invariants(filenameirrsecdata,filenameirrsec1,filenameirrsec2,filenameirrsec3,filenameirrsec4,filenameirrsec5,NBlengths,Deg,preword,prefirrsec)

# filename_deg = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/*GROUP_NAME*"/"*GROUP_NAME*_text_deg.jl";
# -------------------------------------------
#
# Generate primary invariants files
#
# -------------------------------------------
filenameprim1 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text1.jl";
filenameprim2 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text2.jl";
filenameprim3 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text3.jl";
filenameprim4 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text4.jl";
filenameprim5 = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text5.jl";
filenameprimdata = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_prim_invariants.jl";
preword = "# Primary invariants for NBody=$NBody"*"and deg=$Deg \n"
NB_prim = countlines(filenameprimdata)

max_exp_prim = generate_invariants(filenameprimdata,filenameprim1,filenameprim2,filenameprim3,filenameprim4,filenameprim5,NBlengths,Deg,preword,prefprim)
# -------------------------------------------
#
# Secondary invariants (relations with irreducible secondaries)
#
# -------------------------------------------
filenamesec = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_relations_invariants.jl";
NB_secondary = countlines(filenamesec);
# -------------------------------------------
#
# Derivatives of secondary invariants (relations with irreducible secondaries)
#
# -------------------------------------------
filenamesec_d = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_relations_invariants_derivatives.jl";

open(filenamesec_d, "w") do f
end

fileI = open(filenamesec)
line = readlines(fileI)
part1,part2 = split(line[1], "=")
repl1 = replace(part1, prefsec, "d"*prefsec)
open(filenamesec_d, "a") do f
    write(f, repl1, " = @SVector zeros($NBlengths) \n")
end
for i=2:length(line)
    if contains(line[i], "*")
        part1,part2 = split(line[i], "=")
        part2_1,part2_2 = split(part2, "*")
        repl1 = replace(part1, prefsec, "d"*prefsec)
        repl2_1 = replace(part2_1, prefirrsec, "d"*prefirrsec)
        repl2_2 = replace(part2_2, prefirrsec, "d"*prefirrsec)
        open(filenamesec_d, "a") do f
            write(f, repl1, " = ", repl2_1, "*", part2_2, "+", part2_1, "*", repl2_2, "\n")
        end
    else
        open(filenamesec_d, "a") do f
            repl1 = replace(line[i], prefsec, "d"*prefsec)
            repl2 = replace(repl1, prefirrsec, "d"*prefirrsec)
            write(f, repl2, "\n")
        end
    end
end
# -------------------------------------------
#
# Generate function with degrees of primary invariants and secondary
#
# -------------------------------------------
Mon_prim, coef_list_prim, deg_prim = generate_inv_mon(filenameprimdata,NBlengths,Deg)

Mon_irrsec, coef_list_irrsec = generate_inv_mon(filenameirrsecdata,NBlengths,Deg)

# generate list of PolyMon for irreducible secondaries
IrrSecMonPol = CPolyMon{NBlengths}[]
for i=1:length(Mon_irrsec)
    push!(IrrSecMonPol,CPolyMon(CMonList([mon_repr(SVector(Mon_irrsec[i]...))]),[coef_list_irrsec[i]]))
end
IrrSecMonPol


SecMonPol = CPolyMon{NBlengths}[]

fileI = open(filenamesec)
line = readlines(fileI)
# first line contains sec invariant =1, we remove it
push!(SecMonPol,constpoly(NBlengths))
for i=2:length(line)
    part1,part2 = split(line[i], "=")
    Part1 = replace(part1, prefsec, "")
    @assert parse(Int64,Part1) == i
    if contains(line[i], "*")
        part2_1,part2_2 = split(part2, "*")
        Part2_1 = replace(part2_1, prefirrsec, "")
        Part2_2 = replace(part2_2, prefirrsec, "")
        int1 = parse(Int64,Part2_1)
        int2 = parse(Int64,Part2_2)

        IrrSec1 = IrrSecMonPol[int1]
        IrrSec2 = IrrSecMonPol[int2]
        push!(SecMonPol,IrrSec1*IrrSec2)
    else
        Part2 = replace(part2, prefirrsec, "")
        int = parse(Int64,Part2)
        push!(SecMonPol,IrrSecMonPol[int])
    end
end
Deg_sec = [];
for i=1:length(SecMonPol)
    push!(Deg_sec,deg_monlist(SecMonPol[i]))
end
Deg_sec


PrimMonPol = CPolyMon{NBlengths}[]
for i=1:length(Mon_prim)
    push!(PrimMonPol,CPolyMon(CMonList([mon_repr(SVector(Mon_prim[i]...))]),[coef_list_prim[i]]))
end
PrimMonPol

Deg_prim = [];
for i=1:NB_prim
    push!(Deg_prim,deg_monlist(PrimMonPol[i]))
end
Deg_prim

open(filename_deg, "a") do f
    write(f, "tdegrees(::Val{5}) =
          (@SVector [");
    for i=1:NB_prim
        deg = Deg_prim[i]
        write(f, "$deg, ")
    end
    write(f, "]),
    (@SVector [");
    for i=1:NB_secondary
        deg = Deg_sec[i]
        write(f, "$deg, ")
    end
    write(f, "])");
end

# -------------------------------------------
#
# Generate function with all invariants
#
# -------------------------------------------
file = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_invariants.jl";

open(file, "w") do f
    write(f, "module NB5I \n\n")
    write(f, "using NBodyIPs.FastPolys \n")
    write(f, "using StaticArrays \n")
    write(f, "using BenchmarkTools: @btime \n\n")
    write(f, "import NBodyIPs.tdegrees \n\n")
    # write(f, "include(\"fastpolys.jl\") \n")
    # write(f, "using FastPolys \n\n\n\n")

    #write the definition of the degrees of the invariants
    degree_inv = open(io->read(io), filename_deg)
    write(f, degree_inv)
    write(f, "\n\n")

    # write the definition of the constant vectors
    prim1 = open(io->read(io), filenameprim1)
    # read(filenameprim1)
    write(f, prim1)
    irrsec1 = open(io->read(io), filenameirrsec1)
    # read(filenameirrsec1)
    write(f, irrsec1)

    # write the definitions of the types
    write(f, "\n")
    prim2 = open(io->read(io), filenameprim2)
    # read(filenameprim2)
    write(f, prim2)
    irrsec2 = open(io->read(io), filenameirrsec2)
    # read(filenameirrsec2)
    write(f, irrsec2)

    write(f, "\n\n")

    # write the name of the function
    write(f, "function invariants_gen(x1::SVector{$NBlengths, T}) where {T}\n")

    # write the precomputed powers of x
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   x$i = x$im.*x1 \n")
    end

    # write the primary invariants
    write(f, "   #------------------------------------------------\n")
    write(f, "   # Primaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n")
    prim3 = open(io->read(io), filenameprim3)
    # read(filenameprim3)
    write(f, prim3)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec3 = open(io->read(io), filenameirrsec3)
    # read(filenameirrsec3)
    write(f, irrsec3)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec)
    # read(filenamesec)
    write(f, sec)

    #write the return part
    write(f, "\n\n")
    write(f, "return (@SVector [")
    for i=1:NB_prim
        write(f, prefprim, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NB_secondary
        write(f, prefsec, "$i,")
    end
    write(f, "])\n end")



    # -------------------------------------------
    #
    # Generate derivatives of the invariants
    #
    # -------------------------------------------
    write(f, "\n\n\n\n")
    write(f, "function invariants_d_gen(x1::SVector{$NBlengths, T}) where {T}\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   x$i = x$im.*x1 \n")
    end
    write(f, "\n   dx1 = @SVector ones($NBlengths)\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   dx$i = $i * x$im \n")
    end

    # write the primary invariants
    write(f, "   #------------------------------------------------\n")
    write(f, "   # Primaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n")
    prim4 = open(io->read(io), filenameprim4)
    # read(filenameprim4)
    write(f, prim4)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec3 = open(io->read(io), filenameirrsec3)
    # read(filenameirrsec3)
    write(f, irrsec3)

    write(f, "\n\n")
    irrsec4 = open(io->read(io), filenameirrsec4)
    # read(filenameirrsec4)
    write(f, irrsec4)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec_d)
    # read(filenamesec_d)
    write(f, sec)

    #write the return part
    write(f, "\n\n")
    write(f, "return (")
    for i=1:NB_prim
        write(f, "d", prefprim, "$i,")
    end
    write(f, "), (")
    for i=1:NB_secondary
        write(f, "d", prefsec, "$i,")
    end
    write(f, ")\n end")

# -------------------------------------------
#
# Generate both invariants and derivatives of the invariants
#
# -------------------------------------------
    write(f, "\n\n\n\n")
    write(f, "function invariants_ed_gen(x1::SVector{$NBlengths, T}) where {T}\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   x$i = x$im.*x1 \n")
    end
    write(f, "\n   dx1 = @SVector ones($NBlengths)\n")
    for i=2:max(max_exp_irrsec,max_exp_prim)
        im = i-1;
        write(f, "   dx$i = $i * x$im \n")
    end

    # write the primary invariants
    write(f, "   #------------------------------------------------\n")
    write(f, "   # Primaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n")
    prim5 = open(io->read(io), filenameprim5)
    # read(filenameprim5)
    write(f, prim5)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec5 = open(io->read(io), filenameirrsec5)
    # read(filenameirrsec5)
    write(f, irrsec5)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec)
    # read(filenamesec)
    write(f, sec)

    write(f, "\n\n")
    sec = open(io->read(io), filenamesec_d)
    # read(filenamesec_d)
    write(f, sec)

    #write the return part
    write(f, "\n\n")
    write(f, "return (@SVector [")
    for i=1:NB_prim
        write(f, prefprim, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NB_secondary
        write(f, prefsec, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NB_prim
        write(f, "d", prefprim, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NB_secondary
        write(f, "d", prefsec, "$i,")
    end
    write(f, "])\n end \n\n")
    write(f, "end")
end



#Remove the temporary files
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text1.jl");
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text2.jl");
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text3.jl");
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text4.jl");

rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text1.jl");
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text2.jl");
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text3.jl");
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text4.jl");
rm(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_text_deg.jl");
