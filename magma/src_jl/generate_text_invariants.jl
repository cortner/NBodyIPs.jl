using Combinatorics, StaticArrays

include("invariants_generator.jl")
include("misc.jl")

# Parameters
#TODO: check that prefix are the same as the ones in generate_invariants.sh
NBody = 5;
Deg = 6;
prefsec = "SEC" #prefix for the secondaries
prefirrsec = "IS" #prefix for the irreducible secondaries
prefprim = "P" #prefix for the primaries
# --------------
NBlengths = Int(NBody*(NBody-1)/2);

# -------------------------------------------
#
# Generate irreducible secondaries
#
# -------------------------------------------
filenameirrsec1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text1.jl";
filenameirrsec2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text2.jl";
filenameirrsec3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text3.jl";
filenameirrsec4 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text4.jl";
filenameirrsec5 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text5.jl";
filenameirrsecdata = "magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_irr_invariants.jl";
preword = "# Irreducible secondaries for NBody=$NBody"*"and deg=$Deg \n"

NB_irrsec = countlines(filenameirrsecdata)

max_exp_irrsec = generate_invariants(filenameirrsecdata,filenameirrsec1,filenameirrsec2,filenameirrsec3,filenameirrsec4,filenameirrsec5,NBlengths,Deg,preword,prefirrsec)


# -------------------------------------------
#
# Generate primary invariants
#
# -------------------------------------------
filenameprim1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text1.jl";
filenameprim2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text2.jl";
filenameprim3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text3.jl";
filenameprim4 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text4.jl";
filenameprim5 = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text5.jl";
filenameprimdata = "magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_prim_invariants.jl";
preword = "# Primary invariants for NBody=$NBody"*"and deg=$Deg \n"
NB_prim = countlines(filenameprimdata)

max_exp_prim = generate_invariants(filenameprimdata,filenameprim1,filenameprim2,filenameprim3,filenameprim4,filenameprim5,NBlengths,Deg,preword,prefprim)
# -------------------------------------------
#
# Secondary invariants (relations with irreducible secondaries)
#
# -------------------------------------------
filenamesec = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_relations_invariants.jl";
NB_secondary = countlines(filenamesec);
# -------------------------------------------
#
# Derivatives of secondary invariants (relations with irreducible secondaries)
#
# -------------------------------------------
filenamesec_d = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_relations_invariants_derivatives.jl";

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
# Generate function with all invariants
#
# -------------------------------------------
file = "magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_invariants.jl";

open(file, "w") do f
    write(f, "using StaticArrays \n")
    write(f, "using BenchmarkTools: @btime \n\n")
    # write(f, "include(\"fastpolys.jl\") \n")
    # write(f, "using FastPolys \n\n\n\n")

    # write the definition of the constant vectors
    prim1 = read(filenameprim1)
    write(f, prim1)
    irrsec1 = read(filenameirrsec1)
    write(f, irrsec1)

    # write the definitions of the types
    write(f, "\n")
    prim2 = read(filenameprim2)
    write(f, prim2)
    irrsec2 = read(filenameirrsec2)
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
    prim3 = read(filenameprim3)
    write(f, prim3)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec3 = read(filenameirrsec3)
    write(f, irrsec3)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = read(filenamesec)
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
    prim4 = read(filenameprim4)
    write(f, prim4)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec3 = read(filenameirrsec3)
    write(f, irrsec3)

    write(f, "\n\n")
    irrsec4 = read(filenameirrsec4)
    write(f, irrsec4)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = read(filenamesec_d)
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
    prim5 = read(filenameprim5)
    write(f, prim5)

    # write the irreducible secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # Irreducible secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    irrsec5 = read(filenameirrsec5)
    write(f, irrsec5)

    # write all the secondary invariants
    write(f, "\n\n\n   #------------------------------------------------\n")
    write(f, "   # All secondaries\n")
    write(f, "   #------------------------------------------------\n")

    write(f, "\n\n")
    sec = read(filenamesec)
    write(f, sec)

    write(f, "\n\n")
    sec = read(filenamesec_d)
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
    write(f, "])\n end")





#Remove the temporary files
rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text1.jl");
rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text2.jl");
rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text3.jl");
rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text4.jl");

rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text1.jl");
rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text2.jl");
rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text3.jl");
rm("magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_prim_text4.jl");

end
