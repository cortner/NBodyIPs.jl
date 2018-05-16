using Combinatorics, StaticArrays

include("invariants_generator.jl")
include("misc.jl")

# Parameters
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
# filename = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text.jl";
filenameirrsec1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text1.jl";
filenameirrsec2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text2.jl";
filenameirrsec3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text3.jl";
filenameirrsec4 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text4.jl";
filenameirrsecdata = "magma/data/NB_$NBody""_deg_$Deg""_irr_invariants.jl";
preword = "# Irreducible secondaries for NBody=$NBody"*"and deg=$Deg \n"

NB_irrsec = countlines(filenameirrsecdata)

max_exp_irrsec = generate_invariants(filenameirrsecdata,filenameirrsec1,filenameirrsec2,filenameirrsec3,filenameirrsec4,NBlengths,Deg,preword,prefirrsec)


# -------------------------------------------
#
# Generate primary invariants
#
# -------------------------------------------
filenameprim1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text1.jl";
filenameprim2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text2.jl";
filenameprim3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text3.jl";
filenameprim4 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text4.jl";
filenameprimdata = "magma/data/NB_$NBody""_deg_$Deg""_prim_invariants.jl";
preword = "# Primary invariants for NBody=$NBody"*"and deg=$Deg \n"
NB_prim = countlines(filenameprimdata)

max_exp_prim = generate_invariants(filenameprimdata,filenameprim1,filenameprim2,filenameprim3,filenameprim4,NBlengths,Deg,preword,prefprim)
# -------------------------------------------
#
# Secondary invariants (relations with irreducible secondaries)
#
# -------------------------------------------
filenamesec = "magma/data/NB_$NBody"*"_deg_$Deg"*"_relations_invariants.jl";
NB_secondary = countlines(filenamesec);
# -------------------------------------------
#
# Generate function with all invariants
#
# -------------------------------------------
file = "magma/data/NB_$NBody"*"_deg_$Deg"*"_invariants.jl";

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
    write(f, "function invariants(x1::SVector{$NBlengths, T}) where {T}\n")

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
    write(f, "function invariants_d(x1::SVector{$NBlengths, T}) where {T}\n")
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
    irrsec4 = read(filenameirrsec4)
    write(f, irrsec4)

    #write the return part
    write(f, "\n\n")
    write(f, "return (")
    for i=1:NB_prim
        write(f, "d", prefprim, "$i,")
    end
    write(f, "), (")
    for i=1:NB_irrsec
        write(f, "d", prefirrsec, "$i,")
    end
    write(f, ")\n end")

end
