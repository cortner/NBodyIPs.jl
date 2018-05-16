using Combinatorics, StaticArrays

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
filenameirrsec1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text1.jl";
filenameirrsec2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text2.jl";
filenameirrsec3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text3.jl";
filenameirrsecdata = "magma/data/NB_$NBody""_deg_$Deg""_irr_invariants.jl";
preword = "# Irreducible secondaries for NBody=$NBody"*"and deg=$Deg \n"
prefirrsec = "IS"
NB_irrsec = countlines(filenameirrsecdata)

max_exp_irrsec = generate_invariants(filenameirrsecdata,filenameirrsec1,filenameirrsec2,filenameirrsec3,NBlengths,Deg,preword,prefirrsec)


# -------------------------------------------
#
# Generate primary invariants
#
# -------------------------------------------
filenameprim1 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text1.jl";
filenameprim2 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text2.jl";
filenameprim3 = "magma/data/NB_$NBody"*"_deg_$Deg"*"_prim_text3.jl";
filenameprimdata = "magma/data/NB_$NBody""_deg_$Deg""_prim_invariants.jl";
preword = "# Primary invariants for NBody=$NBody"*"and deg=$Deg \n"
prefprim = "P"
NB_prim = countlines(filenameprimdata)

max_exp_prim = generate_invariants(filenameprimdata,filenameprim1,filenameprim2,filenameprim3,NBlengths,Deg,preword,prefprim)


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

    #write the return part
    write(f, "\n\n")
    write(f, "return (@SVector [")
    for i=1:NB_prim
        write(f, prefprim, "$i,")
    end
    write(f, "]), (@SVector [")
    for i=1:NB_irrsec
        write(f, prefirrsec, "$i,")
    end
    write(f, "])\n end")


end
