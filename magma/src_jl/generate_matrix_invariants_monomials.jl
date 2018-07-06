# This code generates the matrix of basis change between the tuples in terms of the invariants and tuples of monomials.
using NBodyIPs
using JLD

include(homedir() * "/.julia/v0.6/NBodyIPs/magma/src_jl/invariants_generator.jl")
include(homedir() * "/.julia/v0.6/NBodyIPs/magma/src_jl/misc.jl")
include(homedir() * "/.julia/v0.6/NBodyIPs/magma/src_jl/inv_monomials.jl")

# Generate monomials with weights: for primaries, irreducible secondaries, and secondaries
NBody = 5;
Deg = 10;
NBlengths = Int(NBody*(NBody-1)/2)

filenameirrsecdata = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_irr_invariants.jl";
filenameprimdata = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_prim_invariants.jl";
filenamesec = homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody"*"_deg_$Deg"*"/NB_$NBody"*"_deg_$Deg"*"_relations_invariants.jl";

prefsec = "SEC" #prefix for the secondaries
prefirrsec = "IS" #prefix for the irreducible secondaries
prefprim = "P" #prefix for the primaries
# for irreducible secondaries


# Generate monomials with weights: for primaries, irreducible secondaries, and secondaries
# for irreducible secondaries
Mon_irrsec, coef_list_irrsec = generate_inv_mon(filenameirrsecdata,NBlengths,Deg)

# generate list of PolyMon for irreducible secondaries
IrrSecMonPol = CPolyMon{NBlengths}[]
for i=1:length(Mon_irrsec)
    push!(IrrSecMonPol,CPolyMon(CMonList([mon_repr(SVector(Mon_irrsec[i]...))]),[coef_list_irrsec[i]]))
end
IrrSecMonPol

# for primaries
Mon_prim, coef_list_prim = generate_inv_mon(filenameprimdata,NBlengths,Deg)

# generate list of PolyMon for primaries
PrimMonPol = CPolyMon{NBlengths}[]
for i=1:length(Mon_prim)
    push!(PrimMonPol,CPolyMon(CMonList([mon_repr(SVector(Mon_prim[i]...))]),[coef_list_prim[i]]))
end
PrimMonPol

# All secondary invariants
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
SecMonPol



M = Int(NBody*(NBody-1)/2)
InvTup = NBodyIPs.gen_tuples(NBody,Deg)

maxpower = zeros(Int,M,1)
powersPrimInv = []
for i=1:2
    maxpower[i] = maximum(InvTup[j][i] for j=1:length(InvTup))
    powersi = []
    push!(powersi,PrimMonPol[i])
    for j=1:maxpower[i]
        @show j
        # push!(powersi,power(PrimMonPol[i],j))
        push!(powersi,PrimMonPol[i]*powersi[j])
    end
    push!(powersPrimInv,powersi)
end
maxpower

InvMonPoly = CPolyMon{M}[]
for i=1:length(InvTup)
    # initialization
    PMonTup = SecMonPol[InvTup[i][end]+1]
    for j=1:M
        if InvTup[i][j] > 0
            PMonTup = PMonTup*powersPrimInv[j][InvTup[i][j]]
            # power(PrimMonPol[j],InvTup[i][j])
        end
    end
    push!(InvMonPoly,PMonTup)
    print("$i ")
end
InvMonPoly

MonBasis = Mon(generate_rep_mon(NBlengths,Deg))
@assert length(InvTup) == length(MonBasis)

M_basis_change = zeros(Float64,length(MonBasis),length(MonBasis))
# express all inv_tuples in terms of monomials

for i=1:length(InvTup)
    InviMon = Mon(InvMonPoly[i])
    InviCoef = Coef(InvMonPoly[i])
    for j=1:length(MonBasis)
        if (MonBasis[j] in InviMon)
            indj = find([MonBasis[j] == InviMon[k] for k=1:length(InviMon)])
            @assert length(indj) == 1
            M_basis_change[i,j] = InviCoef[indj[1]]
        end
    end
end
save(homedir() * "/.julia/v0.6/NBodyIPs/magma/data/NB_$NBody""_deg_$Deg""/NB_$NBody"*"_deg_$Deg"*"_basis_change.jld", "Mbasischange", M_basis_change, "NBody", NBody, "Deg", Deg)


# M_basis_change
#
# cond(M_basis_change)
#
# M_basis_change^(-1)
