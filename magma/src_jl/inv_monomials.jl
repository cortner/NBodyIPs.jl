using Combinatorics, StaticArrays, NBodyIPs
include("../../src/misc.jl")
include("invariants_generator.jl")
include("../../src/polynomials.jl")

import Base: length, +, *

# Tuple : like SVector, but not a vector i.e. cant +, *, etc
# know length, immutable : SVector
# know length, mutable : MVector
# don't know length : Vector

abstract type AbstractMonomial{M} end
# abstract type CompactMonomial{M} <: AbstractMonomial{M} end
abstract type PolyMonomial{M} end

# struct Monomial{M} <: CompactMonomial{M}
#     a::SVector{M, Int}
# end

# struct CMonomial{M} <: CompactMonomial{M}
#     a::SVector{M, Int}
# end

struct Monlist{M} <: AbstractMonomial{M}
    A::Vector{SVector{M, Int}}
end

struct CMonlist{M} <: AbstractMonomial{M}
    A::Vector{SVector{M, Int}}
end
Mon(m::AbstractMonomial{M}) where {M} = m.A

struct PolyMon{M} <: PolyMonomial{M}
   mon::Monlist{M} #monomials
   coef::Vector{Int} #coefficients
end

Mon(m::PolyMon{M}) where {M} = (m.mon).A
Coef(m::PolyMon{M}) where {M} = m.coef

struct CPolyMon{M} <: PolyMonomial{M}
   mon::CMonlist{M} #monomials
   coef::Vector{Int} #coefficients
end
Mon(m::CPolyMon{M}) where {M} = (m.mon).A
Coef(m::CPolyMon{M}) where {M} = m.coef


dim(m::AbstractMonomial{M}) where {M} = M


# lengths redifinitions
length(m::Monlist) = length(m.A)
length(m::CMonlist) = length(m.A)

length(m::PolyMon) = length(m.coef)
length(m::CPolyMon) = length(m.coef)


edges2bo(M::Integer) = (M <= 0) ? 1 : round(Int, 0.5 + sqrt(0.25 + 2 * M))

# simplex_permutations(x::Monomial{M}) =
#    simplex_permutations(Val(edges2bo(M)), x.a)


# Determine the representant of a monomial
function mon_repr(monomial::SVector{M, Int}) where {M}
    perms_unique = unique(simplex_permutations(monomial))
    perms_unique_sort = sort(perms_unique, lt=lexless, rev=true)
    return SVector(perms_unique_sort[1]...)
end


# Check duplicates in a list of monomials
function check_dupl(monomiallist::AbstractMonomial)
    mon_list = Mon(monomiallist)
    @show mon_list
    mon_list_out = Monlist([])

    for i=1:length(mon_list)
        if !(mon_list[i] in mon_list_out)
            @show mon_list[i]
            push!(mon_list_out,mon_list[i])
        end
    end
    @show mon_list_out
    return Monlist(mon_list_out)
end

# Check duplicates in a list of monomials and adds the coefficients when duplicates
function check_dupl_add(monomiallist::PolyMonomial)
    mon_list = Mon(monomiallist)
    coef_list = Coef(monomiallist)
    @assert length(mon_list) == length(coef_list)
    mon_list_out = []
    mon_coef_out = []

    for i=1:length(mon_list)
        @show mon_list
        if !(mon_list[i] in mon_list_out)
            @show mon_list[i]
            push!(mon_list_out,mon_list[i])
            push!(mon_coef_out,coef_list[i])
        else
            ind = find([mon_list[i] == mon_list_out[j] for j=1:length(mon_list_out)])
            mon_coef_out[ind[1]] += coef_list[i]
        end
    end
    return PolyMon(Monlist(mon_list_out),mon_coef_out)
end





# Compute the compact form of a list of monomials
function compact(monomiallist::PolyMonomial)
    # first simplification - remove duplicate monomials
    monomiallist_in = check_dupl_add(monomiallist)
    mon_list = Mon(monomiallist_in)
    coef_list = Coef(monomiallist_in)
    @assert length(mon_list) == length(coef_list)

    mon_list_out = []
    coef_list_out = []

    for i=1:length(mon_list)
        if !(monomial_repr(mon_list[i]) in mon_list_out)
            push!(mon_list_out,mon_repr(mon_list[i]))
            push!(coef_list_out,coef_list_in[i])
        else
            ind = find([mon_repr(mon_list[i]) == mon_repr(mon_list_out[j]) for j=1:length(mon_list_out)])
            @assert coef_list_in[i] == coef_list_out[ind[1]]
        end
    end
    return PolyMonomial(mon_list_out,coef_list_out)
end

Monlist([])

M = 6
mon = SVector(1,2,0,0,0,0)

length(mon)

length(Monlist([mon, mon]))

mon2 = Monlist([mon, mon])

length(mon2.A)

mon_repr(mon)

Mono = PolyMon(Monlist([mon, mon]),Vector([1, 2]))
length(Mono)

mon_repr(mon)

check_dupl_add(Mono::PolyMonomial)
check_dupl(Monlist([mon, mon]))

compact(Mono)




# Compute the compact form of a list of monomials
function compact_form_mon(mon_list)
    # first simplification - remove duplicate monomials
    mon_list_in = check_dupl(mon_list)
    mon_list_out = []

    for i=1:length(mon_list_in)
        if !(monomial_repr(mon_list_in[i]) in mon_list_out)
            push!(mon_list_out,monomial_repr(mon_list_in[i]))
        end
    end
    return mon_list_out
end

# Compute the expanded version of a compact list of monomials
function expanded_form_mon(mon_list,coef_list)
    @assert length(mon_list) == length(coef_list)
    # first simplification - remove duplicate monomials
    mon_list_in, coef_list_in = check_dupl_add(mon_list,coef_list)

    mon_list_out = []
    coef_list_out = []

    for i=1:length(mon_list)
        USP = unique(simplex_permutations(SVector(mon_list_in[i]...)))
        append!(mon_list_out,USP)
        # TODO: more efficient implementation???
        for j=1:length(USP)
            push!(coef_list_out,coef_list_in[i])
        end
    end
    return mon_list_out,coef_list_out
end

function prod_mon_comp(mon1_list,coef1_list,mon2_list,coef2_list)
    mon1_list_ex,coef1_list_ex = expanded_form_mon(mon1_list,coef1_list)
    mon2_list_ex,coef2_list_ex = expanded_form_mon(mon2_list,coef2_list)
    mon_list_out,coef_list_out = prod_mon(mon1_list_ex,coef1_list_ex,mon2_list_ex,coef2_list_ex)
    return compact_form_mon(mon_list_out,coef_list_out)
end

# compute the product of two list of monomials
function prod_mon(mon1_list,coef1_list,mon2_list,coef2_list)
    mon_list_out = []
    coef_list_out = []
    for i=1:length(mon1_list)
        for j=1:length(mon2_list)
            push!(mon_list_out,mon1_list[i]+mon2_list[j])
            push!(coef_list_out,coef1_list[i]*coef2_list[j])
        end
    end
    return check_dupl_add(mon_list_out,coef_list_out)
end

# compute the power of a function given by monomials
# function power(mon_list,coef_list,::Val{1})
#     return check_dupl_add(mon_list,coef_list)
# end
#
# function power(mon_list,coef_list,::Val{2})
#     return prod_mon(mon_list,coef_list,mon_list,coef_list)
# end

# compute the power of a function given by monomials (not compacted)
function power(mon_list,coef_list,p)
    if p == 1
         return check_dupl_add(mon_list,coef_list)
    elseif p == 2
         mon_list_out,coef_list_out = prod_mon(mon_list,coef_list,mon_list,coef_list)
         return check_dupl_add(mon_list_out,coef_list_out)
    else
        mon_list_temp, coef_list_temp = power(mon_list,coef_list,p-1)
        mon_list_out,coef_list_out = prod_mon(mon_list,coef_list,mon_list_temp, coef_list_temp)
        return check_dupl_add(mon_list_out,coef_list_out)
    end
end


multisets(k, n) = map(A -> [sum(A .== i) for i in 1:n],
                      with_replacement_combinations(1:n, k))

# generate all monomial representants up to degree d
function generate_rep_mon(NBlengths,d)
    mon_list_out = []
    for deg = 1:d
        mon_list_out_deg = []
        Decomp_deg = multisets(deg,NBlengths)
        mon_list_in_deg = compact_form_mon(Decomp_deg)
        for i=1:length(mon_list_in_deg)
            if !(monomial_repr(mon_list_in_deg[i]) in mon_list_out_deg)
                push!(mon_list_out_deg,monomial_repr(mon_list_in_deg[i]))
            end
        end
        append!(mon_list_out,mon_list_out_deg)
    end
    return mon_list_out
end


# function invariant_2_monomials(inv_tuple,Invariants)
#
#
# end
