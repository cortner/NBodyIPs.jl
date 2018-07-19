using Combinatorics, StaticArrays
# , NBodyIPs
include("../../src/misc.jl")
include("invariants_generator.jl")
# include("../../src/polynomials.jl")

import Base: length, *

# Tuple : like SVector, but not a vector i.e. cant +, *, etc
# know length, immutable : SVector
# know length, mutable : MVector
# don't know length : Vector

abstract type AbstractMonList{M} end
# abstract type CompactMonomial{M} <: AbstractMonList{M} end
abstract type PolyMonomial{M} end

# struct Monomial{M} <: CompactMonomial{M}
#     a::SVector{M, Int}
# end

# struct CMonomial{M} <: CompactMonomial{M}
#     a::SVector{M, Int}
# end

struct MonList{M} <: AbstractMonList{M}
    A::Vector{SVector{M, Int}}
end

struct CMonList{M} <: AbstractMonList{M}
    A::Vector{SVector{M, Int}}
end
Mon(m::AbstractMonList{M}) where {M} = m.A

struct PolyMon{M} <: PolyMonomial{M}
   mon::MonList{M} #monomials
   coef::Vector{Int} #coefficients
end

Mon(m::PolyMon{M}) where {M} = (m.mon).A
Coef(m::PolyMon{M}) where {M} = m.coef

struct CPolyMon{M} <: PolyMonomial{M}
   mon::CMonList{M} #monomials
   coef::Vector{Int} #coefficients
end
Mon(m::CPolyMon{M}) where {M} = (m.mon).A
Coef(m::CPolyMon{M}) where {M} = m.coef


dim(m::AbstractMonList{M}) where {M} = M


# lengths redifinitions
length(m::MonList) = length(m.A)
length(m::CMonList) = length(m.A)

length(m::PolyMon) = length(m.coef)
length(m::CPolyMon) = length(m.coef)


edges2bo(M::Integer) = (M <= 0) ? 1 : round(Int, 0.5 + sqrt(0.25 + 2 * M))

# defining constant polynomial
constpoly(M::Integer) = CPolyMon(CMonList([@SVector zeros(Int,M)]),[1])

# Determine the representant of a monomial
function mon_repr(monomial::SVector{M, Int}) where {M}
    perms_unique = unique(simplex_permutations(monomial))
    perms_unique_sort = sort(perms_unique, lt=lexless, rev=true)
    return SVector(perms_unique_sort[1]...)
end


# Check duplicates in a list of monomials and remove them
function remove_dupl(monomiallist::MonList{M}) where {M}
    mon_list = Mon(monomiallist)
    mon_list_out = SVector{M, Int}[]

    for i=1:length(mon_list)
        if !(mon_list[i] in mon_list_out)
            push!(mon_list_out,mon_list[i])
        end
    end
    return MonList(mon_list_out)
end


# Remove duplicates in a list of monomials and adds the coefficients when duplicates
function remove_dupl(monomiallist::PolyMon{M}) where {M}
    mon_list = Mon(monomiallist)
    coef_list = Coef(monomiallist)
    @assert length(mon_list) == length(coef_list)
    mon_list_out = SVector{M, Int}[]
    mon_coef_out = Int[]

    for i=1:length(mon_list)
        if !(mon_list[i] in mon_list_out)
            push!(mon_list_out,mon_list[i])
            push!(mon_coef_out,coef_list[i])
        else
            ind = find([mon_list[i] == mon_list_out[j] for j=1:length(mon_list_out)])
            mon_coef_out[ind[1]] += coef_list[i]
        end
    end
    return PolyMon(MonList(mon_list_out),mon_coef_out)
end


# Compute the compact form of a list of monomials
function compact(monomiallist::MonList{M}) where {M}
    # first simplification - remove duplicate monomials
    monomiallist_in = remove_dupl(monomiallist)
    mon_list = Mon(monomiallist_in)
    mon_list_out = SVector{M, Int}[]
    for i=1:length(mon_list)
        if !(mon_repr(mon_list[i]) in mon_list_out)
            push!(mon_list_out,mon_repr(mon_list[i]))
        end
    end
    return CMonList(mon_list_out)
end

# Compute the compact form of a list of monomials
function compact(monomiallist::PolyMon{M}) where {M}
    # first simplification - remove duplicate monomials
    monomiallist_in = remove_dupl(monomiallist)
    mon_list = Mon(monomiallist_in)
    coef_list = Coef(monomiallist_in)
    @assert length(mon_list) == length(coef_list)

    mon_list_out = SVector{M, Int}[]
    coef_list_out = Int[]

    for i=1:length(mon_list)
        if !(mon_repr(mon_list[i]) in mon_list_out)
            push!(mon_list_out,mon_repr(mon_list[i]))
            push!(coef_list_out,coef_list[i])
        else
            # ind = find([mon_repr(mon_list[i]) == mon_repr(mon_list_out[j]) for j=1:length(mon_list_out)])
            # @assert coef_list[i] == coef_list_out[ind[1]]
        end
    end
    return CPolyMon(CMonList(mon_list_out),coef_list_out)
end


# Compute the expanded version of a compact list of monomials or polymon
function expand_mon(monomiallist::CPolyMon{M}) where {M}
    mon_list = Mon(monomiallist)
    coef_list = Coef(monomiallist)
    @assert length(mon_list) == length(coef_list)

    mon_list_out = SVector{M, Int}[]
    coef_list_out = Int[]

    for i=1:length(mon_list)
        USP = unique(simplex_permutations(mon_list[i]))
        append!(mon_list_out,USP)
        append!(coef_list_out,collect(Iterators.repeated(coef_list[i],length(USP))))
    end
    return PolyMon(MonList(mon_list_out),coef_list_out)
end


# Compute the expanded version of a compact list of monomials or polymon
function deg_monlist(monomiallist::CPolyMon{M}) where {M}
    mon = Mon(monomiallist)[1]
    return sum(mon)
end


function expand_mon(monomiallist::CMonList{M}) where {M}
    mon_list = Mon(monomiallist)
    mon_list_out = SVector{M, Int}[]

    for i=1:length(mon_list)
        USP = unique(simplex_permutations(mon_list[i]))
        append!(mon_list_out,USP)
    end
    return MonList(mon_list_out)
end


# compute the product of 2 lists of monomials
function *(PMon1::CPolyMon{M},PMon2::CPolyMon{M}) where {M}
    # Exp_PMon1 = expand_mon(PMon1)
    Exp_PMon2 = expand_mon(PMon2)
    mon1_list = Mon(PMon1)
    coef1_list = Coef(PMon1)
    mon2_list = Mon(Exp_PMon2)
    coef2_list = Coef(Exp_PMon2)

    mon_list_out = SVector{M, Int}[]
    coef_list_out = Int[]
    for i=1:length(mon1_list)
        for j=1:length(mon2_list)
            push!(mon_list_out,mon1_list[i]+mon2_list[j])
            push!(coef_list_out,coef1_list[i]*coef2_list[j])
        end
    end
    return compact(PolyMon(MonList(mon_list_out),coef_list_out))
end

# # compute the product of 2 lists of monomials
# function *(PMon1::CPolyMon{M},PMon2::CPolyMon{M}) where {M}
#     Exp_PMon1 = expand_mon(PMon1)
#     Exp_PMon2 = expand_mon(PMon2)
#     mon1_list = Mon(Exp_PMon1)
#     coef1_list = Coef(Exp_PMon1)
#     mon2_list = Mon(Exp_PMon2)
#     coef2_list = Coef(Exp_PMon2)
#
#     mon_list_out = SVector{M, Int}[]
#     coef_list_out = Int[]
#     for i=1:length(mon1_list)
#         for j=1:length(mon2_list)
#             push!(mon_list_out,mon1_list[i]+mon2_list[j])
#             push!(coef_list_out,coef1_list[i]*coef2_list[j])
#         end
#     end
#     return compact(PolyMon(MonList(mon_list_out),coef_list_out))
# end

function *(Mon1::PolyMon{M},PMon2::CPolyMon{M}) where {M}
    PMon1 = compact(Mon1)
    return *(PMon1,PMon2)
end


function *(PMon1::CPolyMon{M},Mon2::PolyMon{M}) where {M}
    return *(Mon2,PMon1)
end


function *(Mon1::PolyMon{M},Mon2::PolyMon{M}) where {M}
    PMon1 = compact(Mon1)
    PMon2 = compact(Mon2)
    return *(PMon1,PMon2)
end



# compute the power of a function given by monomials
function power(PMon::CPolyMon{M},n::Int) where {M}
    if n == 0
        return constpoly(M)
    elseif n == 1
        return PMon
    else
        return PMon*(power(PMon,n-1))
    end
end

power(Mon::PolyMon{M},n::Int) where {M} = power(compact(Mon),n)



multisets(k, n) = map(A -> [sum(A .== i) for i in 1:n],
                      with_replacement_combinations(1:n, k))

# generate all monomial representants up to degree d
function generate_rep_mon(M,d)
    mon_list_out = SVector{M, Int}[]
    for deg = 1:d
        mon_list_out_deg = SVector{M, Int}[]
        Decomp_deg = multisets(deg,M)
        for i=1:length(Decomp_deg)
            monrep = mon_repr(SVector(Decomp_deg[i]...))
            if !(monrep in mon_list_out_deg)
                push!(mon_list_out_deg,monrep)
            end
        end
        append!(mon_list_out,mon_list_out_deg)
    end
    return MonList(mon_list_out)
end


# # Small tests
# M = 6
# mon = SVector(1,2,0,0,0,0)
# length(mon)
# mon_repr(mon)
#
# listmon = MonList([mon, mon])
# length(listmon)
# remove_dupl(listmon)
# clistmon = compact(listmon)
# expand_mon(clistmon)
#
# Monopoly = PolyMon(MonList([mon, mon]),Vector([1, 2]))
# length(Monopoly)
# remove_dupl(Monopoly)
# cMonopoly = compact(Monopoly)
# expand_mon(cMonopoly)
# cMonopoly*cMonopoly
# cMonopoly*Monopoly
# Monopoly*cMonopoly
# Monopoly*Monopoly
#
# monzero = SVector(0,0,0,0,0,0)
# mon_repr(monzero)
# M0P = CPolyMon(CMonList([monzero, monzero]),Vector([1, 2]))
# M0P*M0P
#
# power(cMonopoly,3)
