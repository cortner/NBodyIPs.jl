using Combinatorics, StaticArrays, NBodyIPs
include("misc.jl")
include("invariants_generator.jl")
include("../../src/polynomials.jl")


# Determine the representant of a monomial
function monomial_repr(monomial)
    perms_unique = unique(simplex_permutations(SVector(monomial...)))
    perms_unique_sort = sort(perms_unique, lt=lexless, rev=true)
    return perms_unique_sort[1]
end


# Check duplicates in a list of monomials and adds the coefficients when duplicates
function check_dupl_add(mon_list,coef_list)
    @assert length(mon_list) == length(coef_list)
    mon_list_out = []
    mon_coef_out = []

    for i=1:length(mon_list)
        if !(mon_list[i] in mon_list_out)
            push!(mon_list_out,mon_list[i])
            push!(mon_coef_out,coef_list[i])
        else
            ind = find([mon_list[i] == mon_list_out[j] for j=1:length(mon_list_out)])
            mon_coef_out[ind[1]] += coef_list[i]
        end
    end
    return mon_list_out,mon_coef_out
end

# Check duplicates in a list of monomials
function check_dupl(mon_list)
    mon_list_out = []

    for i=1:length(mon_list)
        if !(mon_list[i] in mon_list_out)
            push!(mon_list_out,mon_list[i])
        end
    end
    return mon_list_out
end


# Compute the compact form of a list of monomials
function compact_form_mon(mon_list,coef_list)
    @assert length(mon_list) == length(coef_list)
    # first simplification - remove duplicate monomials
    mon_list_in, coef_list_in = check_dupl_add(mon_list,coef_list)

    mon_list_out = []
    coef_list_out = []

    for i=1:length(mon_list_in)
        if !(monomial_repr(mon_list_in[i]) in mon_list_out)
            push!(mon_list_out,monomial_repr(mon_list_in[i]))
            push!(coef_list_out,coef_list_in[i])
        else
            ind = find([monomial_repr(mon_list[i]) == monomial_repr(mon_list_out[j]) for j=1:length(mon_list_out)])
            @assert coef_list_in[i] == coef_list_out[ind[1]]
        end
    end
    return mon_list_out,coef_list_out
end

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


function invariant_2_monomials(inv_tuple,Invariants)
    

end

Mon = (0,0,1,0,0,0)




USP = unique(simplex_permutations(SVector(Mon...)))

check_dupl_add([USP; USP],[1,2,3,4,5,6,1,2,3,4,5,6])

compact_form_mon(USP)

prod_mon(USP,[1,1,1,1,1,1], USP,[1,1,1,1,1,1])


compact_form_mon([USP; USP],[1,1,1,1,1,1,2,2,2,2,2,2])

expanded_form_mon([Mon],[1])

NBodyIPs.gen_tuples(5,5)

generate_rep_mon(10,5)

mon_list,coef_list = power(USP,[1,1,1,1,1,1],3)

compact_form_mon(mon_list,coef_list)
