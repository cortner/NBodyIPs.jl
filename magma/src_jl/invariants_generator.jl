using Combinatorics, StaticArrays
include("misc.jl")


function generate_monomials(filename,NBlengths,Deg=10)
    # Deg is the max. degree
    NB_inv = countlines(filename); #nb of invariants

    Monomials = []; #indices with exponents
    Monomials_simple = []; #indices without exponents (all 1)

    file = open(filename)
    line = readlines(file)

    for i=1:length(line)
        Mono_temp = 0*Vector{Int64}(NBlengths)
        Mono_sim_temp = 0*Vector{Int64}(NBlengths)
        for j=1:NBlengths
            if contains(line[i], "x[$j]")
                if contains(line[i], "x[$j]^")
                    for k=1:Deg
                        if contains(line[i], "x[$j]^$k")
                            Mono_temp[j] = k
                            Mono_sim_temp[j] = 1
                        end
                    end
                else
                    Mono_temp[j] = 1
                    Mono_sim_temp[j] = 1
                end
            end
        end
        push!(Monomials,Mono_temp);
        push!(Monomials_simple,Mono_sim_temp);
    end
    return NB_inv,Monomials,Monomials_simple
end

function perm_2_indice(Perms)
   L = length(Perms);
   deg = length(find(Perms[1]));
   M = length(Perms[1]);

   Indices = Array{Int64, 2}(L,deg);
   # zeros(L,deg);

      for j=1:L
         Ind_temp = sortperm(Perms[j],rev=true)
         Indices[j,:] = Ind_temp[1:deg]
      end
   return Indices
end


perm_2_indice([SVector(1,2,3,0,0,0)])

function monomial_2_vec_exp(monomial)
    perm_deg = length(find(monomial))
    Vec_ind = perm_2_indice(unique(simplex_permutations(monomial)))
    exponent = sort(monomial,rev=true)[1:perm_deg]
    return exponent, Vec_ind
end

monomial_2_vec_exp(SVector(4,1,1,0,0,0),SVector(1,1,1,0,0,0))


function vec_exp_2_file(filename1,filename2,filename3,exponent,Vec_ind,prefix,number)
   nb_vec = size(Vec_ind,2);
    # filename1: write the definition of the vectors
    # filename2: write the definition of the types containing the vectors
    # filename3: write what goes inside the function
    open(filename1, "a") do f
       for s=1:nb_vec
          write(f, "const ")
          write(f, prefix, "$number", "_$s", " = (" )
          for j = 1:size(Vec_ind,1)
              k = Vec_ind[j,s]
              write(f, "$k,")
          end
          write(f, ") \n")
       end
    end
    open(filename2, "a") do f
       write(f, "const ")
       write(f, prefix, "$number", " = Val((")
       for s=1:nb_vec
          write(f, prefix, "$number", "_$s,")
       end
       write(f, ")) \n")
    end
    open(filename3, "a") do f
       write(f, prefix,  "$number", " = fpoly((")
       for i=1:length(exponent)
          expi = exponent[i];
          write(f, "x$expi,")
       end
       write(f, ") , Main.", prefix, "$number", ") \n")
    end
end

22

function monomial_2_file(filename1,filename2,filename3,monomial,prefix,number)
   exponent, Vec_ind = monomial_2_vec_exp(monomial)
   vec_exp_2_file(filename1,filename2,filename3,exponent,Vec_ind,prefix,number)
end

monomial_2_file("test_mono1","test_mono2","test_mono3",SVector(1,0,0,0,0,0),"IS",3)

22



function generate_invariants(filenamedata,filename,NBlengths,Deg,preword,prefix)
    (NB_inv,Monomials,Monomials_simple) = generate_monomials(filenamedata,NBlengths,Deg)
    filename1 = filename*"1"
    filename2 = filename*"2"
    filename3 = filename*"3"
    open(filename, "w") do f
       write(f, preword )
    end
    open(filename1, "w") do f
      # resets the file
    end
    open(filename2, "w") do f
      # resets the file
    end
    open(filename3, "w") do f
      # resets the file
    end
    for j=1:NB_inv
       monomial = SVector(Monomials[j]...);
       monomial_2_file(filename1,filename2,filename3,monomial,prefix,j)
       # TODO: keep track of matrices
    end




end

# generate_invariants(filenamedata,filename,NBlengths,Deg,preword,prefix)

22


# function generate_invariants(filenamedata,filename,NBlengths,Deg,preword,preff)
#     (NB_inv,Monomials,Monomials_simple) = generate_monomials(filenamedata,NBlengths,Deg)
#     open(filename, "w") do f
#        write(f, preword )
#     end
#     for j=1:NB_inv
#        pref = preff*"$j";
#        Perms = SVector(Monomials[j]...)
#        Permsref = SVector(Monomials_simple[j]...)
#        perm_deg = sum(Permsref)
#        nb_diff_coef = length(unique(Perms))-1
#        Perm_inv = sort(Perms,rev=true)
#        # TODO: keep track of matrices
#
#        # different case with respect to the number of variables in the monomial
#        if (perm_deg==1)
#            # only 1 variable -> sum of x to some power
#            open(filename, "a") do f
#               exponent = Perm_inv[1];
#               write(f, pref)
#               write(f, " = ")
#               write(f, "sum(x$exponent)\n\n")
#            end
#
#        elseif (perm_deg==2)
#            A = _generate_matrix_Deg_2(NBlengths,Perms)
#            exp1 = Perm_inv[1];
#            exp2 = Perm_inv[2];
#            open(filename, "a") do f
#               write(f, preff,"A$j=")
#           end
#               printmatrix(filename,A)
#           open(filename, "a") do f
#               write(f,"\n")
#               write(f, pref)
#               write(f, " = ")
#               write(f, "x$exp1","'*",preff,"A$j","*"," x$exp2")
#               write(f, "\n\n")
#            end
#
#        else # (perm_deg>2)
#            for s=1:perm_deg
#                Ind = perm_2_indice(unique(simplex_permutations(Perms)))
#                open(filename, "a") do f
#                    # First vector
#                    write(f, pref, "_$s = ")
#                    write(f, "@SVector [")
#                    for j = 1:size(Ind,1)
#                        k = Ind[j,s]
#                        # write(f, "x$exp[$k],")
#                        write(f, "$k,")
#                    end
#                    write(f, "] \n")
#                end
#            end
#            if (perm_deg == 3)
#                exp1 = Perm_inv[1]
#                exp2 = Perm_inv[2]
#                exp3 = Perm_inv[3]
#                open(filename, "a") do f
#                    write(f, pref, " = ")
#                    write(f, "dot(x$exp1[", pref, "_1].* x$exp2[", pref, "_2],x$exp3[", pref, "_3])\n \n")
#                end
#            elseif (perm_deg == 4)
#                exp1 = Perm_inv[1]
#                exp2 = Perm_inv[2]
#                exp3 = Perm_inv[3]
#                exp4 = Perm_inv[4]
#                open(filename, "a") do f
#                    write(f, pref, " = ")
#                    write(f, "dot(x$exp1[", pref, "_1].* x$exp2[", pref, "_2],x$exp3[", pref, "_3].* x$exp4[", pref, "_4])\n \n")
#                end
#            elseif (perm_deg == 5)
#                exp1 = Perm_inv[1]
#                exp2 = Perm_inv[2]
#                exp3 = Perm_inv[3]
#                exp4 = Perm_inv[4]
#                exp5 = Perm_inv[5]
#                open(filename, "a") do f
#                    write(f, pref, " = ")
#                    write(f, "dot(x$exp1[", pref, "_1].* x$exp2[", pref, "_2],x$exp3[", pref, "_3].* x$exp4[", pref, "_4].* x$exp5[", pref ,"_5])\n \n")
#                end
#            else
#                error("not implemented (perm_deg larger than 5)")
#            end
#        end
#     end
#     return 0
#
# end

# function printmatrix(filename,A)
#     M = size(A,1);
#     N = size(A,2);
#     open(filename, "a") do f
#        write(f, "[")
#        for m=1:(M-1)
#            for n=1:N
#                element = A[m,n]
#                 write(f, "$element")
#                 write(f, " ")
#             end
#             write(f, "; ")
#         end
#             for n=1:N
#             element = A[M,n]
#              write(f, "$element")
#              write(f, " ")
#          end
#         write(f, "]\n")
#     end
# end

# function _generate_matrix_Deg_2(NBlengths,Perms)
#     Ind = perm_2_indice(unique(simplex_permutations(Perms)))
#     A = zeros(Int64, NBlengths,NBlengths);
#     for i=1:size(Ind, 1)
#         A[Ind[i,1],Ind[i,2]] = 1;
#     end
#     return A
# end
