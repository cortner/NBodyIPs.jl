using Combinatorics, StaticArrays
include("misc.jl")


# function generate_monomials(filename,NBvar)
#    #generate leading monomials as a tuple from a choosen file
#     NB_inv = countlines(filename); #nb of invariants
#
#     Monomials = []; #indices with exponents
#     Monomials_simple = []; #indices without exponents (all 1)
#     deg = [];
#
#     file = open(filename)
#     line = readlines(file)
#
#     for i=1:length(line)
#         Mono_temp = 0*Vector{Int64}(NBvar)
#         Mono_sim_temp = 0*Vector{Int64}(NBvar)
#         for j=1:NBvar
#             if contains(line[i], "x[$j]")
#                 if contains(line[i], "x[$j]^")
#                         Mono_sim_temp[j] = 1
#                         if contains(line[i], "x[$j]^$k")
#                             Mono_temp[j] = k
#
#                         end
#                     end
#                 else
#                     Mono_temp[j] = 1
#                     Mono_sim_temp[j] = 1
#                 end
#             end
#         end
#         push!(deg,sum(Mono_temp));
#         push!(Monomials,Mono_temp);
#         push!(Monomials_simple,Mono_sim_temp);
#     end
#     return NB_inv,Monomials,Monomials_simple,deg
# end


# function generate_inv_mon(filename,NBvar,Deg=10)
#     NB_inv,Monomials,Monomials_simple,deg = generate_monomials(filename,NBvar,Deg)
#
#     coef_list = []
#     for j=1:NB_inv
#         push!(coef_list,1)
#     end
#     return Monomials, coef_list, deg
# end

function perm_2_indice(Perms)
   # convert a list of tuples into an array of indices with non-zero entries in reverse order. All the tuples should have the same number of nonzero entries.
   L = length(Perms);
   deg = length(find(Perms[1]));
   M = length(Perms[1]);

   Indices = Array{Int64, 2}(L,deg);
   for j=1:L
      Ind_temp = sortperm(Perms[j],rev=true)
      Indices[j,:] = Ind_temp[1:deg]
   end
   return Indices
end


function monomial_2_vec_exp(monomial)
   #converts a tuple representing a monomial into an array containing the indices with non-zero entries in simplex_permutations of the monomial.
   # The exponent of each column is given by "exponent".
    perm_deg = length(find(monomial))
    Vec_ind = perm_2_indice(unique(simplex_permutations(monomial)))
    exponent = sort(monomial,rev=true)[1:perm_deg]
    return exponent, Vec_ind
end



function vec_exp_2_file(filename1,filename2,filename3,filename4,filename5,exponent,Vec_ind,prefix,number)
   # generate different files with parts of the invariants in them.
   nb_vec = size(Vec_ind,2);
    # filename1: write the definition of the vectors
    # filename2: write the definition of the types containing the vectors
    # filename3: write what goes inside the function
    # filename4: write what goes inside the function for the derivatives
    # filename5: write what goes inside the function for both evaluation and derivatives
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
       write(f, "\n")
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
       write(f, ") , NB5I.", prefix, "$number", ") \n")
    end
    open(filename4, "a") do f
       write(f, "d", prefix,  "$number", " = fpoly_d((")
       for i=1:length(exponent)
          expi = exponent[i];
          write(f, "x$expi,")
       end
       write(f, "),(")
       for i=1:length(exponent)
          expi = exponent[i];
          write(f, "dx$expi,")
       end
       write(f, ") , NB5I.", prefix, "$number", ") \n")
    end
    open(filename5, "a") do f
      write(f, prefix,  "$number", ", d", prefix,  "$number", " = fpoly_ed((")
      for i=1:length(exponent)
          expi = exponent[i];
          write(f, "x$expi,")
      end
      write(f, "),(")
      for i=1:length(exponent)
          expi = exponent[i];
          write(f, "dx$expi,")
      end
      write(f, ") , NB5I.", prefix, "$number", ") \n")
    end
end


function monomial_2_file(filename1,filename2,filename3,filename4,filename5,monomial,prefix,number)
   exponent, Vec_ind = monomial_2_vec_exp(monomial)
   vec_exp_2_file(filename1,filename2,filename3,filename4,filename5,exponent,Vec_ind,prefix,number)
   return exponent
end


function generate_invariants(filenamedata,filename1,filename2,filename3,filename4,filename5,NBvar,Deg,preword,prefix)
    max_exp = 1;
    (NB_inv,Monomials,Monomials_simple) = generate_monomials(filenamedata,NBvar,Deg)
    open(filename1, "w") do f
      write(f, preword, " # : definitions at the beginning of the file \n")
    end
    open(filename2, "w") do f
      write(f, preword, " # : definitions of the types at the beginning of the file \n")
    end
    open(filename3, "w") do f
      write(f, preword, " # : what goes in the function for the evaluation \n")
    end
    open(filename4, "w") do f
      write(f, preword, " # : what goes in the function for the derivatives \n")
    end
    open(filename5, "w") do f
      write(f, preword, " # : what goes in the function for the evaluation and derivatives \n")
    end
    for j=1:NB_inv
       monomial = SVector(Monomials[j]...);
       exponent = monomial_2_file(filename1,filename2,filename3,filename4,filename5,monomial,prefix,j)
       max_exp_temp = maximum(exponent)
       if max_exp_temp > max_exp
          max_exp = max_exp_temp
       end
       # TODO: keep track of matrices
    end

    return max_exp


end
