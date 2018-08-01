# Functions used to generate the function for the invariants and their derivatives



function generate_filename(filename)
   return filename*"1.jl",filename*"2.jl",filename*"3.jl",filename*"4.jl",filename*"5.jl"
end

function perm_2_indice(Perms)
   # convert a list of tuples (monomials) into an array of indices with non-zero entries in reverse order. All the tuples should have the same number of nonzero entries.
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



function vec_exp_2_file(filename,exponent,Vec_ind,prefix,number)
   filename1,filename2,filename3,filename4,filename5 = generate_filename(filename)
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


function monomial_2_file(filename,monomial,prefix,number)
   # writing in the different files for one given monomial
   exponent, Vec_ind = monomial_2_vec_exp(monomial)
   vec_exp_2_file(filename,exponent,Vec_ind,prefix,number)
   return exponent
end


function generate_invariants(filenamedata,filename,preword,prefix,Monomials)
   # writing in the different files for a vector of monomials
   filename1,filename2,filename3,filename4,filename5 = generate_filename(filename)

   NB_inv = length(Monomials)
    max_exp = 1;

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
       exponent = monomial_2_file(filename,monomial,prefix,j)
       max_exp_temp = maximum(exponent)
       if max_exp_temp > max_exp
          max_exp = max_exp_temp
       end
       # TODO: keep track of matrices
    end

    return max_exp


end
