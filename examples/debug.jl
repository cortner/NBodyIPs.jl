
using JuLIP, NBodyIPs


B = NBodyIPs.polys_fourbody(1) |> display
B = NBodyIPs.polys_fourbody2(4) |> display

#some tests to count the number of monomials in polys_fourbody2
Deg_max = 2:8

Check_nb_basis_fct = zeros(length(Deg_max),2)

for (ideg,deg_max) in enumerate(Deg_max)
   theory_nb_basis_fct = 0; #do not count [0,0,0,0,0,0]
   for i=1:deg_max
      for m = 1:6
         for α in collect(partitions(i, m))
            # any terms not included get zeros appended
            append!(α, zeros(Int, 6 - length(α)))
            #display(unique(permutations(α)))
            theory_nb_basis_fct += length(unique(permutations(α)))
         end
      end
   end
   Check_nb_basis_fct[ideg,1] = theory_nb_basis_fct
   Check_nb_basis_fct[ideg,2] = sum(length(NBodyIPs.polys_fourbody2(deg_max)[i]) for i=1:length(NBodyIPs.polys_fourbody2(deg_max)))
end

#the numbers should correspond.
Check_nb_basis_fct |> display

deg_max = 2
length(NBodyIPs.polys_fourbody2(deg_max)) |>display
for i=1:length(NBodyIPs.polys_fourbody2(deg_max))
   display(NBodyIPs.polys_fourbody2(deg_max)[i][1]')
end

for i=1:length(NBodyIPs.polys_fourbody(deg_max))
   display(NBodyIPs.polys_fourbody(deg_max)[i][1]')
end
