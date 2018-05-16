

fileoutput = "magma/data/test.jl";

fileinput = "magma/data/NB_5_deg_8_relations_invariants.jl";


open(fileoutput, "w") do f

end

fileI = open(fileinput)
line = readlines(fileI)
for i=1:length(line)
    if contains(line[i], "*")
        part1,part2 = split(line[i], "=")
        part2_1,part2_2 = split(part2, "*")
        repl1 = replace(part1, "SEC", "dSEC")
        repl2_1 = replace(part2_1, "IS", "dIS")
        repl2_2 = replace(part2_2, "IS", "dIS")
        open(fileoutput, "a") do f
            write(f, repl1, " = ", repl2_1, "*", part2_2, "+", part2_1, "*", repl2_2, "\n")
        end
    else
        open(fileoutput, "a") do f
            repl1 = replace(line[i], "SEC", "dSEC")
            repl2 = replace(repl1, "IS", "dIS")
            write(f, repl2, "\n")
        end
    end
end
