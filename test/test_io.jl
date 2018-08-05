using Base.Test
using NBodyIPs
using NBodyIPs.Polys


println("Check (De-)Dictionisation of `BLDictionary`")
D = BLDictionary("r -> (2.0/r)^3", "(:cos, 6.0, 9.0)")
Ds = Dict(D)
D1 = BLDictionary(Ds)
println(@test D1 == D)


println("generate some basis functions")
rcuts = [9.2, 6.2, 4.5]   # 9.2
TRANSFORM = "r -> (2.9/r)^3"
CUTOFF = ["(:cos, $(0.66*rcut), $(rcut))" for rcut in rcuts]
B1 = [NBPoly(1.0)]
B2 = bl_basis(2, TRANSFORM, CUTOFF[1], 12)
B3 = bl_basis(3, TRANSFORM, CUTOFF[2], 10)
B4 = bl_basis(4, TRANSFORM, CUTOFF[3], 6)
B = [B1; B2; B3; B4]
c = rand(length(B))
IP = NBodyIP(B, c)
println("Check (De-)Dictionisation of `NBodyIP`")
D_IP = Dict(IP)
IP1 = NBodyIP(D_IP)
println(@test IP1 == IP)

# store it as a JSON file and load it again
fname = tempname() * ".json"
save_ip(fname, IP)
IP2 = load_ip(fname)
println(@test IP2 == IP)
run(`rm $fname`)
