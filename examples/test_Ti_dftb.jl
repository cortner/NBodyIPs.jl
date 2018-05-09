
using NBodyIPs, JuLIP

datafile = "~/Dropbox/PIBmat/Ti_DFTB_Data/Ti_N54_vartemp_virials_ASE.xyz"
data = NBodyIPs.Data.read(datafile, index=":")
shuffle!(data)
train_data = data
test_data = data[1:10:end]

# some notes on the orders of magnitude
# E per atom ~ 6.0 eV, std ca 0.03
# mean force on each atom is ~ 1 eV / A
# but range of forces is between 1 and 12 eV / A

# construct a sequence of bases
# -----------------------------

println("generating basis functions")

BASES = []

B1 = [NBody(1.0)]

rcut2 = 9.2
rcut3 = 6.5

TRANSFORM = (@analytic r -> (2.9/r)^3)
CUTOFF2 = (:cos, 0.66*rcut2, rcut2)

D2 = Dictionary(TRANSFORM, CUTOFF2)
B2 = gen_basis(2, D2, 14)

for deg in [10,12,14]
   B = [B1; gen_basis(2, D2, deg)]
   push!(BASES, (B, D2, "2 / $(length(B)) / $rcut2"))
end

println("   first a few 3-body bases ...")
for (deg, rcut) in zip( [6, 8, 10],
                        [rcut3, rcut3, rcut3, rcut3] )
   CUTOFF3 = (:cos, 0.66*rcut, rcut)
   D = Dictionary(TRANSFORM, CUTOFF3)
   B3 = gen_basis(3, D, deg)
   B = [B1; B2; B3]
   push!(BASES, (B, D, "2+3 / $(length(B)) / $rcut2+$rcut"))
end

# D3 = Dictionary(TRANSFORM, (:cos, 0.66*rcut3, rcut3))
# B3 = gen_basis(3, D3, 8)
# for (deg, rcut) in zip([4, 6, 8],
#                        [4.5, 4.5, 4.5, 4.5])
#    CUTOFF4 = (:cos, 0.66*rcut, rcut)
#    D = Dictionary(TRANSFORM, CUTOFF4)
#    B4 = gen_basis(4, D, deg)
#    B = [B1; B2; B3; B4]
#    push!(BASES, (B, D, "2+3+4 / $(length(B)) / $rcut2...$rcut") )
# end

rmsE = Float64[]
rmsF = Float64[]
maeE = Float64[]
maeF = Float64[]

length_bases = [length(B[1]) for B in BASES]
@show length_bases

println("For each generated basis test the fit: " )
for (B, D, description) in BASES
   println(" - ", description)
   # standard least squares (see NBodyIPs/src/fitting.jl)
   # nforces = number of (randomly chosen) forces per configuration added
   #           to the LSQ problem
   ndata = min(length(train_data), 2 * length(B))
   c = regression(B, train_data[1:ndata], nforces = Inf, stab = 0.0)
   @show norm(c, Inf)
   # construct an IP from the the basis and the weights
   IP = NBodyIP(B, c)
   # check error => the normalisation is w.r.t. natoms, not a genuine
   # relative error; we can discuss
   rE, rF, mE, mF = fiterrors(IP, test_data)
   push!(rmsE, rE); push!(rmsF, rF)
   push!(maeE, mE); push!(maeF, mF)
   println("   E-rms, E-mae on testset = ", rE, ", ", mE)
   println("   F-rms, F-mae on testset = ", rF, ", ", mF)
end


# using JLD
# JLD.save("backup.jld", "desc", [B[3] for B in BASES],
#          "rmsE", rmsE, "rmsF", rmsF, "maeE", maeE, "maeF", maeF)

using DataFrames
df = DataFrame(:desc => [B[3] for B in BASES])
df[Symbol("rms-E")] = rmsE
df[Symbol("rms-F")] = rmsF
df[Symbol("mae-E")] = maeE
df[Symbol("mae-F")] = maeF
println("Energy and Force Errors for Ti-DFTB Database: ")
println(df)
