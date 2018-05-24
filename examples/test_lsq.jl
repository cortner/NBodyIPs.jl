
using NBodyIPs, JuLIP, Base.Test

datafile = homedir() * "/Dropbox/PIBmat/Ti_DFTB_Data/Ti_N54_vartemp_virials_ASE.xyz"
data = NBodyIPs.Data.read(datafile, index=":")
train_data = data
test_data = data[1:5:end]

println("generating basis functions")

B1 = [NBody(1.0)]
TRANSFORM = (@analytic r -> (2.9/r)^3)

rcut2 = 9.2
CUTOFF2 = (:cos, 0.66*rcut2, rcut2)
D2 = Dictionary(TRANSFORM, CUTOFF2)
B2 = gen_basis(2, D2, 8)

rcut3 = 6.5
CUTOFF3 = (:cos, 0.66*rcut3, rcut3)
D3 = Dictionary(TRANSFORM, CUTOFF3)
B3 = gen_basis(3, D3, 8)

rcut4 = 1.5 * rnn(:Ti)
CUTOFF4 = (:cos, 0.66*rcut4, rcut4)
D4 = Dictionary(TRANSFORM, CUTOFF4)
B4 = gen_basis(4, D4, 8)

# rcut5 = 1.5 * rnn(:Ti)
# CUTOFF5 = (:cos, 0.66*rcut5, rcut5)
# D5 = Dictionary(TRANSFORM, CUTOFF5)
# B5 = gen_basis(5, D5, 6)


B = [B1; B2; B3; B4]

Ψold, Yold, _ = NBodyIPs.assemble_lsq_old(B, data[1:10], nforces = Inf)
Ψ, Y, _ = NBodyIPs.assemble_lsq(B, data[1:10], nforces = Inf)
(@test Yold ≈ Y) |> println
(@test Ψold ≈ Ψ) |> println
