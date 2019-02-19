

using JuLIP, NBodyIPs, LinearAlgebra 

path = homedir() * "/Dropbox/PIBmat/Ti_DFTB_Data/Ti_N54_vartemp.xyz"

# Loading the training data
data = NBodyIPs.Data.read(path)
data_sm = data[1:10]

println("generating basis functions...")
r0 = rnn(:W)
rcuts = [5.2, 5.2, 4.5]  # 2-body, 3-body, 4-body
degrees = [10, 8, 6]
TRANSFORM = "@analytic r -> ($r0/r)^3"

# generate everything else
CUTOFFS = ["(:cos, $(0.66*rcut), $rcut)" for rcut in rcuts]
dicts = Dictionary.(TRANSFORM, CUTOFFS)
B = vcat(NBody(1.0), gen_basis.([2,3,4], dicts, degrees, purify=false)...)
@show length(B)

B3 = B[ findall(bodyorder.(B) .< 4) ]
# normalize_basis!(B3, data_sm)

Φ, _, _ = NBodyIPs.assemble_lsq(B, data_sm; verbose=true, nforces=10)

Q, R = qr(Φ)
A = R' * R



@show norm(A - Φ'*Φ, Inf)

display(A[1:13,1:13])

idx = 12  # this should be a 2-body
idx = 20  # this should be a pure 3-body
interesting_b = B3[idx]
# TODO how do we figure out that this is almost in the 2-body set
#  Φ[:,idx] almost in the span of Φ[:, 1:11]
Φ2 = Φ[:,1:11]
Q2, R2 = qr(Φ2)
for idx = 12:length(B3)
   c = Q2' * Φ[:,idx]
   err = Q2 * c - Φ[:,idx]
   @show B3[idx].t
   @show norm(err, Inf) / norm(Φ[:,idx], Inf)
end

# # [3] do the regression and construct and IP
# c = regression(B, train_data, nforces = Inf, stab = 0.0)
# @show norm(c, Inf)
# IP = NBodyIP(B, c)
# # check error => the normalisation is w.r.t. natoms, not a genuine
# # relative error; we can discuss
# rE, rF, mE, mF = fiterrors(IP, test_data)
# println("   E-rms, E-mae on testset = ", rE, ", ", mE)
# println("   F-rms, F-mae on testset = ", rF, ", ", mF)
#
# # [4] save the IP
# NBodyIPs.IO.write(@__DIR__() * "/ip.jld", IP)
