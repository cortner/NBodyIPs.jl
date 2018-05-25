

using JuLIP, NBodyIPs

path = homedir() * "/Dropbox/PIBmat/Ti_DFTB_Data/Ti_N54_vartemp.xyz"

# Loading the training data
data = NBodyIPs.Data.read(path)
data_sm = data[1:10]

println("generating basis functions...")
r0 = rnn(:W)
rcuts = [7.2, 5.2, 4.5]  # 2-body, 3-body, 4-body
degrees = [10, 8, 6]
TRANSFORM = "@analytic r -> ($r0/r)^3"

# generate everything else
CUTOFFS = ["(:cos, $(0.66*rcut), $rcut)" for rcut in rcuts]
dicts = Dictionary.(TRANSFORM, CUTOFFS)
B = vcat(NBody(1.0), gen_basis.([2,3,4], dicts, degrees, purify=false)...)
@show length(B)

B3 = B[ find(bodyorder.(B) .< 4) ]
# normalize_basis!(B3, data_sm)

Ψ, Y, _ = NBodyIPs.assemble_lsq(B, data_sm; verbose=true, nforces=3)

Q, R = qr(Ψ)
A = R' * R

n = size(Q,1)
m = size(Q,2)

# w_vec = (m/n)./diag(Q * Q')
w_vec = (m/n) * [ 1/norm(Q[n,:])^2 for n = 1:size(Q,1) ]
@show w_vec
@show minimum(w_vec)
@show maximum(w_vec)

Wopt = (size(Q,2)/size(Q,1)) * Diagonal([ 1/(norm(Q[n,:])^2) for n = 1:size(Q,1) ])



# Wopt = diagm(w_vec)

cstab=0

Aopt = Q' * (Wopt * Q) + cstab * eye(size(R, 1))
bopt = Q' * (Wopt * Y)
copt = R \ (Aopt \ bopt)

sqrtW = sqrt.(W)
sqrtWopt = sqrt(size(Q,2)/size(Q,1)) * Diagonal([ 1/sqrt(norm(Q[n,:])^2) for n = 1:size(Q,1) ])
tildeΨ = sqrtWopt * Ψ
tildeQ, tildeR = qr(tildeΨ)
tildeY = sqrtWopt * Y

A = (1 + cstab) * tildeR' * tildeR
b = tildeR' * tildeQ' * tildeY
cnew = A \ b

norm(copt-cnew,Inf)

copt2 = ((R' * Aopt) * R) \ (R' * bopt)
copt2b = R' \ (((R' * Aopt) * R) \ bopt)
copt3 =  (Aopt * R) \  bopt

A = R' * (Aopt * R)
b = R' * bopt
c = A \ b


zopt = Ψ * copt - Y
maximum(abs.(zopt))
@show norm(zopt.^2)
rmsopt = sqrt(dot( zopt, zopt) / length(Y))
rmsoptw = sqrt(dot( Wopt * zopt, zopt) / length(Y))

println("naive rms error on training set with Wopt: ", rmsoptw)



regression(B, data_sm;
                    verbose = true,
                    nforces=3, usestress=false,
                    stabstyle=:basis, cstab=0.,
                    weights=:CM,
                    regulariser = nothing)

Wnaive = eye(length(w_vec),length(w_vec));

Anaive = Q' * (Wnaive * Q) + cstab * eye(size(R, 1))
bnaive = Q' * (Wnaive * Y)
cnaive = R \ (Anaive \ bnaive)

znaive = Ψ * cnaive - Y
@show norm(znaive.^2)
maximum(abs.(znaive))
rmsnaive = sqrt(dot(Wnaive * znaive, znaive) / length(Y))
println("naive rms error on training set with W=Id: ", rmsnaive)


# @show vecnorm(A - Φ'*Φ, Inf)
#
# display(A[1:13,1:13])
#
# idx = 12  # this should be a 2-body
# idx = 20  # this should be a pure 3-body
# interesting_b = B3[idx]
# # TODO how do we figure out that this is almost in the 2-body set
# #  Φ[:,idx] almost in the span of Φ[:, 1:11]
# Φ2 = Φ[:,1:11]
# Q2, R2 = qr(Φ2)
# for idx = 12:length(B3)
#    c = Q2' * Φ[:,idx]
#    err = Q2 * c - Φ[:,idx]
#    @show B3[idx].t
#    @show norm(err, Inf) / norm(Φ[:,idx], Inf)
# end

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
