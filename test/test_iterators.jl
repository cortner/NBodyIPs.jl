# This file contains some code that is chared across several tests

using NeighbourLists, NBodyIPs, JuLIP, StaticArrays, Test,
      ForwardDiff, Test
using NBodyIPs: NBodyFunction
import NBodyIPs: descriptor
import JuLIP.Potentials: evaluate, cutoff, evaluate_d!

# --------- MANY BODY CODE THAT IS SHARED ACROSS TESTS ------------

@pot struct TestBL{N} <: NBodyFunction{N, NBodyIPs.NullSiteDesc}
   r0::Float64
   rcut::Float64
   vN::Val{N}
end

function edgelengths1(Rs, J)
   r = [ norm(Rs[j])         for j in J ]
   S = [ Rs[j] / norm(Rs[j]) for j in J ]
   for i = 1:length(J)-1, j = i+1:length(J)
      push!(r, norm(Rs[J[j]]-Rs[J[i]]))
      push!(S, (Rs[J[i]]-Rs[J[j]]) / norm(Rs[J[i]]-Rs[J[j]]))
   end
   return r, S
end

# function edgelengths(Rs, J)
#    r = NBodyIPs.edge_lengths(Rs, J)

descriptor(::TestBL) = NBodyIPs.NullSiteDesc()

evaluate(V::TestBL{N}, Rs, ii, J) where {N} =
      fnbody(edgelengths1(Rs, J)[1], V.r0, V.rcut) / N

function evaluate_d!(dVsite, V::TestBL{N}, Rs, ii, J) where {N}
   r, S = edgelengths1(Rs, J)
   dV = (fnbody_d(r, V.r0, V.rcut) / N) .* S
   idx = 0
   for i = 1:length(J)
      idx += 1
      dVsite[J[i]] += dV[idx]
   end
   for i = 1:length(J)-1, j = i+1:length(J)
      idx += 1
      dVsite[J[i]] += dV[idx]
      dVsite[J[j]] -= dV[idx]
   end
   return dVsite
end

evaluate(V::TestBL, rs::Vector{Float64}) = fnbody(rs, V.r0, V.rcut)

cutoff(V::TestBL) = V.rcut

sqrt1p(x) = expm1( 0.5 * log1p( x ) )   # sqrt(1 + x) - 1
sqrt1p_d(x) = 0.5 / sqrt(1+x)

function fnbody(r, r0, rcut)
   E = sum( exp.(1.0 .- r./r0) )
   F = prod( (r./rcut - 1.0) .* (r .< rcut) )
   return sqrt1p(E) * F^2
end

function fnbody_d(r, r0, rcut)
   E = sum( exp.(1.0 .- r./r0) )
   F = prod( (r./rcut - 1.0) .* (r .< rcut) )
   ∇f = (- sqrt1p_d(E) * F^2 / r0) .* (exp.(1.0 .- r./r0))
   ∇f += (sqrt1p(E) * 2.0 / rcut * F^2) .* (r./rcut .- 1.0 - 1e-14).^(-1)
   return ∇f
end

# Generate a MODEL N-Body function
gen_fnbody(N; r0=1.0, rcut = N < 5 ? 2.5 * r0 : 1.7 * r0) =
   TestBL(r0, rcut, Val(N))

println("Checking that `fnbody` is correct...")
fnbody_ad(rs, r0, rcut) = ForwardDiff.gradient( t -> fnbody(t, r0, rcut), rs )
for n = 1:10
   print(".")
   r = SVector(rand(n)...)
   @test fnbody_d(r, 1.0, 2.0) ≈ fnbody_ad(r, 1.0, 2.0)
end
println()




function naive_n_body(X::Vector{JVec{T}}, V::NBodyFunction{N}) where {T, N}
   nX = length(X)
   E = 0.0
   cnt = 0
   fact_N = factorial(N)
   start = CartesianIndex(ntuple(_->1, N))
   stop = CartesianIndex(ntuple(_->nX, N))
   for j in CartesianIndices(start, stop)
      s = zeros(T, (N*(N-1))÷2); n = 0
      for a = 1:N-1, b = a+1:N
         n += 1
         s[n] = norm(X[j[a]] - X[j[b]])
      end
      if 0 < minimum(s) && maximum(s) <= V.rcut
         # each of these terms occur factorial(N) times, i.e. for
         # every permutation of the indices!
         E += V(s) / fact_N
      end
   end
   return E
end


println("--------------------------------------")
println("    Testing EdgeLength Implementation")
println("--------------------------------------")
for N in [2,2,3,3,4,4,5]
   Rs = rand(JVecF, 20)
   J = unique(rand(1:20, 30))
   if length(J) < N
      continue
   end
   J = SVector(J[1:N]...)
   r1, S1 = edgelengths1(Rs, J)
   r = NBodyIPs.edge_lengths(Rs, J)
   (@test r1 ≈ Vector(r)) |> print; print(" ")
   # (@test S1 ≈ Vector(S)) |> print; print(" ")
end
println()

println("--------------------------------------")
println("    Testing NBodyIterator")
println("--------------------------------------")
println("Check that the energy is consistent with a naive implementation")
MM = [2,2,3,3,3,4,4,5]  # body orders
println("   N     Nat    =>   |Emr-Enaive|")
for N in MM
   # create N-body
   VN = gen_fnbody(N; r0 = rnn(:Cu))
   if N < 5
      # create a not-too-large copper cell
      at = bulk(:Cu, cubic=true, pbc = false) * 2
   else
      at = bulk(:Cu, cubic=true, pbc = false) * (1,2,2)
   end

   rattle!(at, 0.1)
   nat = length(at)

   # assemble energy via neighbourlist and map-reduce
   Emr = energy(VN, at)
   # assemble energy naively
   Enaive = naive_n_body(positions(at), VN)

   print("   $N      $nat    =>   $(abs(Emr - Enaive)): ")
   (@test Emr ≈ Enaive) |> println
end


println("----------------------------------------------")
println("Finite-difference tests for forces and stress")
MM = [2,3,4,5] # [2,3,4,5]  # body orders
for N in MM
   # create N-body
   VN = gen_fnbody(N; r0 = rnn(:Cu))
   # create a not-too-large copper cell
   at = bulk(:Cu, cubic=true, pbc = true)
   rattle!(at, 0.03)
   set_constraint!(at, VariableCell(at))
   set_calculator!(at, VN)

   nat = length(positions(at))
   dE = gradient(at)
   dEh = nothing
   E = energy(at)
   println("[N = $N]")
   @printf("    h    | error \n")
   @printf("---------|----------- \n")
   x = dofs(at)
   errors = []
   for p = 2:11
      h = 0.1^p
      dEh = copy(dE)
      for n = 1:length(x)
         x[n] += h
         set_dofs!(at, x)
         dEh[n] = (energy(at) - E) / h
         x[n] -= h
      end
      push!(errors, vecnorm(dE - dEh, Inf))
      @printf(" %1.1e | %4.2e  \n", h, errors[end])
   end
   (@test minimum(errors) <= 1e-3 * maximum(errors)) |> println
end
