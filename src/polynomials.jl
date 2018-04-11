

using MacroTools,  Combinatorics, Calculus, ForwardDiff, StaticArrays

import Base.parse

export nbody_tuples, nbody_alltuples, nbody_onlytuples

const CRg = CartesianRange
const CInd = CartesianIndex
const Tup{M} = NTuple{M, Int}
const VecTup{M} = Vector{NTuple{M, Int}}

using FunctionWrappers: FunctionWrapper
const FWrap{N, T} = FunctionWrapper{T, Tuple{SVector{N,T}}}
const GWrap{N, T} = FunctionWrapper{SVector{N,T}, Tuple{SVector{N,T}}}

dim(::FWrap{DIM}) where DIM = DIM
dim(::GWrap{DIM}) where DIM = DIM
bodyorder(::FWrap{DIM}) where DIM = ceil(Int, sqrt(2*DIM))    # DIM = N * (N-1) / 2
bodyorder(::GWrap{DIM}) where DIM = ceil(Int, sqrt(2*DIM))



"""
`vec2ind(exorstr)`

take an expression, such as,
```
:( x[1]^2 * x[2]^3 * x[3] + x[1] * x[2] * x[3] + . . . )
```
and convert it to
```
:( x1^2 * x2^3 * x3 + x1 * x2 * x3 + . . . )
```
"""
function vec2ind(ex, dim, sym)
   str = string(ex)
   for n = 1:dim
      str = replace(str, "$sym[$n]", "$sym$n")
   end
   return parse(str)
end

"""
`ind2vec(exorstr)`

take an expression, such as,
```
:( x1^2 * x2^3 * x3 + x1 * x2 * x3 + . . . )
```
and convert it to
```
:( x[1]^2 * x[2]^3 * x[3] + x[1] * x[2] * x[3] + . . . )
```
"""
function ind2vec(ex, dim, sym)
   str = string(ex)
   for n = dim:-1:1
      str = replace(str, "$sym$n", "$sym[$n]")
   end
   return parse(str)
end


include("polynomials_legacy.jl")


# ================== hacked-together basis generation  ============
# ================== up to 4-body terms                ============

# 4-body = 4-simplex has 4 corners (atom positions) and 6 edges  rᵢⱼ
#
# edge index : A[1]  A[2]  A[3]  A[4]  A[5]  A[6]
# edge length: r12   r13   r14   r23   r24   r34
# (where rij = |xᵢ - xⱼ| with xᵢ the corner positions)

const b4_e_inds = [0 1 2 3
                   1 0 4 5
                   2 4 0 6
                   3 5 6 0]

const b4_e_proj = [1 1 1 0 0 0
                   1 0 0 1 1 0
                   0 1 0 1 0 1
                   0 0 1 0 1 1]

# 5-body = 5-simplex has 5 corners (atom positions) and 10 edges  rᵢⱼ
#
# edge index : A[1]  A[2]  A[3]  A[4]  A[5]  A[6] A[7]  A[8]  A[9]  A[10]
# edge length: r12   r13   r14   r15   r23   r24  r25   r34   r35   r45
# (where rij = |xᵢ - xⱼ| with xᵢ the corner positions)

const b5_e_inds = [0 1 2 3 4
                   1 0 5 6 7
                   2 5 0 8 9
                   3 6 8 0 10
                   4 7 9 10 0]

const b5_e_proj = [1 1 1 1 0 0 0 0 0 0
                   1 0 0 0 1 1 1 0 0 0
                   0 1 0 0 1 0 0 1 1 0
                   0 0 1 0 0 1 0 1 0 1
                   0 0 0 1 0 0 1 0 1 1]

const πb3 = collect(permutations(1:3))

"""
convert a permutation of simplex corners into a permutation of
simplex edges (for 4 body)
"""
S4_to_S6(π::Vector{Int}, b4_e_inds=NBodyIPs.b4_e_inds) = Int[
   b4_e_inds[π[1], π[2]], b4_e_inds[π[1], π[3]], b4_e_inds[π[1], π[4]],
   b4_e_inds[π[2], π[3]], b4_e_inds[π[2], π[4]], b4_e_inds[π[3], π[4]] ]

"""
convert a permutation of simplex corners into a permutation of
simplex edges (for 5 body)
"""
S5_to_S10(π::Vector{Int}, b5_e_inds=NBodyIPs.b5_e_inds) = Int[
   b5_e_inds[π[1], π[2]], b5_e_inds[π[1], π[3]], b5_e_inds[π[1], π[4]],
   b5_e_inds[π[1], π[5]], b5_e_inds[π[2], π[3]], b5_e_inds[π[2], π[4]],
   b5_e_inds[π[2], π[5]], b5_e_inds[π[3], π[4]], b5_e_inds[π[3], π[5]],
   b5_e_inds[π[4], π[5]] ]


"""
generate all permutations of A that correspond to permutations of corners,
then keep only the unique ones so that we don't double-count.
"""
simplex_permutations(::Val{5}, A) =
   unique(  [ A[S5_to_S10(πX)]
              for πX in permutations(1:5) ]  )
simplex_permutations(::Val{4}, A) =
   unique(  [ A[S4_to_S6(πX)]
              for πX in permutations(1:4) ]  )
simplex_permutations(::Val{3}, A) = unique( [ A[π] for π in πb3 ] )
simplex_permutations(::Val{2}, A) = [ A ]



# len = length of dictionary, not including f(x) = 1
nbody_alltuples(N::Integer, len::Integer) =
      nbody_alltuples(Val(N), Val((N*(N-1))÷2), len)

Base.ntuple(n::Integer, ::Val{M}) where {M} = ntuple(_->n, M)

function nbody_alltuples(vN::Val{N}, vM::Val{M}, len::Integer) where {N, M}
   # representation of the basis functions and Lists of Tuples
   basis = VecTup{M}[]
   # store all the tuples that are already in a basis function
   alldone = Tup{M}[ ntuple(0, vM) ]
   # for i1 = 0:len, i2 = 0:len, ..., iM = 0:len
   #   (len+1)^M
   for I in CRg(CInd(ntuple(0, vM)), CInd(ntuple(len, vM)))
      A = I.I
      if sum(A) > len
         continue
      end
      if !(A ∈ alldone)
         P = simplex_permutations(vN, A)
         append!(alldone, P)
         push!(basis, P)
      end
	end
   return basis
end

# #Select only Nbody tuples
nbody_onlytuples(N::Integer, len::Integer) =
      nbody_onlytuples(Val(N), Val((N*(N-1))÷2), len)

function nbody_onlytuples(vN::Val{N}, vM::Val{M}, len::Integer) where {N, M}
   # representation of the basis functions and Lists of Tuples
   basis = VecTup{M}[]
   # store all the tuples that are already in a basis function
   alldone = Tup{M}[ ntuple(0, vM) ]
   # for i1 = 0:len, i2 = 0:len, ..., iM = 0:len
   #   (len+1)^M
   for I in CRg(CInd(ntuple(0, vM)), CInd(ntuple(len, vM)))
      A = I.I
      if sum(A) > len
         continue
      end
      if !(A ∈ alldone)
         P = simplex_permutations(vN, A)
         append!(alldone, P)
           if trueNbody(A)
            push!(basis, P)
           end
      end
	end
   return basis
end

 function trueNbody(A)
    if length(A)==6
       S = sum((A.*b4_e_proj[i,:]) == [0,0,0,0,0,0] for i=1:4)
       return (S == 0)
    elseif length(A)==10
       S = sum((A.*b5_e_proj[i,:]) == [0,0,0,0,0,0] for i=1:5)
       return (S == 0)
     else
        error("pb: not implemented in trueNbody")
     end
 end

trueNbody((1,1,1,1,1,1,1,1,1,1))


nbody_tuples(N::Integer, len::Integer) =
      nbody_tuples(Val(N), Val((N*(N-1))÷2), len)


nbody_tuples(::Val{2}, ::Val{1}, len::Integer) =
      [ [ntuple(i, Val(1))] for i = 1:len ]

function nbody_tuples(::Val{3}, ::Val{3}, len::Integer)
   B = VecTup{3}[]
	for i in 2:len, m in 2:3, α in collect(partitions(i, m))
      append!(α, zeros(Int, 3 - length(α)))
      push!(B, simplex_permutations(Val(3), tuple(α...)))
   end
   return B
end

function nbody_tuples(vN::Val{4}, vM::Val{6}, len::Integer)
   # representation of the basis functions
   basis = VecTup{6}[]

   alldone = Tup{6}[]
   # add the strange deg-2 terms
   for i = 1:len, j = i:len
      if i + j > len
         continue
      end
      α = zeros(Int, 6)
      α[b4_e_inds[1,2]] = i
      α[b4_e_inds[3,4]] = j
      αt = tuple(α...)
      if !(αt ∈ alldone)
         P = simplex_permutations(Val(4), αt)
         append!(alldone, P)
         push!(basis, P)
      end
   end

   # add the strange deg-3 terms
   alldone = Tup{6}[]
   for i1 = 1:len, i2 = 1:len, i3 = 1:len
      if i1+i2+i3 > len
         continue
      end

      # type 1
      α = zeros(Int, 6)
      α[b4_e_inds[1,2]] = i1
      α[b4_e_inds[2,3]] = i2
      α[b4_e_inds[3,4]] = i3
      αt = tuple(α...)
      if !(αt ∈ alldone)
         A = simplex_permutations(Val(4), αt)
         append!(alldone, A)
         push!(basis, A)
      end

      # type 2
      α = zeros(Int, 6)
      α[b4_e_inds[1,2]] = i1
      α[b4_e_inds[1,3]] = i2
      α[b4_e_inds[1,4]] = i3
      αt = tuple(α...)
      if !(αt ∈ alldone)
         A = simplex_permutations(Val(4), αt)
         append!(alldone, A)
         push!(basis, A)
      end
   end

   # partitions(deg, m) -> all tuples of length m summing to deg
   #
   # generate all partitions of `i` from m integers where
   # dim_lo=4 ≦ m ≦ dim_hi=6. Fewer terms means it is an (N-1)-body term.
   # More are not allowed since that would involve at least (N+1)-bodies;
   #    `i` runs from dim_lo to length(dict); this gives all possible
   #    ordered tuples ⇔ multi-variate polynomials of sum-degree between
   #    dim_lo and length(dict)

	for      i in 4:len,     # (sum of tuple is between 4 and dict_len)
            m = 4:6,   # (length of tuple is between 4 and 6)
            α in collect(partitions(i, m))
      # any terms not included get zeros appended
      append!(α, zeros(Int, 6 - length(α)))
      # store which tuples we've already covered
      alldone = Tup{6}[]
      # look at all permutations of α that actually modify α
      for A in permutations(α)
         At = tuple(A...)
         if !(At ∈ alldone)   # (not yet encountered)
            # not need to generate 4! = 24 (unique) permutations of A
            # that correspond to permutations of the corners
            P = simplex_permutations(Val(4), At)
            # then add all of these to `alldone` so that we don't generate that
            # basis function a second time!
            append!(alldone, P)
            # the vector of tuples P represents a basis function
            push!(basis, P)
         end
      end
	end
   return basis
end


# ============= Generate Functions from Tuples ============

function parse(A::VecTup{M}, dict, sym; simplify=true, wrap=false) where {M}

   # add a `1` to the dictionary which corresponds to zero-entries in the
   # tuples
   dict = ["1"; dict]

   # generate a string representing the function
   fstr = ""
   for α ∈ A
      fstr = fstr * " + " * replace(dict[α[1]+1], sym, "$(sym)1")
      for i = 2:length(α)
         fstr = fstr * " * " * replace(dict[α[i]+1], sym, "$(sym)$i")
      end
   end
   fstr = fstr[4:end] # remove the first occurance of +

   # generate an expression function
   fex = parse(fstr)
   if simplify
      fex = Calculus.simplify(fex)
   end

   # generate a function
   s = Symbol(sym)
   f = eval(:($s -> $(ind2vec(fex, M, sym))))
   df = x -> ForwardDiff.gradient(f, x)

   if wrap
      f = FWrap{M, Float64}(f)
      df = GWrap{M, Float64}(df)
   end

   return fex, f, df
end

function parse(B::Vector{VecTup{M}}, dict, sym; wrap = false, kwargs...) where M
   exs = Expr[]
   if wrap
      fs = FWrap{M, Float64}[]
      dfs = GWrap{M, Float64}[]
   else
      fs = Function[]
      dfs = Function[]
   end
   for b in B
      ex, f, df = parse(b, dict, sym; kwargs...)
      push!(exs, ex)
      push!(fs, f)
      push!(dfs, df)
   end
   return exs, fs, dfs
end
