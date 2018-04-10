

using MacroTools,  Combinatorics, FunctionWrappers, Calculus, ForwardDiff,
      StaticArrays

using FunctionWrappers: FunctionWrapper
const FWrap{N, T} = FunctionWrapper{Float64, Tuple{SVector{N,T}}}
const GWrap{N, T} = FunctionWrapper{SVector{N,T}, Tuple{SVector{N,T}}}

export psym_polys, psym_polys_tot, dict

const CRg = CartesianRange
const CInd = CartesianIndex


dict(::Val{:poly}, n) =
    ["r^$i" for i = 1:n], "r"

dict(v::Val{:poly}, n, rcut) = dict(v, n)

dict(::Val{:poly1}, n, rcut) =
    ["x^$i * (x^(-1) - $(1/rcut) + $(1/rcut^2) * (x-$rcut))" for i = 0:n-1], "x"

dict(::Val{:poly2}, n, rcut) =
    ["(x*$(1/rcut)-1.0)^$(2+i)" for i = 0:n-1], "x"

dict(::Val{:inv1}, n, rcut) =
    ["x^$(-i) - $(rcut^(-i)) + $(i * rcut^(-i-1)) * (x - $rcut)" for i = 1:n], "x"

dict(::Val{:inv2}, n, rcut) =
    ["x^$(-i) * (x*$(1/rcut)-1.0)^2" for i = 0:n-1], "x"

dict(::Val{:exp1}, n, rcut) =
    ["exp(-$i * x) - $(exp(-i*rcut)) + $(i * exp(-i*rcut)) * (x-$rcut)" for i = 1:n], "x"

dict(sym, args...) = dict(Val(sym), args...)

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


"""
return a list of all the unique permutations of a vector alpha
"""
function uniqueperms(alpha::Vector{Int64})
   size = length(alpha)
   sort!(alpha)
   isFinished = 0
   perms = Vector{Int64}[]
   it = 0
   while isFinished == 0
      push!(perms, copy(alpha))
      k = -1
      for i = (size-1):-1:1
         if alpha[i] < alpha[i+1]
            k = i
            break
         end
      end
      if k == -1
         isFinished = 1
         break
      else
         ceilIndex = k + 1
         for i = k + 2:size
            if alpha[ceilIndex] > alpha[i] && alpha[i] > alpha[k]
               ceilIndex = i
            end
         end
         t = alpha[ceilIndex]
         alpha[ceilIndex] = alpha[k]
         alpha[k] = t
      end
      temp = zeros(Int64,size-k)
      for i = 1:(size-k)
         temp[i] = alpha[i + k]
      end
      sort!(temp)
      for i = 1:(size-k)
         alpha[i + k] = temp[i]
      end
   end
   return perms
end

"""
`psym_monomial(alpha, dict, var) -> Expr, Function`

Function to construct the symbol expression and a wrapped function for the
monomial symmetric polynomial corresponding to the vector alpha. We replace
`x[i]^1` by `x[i]`, which may be pointless/slow. Do not replace `x[i]^0` by `1`,
`:(1)` is not treated as an expression. Check that the letter used for the
variable is not used elsewhere (e.g x in exp) as the variable is automatically
replaced.

### Example
```
ex, f = psym_monomial([0,1,1], ["y^0","y^1","y^2"], "y")
f([1., 2., 3.])
```
"""
function psym_monomial(alpha, dict, sym; simplify = true)
   dict = ["1"; dict]
   dim = length(alpha)
   ex = ""
	for p in uniqueperms(alpha)
      ext = ""
      for i = 1:dim
			ext = "$ext*" * replace(dict[p[i]+1], sym, "$(sym)$i")
		end
      ex = ex * "+$(ext[2:end])"
	end
	ex = parse(ex[2:end])

   if simplify
      exs = Calculus.simplify(ex)
   else
      exs = ex
   end
   s = Symbol(sym)
   f = eval(:($s -> $(ind2vec(exs, dim, sym))))
   df = x -> ForwardDiff.gradient(f, x)
	return ex, f, df
end

psym_monomial(alpha, t::Tuple; kwargs...) =
      psym_monomial(alpha, t[1], t[2]; kwargs...)


"""
`psym_polys(dim, dict, sym):`

collect all permutation invariant polynomials up to degree `deg`
in dimension dim. This assembles up to

### Example
```
exs, fs = PermPolys(3, ["y^0","y^1","y^2"], "y")
fs[1]([1., 2., 3.])
```
"""
function psym_polys(dim::Integer, dict, sym; simplify = true)
	polys_ex = Expr[]
	polys_f = Function[]
	polys_df = Function[]
	for i in 1:length(dict)
		for m = 1:dim, alpha in collect(partitions(i, m))
         append!(alpha, zeros(Int, dim - length(alpha)))
         mex, mf, mdf = psym_monomial(alpha, dict, sym; simplify=simplify)
         push!(polys_ex, mex)
         push!(polys_f, mf)
         push!(polys_df, mdf)
      end
	end
	return polys_ex, polys_f, polys_df
end


"""
which dimensionality corresponds to a body-order
"""
nbody_dim(bo::Integer) = (bo * (bo-1)) ÷ 2

"""
`psym_polys_nbody(bo::Integer, dict, sym)`

* `bo` : body order
* `dict` : 1D dictionary
* `sym` : symbol used in the dictionary
"""
function psym_polys_nbody(N::Integer, dict, sym; simplify = true)
	polys_ex = Expr[]
	polys_f = Function[]
	polys_df = Function[]
   # get the lower and upper dimensionality for genuine N-body terms
   dim_lo = nbody_dim(N-1)+1
   dim_hi = nbody_dim(N)
   if length(dict) < dim_lo
      warn("the length of the dictionary is too short for $N-body terms")
   end
	for i in dim_lo:length(dict)
		for m = dim_lo:dim_hi, alpha in collect(partitions(i, m))
         append!(alpha, zeros(Int, dim_hi - length(alpha)))
         mex, mf, mdf = psym_monomial(alpha, dict, sym; simplify=simplify)
         push!(polys_ex, mex)
         push!(polys_f, mf)
         push!(polys_df, mdf)
      end
	end
	return polys_ex, polys_f, polys_df
end



# TODO: psym_polys_tot   has been neglected a bit, needs to be updated!
"""
`psym_polys_tot(dim::Integer, dict, sym)

collects all permutation invariant polynomials
based on all the one-variable Basis_fct in dimension dim.

### Example
```
psym_polys_tot(3, ["y^0","y^1","y^2"], "y")
```
"""
function psym_polys_tot(dim::Integer, dict, sym; simplify = false)
	polys_ex = [psym_monomial([0], dict, sym)[1]]
	polys_f = [psym_monomial([0], dict, sym)[2]]
	for i in 1:(dim*(length(dict)))
		for alpha in collect(partitions(i))
			if maximum(alpha)<=(length(dict)-1)
	         if length(alpha) == dim
	            push!(polys_ex, psym_monomial(alpha, dict, sym)[1])
		         push!(polys_f, psym_monomial(alpha, dict, sym)[2])
	         elseif length(alpha) < dim
	            add = zeros(Int64, dim - length(alpha))
	            append!(alpha, add)
	            push!(polys_ex, psym_monomial(alpha, dict, sym)[1])
               push!(polys_f, psym_monomial(alpha, dict, sym)[2])
	         end
			end
      end
	end
	return polys_ex, polys_f
end


# ================== hacked-together four-body terms ============

# 4-body = 4-simplex has 4 corners (atom positions) and 6 edges  rᵢⱼ
#
# edge index : A[1]  A[2]  A[3]  A[4]  A[5]  A[6]
# edge length: r12   r13   r14   r23   r24   r34
# (where rij = |xᵢ - xⱼ| with xᵢ the corner positions)

const b4_e_inds = [0 1 2 3
                   1 0 4 5
                   2 4 0 6
                   3 5 6 0]

"""
convert a permutation of simplex corners into a permutation of
simplex edges
"""
S4_to_S6(π::Vector{Int}, b4_e_inds=NBodyIPs.b4_e_inds) = Int[
   b4_e_inds[π[1], π[2]], b4_e_inds[π[1], π[3]], b4_e_inds[π[1], π[4]],
   b4_e_inds[π[2], π[3]], b4_e_inds[π[2], π[4]], b4_e_inds[π[3], π[4]] ]

"""
generate all permutations of A that correspond to permutations of corners,
then keep only the unique ones so that we don't double-count.
"""
fourbody_permutations(A::Vector{Int}) =
   unique(  [ A[S4_to_S6(πX)]
              for πX in permutations(1:4) ]  )


function polys_fourbody(dict_len::Integer)
   basis = Vector{Vector{Int}}[]  # representation of the basis functions
   # get the lower and upper dimensionality for genuine N-body terms
   N = 4
   dim_lo = nbody_dim(N-1)+1   # 4
   dim_hi = nbody_dim(N)       # 6

   alldone = Vector{Int}[]

   # add the strange deg-2 terms
   for i = 1:dict_len, j = i:dict_len
      if i + j > dict_len
         continue
      end
      α = zeros(Int, 6)
      α[b4_e_inds[1,2]] = i
      α[b4_e_inds[3,4]] = j
      if !(α ∈ alldone)
         A = fourbody_permutations(α)
         append!(alldone, A)
         push!(basis, A)
      end
   end

   # add the strange deg-3 terms
   for i1 = 1:dict_len, i2 = 1:dict_len, i3 = 1:dict_len
      if i1+i2+i3 > dict_len
         continue
      end
      α = zeros(Int, 6)
      α[b4_e_inds[1,2]] = i1
      α[b4_e_inds[2,3]] = i2
      α[b4_e_inds[3,4]] = i3
      if !(α ∈ alldone)
         A = fourbody_permutations(α)
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
	for      i in 4:dict_len,     # (sum of tuple is between 4 and dict_len)
            m = dim_lo:dim_hi,   # (length of tuple is between 4 and 6)
            α in collect(partitions(i, m))
      # any terms not included get zeros appended
      append!(α, zeros(Int, dim_hi - length(α)))
      # store which tuples we've already covered
      alldone = Vector{Int}[]
      # look at all permutations of α that actually modify α
      for A in permutations(α)
         if !(A ∈ alldone)   # (not yet encountered)
            # not need to generate 4! = 24 (unique) permutations of A
            # that correspond to permutations of the corners
            P = fourbody_permutations(A)
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

Base.Vector(I::CartesianIndex) = Int[I.I...]

# len = length of dictionary, not including f(x) = 1
function polys_fourbody2(len::Integer)
   # representation of the basis functions and Lists of Tuples
   basis = Vector{Vector{Int}}[]
   # store all the tuples that are already in a basis function
   alldone = Vector{Int}[ [0,0,0,0,0,0] ]
   # for i1 = 0:len, i2 = 0:len, ..., i6 = 0:len
   #   (len+1)^6
   for I in CRg(CInd{6}(0,0,0,0,0,0), CInd{6}(ntuple(_->len, 6)))
      A = Vector(I)
      if sum(A) > len
         continue
      end
      if !(A ∈ alldone)
         P = fourbody_permutations(A)  # P::Vector{Vector{Int}}
         append!(alldone, P)
         push!(basis, P)
      end
	end
   return basis
end


function gen_fun(A::Vector{Vector{Int}}, dict, sym; simplify=true)
   dict = ["1"; dict]
   dim = length(A[1])

   # generate a string
   fstr = ""
   for α ∈ A
      α .+= 1
      fstr = fstr * " + " * replace(dict[α[1]], sym, "$(sym)1")
      for i = 2:length(α)
         fstr = fstr * " * " * replace(dict[α[i]], sym, "$(sym)$i")
      end
   end

   # generate an expression function
   fex = parse(fstr)
   if simplify
      fex = Calculus.simplify(fex)
   end

   # generate a function
   s = Symbol(sym)
   f = eval(:($s -> $(ind2vec(fex, dim, sym))))
   df = x -> ForwardDiff.gradient(f, x)

   return fex, f, df
end


function polys_fourbody(dict, sym; simplify=true)
   A = polys_fourbody(length(dict))
   return gen_fun(A, dict, sym; simplify=simplify)
end
