

using MacroTools,  Combinatorics, FunctionWrappers, Calculus, ForwardDiff,
      StaticArrays

using FunctionWrappers: FunctionWrapper
const FWrap{N, T} = FunctionWrapper{Float64, Tuple{SVector{N,T}}}
const GWrap{N, T} = FunctionWrapper{SVector{N,T}, Tuple{SVector{N,T}}}

export psym_polys, psym_polys_tot, dict

# @simple_rule (x^0) 1   # REVISIT Espresso

dict(::Val{:poly}, n) =
    ["r^$i" for i = 0:n-1], "r"

dict(v::Val{:poly}, n, rcut) = dict(v, n)

dict(::Val{:poly1}, n, rcut) =
    ["x^$i * (x^(-1)-$(rcut^(-1))+$(rcut^(-2))*(x-$rcut)" for i = 0:n-1], "x"

dict(::Val{:poly2}, n, rcut) =
    ["(x*$(1/rcut)-1.0)^$(2+i)" for i = 0:n-1], "x"

dict(::Val{:inv1}, n, rcut) =
    ["x^$(-i) - $(rcut^(-i)) + $(i * rcut^(-i-1)) * (x - $rcut)" for i = 1:n], "x"

dict(::Val{:inv2}, n, rcut) =
    ["x^$(-i) * (x*$(1/rcut)-1.0)^2" for i = 0:n-1], "x"

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
   for n = 1:dim
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
