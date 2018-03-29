

using MacroTools, Espresso, Combinatorics, FunctionWrappers

using FunctionWrappers: FunctionWrapper
const F64fun = FunctionWrapper{Float64, Tuple{AbstractVector{Float64}}}

export psym_polys, psym_polys_tot, dict

@simple_rule (x^0) 1

dict(::Val{:poly}, n) =
    ["r^$i" for i = 0:n-1], "r"

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
function psym_monomial(alpha, dict = ["y^0","y^1","y^2"], sym = "y";
                        simplify = false)
	perms = uniqueperms(alpha)
	for j in 1:length(perms)
		ext = :()
		ext = replace(dict[perms[j][1]+1], sym, "(x[1])")
		for i = 2:length(alpha)
			ext = "$ext*" * replace(dict[perms[j][i]+1], sym, "(x[$i])")
		end
		if j == 1
			ex = "$ext"
		else
			ex = "$ex+$ext"
		end
	end
	ex = parse(ex)

    if simplify
        # two simplification passes (still doesn't catch everything)
        ex = Espresso.simplify(ex)
        ex = Espresso.simplify(ex)
    end
	f = F64fun(eval(:(x -> $(ex))))

	return [ex, f]
end

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
function psym_polys(dim::Integer, dict, sym; simplify = false)
    deg = length(dict)-1
    display(deg)
    mex, mf = psym_monomial([0], dict, sym; simplify=simplify)
	polys_ex = [mex]
	polys_f = [mf]
	for i in 1:deg
		for alpha in collect(partitions(i))
            if length(alpha) == dim
                mex, mf = psym_monomial(alpha, dict, sym; simplify=simplify)
                push!(polys_ex, mex)
				push!(polys_f, mf)
            elseif length(alpha) < dim
                add = zeros(Int64, dim - length(alpha))
                append!(alpha, add)
                mex, mf = psym_monomial(alpha, dict, sym; simplify=simplify)
                push!(polys_ex, mex)
				push!(polys_f, mf)
            end
        end
	end
	return polys_ex, polys_f
end

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
