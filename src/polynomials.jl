

# module Polynomials


using MacroTools
using Combinatorics
import FunctionWrappers
import FunctionWrappers: FunctionWrapper
F64fun = FunctionWrapper{Float64, Tuple{AbstractVector{Float64}}}

export MonSymPol, PermPolys, uniqueperms, Fctbasisgen, PermPolys_all



struct PermPolys
	polysSymbol::Array{Expr}
	polys::Array{F64fun}
end

struct PermPolys_all
	polysSymbol::Array{Expr}
	polys::Array{F64fun}
end

struct MVPolys
	polysSymbol::Array{Expr}
	polys::Array{F64fun}
end


#Returns a matrix of one-variable basis functions
function Fctbasisgen(Nbparticles,Deg, BasisFct = "x^j")
	BasisFcts = Array{String, 2}(length(Deg), 1)
	for (ik,Degk) in enumerate(Deg)
		BasisFcts[ik] = replace(BasisFct, "j", "$Degk")
	end
	Fcts = Array{String, 2}(length(BasisFcts), Nbparticles)
	for i=1:length(BasisFcts)
		for j=1:Nbparticles
			exout = replace(BasisFcts[i], "x", "(x[$j])")
			Fcts[i,j] = exout
		end
	end
	return Fcts
end
#Example
#BF = Fctbasisgen(3, collect(0:1:1), "x^j")


# Function which returns a list of all the unique permutations of a vector alpha.
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

# Function to construct the symbol expression and a wrapped function for the
# monomial symmetric polynomial corresponding to the vector alpha.
# We replace x[i]^1 by x[i], which may be pointless/slow.
# Do not replace x[i]^0 by 1, :(1) is not treated as an expression.
function MonSymPol(alpha, BasisFcts = ["y^0","y^1","y^2"])
	perms = uniqueperms(alpha)
	for j in 1:length(perms)
		ext = :()

		ext = replace(BasisFcts[perms[j][1]+1], "y", "(x[1])")

		for i = 2:length(alpha)
			ext = "$ext*" * replace(BasisFcts[perms[j][i]+1], "y", "(x[$i])")
		end
		if j == 1
			ex = "$ext"
		else
			ex = "$ex+$ext"
		end
	end
	ex = parse(ex)
	ex2 = ex
	ex2 = :(x -> $(ex2))
	f = F64fun(eval(ex2))
	return [ex, f]
end
#Example
#MSP = MonSymPol([0,1,1], ["x[1]", "x[2]", "x[3]"])
#MSP[2]([1., 2., 3.])

function MonSymFct(alpha, Basis_fcts = ["(x[1])" "(x[2])" ; "((x[1])^2)" "((x[2])^2)"])
	ext = :()
	perms = uniqueperms(alpha)
	for i in 1:length(perms)
		exttemp = :()
		for j in 1:length(alpha)
			if j==1
				exttemp = Basis_fcts[perms[i][j],j];
			else
				exttemp = "$exttemp*" * Basis_fcts[perms[i][j],j]
			end
		end
		if i==1
			ext = "$exttemp";
		else
			ext = "$ext+$exttemp";
		end
	end
	ex = "$ext"
	ex = parse(ex)
	ex2 = ex
	ex2 = :(x -> $(ex2))
	display(ex2)
	f = F64fun(eval(ex2))
	return [ex, f]
end
#Example
# Symfct = MonSymFct([1, 2, 2], BF)
# SF = Symfct[2]
# SF([1., 2., 3.])





function PermFct(Basis_fct::Array{String,2},maxdeg::Int64)
	deg = size(Basis_fct,1);
	display(deg)
	dim = size(Basis_fct,2);

	polys1 = []
	polys2 = []
	for i in 1:maxdeg #(deg*dim)
		for alpha in collect(partitions(i))
            if (length(alpha) == dim) && (maximum(alpha)<=deg)
                push!(polys1, MonSymFct(alpha,Basis_fct)[1])
				push!(polys2, MonSymFct(alpha,Basis_fct)[2])

			end
        end
	end
	return PermFct(polys1, polys2)
end


# Function which collects all permutation invariant polynomials up to degree deg
# in dimension dim.
function PermPolys(deg::Int64, dim::Int64, Basis_fct = ["x^0","x^1","x^2"])
	polys1 = [MonSymPol([0],Basis_fct)[1]]
	polys2 = [MonSymPol([0],Basis_fct)[2]]
	for i in 1:deg
		for alpha in collect(partitions(i))
            if length(alpha) == dim
                push!(polys1, MonSymPol(alpha,Basis_fct)[1])
				push!(polys2, MonSymPol(alpha,Basis_fct)[2])
            elseif length(alpha) < dim
                add = zeros(Int64, dim - length(alpha))
                append!(alpha, add)
                push!(polys1, MonSymPol(alpha,Basis_fct)[1])
				push!(polys2, MonSymPol(alpha,Basis_fct)[2])
            end
        end

	end
	return PermPolys(polys1, polys2)
end
#Example
#PermP = PermPolys(1,3,["x[1]","x[2]","x[3]"])
#PermPP = PermP.polys
#PermPP[1]([1., 2., 3.])

# Function which collects all permutation invariant polynomials
#based on all the one-variable Basis_fct
# in dimension dim.
function PermPolys_all(dim::Int64, Basis_fct = ["x^0","x^1","x^2"])
	polys1 = [MonSymPol([0],Basis_fct)[1]]
	polys2 = [MonSymPol([0],Basis_fct)[2]]
	for i in 1:(dim*(length(Basis_fct)))
		for alpha in collect(partitions(i))
			if maximum(alpha)<=(length(Basis_fct)-1)
	            if length(alpha) == dim
	                push!(polys1, MonSymPol(alpha,Basis_fct)[1])
					push!(polys2, MonSymPol(alpha,Basis_fct)[2])
	            elseif length(alpha) < dim
	                add = zeros(Int64, dim - length(alpha))
	                append!(alpha, add)
	                push!(polys1, MonSymPol(alpha,Basis_fct)[1])
					push!(polys2, MonSymPol(alpha,Basis_fct)[2])
	            end
			end
        end
	end
	return PermPolys_all(polys1, polys2)
end




# end
