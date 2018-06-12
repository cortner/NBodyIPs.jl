
module OrthPolys

# TODO: Clenshaw Algorithm 

function nonorth!(p::AbstractVector{T}, x::T)
   p[1] = one(T)
   for n = 2:length(p)
      p[n] = p[n-1] * x
   end
   return p
end

function legendre!(p::AbstractVector{T}, x::T)
   p[1] = one(T)
   p[2] = x
   for n = 3:length(p)
      # n P_n(x) = (2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x)
      p[n] = (2-1/n) * x * p[n-1] - (1-1/n) * p[n-2]
   end
   return p
end

function chebyshev!(p::AbstractVector{T}, x::T)
   p[1] = one(T)
   p[2] = x
   for n = 3:length(p)
      # T_{n}(x) = 2xT_{n-1}(x) - T_{n-2}(x)
      p[n] = 2 * x * p[n-1] - p[n-2]
   end
   return p
end


end
