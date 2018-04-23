#
# using XGrad
#
# exf = :( Q = sum(r) )
# xdiff(exf; r=rand(3))

import Base.count

function _count(Rg, f)
   c = 0.0
   for I in Rg
      if length(unique(I.I)) == length(I)
         c += prod(f[i] for i in I.I)
      end
   end
   return c
end

function count(n, N, p, d)
   x = collect(-N:N)
   o = ones(2*N+1)
   if d == 1
      X = x
   elseif d == 2
      X = [ (x * o')[:] (o * x')[:] ]
   elseif d == 3
      X = [ kron(x, o, o)[:] kron(o, x, o)[:] kron(o, o, x)[:]]
   else
      error("d must be 1, 2, or 3")
   end
   f = setdiff( sqrt.(sum(abs2, X, 2)),  0.0).^(-p)
   M = length(f)

   return  _count(
      CartesianRange(
         CartesianIndex( ntuple(_->1, n) ),
         CartesianIndex( ntuple(_->M, n) )
      ),
      f
   )
end

Q = [ count(n, 7, 8.0, 3)  for n = 1:4 ]
