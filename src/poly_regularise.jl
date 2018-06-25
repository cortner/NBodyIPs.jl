
using Roots: find_zero, Bisection, Order2, FalsePosition
using Sobol: SobolSeq


"""
`regularise_2b(B, r0, r1; creg = 1e-2, Nquad = 20)`

construct a regularising stabilising matrix that acts only on 2-body
terms.

* `B` : basis
* `r0, r1` : upper and lower bound over which to integrate
* `creg` : multiplier
* `Nquad` : number of quadrature / sample points

```
   I = ∑_r h | ∑_j c_j ϕ_j''(r) |^2
     = ∑_{i,j} c_i c_j  ∑_r h ϕ_i''(r) ϕ_j''(r)
     = c' * M * c
   M_ij = ∑_r h ϕ_i''(r) ϕ_j''(r) = h * Φ'' * Φ''
   Φ''_ri = ϕ_i''(r)
```
"""
function regularise_2b(B::Vector, r0::Number, r1::Number;
                       creg = 1e-2, Nquad = 20)
   I2 = find(bodyorder.(B) .== 2)
   rr = linspace(r0, r1, Nquad)
   Φ = zeros(Nquad, length(B))
   h = (r1 - r0) / (Nquad-1)
   for (ib, b) in zip(I2, B[I2]), (iq, r) in enumerate(rr)
      Φ[iq, ib] = evaluate_dd(b, r) * sqrt(creg * h)
   end
   return Φ
end

function inv_transform2(x::T, r0::T, r1::T, transform) where {T}
   # Secant Bisection Method (transform is monotone)
   t0 = transform(r0) - x
   t1 = transform(r1) - x
   while abs(r1 - r0) > 1e-7
      r = (r1 * t0 - r0 * t1) / (t0 - t1)
      t = transform(r) - x
      if abs(t) < 1e-8
         return r
      end
      if t*t0 > 0
         r0 = r
         t0 = t
      else
         r1 = r
         t1 = t
      end
   end
   return r
end

function inv_transform(x, r0, r1, transform)
   f = r -> transform(r) - x
   r_a = find_zero(f, (r0, r1), FalsePosition(), xatol=1e-2)
   r_b = find_zero(f, r_a, Order2(), xatol = 1e-7)
end

function is_simplex(::Val{3}, r)
   return ((r[1] <= r[2]+r[3]) &&
           (r[2] <= r[3]+r[1]) &&
           (r[3] <= r[1]+r[2]))
end

function is_simplex(::Val{4}, r)
end


"""
* N : body-order
* D : dictionary
* r0, r1 : upper and lower bound on the bond-lengths
"""
function _sobol_seq(npoints, N, D, r0, r1)
   # compute the boundary in transformed coordinates
   x0 = D.transform(r1)
   x1 = D.transform(r0)
   # define the inverse transform
   inv_t = x -> inv_transform(x, r0, r1, D.transform, D.transform_d)
   # call the inner sobol function (function barrier)
   return _sobol_inner(Val((N*(N-1)) ÷ 2), npoints, inv_t, x0, x1)
end

function _sobol_inner(::Val{DIM}, npoints, inv_t, x0::T, x1::T ) where {DIM, T}
   # upper and lower bounds on the transformed variables
   @assert x0 < x1
   # Sobol sequence in the [x0, x1]^dim hypercube
   s = SobolSeq(DIM, x0, x1)
   # temporary storage
   t = zero(MVector{DIM, T})
   # output storage
   X = SVector{DIM, T}[]
   # generate points
   while length(X) < npoints
      next!(s, t)
      if is_simplex(inv_t.(t))
         push!(X, SVector(t))
      end
      print("*")
   end
   return X
end

# """
# regularise(N, B,
#
# * `N::Integer` : body order
# * `B::Vector` : basis
# """
# function regularise(N::Integer, B::Vector,
#
#    # get the indices of the N-body basis functions
#    Ib = find(bodyorder.(B) .== N)
#    # assume all have the same dictionary
#    # if not, then this regularisation is not valid
#    @assert all(b.D == B[Ib[1]].D  for b in B[Ib])
#
#    # construct a low discrepancy sequence
#
#
# end
