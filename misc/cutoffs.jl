using Plots, NBodyIPs, JuLIP
using NBodyIPs: Cutoff
using BenchmarkTools

##
r0, rnn, rc = 0.0, 1.0, 2+√2
xx = range(r0-0.5, rc+0.5, length=300)

ξ1 = @analytic r -> (r - rnn)/(rc-r0) * ((r-r0)/(rc-rnn) + (r-rc)/(r0-rnn))
ξ2 = @analytic r -> (r-rnn) * exp(-r)

plot(xx, ξ1.(xx))
scatter!([r0, rnn, rc], [-1,0,1], ms=10)

##

using LinearAlgebra, Roots
r0, rnn, rc = 0.75, 1.0, 2.1
xx = range(r0-0.5, rc+0.5, length=300)

# ξ(r) = A (1/(1+r) - 1/2) + B (r-1)
c = qr( [1/(r0+1)-1/2  r0-1; 1/(rc+1)-1/2  rc-1] ) \ [-1.0, 1.0]
ξ2 = r -> c[1]*(1/(r+1)-1/2) + c[2] * (r - 1)

# ξ(r) = C(exp( λ (r/r0 - 1) ) - 1)
gg = λ -> exp( λ*(rc/rnn - 1) ) + exp( λ*(r0/rnn - 1) ) - 2
@show λopt = find_zero(gg, -2.5)
C = 1 / (exp( λopt*(rc/rnn - 1) ) - 1)
ξ3  = r -> C * (exp( λopt*(r/rnn - 1) ) - 1)

# cubic, minimising curvature
#  (r-rnn)/(rc-rnn) * (1 + A (r-rc) + B (r-rc)^2)
B = (7*(r0^2 + r0*rc + rc^2 - 3*(r0 + rc)*rnn + 3*rnn^2)) / (
     2*(r0 - rc)*(r0 - rnn)*(2*r0^2 + 3*r0*rc + 2*rc^2 -
     7*(r0 + rc)*rnn + 7*rnn^2))
A = 2/(-r0 + rc) + 3/(4*(-r0 + rnn)) + (7*(rc - rnn)) / (
    4*(2*r0^2 + 3*r0*rc + 2*rc^2 - 7*(r0 + rc)*rnn + 7*rnn^2))
ξ4 = r -> (r-rnn)/(rc-rnn) * (1 + A*(r-rc) + B*(r-rc)^2)
χ = x -> (x+1)^2 * (x-1)^2


plot(xx, ξ1.(xx))
plot!(xx, ξ2.(xx))
plot!(xx, ξ3.(xx))
plot!(xx, ξ4.(xx))
scatter!([r0, rnn, rc], [-1,0,1], ms=10)

##
χ = x -> (x+1)^2 * (x-1)^2
penv2s = r -> χ(ξ3(r)) * (r0 < r < rc)
plot(xx, penv2s)
scatter!([r0, rc], [0,0])

##
cutcos = Cutoff((:cos2s, 0.7, 0.9, 2.3, 3.0))
cutpenv = Cutoff((:penv2s, 2, 0.7, 1.0, 3.0))

r0, rnn, rc = 0.7, 1.0, 3.0
ξ = @analytic r -> (r - rnn)/(rc-r0) * ((r-r0)/(rc-rnn) + (r-rc)/(r0-rnn))


env = @analytic x -> (x-1)^2 * (x+1)^2
f = r -> env(ξ(r)) * (r0 < r < rc)


xx = range(0, 3.5, length=300)
P1 = plot(xx, cutcos.f.(xx), label = "cos2s")
plot!(xx, cutpenv.f.(xx), label = "penv2s")



## one-sided cutoffs
cosenv(x) = (1 + cos(pi*x))/2
p2env(x) = (1-x)^2
p2envsym(x) = (1-x)^2*(1+x)^2
p3env(x) = (1-x)^3
p3envsym(x) = (1-x)^3*(1+x)^3

xx = range(0, 1, length=100)
P1 = plot(xx, cosenv.(xx), label = "cosenv")
plot!(xx, p2env.(xx), label = "p2env")
plot!(xx, p2envsym.(xx), label = "p2envsym")
plot!(xx, p3env.(xx), label = "p3env")
plot!(xx, p3envsym.(xx), label = "p3envsym")

## two-sided cutoffs

xx = range(-1, 1, length=100)
P2 = plot(xx, cosenv.(xx), label = "cosenv")
plot!(xx, p2envsym.(xx), label = "p2envsym")
plot!(xx, p3envsym.(xx), label = "p3envsym")
plot!(xx, cosenv.(xx).^2, label = "cosenv^2")

plot(P1, P2)

x = 1.456

function runN(f, x, N)
   s = 0.0
   for n = 1:N
      x += 1e-9
      s += f(x)
   end
   return s
end

penv = Cutoff((:penv2s, 2, 0.7, 1.0, 2.5))
cos2s = Cutoff((:cos2s, 0.7, 0.9, 1.7, 2.5))

@btime(runN($(penv.f), $x, 1_000))
@btime(runN($(cos2s.f), $x, 1_000))



@btime runN($cosenv, $x, 1_000)
@btime runN($p2env, $x, 1_000)
@btime runN($p2envsym, $x, 1_000)
@btime runN($p3env, $x, 1_000)
@btime runN($p3envsym, $x, 1_000)
