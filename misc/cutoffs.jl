using Plots

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

x = 0.456

function runN(f, x, N)
   s = 0.0
   for n = 1:N
      x += 1e-9
      s += f(x)
   end
   return s
end

@btime runN($cosenv, $x, 1_000)
@btime runN($p2env, $x, 1_000)
@btime runN($p2envsym, $x, 1_000)
@btime runN($p3env, $x, 1_000)
@btime runN($p3envsym, $x, 1_000)
