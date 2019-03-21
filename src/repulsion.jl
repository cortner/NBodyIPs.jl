
module Repulsion

import NBodyIPs: NBodyFunction, descriptor, evaluate_I, evaluate_I_ed,
                 inv_transform, fcut

import JuLIP: AbstractCalculator, energy, forces, virial
import JuLIP.Potentials: @pot, evaluate, evaluate_d, PairPotential, @D, cutoff,
                         @analytic

struct RepulsiveCore{TV1, TV2, DT} <: NBodyFunction{2, DT}
   Vout::TV1         # the outer pair potential
   Vin::TV2          # the inner (repulsive) pair potential
   ri::Float64       # the interface between the two
   D::DT             # shouldn't be needed, but can't get it to work otherwise...
end

@pot RepulsiveCore

cutoff(V::RepulsiveCore) = cutoff(V.Vout)
descriptor(V::RepulsiveCore) = V.D

function evaluate_I(V::RepulsiveCore, II)
   I1, I2 = II
   @assert length(I1) == 1
   r = inv_transform(V.D.transform, I1[1])
   @show r, V.ri
   @show V.Vout(r), V.Vin(r)
   @show evaluate_I(V.Vout, II)
   @show I1, I2
   if r > V.ri
      return evaluate_I(V.Vout, II)
   end
   return V.Vin(r) / fcut(V.Vout.D, r)
end

function evaluate_I_ed(V::RepulsiveCore, II)
   I1, I2, dI1, dI2 = II
   @assert length(I1) == 1
   r = I1[1]
   if r > V.ri
      evaluate_I_ed(V.Vout, II)
   end
   dV = @D V.Vin(r)
   return dV * dI1[1]
end


function RepulsiveCore(Vout::NBodyFunction{2}, ri)
   v = Vout(ri)
   dv = @D Vout(ri)
   if dv >= 0.0
      @error("The slope `Vout'(ri)` is required to be negative")
   end
   if dv > -1.0
      @warn("""The slope `Vout'(ri) = $dv` doesn't seem steep enough to attach a
               repulsive core. Proceed at your own risk.""")
   end
   if v <= 0.0
      @error("it is required that `Vout(ri) > 0`.")
   end
   if v <= 1.0
      @warn("""Ideally the repulsive core should not be attached at small
               values of `Vout(ri) = $v`. Proceed at your own risk.""")
   end
   # B e^{-A (r/ri-1)} * ri/r
   #    => B = Vout(ri)
   # dv = - A*B/ri e^{-A (r/ri-1)} * ri/r - B*ri*e^{...} / r^2
   #    = - A/ri * v - 1/ri * v = - (1+A)/ri * v
   #    => -(1+A)/ri * v = dv
   #    => 1+A = - ri dv / v
   #    => A = -1 - ri dv / v
   Vin = let A = -1 - ri * dv / v, B = v, ri = ri
      @analytic r -> B * exp( - A * (r/ri-1) ) * ri/r
   end
   @show v, dv
   @show Vin(ri), (@D Vin(ri))
   # construct the piecewise potential
   return RepulsiveCore(Vout, Vin, ri, descriptor(Vout))
end

end

# EXAMPLE CODE => TURN THiS INTO A TEST!
# ## try out the repulsive potential
# Vfit = IPall.components[2]
# Vrep = NBodyIPs.Repulsion.RepulsiveCore(Vfit, 2.1)
#
# rp = range(0.3*r0, 9.0, length=500)
# plot(rp, Vfit.(rp), lw=2, label ="fit")
# plot!(rp, Vrep.Vin.(rp), lw=2, label ="inner", ylims = [-2.0, 10.0] )
# plot!(rp, Vrep.(rp), label = "combined")
