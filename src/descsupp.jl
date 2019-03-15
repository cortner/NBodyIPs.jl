import JuLIP

# TODO: get rid of this file altogether

using JuLIP.Potentials: @analytic,
                        AnalyticFunction,
                        F64fun

using Roots

import Base: Dict,
             ==

# -------------- Space Tranformations ---------------



# -------------- Cut-off Mechanisms ---------------

fcut_analyse(descr::String) = fcut_analyse(eval(Meta.parse(descr)))
fcut_analyse(args...) = fcut_analyse(args)
fcut_analyse(C::NBodyIPs.Cutoffs.NBCutoff) = C

function fcut_analyse(args::Tuple)
   sym = args[1]::Symbol
   args = args[2:end]
   if Symbol(sym) == :cos
      return CosCut(args...)
   elseif Symbol(sym) == :square
      return PolyCut(2, args[1])
   elseif Symbol(sym) == :cos2s
      return CosCut2s(args...)
   elseif Symbol(sym) == :penv
      return PolyCutSym(args...)
   elseif Symbol(sym) == :penv2s
      return PolyCut2s(args...)

   # elseif Symbol(sym) == :sw
   #    return let L=args[1], rcut=args[2]
   #       r -> cutsw(r, rcut, L), r -> cutsw_d(r, rcut, L), rcut
   #    end
   # elseif Symbol(sym) == :spline
   #    return let rc1=args[1], rc2=args[2]
   #       r -> cutsp(r, rc1, rc2), r -> cutsp_d(r, rc1, rc2), rc2
   #    end
   # elseif Symbol(sym) == :s2rat
   #    return let rnn=args[1], rin = args[2], rcut = args[3], p = args[4]
   #       f = @analytic r -> ( ((rnn/r)^p - (rnn/rcut)^p)^2 * ((rnn/r)^p - (rnn/rin)^p)^2 )
   #       r -> f.f(r) * (rin < r < rcut), r -> f.f_d(r) * (rin < r < rcut), rcut
   #    end

   else
      error("fcut_analyse: unknown symbol $(sym) for fcut.")
   end
end
