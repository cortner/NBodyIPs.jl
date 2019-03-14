import JuLIP

using JuLIP.Potentials: @analytic,
                        AnalyticFunction,
                        F64fun

using Roots

import Base: Dict,
             ==

# -------------- Space Tranformations ---------------


==(T1::SpaceTransform, T2::SpaceTransform) = (T1.id == T2.id)

function SpaceTransform(strans::String; fwrap = true)
   strans0 = strans
   # if @analytic is a substring then we don't do anything
   if !occursin(r"@analytic", strans)
      # but if not, then we next look for ->
      if !occursin(r"->", strans)
         # if -> is not a substring then we assume that strans is of the form
         # "(r0/r)^4" or similar i.e. explicitly uses r as the variable.
         strans = "@analytic r -> " * strans
      else
         # @analytic is not a substring but -> is a substring. e.g.
         # r -> (r0/r)^3 or s -> (r0/s)^3. we add @analytic to compile it
         strans = "@analytic " * strans
      end
   end
   ftrans = eval(Meta.parse(strans))
   if fwrap
      return SpaceTransform(strans0, F64fun(ftrans.f), F64fun(ftrans.f_d))
   end
   return SpaceTransform(strans0, ftrans.f, ftrans.f_d)
end

Dict(t::SpaceTransform) = Dict( "__id__" => "SpaceTransform",
                                "defn" => t.id )

SpaceTransform(D::Dict) = SpaceTransform(D["defn"])

hash(::BASIS, t::SpaceTransform) = hash((SpaceTransform, t.id))


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
