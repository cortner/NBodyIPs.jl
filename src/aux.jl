import JuLIP

using JuLIP.Potentials: @analytic,
                        cutsw,
                        cutsw_d,
                        coscut,
                        coscut_d,
                        AnalyticFunction

const cutsp = JuLIP.Potentials.fcut
const cutsp_d = JuLIP.Potentials.fcut_d

import Base: Dict,
             ==

# -------------- Space Tranformations ---------------

struct SpaceTransform{FT, FDT}
   id::String
   f::FT
   f_d::FDT
end

==(T1::SpaceTransform, T2::SpaceTransform) = (T1.id == T2.id)

function SpaceTransform(strans::String)
   strans0 = strans
   # if @analytic is a substring then we don't do anything
   if !ismatch(r"@analytic", strans)
      # but if not, then we next look for ->
      if !ismatch(r"->", strans)
         # if -> is not a substring then we assume that strans is of the form
         # "(r0/r)^4" or similar i.e. explicitly uses r as the variable.
         strans = "@analytic r -> " * strans
      else
         # @analytic is not a substring but -> is a substring. e.g.
         # r -> (r0/r)^3 or s -> (r0/s)^3. we add @analytic to compile it
         strans = "@analytic " * strans
      end
   end
   ftrans = eval(parse(strans))
   return SpaceTransform(strans0, ftrans.f, ftrans.f_d)
end

Dict(t::SpaceTransform) = Dict( "__id__" => "SpaceTransform",
                                "defn" => t.id )

SpaceTransform(D::Dict) = SpaceTransform(D["defn"])




# -------------- Cut-off Mechanisms ---------------

struct Cutoff{FT, DFT}
   sym::Symbol
   params::Vector{Float64}
   f::FT
   f_d::DFT
   rcut::Float64
end

==(C1::Cutoff, C2::Cutoff) = (C1.sym == C2.sym) && (C1.params == C2.params)

Cutoff(args...) = Cutoff(args)

function Cutoff(args::Tuple)
   f, f_d, rcut = fcut_analyse(args)
   return Cutoff(args[1], [args[2:end]...], f, f_d, rcut)
end

function fcut_analyse(args::Tuple)
   sym = args[1]::Symbol
   args = args[2:end]
   if Symbol(sym) == :cos
      rc1, rc2 = args
      return let rc1=args[1], rc2=args[2]
         r -> coscut(r, rc1, rc2), r -> coscut_d(r, rc1, rc2), rc2
      end

   elseif Symbol(sym) == :sw
      return let L=args[1], rcut=args[2]
         r -> cutsw(r, rcut, L), r -> cutsw_d(r, rcut, L), rcut
      end

   elseif Symbol(sym) == :spline
      return let rc1=args[1], rc2=args[2]
         r -> cutsp(r, rc1, rc2), r -> cutsp_d(r, rc1, rc2), rc2
      end

   elseif Symbol(sym) == :square
      return let rcut=args[1]
         F = (@analytic r -> (r - rcut)^2)
         F.f, F.f_d, rcut
      end

   elseif Symbol(sym) == :s2rat
      return let rnn=args[1], rin = args[2], rcut = args[3], p = args[4]
         f = @analytic r -> ( ((rnn/r)^p - (rnn/rcut)^p)^2 * ((rnn/r)^p - (rnn/rin)^p)^2 )
         r -> f.f(r) * (rin < r < rcut), r -> f.f_d(r) * (rin < r < rcut), rcut
      end

   elseif Symbol(sym) == :cos2s
      return let ri1 = args[1], ri2 = args[2], ro1 = args[3], ro2 = args[4]
         (  r -> (1-coscut(r, ri1, ri2)) * coscut(r, ro1, ro2),
            r -> (- coscut_d(r, ri1, ri2) * coscut(r, ro1, ro2)
                  + (1-coscut(r, ri1, ri2)) * coscut_d(r, ro1, ro2)),
            ro2 )
         end
   else
      error("Dictionary: unknown symbol $(sym) for fcut.")
   end
end


Dict(fcut::Cutoff) = Dict( "__id__" => "Cutoff",
                           "sym" => String(fcut.sym),
                           "params" => fcut.params )

Cutoff(D::Dict) = Cutoff(Symbol(D["sym"]), D["params"]...)
