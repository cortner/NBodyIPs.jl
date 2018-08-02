
using JuLIP.Potentials: @analytic,
                        cutsw,
                        cutsw_d,
                        coscut,
                        coscut_d

const cutsp = JuLIP.Potentials.fcut
const cutsp_d = JuLIP.Potentials.fcut_d


import Base: Dict,
             ==

# -------------- Space Tranformations ---------------

struct SpaceTransform{FT, FDT}
   id::String
   f::FT
   df::FDT
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
   ftrans = eval(parse(ftrans_analyse(strans)))
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
   df::DFT
   rcut::Float64
end

==(C1::Cutoff, C2::Cutoff) = (C1.sym == C2.sym) && (C1.params == C2.params)

Cutoff(args...) = Cutoff(args)

function Cutoff(args::Tuple)
   fcut, rcut = fcut_analyse(args)
   return Cutoff(args[1], [args[2:end]...], fcut.f, fcut.f_d, rcut)
end

function fcut_analyse(args::Tuple)
   sym = args[1]::Symbol
   args = args[2:end]
   if Symbol(sym) == :cos
      rc1, rc2 = args
      return let rc1=rc1, rc2=rc2
                  AnalyticFunction( r -> coscut(r, rc1, rc2),
                                    r -> coscut_d(r, rc1, rc2),
                                    nothing ); end, rc2

   elseif Symbol(sym) == :sw
      L, rcut = args
      return let L=L, rcut=rcut
                  AnalyticFunction( r -> cutsw(r, rcut, L),
                                    r -> cutsw_d(r, rcut, L),
                                    nothing ); end, rcut

   elseif Symbol(sym) == :spline
      rc1, rc2 = args
      return let rc1=rc1, rc2=rc2
                  AnalyticFunction( r -> cutsp(r, rc1, rc2),
                                    r -> cutsp_d(r, rc1, rc2),
                                    nothing );end , rc2

   elseif Symbol(sym) == :square
      rcut = args[1]
      return let rcut=rcut; (@analytic r -> (r - rcut)^2); end, rcut

   elseif Symbol(sym) == :s2rat
      return let rnn=args[1], rin = args[2], rcut = args[3], p = args[4]
         f = @analytic r -> ( ((rnn/r)^p - (rnn/rcut)^p)^2 * ((rnn/r)^p - (rnn/rin)^p)^2 )
         AnalyticFunction(r -> f.f(r) * (rin < r < rcut),
                          r -> f.f_d(r) * (rin < r < rcut),
                          nothing) end, args[2]

   elseif Symbol(sym) == :cos2s
      return let ri1 = args[1], ri2 = args[2], ro1 = args[3], ro2 = args[4]
            AnalyticFunction(
                  r -> (1-coscut(r, ri1, ri2)) * coscut(r, ro1, ro2),
                  r -> (- coscut_d(r, ri1, ri2) * coscut(r, ro1, ro2)
                        + (1-coscut(r, ri1, ri2)) * coscut_d(r, ro1, ro2)),
                  nothing) end, args[4]
   else
      error("Dictionary: unknown symbol $(sym) for fcut.")
   end
end


Dict(fcut::Cutoff) = Dict( "__id__" => "Cutoff",
                           "sym" => String(fcut.sym),
                           "params" => fcut.params )

Cutoff(D::Dict) = Cutoff(Symbol(D["sym"]), D["params"]...)
