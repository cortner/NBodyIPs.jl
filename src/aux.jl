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

import JuLIP: cutoff

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

cutoff(C::Cutoff) = C.rcut

==(C1::Cutoff, C2::Cutoff) = (C1.sym == C2.sym) && (C1.params == C2.params)

Cutoff(descr::String) = Cutoff(eval(parse(descr)))

Cutoff(args...) = Cutoff(args)

function Cutoff(args::Tuple)
   f, f_d, rcut = fcut_analyse(args)
   return Cutoff(args[1], [args[2:end]...], f, f_d, rcut)
end

function fcut_analyse(args::Tuple)
   sym = args[1]::Symbol
   args = args[2:end]
   if Symbol(sym) == :cos
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


fcut(C::Cutoff, r::Number) = C.f(r)
fcut_d(C::Cutoff, r::Number) = C.f_d(r)


# ---------------------------------------------------------------
# NOT TECHNICALLY MONOMIALS, BUT VERY SIMILAR OBJECT:
#   TODO: combine fcut_d and monomial_d into a single function
# ---------------------------------------------------------------

fcut(C::Cutoff, r::SVector) = prod(C.f, r)

fcut_d(C::Cutoff, r::SVector{1}) = C.f(r[1]), C.f_d(r[1])

function fcut_d(C::Cutoff, r::SVector{3})
   f = C.f.(r)
   f_d = C.f_d.(r)
   f12 = f[1]*f[2]
   f23 = f[2]*f[3]
   return f12*f[3], @SVector [f_d[1]*f23, f[1]*f_d[2]*f[3], f12*f_d[3]]
end

function fcut_d(C::Cutoff, r::SVector{6})
   # f = C.f.(r)
   # f_d = C.f_d.(r)
   f1 = C.f(r[1])
   f2 = C.f(r[2])
   f3 = C.f(r[3])
   f4 = C.f(r[4])
   f5 = C.f(r[5])
   f6 = C.f(r[6])
   f12 = f1*f2
   f123 = f12*f3
   f56 = f5*f6
   f456 = f4*f56
   f3456 = f3*f456
   f1234 = f123*f4
   return f123*f456,
          @SVector [C.f_d(r[1])*f2*f3456,
                    C.f_d(r[2])*f1*f3456,
                    C.f_d(r[3])*f12*f456,
                    C.f_d(r[4])*f123*f56,
                    C.f_d(r[5])*f1234*f6,
                    C.f_d(r[6])*f1234*f5]
end



# @generated function fcut_d(D::Cutoff, r::SVector{M, T}) where {M, T}
#    exprs = Expr[]
#
#    # evaluate the scalar monomials
#    #  f = SVector{...}( x[1]^α[1], ...)
#    # df = SVector{...}( α[1] * x[1]^(α[1]-1), ...)
#    ex_f =   "f = @SVector $T["
#    ex_df = "df = @SVector $T["
#    for i = 1:M
#       ex_f  *=   "fcut(D, r[$i]), "
#       ex_df *= "fcut_d(D, r[$i]), "
#    end
#    ex_f  =  ex_f[1:end-2] * "]"
#    ex_df = ex_df[1:end-2] * "]"
#    push_str!(exprs, ex_f)
#    push_str!(exprs, ex_df)
#
#    # evaluate the multi-variate monomial
#    push_str!(exprs, "fc = 1.0")
#    for i = 1:M
#       push_str!(exprs, "fc *= f[$i]")
#    end
#
#    for i = 1:M
#       push_str!(exprs, "fc_d_$i = 1.0")
#       for j = 1:M
#          if i == j
#             push_str!(exprs, "fc_d_$i *= df[$j]")
#          else
#             push_str!(exprs, "fc_d_$i *= f[$j]")
#          end
#       end
#    end
#
#    # evaluate the derivative
#    ex_dm = "fc_d = @SVector $T["
#    for i = 1:M
#       ex_dm *= "fc_d_$i, "
#    end
#    ex_dm = ex_dm[1:end-2] * "]"
#    push_str!(exprs, ex_dm)
#
#    quote
#       $(Expr(:meta, :inline))
#       @inbounds $(Expr(:block, exprs...))
#       return fc, fc_d
#    end
# end
