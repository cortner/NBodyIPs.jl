

# fc = fcut($desc, $rθ)

using NBodyIPs, StaticArrays, BenchmarkTools, JuLIP, Test, Profile,
      LinearAlgebra
using NBodyIPs.Polys: NBPoly, StNBPoly
using NBodyIPs: evaluate


TRANSFORM = "r -> (2.9/r)^3"
rcut = 7.5
desc = BondLengthDesc(TRANSFORM, (:cos, (0.66*rcut), rcut))
r = 9.5 * (@SVector rand(6))
x = (2.9./r).^3
@btime NBodyIPs.fcut($desc, $r)
@btime NBodyIPs.invariants($desc, $x)

@which NBodyIPs.fcut(desc, r)

@btime NBodyIPs.fcut($(desc.cutoff), $x)
@which NBodyIPs.fcut(desc.cutoff, x)
@btime NBodyIPs.fcut($(desc.cutoff), 4.567)
@which NBodyIPs.fcut(desc.cutoff, 4.567)

fff = desc.cutoff.f
@btime $fff(4.567)

f1, df1, rc = let rc1=3.456, rc2=7.5
   r -> JuLIP.Potentials.coscut(r, rc1, rc2), r -> JuLIP.Potentials.coscut_d(r, rc1, rc2), rc2
end
@btime $f1(4.567)

f1w = JuLIP.Potentials.F64fun(f1)
@btime $f1w(4.567)
@btime $(fff.obj.x)(4.567)
fffw = JuLIP.Potentials.F64fun(fff.obj.x)
@btime $fffw(4.567)

cutoff1 = NBodyIPs.Cutoff(desc.cutoff.sym,
                 desc.cutoff.params,
                 JuLIP.Potentials.F64fun(f1),
                 JuLIP.Potentials.F64fun(df1),
                 desc.cutoff.rcut)
@btime $(cutoff1.f)(4+rand())
@btime NBodyIPs.fcut($cutoff1, 4+rand())
@btime NBodyIPs.fcut($(desc.cutoff), 4+rand())
@btime $(desc.cutoff.f.obj.x)(4+rand())
@btime $(JuLIP.Potentials.F64fun(desc.cutoff.f.obj.x))(4+rand())

x = 4 .+ (@SVector rand(6)) .* 2
@btime NBodyIPs.fcut_d($cutoff1, $x)
@btime NBodyIPs.fcut_d_new($cutoff1, $x)
@btime NBodyIPs.fcut_d_new2($cutoff1, $x)

for n = 1:100
   x = 4 .+ (@SVector rand(6)) .* 2
   println(NBodyIPs.fcut(cutoff1, x), " : ",
          norm(NBodyIPs.fcut_d(cutoff1, x)[2] .-
               NBodyIPs.fcut_d_old(cutoff1, x)[2]))
end


##
# # 3-body potential
# B3 = nbpolys(3, desc, 8)
# c = rand(length(B3))
# V3 = NBPoly(B3, c, desc)
# V3sp = StNBPoly(V3)
#
# (R, ii, J) = (SArray{Tuple{3},Float64,1,3}[[1.58988, 1.58644, 1.58548], [3.16204, -0.00261097, -0.00218217], [-4.73656, 1.58054, 1.57758], [4.74344, 1.58054, 1.57758], [-3.15081, 0.00344743, 0.000898426], [-1.5723, 1.58215, 1.57658], [0.00678987, 3.16114, 0.00410215], [1.58579, -4.74039, 1.58146], [1.58579, 4.73961, 1.58146], [3.1617, 3.15795, 0.00452045], [-3.15924, 3.1666, 0.00571087], [-1.57802, -4.74068, 1.5814], [-1.57802, 4.73932, 1.5814], [0.00726609, -3.16374, 0.00326648], [1.5791, -1.57564, 1.57527], [3.15877, -3.15342, 0.00183828], [-4.73569, -1.58259, 1.58014], [4.74431, -1.58259, 1.58014], [-3.1615, -3.16062, 0.00499946], [-1.57353, -1.57909, 1.57556], [0.00953652, 0.00570647, 3.1631], [1.58575, 1.58157, -4.73906], [1.58575, 1.58157, 4.74094], [3.16764, 0.00644966, 3.16233], [-3.1569, -0.00360474, 3.15531], [-1.57745, 1.58326, -4.7442], [-1.57745, 1.58326, 4.7358], [-0.00115697, 3.16279, 3.15497], [3.1669, 3.15745, 3.16212], [-3.15976, 3.16027, 3.15724], [0.000180625, -3.158, 3.15627], [1.5873, -1.57522, -4.73615], [1.5873, -1.57522, 4.74385], [3.16105, -3.15639, 3.15749], [-3.16145, -3.15836, 3.16352], [-1.58108, -1.57531, -4.74534], [-1.58108, -1.57531, 4.73466], [0.00234916, -0.00279929, -3.15585], [1.58104, 1.5818, -1.57899], [3.16084, 0.00199752, -3.15851], [-4.73986, 1.5849, -1.58027], [4.74014, 1.5849, -1.58027], [-3.15168, 0.00143851, -3.15614], [-1.57103, 1.58238, -1.58417], [0.00348743, 3.1621, -3.15689], [1.5891, -4.74034, -1.57682], [1.5891, 4.73966, -1.57682], [3.16149, 3.16205, -3.15535], [-3.15445, 3.16553, -3.16449], [-1.57679, -4.73654, -1.585], [-1.57679, 4.74346, -1.585], [0.00439256, -3.15483, -3.15568], [1.58515, -1.58148, -1.57721], [3.16848, -3.16355, -3.15797], [-4.73615, -1.57912, -1.58197], [4.74385, -1.57912, -1.58197], [-3.15364, -3.15896, -3.16494], [-1.5789, -1.58257, -1.58277]], 1, [57, 58])
# J = SVector(J...)
#
# @btime NBodyIPs.evaluate($V3sp, $desc, $R, $ii, $J)
#
# @code_warntype NBodyIPs.evaluate(V3sp, desc, R, ii, J)
#
# # function evaluate(V::NBodyFunction{N},
# #                   desc::NBSiteDescriptor,
# #                   Rs::AbstractVector{JVec{T}},
# #                   i::Int,
# #                   J::SVector{K, Int}) where {N, T, K}
# #    # get the physical descriptor: bond-lengths (+ bond-angles)
# #    rθ = ricoords(desc, Rs, J)
# #    # check whether to skip this N-body term?
# #    skip_simplex(desc, rθ) && return zero(T)
# #    # compute the cutoff (and skip this site if the cutoff is zero)
# #    fc = fcut(desc, rθ)
# #    # @btime fcut($desc, $rθ)
# #    # @show typeof(fc), fc
# #    fc == 0 && return zero(T)
# #    # compute the invariants (this also applies the transform)
# #    II = invariants(desc, rθ)
# #    # @btime invariants($desc, $rθ)
# #    # @show typeof(II), II
# #    # exit()
# #    # evaluate the inner potential function (e.g. polynomial)
# #    return evaluate_I(V, II) * fc
# # end
#
# # @btime NBodyIPs.ricoords($desc, $R, $J) # ok
# rθ = NBodyIPs.ricoords(desc, R, J)
# # @btime NBodyIPs.skip_simplex($desc, $rθ)  # ok
# fc = NBodyIPs.fcut(desc, rθ)
# @btime NBodyIPs.fcut($desc, $rθ)
# II = NBodyIPs.invariants(desc, rθ)
# @btime NBodyIPs.invariants($desc, $rθ)
# @btime NBodyIPs.evaluate_I($V3sp, $II)
#
# x = @SVector rand(3)
# @btime NBodyIPs.invariants($desc, $x)
