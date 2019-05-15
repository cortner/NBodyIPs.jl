

module Experimental

using NBodyIPs: NBodyFunction, BondLengthDesc, ClusterBLDesc,
                descriptor, bodyorder, NBodyIP, transform, cutoff,
                fast
using NBodyIPs.Polys: NBPoly, StNBPoly
using NBodyIPs.EnvIPs: EnvIP, EnvPoly
using JuLIP: AbstractCalculator

function faster_blpot(IP)
   components = AbstractCalculator[]
   for V in IP.components
      modify = false
      if (V isa NBPoly)
         if (bodyorder(V) >= 2) && (descriptor(V) isa BondLengthDesc)
            modify=true
         end
      end
      if modify
         desc_V = descriptor(V)
         desc = ClusterBLDesc(transform(desc_V), desc_V.cutoff)
         Vcl = NBPoly(V.t, V.c * bodyorder(V), desc)
         Vff = fast(Vcl)
      elseif V isa NBPoly # NBPoly but not a 3B or higher BL
         Vff = fast(V)
      else
         Vff = V
      end
      push!(components, Vff)
   end
   return NBodyIP(components)
end


function faster_envpot(IP)
   @warn("""Warning - this is an experimental code; I'm not checking
            compatibility of the individual components!""")
   components = AbstractCalculator[]
   envs = EnvIP[]
   for V in IP.components
      if !(V isa EnvIP)
         push!(components, V)
      else
         push!(envs, V)
      end
   end
   envdegs = [V.t for V in envs]           # env degrees
   @assert maximum(envdegs) <= 4    # TODO: why?!?!?
   bos = [bodyorder(V.Vr) for V in envs]   #
   for bo = 2:4  # loop over body-orders
      ibo = findall(bos .== bo)  # find the relevant subset
      isempty(ibo) && continue
      # check that _all_ envdegs are present
      @assert sort(envdegs[ibo]) == 0:(length(ibo)-1)
      Vrs = Vector{Any}(undef, length(ibo))
      for i in ibo
         Vrs[envs[i].t+1] = envs[i].Vr
      end
      # check that all Vrs have the same type
      # and also make them fast!
      Vrs = [fast(V) for V in Vrs]
      @assert isconcretetype(eltype(Vrs))
      # construct the combined EnvPoly
      Vn = envs[ibo[1]].Vn
      str_Vn = envs[ibo[1]].str_Vn
      push!(components, EnvPoly(Vrs, Vn, str_Vn))
   end
   return NBodyIP(components)
end


end
