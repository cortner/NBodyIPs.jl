

module Experimental

using NBodyIPs: NBodyFunction, BondLengthDesc, ClusterBLDesc,
                descriptor, bodyorder, NBodyIP, transform, cutoff
using NBodyIPs.Polys: NBPoly, StNBPoly

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




function veryfast_envpot()

end

end
