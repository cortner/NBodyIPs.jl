
using JSON

struct XJld2 end
struct XJson end
struct XJld end


Dict(IP::NBodyIP) = Dict("__id__" => "NBodyIP",
                         "components" => Dict.(IP.components))

_decode_dict(D::Dict) = convert(Val(Symbol(D["__id__"])), D)

NBodyIP(D::Dict) = NBodyIP(_decode_dict.(D["components"]))

Base.convert(::Val{:NBodyIP}, D::Dict) = NBodyIP(D)


function _checkextension(fname)
   if fname[end-3:end] == "jld2"
      return XJld2()
   elseif fname[end-2:end] == "jld"
      return XJld()
   elseif fname[end-3:end] == "json"
      return XJson()
   end
   error("unknown extension for `$fname`")
end

save_ip(fname::AbstractString, IP) =
   save_ip(_checkextension(fname), fname, IP)

load_ip(fname::AbstractString) =
   load_ip(_checkextension(fname), fname)

save_ip(::Union{XJld2, XJld}, fname, IP) =
   save(fname, "IP", Dict(IP))

function load_ip(::Union{XJld2, XJld}, fname)
   warn("""`load_ip` and `save_ip` do not directly support jld and jld2. In order
           to load an `NBodyIP` stored in one of those formats, please use `FileIO`,
           load the IP as a `ipdict::Dict` and then decode it using
           `NBodyIP(ipdict)`. To save an `ip::NBodyIP` as a jld or jld2 file,
           convert it to a `Dict` using `Dict(ip)` and then save it using
           `FileIO.save`.""")
end

function save_ip(::XJson, fname, IP)
   f = open(fname, "w")
   print(f, JSON.json(Dict(IP)))
   close(f)
end

function load_ip(::XJson, fname)
   IPj = JSON.parsefile(fname)
   @assert IPj["__id__"] == "NBodyIP"
   return NBodyIP(IPj)
end
