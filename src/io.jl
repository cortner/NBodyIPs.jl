
using JSON

struct XJld2 end
struct XJson end
struct XJld end



_decode_dict(D::Dict) = convert(Val(Symbol(D["__id__"])), D)

# ------------ writing and reading NBodyIP ------------

NBodyIP(D::Dict) = NBodyIP(_decode_dict.(D["components"]))

Dict(IP::NBodyIP) = Dict("__id__" => "NBodyIP",
                         "components" => Dict.(IP.components))

Base.convert(::Val{:NBodyIP}, D::Dict) = NBodyIP(D)

# ------------ writing and reading nothing ------------

Dict(::Nothing) = Dict("__id__" => "Void")

Base.convert(::Val{:Nothing}, D::Dict) = nothing

# ------------ FileIO -------------------

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

save_ip_args(IP) = IP, Dict{String, Any}()
save_ip_args(IP, info) = IP, info

"""
`save_ip(fname, IP, info)`

where `IP` must be an `NBodyIP` and `info` a `Dict`
"""
save_ip(fname::AbstractString, args...) =
   save_ip(_checkextension(fname), fname, save_ip_args(args...)...)

"""
`load_ip(fname) -> IP, info`

where `IP`is an `NBodyIP` and `info` a `Dict`
"""
load_ip(fname::AbstractString) =
   load_ip(_checkextension(fname), fname)

save_ip(::Union{XJld2, XJld}, args...) = err_jld()
load_ip(::Union{XJld2, XJld}, args...) = err_jld()

err_jld() =
   @warn("""`load_ip` and `save_ip` do not directly support jld and jld2. In order
           to load an `NBodyIP` stored in one of those formats, please use `FileIO`,
           load the IP as a `ipdict::Dict` and then decode it using
           `NBodyIP(ipdict)`. To save an `ip::NBodyIP` as a jld or jld2 file,
           convert it to a `Dict` using `Dict(ip)` and then save it using
           `FileIO.save`.""")



function save_ip(::XJson, fname, IP, info)
   f = open(fname, "w")
   print(f, JSON.json(Dict("IP" => Dict(IP), "info" => info)))
   close(f)
end

function load_ip(::XJson, fname)
   IPj = JSON.parsefile(fname)
   if haskey(IPj, "__id__")
      info("`load_ip`: Old IP filetype")
      @assert IPj["__id__"] == "NBodyIP"
      return NBodyIP(IPj), IPj["info"], Dict{String, Any}()
   end

   # new filetype
   @assert IPj["IP"]["__id__"] == "NBodyIP"
   IP = NBodyIP(IPj["IP"])
   if haskey(IPj, "info")
      info = IPj["info"]
   else
      info = Dict{String,Any}()
   end
   return IP, info
end
