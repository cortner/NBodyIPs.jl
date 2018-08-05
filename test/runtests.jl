using NBodyIPs
using JuLIP, Base.Test, StaticArrays

# TODO: * monomials testset
#       * fast_polys testset

@testset "NBodyIPs" begin
   @testset "Iterators" begin include("test_iterators.jl") end
   # @testset "Invariants" begin include("test_invariants.jl") end
   # @testset "Bond-Length Polys" begin include("test_blpolys.jl") end
   # @testset "Static Bond-Length Polys" begin include("test_stblnbody.jl") end
   # @testset "NBodyIP IO" begin include("test_io.jl") end
   # @testset "EnvironmentIPs" begin include("test_environ.jl") end
end
