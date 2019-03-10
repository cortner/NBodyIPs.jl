using NBodyIPs
using JuLIP, Test, StaticArrays

# TODO: * monomials testset
#       * fast_polys testset
#       * OneBody

@testset "NBodyIPs" begin
   @testset "Iterators" begin include("test_iterators.jl") end
   @testset "Invariants" begin include("test_invariants.jl") end
   @testset "Cutoffs" begin include("test_cutoffs.jl") end
   @testset "Polys" begin include("test_polys.jl") end
   @testset "Static Polys" begin include("test_stpolys.jl") end
   @testset "NBodyIP IO" begin include("test_io.jl") end
   @testset "EnvironmentIPs" begin include("test_environ.jl") end
end
