using StellaratorGrids
using Test
using VMEC

@testset "StellaratorGrids.jl" begin
    include("read_wall_test.jl")
    include("expand_test.jl")
end
