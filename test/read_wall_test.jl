@testset "read Stellarator grid tests" begin
    θarray = range(0, π/2, 10)
    fw = read_flare_wall(joinpath(@__DIR__, "test_wall.txt"), θarray = θarray, units="m")
    ζmin = fw.ζs[1]
    ζmax = fw.ζs[end]
    θmin = fw.θs[1,1] 
    θmax = fw.θs[end, 1]
    atol = 1.0E-6
    @testset "Test read values" begin
        @test fw.nζ == 4
        @test fw.nθ == 10
        @test fw.nfp == 3
        @test size(fw.R) == size(fw.Rp)
        @test size(fw.R) == size(fw.Zp)
        @test size(fw.R) == size(fw.Z)
        @test size(fw.R) == size(fw.θs)
        @test length(fw.ζs) == 4
        @test isapprox(ζmin, 0.0, atol=atol)   
        @test isapprox(ζmax, π/8, atol=atol)   
        @test isapprox(θmin, 0.0, atol=atol)   
        @test isapprox(θmax, π/2, atol=atol)   
    end
        
    @testset "test resample" begin
        ζsnew = range(π/16, π/8, 8)
        θsnew = range(π/5, 2π/5, 5)
        fw_new = resample_flare_wall(fw, θsnew, ζsnew, units="m")
        @test fw_new.nζ == 8
        @test fw_new.nθ == 5
        @test fw_new.nfp == 3
        @test size(fw_new.R) == size(fw_new.Rp)
        @test size(fw_new.R) == size(fw_new.Zp)
        @test size(fw_new.R) == size(fw_new.Z)
        @test size(fw_new.R) == size(fw_new.θs)
        @test length(fw_new.ζs) == 8
        @test isapprox(fw_new.ζs[1], π/16, atol=atol)   
        @test isapprox(fw_new.ζs[end], π/8, atol=atol)   
        @test isapprox(fw_new.θs[1,1], π/5, atol=atol)   
        @test isapprox(fw_new.θs[end,1], 2π/5, atol=atol)   
    end
end
