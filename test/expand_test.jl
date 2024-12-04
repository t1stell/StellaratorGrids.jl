@testset "read Stellarator grid tests" begin
    θarray = range(0, π/2, 10)
    fw = read_flare_wall(joinpath(@__DIR__, "test_wall.txt"), θarray = θarray, units="m")
    ζrange = range(fw.ζs[1], fw.ζs[end], length(fw.ζs))
    ζmin = fw.ζs[1]
    ζmax = fw.ζs[end]
    θmin = fw.θs[1,1] 
    θmax = fw.θs[end, 1]
    atol = 1.0E-10

    vmec = readVmecWout("wout_circular_tokamak.nc")
    vmecsurf = VmecSurface(1.0, vmec)
    
    function check_array(a1, a2)
        @test length(a1) == length(a2)
        for i in 1:length(a1)
            @test isapprox(a1[i], a2[i], atol=atol)
        end
    end
    function check_matrix(m1, m2)
        s1 = size(m1)
        s2 = size(m2)
        @test s1 == s2
        for i in 1:s1[1]
            check_array(m1[i,:], m2[i,:])
        end
    end
            


    @testset "basic expansion" begin
        expw = expanded_wall_simple(vmecsurf, ζrange, θarray, 0.05, uniform=true, wrapθ=false)
        check_array(fw.ζs, expw.ζs)
        check_matrix(fw.θs, expw.θs)
        check_matrix(fw.Rp, expw.Rp)
        check_matrix(fw.Zp, expw.Zp)
        for ζ in range(π/16, π/8, 5)
            for θ in range(π/5, 2π/5, 7)
                @test isapprox(fw.R(θ, ζ), expw.R(θ, ζ))
                @test isapprox(fw.Z(θ, ζ), expw.Z(θ, ζ))
            end
        end
    end
end
