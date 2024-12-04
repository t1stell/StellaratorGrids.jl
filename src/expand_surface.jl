function find_ζhat(vmecSurface::VmecSurface, ζ::Float64, θ::Float64, Δ::Float64;
        min_dψ = nothing)

    function zeta_diff(ζhat::Float64)
        vmecCoords = VmecCoordinates(vmecSurface.s, θ, ζhat)
        cartCoords = CartesianFromVmec()(vmecCoords, vmecSurface)
        basis = basis_vectors(Contravariant(), CartesianFromVmec(), vmecCoords, vmecSurface)
        if min_dψ == nothing
            n = Δ*basis[:,1]/abs(basis, 1)
        else
            n = Δ*basis[:,1]*min_dψ/abs(basis,1)^2
        end
        newCoords = cartCoords + n
        newζ = atan(newCoords[2], newCoords[1])
        return newζ - ζ
    end
    return mod(find_zero(zeta_diff, ζ),2π)
end
"""
function expand_surface(vmecSurface::VmecSurface, ζ::Float64, Δ::Float64;
                       θres = 100)
  θRange = range(0, 2π, θres)
  return expand_surface(vmecSurface, ζ, Δ, θRange)
end

function expand_surface(vmecSurface::VmecSurface, ζ::Float64, Δ::Float64,
                        θRange = Union{T, AbstractRange}) where {T}
  r = Vector{Float64}(undef, length(θRange))
  z = similar(r)


  for (i,θ) in enumerate(θRange)
    ζhat = find_ζhat(vmecSurface, ζ, θ, Δ)
    vmecCoords = VmecCoordinates(vmecSurface.s, θ, ζhat)
    cartCoords = CartesianFromVmec()(vmecCoords, vmecSurface)
    basis = basis_vectors(Contravariant(), CartesianFromVmec(), vmecCoords, vmecSurface)
    n = Δ*basis[:,1]/abs(basis, 1)
    newCoords = cartCoords + n
    r[i] = sqrt(newCoords[1]^2 + newCoords[2]^2)
    z[i] = newCoords[3]

  end
  return (r, z)
end
"""
function expand_surface(vmecSurface::VmecSurface, ζ::Float64, Δ::Float64; 
        min_dψ = nothing,
        θres = 100)
    θRange = range(0,2π, θres)
    return expand_surface(vmecSurface, ζ, Δ, θRange, min_dψ=min_dψ)
end


function expand_surface(vmecSurface::VmecSurface, ζ::Float64, Δ::Float64, 
        θRange::Union{T, AbstractRange};
        min_dψ=nothing, α=1) where {T}
    r = Vector{Float64}(undef, length(θRange))
    z = similar(r)
    for (i,θ) in enumerate(θRange)
        ζhat = find_ζhat(vmecSurface, ζ, θ, Δ, min_dψ=min_dψ)
        vmecCoords = VmecCoordinates(vmecSurface.s, θ, ζhat)
        cartCoords = CartesianFromVmec()(vmecCoords, vmecSurface)
        basis = basis_vectors(Contravariant(), CartesianFromVmec(), vmecCoords, vmecSurface)
        norm = abs(basis, 1)
        if min_dψ == nothing
            n = Δ*basis[:,1]/norm
        else
            n = Δ*(basis[:,1]/norm)*(min_dψ/norm)^α 
        end
        newCoords = cartCoords + n
        r[i] = sqrt(newCoords[1]^2 + newCoords[2]^2)
        z[i] = newCoords[3]
    end
    return (r, z)
end

function expand_surface_2d(vmecSurface::VmecSurface, ζ::Float64, Δ::Float64;
        θres = 400,  min_dψ=nothing, α=1) 
    θRange = range(0, 2π, θres)
    r = Vector{Float64}(undef, length(θRange))
    z = similar(r)
    for (i, θ) in enumerate(θRange)
        vmecCoords = VmecCoordinates(vmecSurface.s, θ, ζ)
        cartCoords = CartesianFromVmec()(vmecCoords, vmecSurface)
        basis = basis_vectors(Contravariant(), CartesianFromVmec(), vmecCoords, vmecSurface)
        dψds = basis[:,1]
        dψds_r = dψds[1] * cos(ζ) + dψds[2] * sin(ζ)
        dψds_proj = [dψds_r, dψds[3]]
        norm = sqrt(dψds_proj[1]^2 + dψds_proj[2]^2)
        if min_dψ == nothing
            n = Δ*dψds_proj/norm
        else
            n = Δ*(dψds_proj/norm)*(min_dψ/norm)^α
        end
        r[i] = sqrt(cartCoords[1]^2 + cartCoords[2]^2) + n[1]
        z[i] = cartCoords[3] + n[2]
    end
    return (r,z)
end



function get_surf(vmecSurface::VmecSurface, ζ::Float64; θres = 100)
    rsurf = Vector{Float64}(undef, θres)
    zsurf = similar(rsurf)
    for (i,θ) in enumerate(range(0, 2π, θres))
        vmecCoords = VmecCoordinates(vmecSurface.s, θ, ζ)
        cartCoords = CartesianFromVmec()(vmecCoords, vmecSurface)
        rsurf[i] = sqrt(cartCoords[1]^2 + cartCoords[2]^2)
        zsurf[i] = cartCoords[3]
    end
    return (rsurf, zsurf)
end

function get_min_dψ(vmecSurface::VmecSurface; θres = 100, ζres = 100)
    min_dψds = 1.0E6
    for θ in range(0, 2π, θres)
        for ζ in range(0, 2π/vmecSurface.nfp, ζres)
            vmecCoords = VmecCoordinates(vmecSurface.s, θ, ζ)
            basis = basis_vectors(Contravariant(), CartesianFromVmec(),vmecCoords, vmecSurface)
            dψds = abs(basis, 1)
            if dψds < min_dψds
                min_dψds = dψds
            end
        end
    end
    return min_dψds
end

function plot_expanded_surface(vmecSurface::VmecSurface, ζ::Float64, 
        Δ::Float64, θres::Int64;
        α=nothing, figName="")
    (rsurf, zsurf) = get_surf(vmecSurface, ζ, θres=θres)
    if α == nothing
        (r, z) = expand_surface(vmecSurface, ζ, Δ, θres = θres)
    else
        min_dψ = get_min_dψ(vmecSurface, θres = θres, ζres = 100)
        (r, z) = expand_surface(vmecSurface, ζ, Δ, min_dψ = min_dψ, θres = θres, α=α)
    end
    #Plots.scalefontsizes(2)
    plot(rsurf, zsurf, aspect_ratio=:equal, legendfontsize=14,label="Orig.", 
         xlabel="R (m)", ylabel="Z (m)")
    plot!(r, z, label="Expand "*string(Δ), legendfontsize=14)
    if length(figName) > 0
        savefig(figName)
    end
end

function plot_expanded_surface(vmecSurface::VmecSurface, ζ::Float64, 
        Δlist::Vector{Float64}, θres::Int64;
        α=nothing, figName="", threed=true)
    (rsurf, zsurf) = get_surf(vmecSurface, ζ, θres=θres)
    plot(rsurf, zsurf, aspect_ratio=:equal, legendfontsize=14,label="Orig.", 
         xlabel="R (m)", ylabel="Z (m)", size=(800,600), linewidth=2)
    min_dψ = get_min_dψ(vmecSurface, θres = θres, ζres = 100)

    for Δ in Δlist
        if α == nothing && threed
            (r, z) = expand_surface(vmecSurface, ζ, Δ, θres = θres)
        elseif threed
            (r, z) = expand_surface(vmecSurface, ζ, Δ, min_dψ = min_dψ, θres = θres)
        elseif α == nothing
            (r, z) = expand_surface_2d(vmecSurface, ζ, Δ, θres = θres)
        else
            (r, z) = expand_surface_2d(vmecSurface, ζ, Δ, min_dψ = min_dψ, θres = θres)
        end
        plot!(r, z, label="Expand "*string(Δ), legendfontsize=14, 
              legend=:outerright, linewidth=2)
    end
    if length(figName) > 0
        savefig(figName)
    end
end

function connecting_wall(rstart::Float64,rend::Float64,zstart::Float64,
        zend::Float64,resolution::Int)
    rdel = (rend - rstart)/(resolution+2)
    zdel = (zend - zstart)/(resolution+2)
    r = range(rstart+rdel,rend-rdel,resolution)
    z = range(zstart+zdel,zend-zdel,resolution)
    return r,z
end


"""
    function expanded_wall(vmecSurface::VmecSurface, ζRange::StepRangeLen, 
                       θRanges::Vector{StepRangeLen}, 
                       ΔDivertor::Float64, ΔWall::Float64)

Computes a wall for output to flare.  This includes both divertor plates and 
a surrounding wall.  The divertor plates, are at distance ΔDivertor
and the wall is at ΔWall. 

The θRange values provide the poloidal ranges for the divertor plates, 
which will also extend across the ζRanges.  The wall will be poloidally
closed in those regions with non divertor regions being at the wall distance


"""
function expanded_wall(vmecSurface::VmecSurface, 
        ζRange::StepRangeLen{T, R, S, L}, 
        θRanges::Vector{StepRangeLen{T, R, S, L}}, 
        ΔDivertor::Float64, ΔWall::Float64;
        wall_res = 20, uniform_div = false, 
        uniform_wall = true, wf = "wall.txt",
        cm = false,
        α_div=1, α_wall=1,
        showplot=nothing, smoothbuff = π/10,
        buffwidth = 10) where {T, R, S, L}



    if !uniform_wall || !uniform_div
        min_dψ = get_min_dψ(vmecSurface, θres = 500, ζres = 500)
    end

    #calculate the number of toroidal points
    count = 0
    np = 0
    for θRange in θRanges
        count += 1
        np += length(θRange)
    end
    np += count*wall_res
    #assume 10 points for each buffer region
    if smoothbuff > 0
        np += 2 * count * buffwidth
    end
    np += 1 #the wraparound point
    nt = length(ζRange)


    rw = nothing
    zw = nothing
    roriginal = nothing
    zoriginal = nothing

    io = open(wf,"w")
    write(io,"#Wall\n")
    s = string(nt)*"\t"*string(np)*"\t"*string(vmecSurface.nfp)*" 0.0 0.0\n"
    write(io, s)

    for (i,ζ) in enumerate(ζRange)

        #plotting
        if showplot == i
            (rsurf, zsurf) = get_surf(vmecSurface, ζ, θres = 200)
            plot(rsurf, zsurf, aspect_ratio=:equal, legend=false,
                 xlabel="R (m)", ylabel="Z (m)", size=(800,600), linewidth=2, 
                 color=:black)
        end



        ζdeg = ζ*180/π #output needs to have these in degrees
        write(io,string(ζdeg)*"\n")
        #for each divertor segment, make the divertor, then make the wall
        for (j, θRange) in enumerate(θRanges)
            #first divertor
            if uniform_div
                (r, z) = expand_surface(vmecSurface, ζ, ΔDivertor, θRange)   
            else
                (r, z) = expand_surface(vmecSurface, ζ, ΔDivertor, θRange, min_dψ = min_dψ, 
                                        α=α_div)   
            end
            if showplot == i
                plot!(r,z, linewidth=2, color=:red)
            end
            #save wraparound point for end
            if j == 1
                roriginal = r[1]
                zoriginal = z[1]
            else
                #make buffer region
                (rb, zb) = connecting_wall(rw[end], r[1], zw[end], z[1], buffwidth)

                if showplot == 1 
                    plot!(rb, zb, linewidth=2, color=:blue)
                end
                write_segment(rb,zb,io, cm=cm)
            end
            write_segment(r,z,io, cm=cm)



            #note there will be issues if buffer regions are larger than gaps...
            #
            #wall start and end
            #get startpoint
            θstart = θRange[end] + smoothbuff
            if j < length(θRanges) #another set exists, so the wall 
                #ends where the next one begins
                θend = θRanges[j+1][1]
            else #this ends where the first one started
                θend = θRanges[1][1]
            end
            θend -= smoothbuff
            #handle wraparounds
            if θend < θstart
                θend += 2π
            end
            θWallRange = range(θstart, θend, wall_res)

            if uniform_wall
                (rw, zw) = expand_surface(vmecSurface, ζ, ΔWall, θWallRange)
            else
                (rw, zw) = expand_surface(vmecSurface, ζ, ΔWall, θWallRange, min_dψ=min_dψ,
                                          α = α_wall)
            end


            #first buffer region
            (rb, zb) = connecting_wall(r[end],rw[1],z[end],zw[1],buffwidth)
            if showplot == i
                plot!(rb, zb, linewidth=2, color=:blue)
            end


            write_segment(rb,zb,io)
            write_segment(rw,zw,io)

            #write to file here

            if showplot == i
                plot!(rw, zw, linewidth=2, color=:green)
            end

            #make final buffer
            if j == length(θRanges)
                (rb, zb) = connecting_wall(rw[end],roriginal,zw[end],
                                           zoriginal,buffwidth)

                if showplot == 1
                    plot!(rb, zb, linewidth=2, color=:blue)
                end
                write_segment(rb,zb,io)
                write(io, string(roriginal)*"\t"*string(zoriginal)*"\n")
            end

        end 

    end  
    if showplot != nothing
        savefig("wallfig.png")
    end
    close(io)

end

function buffer_forward!(v::Vector{Float64},n::Int,width::Int)
    for j in 1:width
        #note we want to start one back since that's when the transition occurs
        #since i = 1 is not a valid argument, this should be ok
        i = n+j-2
        #check that we haven't wrapped around, note skip the end, we'll set 0,2π equiv in post
        if i>=length(v)
            i = i+1 - length(v)
        end
        #check if we're not supposed to adjust this value anymore
        #if it's not zero, it means we're back at a wall point
        if v[i] != 0
            continue
        end
        #set the value
        v[i] = (width - float(j))/width
    end
end

function buffer_reverse!(v::Vector{Float64},n::Int,width::Int)
    for j in 1:width
        i = n -j 
        #check wraparound
        if i < 0
            i = i + length(v)
        end
        #calculate value first
        val = (width - float(j))/width
        #if this else condition is satisfied, we've gone as far as we should
        val < v[i] ? v[i] == val : return
    end
end


"""
    function expanded_wall_auto(vmecSurface::VmecSurface, 
                       ζRange::StepRangeLen, 
                       θRanges::StepRangeLen, 
                       target_dψ::Float64,
                       ΔDivertor::Float64, ΔWall::Float64;
                       uniform_div = false, uniform_wall=true,
                       wf = "wall.txt", α_div = 1, α_wall = 1,
                       showplot=nothing, buffwidth=10)

Computes a surface that is expanded from a vmecSurface object and outputs 
a wall file with the suitable input format for FLARE. Through the use of a
target dψ/ds value, the code will attempt to place divertors only at suitable
poloidal locations, placing walls at other locations.

# Arguments
 - `vmecSurface::VmecSurface`: the VMEC surface file to use as the expansion base
 - `ζRange::StepRangeLen`: The range of toroidal angles to expand the surface at
 - `θRange::StepRangeLen`: The range of poloidal angles (calculated on the VMEC surface) to expand outwards at
 - `target_dψ::Float64`: Threshold value to distinguish between divertor and wall sections (see below)
 - `ΔDivertor::Float64`: Distance in meters for the divertor (which may be scaled depending on optional arguments)
 - `ΔWall::Float64`: Same as above, but for the wall

# Optional Arguments
 - `uniform_div = false`: If false, divertor will have an adjustable length that scales with the α_div parameter. If true, it will be uniform distance from the VmecSurface
 - `uniform_wall = true`: Same as above, but for the wall
 - `wf="wall.txt"`: Output file for the resulting surface
 - `α_div=1`: α value to be used if divertor is non-uniform (see below)
 - `α_wall=1`: Same as above, but for the wall
 - `showplot = nothing`: For diagnostic purposes, set to an integer it will save a plot to "wallfig.png" of the toroidal _index_ 
 - `buffwidth = 10`: Sets the number of poloidal index values between wall-divertor transition and provides a linear interpolation between those points

The main idea of the code is to expand a surface, creating both a wall and a 
divertor, and making use of dψ/ds to automatically distinguish whether a given
poloidal value will be a wall or divertor value. 

At each poloidal point the value of dψ/ds is calculate.  Values lower than 
`target_dψ` are designated at divertors, while larger values are walls. It
is possible to set divertors or walls everywhere by making `target_dψ` very
large or very small.

If `uniform_div` or `uniform_wall` is false, then the value to expand is locally
calculated using

Δ' = Δ * (|dψ/ds_min|/|dψ/ds|)^α

where Δ and α are user defined, dψ/ds_min is the minimum value of dψ/ds 
on the surface and |dψ/ds| is the local value

"""
function expanded_wall_auto(vmecSurface::VmecSurface, 
        ζRange::StepRangeLen{T, RR, S, L}, 
        θRange::StepRangeLen{T, RR, S, L},
        target_dψ::Float64,
        ΔDivertor::Float64, ΔWall::Float64;
        uniform_div = false, 
        uniform_wall = true, wf = nothing,
        α_div = 1, α_wall = 1,
        units = "m",
        wrapθ = true,
        showplot=nothing, buffwidth = 10) where {T, RR, S, L}


    if !uniform_wall || !uniform_div
        min_dψ = get_min_dψ(vmecSurface, θres = 500, ζres = 500)
    end
    nt = length(ζRange)
    np = length(θRange)

    if units == "cm"
        mult = 100
    else
        mult=1
    end

    #io = open(wf,"w")
    #write(io,"#Wall\n")
    #s = string(nt)*"\t"*string(np)*"\t"*string(vmecSurface.nfp)*" 0.0 0.0\n"
    #write(io, s)

    dψ = Vector{Float64}(undef, np)
    rwall = similar(dψ)
    zwall = similar(dψ)
    Rp = Array{Float64}(undef, np, nt)
    Zp = similar(Rp)
    θs = similar(Rp)
    ζs = Vector{Float64}(undef, nt)

    for (iζ,ζ) in enumerate(ζRange)
        ζdeg = ζ*180/π #output needs to have these in degrees
        ζs[iζ] = ζ
        #write(io,string(ζdeg)*"\n")

        #plotting
        if showplot == iζ
            (rsurf, zsurf) = get_surf(vmecSurface, ζ, θres = 200)
            plot(rsurf, zsurf, aspect_ratio=:equal, legend=false,
                 xlabel="R (m)", ylabel="Z (m)", size=(800,600), linewidth=2, 
                 color=:black)
        end

        #Pass one: label all points as wall (0) or divertor (1)
        for (i,θ) in enumerate(θRange)
            vmecCoords = VmecCoordinates(vmecSurface.s, θ, ζ)
            basis = basis_vectors(Contravariant(), CartesianFromVmec(),vmecCoords, vmecSurface)
            dψds = abs(basis, 1)
            dψds > target_dψ ? dψ[i] = 0 : dψ[i] = 1
        end

        #Pass two: label the buffer points
        #This will be a logical mess but let's get something that
        #works first and make it clean later
        #basic idea is to find the transition points and count forward or backward
        #and watch out for intersects, and properly handle the wraparound
        for (i,θ) in enumerate(θRange)
            if i == 1
                continue
            end
            #last and first point must be the same
            if i == length(θRange) && wrapθ
                dψ[i] = dψ[1]
                continue
            end
            if dψ[i-1] == 1 && dψ[i] == 0
                #extract this for readability
                buffer_forward!(dψ, i, buffwidth)
            elseif dψ[i-1] == 0 && dψ[i] == 1
                buffer_reverse!(dψ, i, buffwidth)
            end
        end

        #Pass 3, assign the values
        for (i,θ) in enumerate(θRange)
            if i == length(θRange) && wrapθ
                rwall[i] = rwall[1]
                zwall[i] = zwall[1]
                continue
            end
            if dψ[i] < 1

                #calculate the wall position
                if uniform_wall
                    (rw, zw) = expand_surface(vmecSurface, ζ, ΔWall, θ)
                else
                    (rw, zw) = expand_surface(vmecSurface, ζ, ΔWall, θ, min_dψ=min_dψ,
                                              α=α_wall)
                end
                rw = rw[1]
                zw = zw[1]
            else
                rw = 0
                zw = 0
            end
            if dψ[i] > 0
                #calculate the div position
                if uniform_div
                    (rd, zd) = expand_surface(vmecSurface, ζ, ΔDivertor, θ)   
                else
                    (rd, zd) = expand_surface(vmecSurface, ζ, ΔDivertor, θ, min_dψ=min_dψ,
                                              α = α_div)   
                end
                rd = rd[1]
                zd = zd[1]
            else
                rd = 0
                zd = 0
            end
            rwall[i] = dψ[i] * rd + (1-dψ[i]) * rw
            zwall[i] = dψ[i] * zd + (1-dψ[i]) * zw
        end
        if showplot == iζ
            plot!(rwall, zwall, color=:red)
        end
        #write_segment(rwall.*mult,zwall.*mult,io)
        Rp[:,iζ] = rwall[:]
        Zp[:,iζ] = zwall[:]
        θs[:,iζ] = collect(θRange)

    end
    
    R, Z = StellaratorGrids.create_wall_splines(ζRange, θRange, Rp, Zp, "none")
    fw = FlareWall(nt, np, Rp, Zp, collect(ζRange), θs, vmecSurface.nfp, R, Z, "none") 

    if showplot != nothing
        savefig("wallfig.png")
    end
    if wf != nothing
        StellaratorGrids.write_wall(fw, wf, units = units)
    end
    return fw
end


"""
    function expanded_wall_simple(vmecSurface::VmecSurface, 
                       ζRange::StepRangeLen, 
                       θRanges::StepRangeLen, 
                       Δ;
                       uniform = false,
                       wf = "wall.txt", α = 1,
                       showplot=nothing)

Simplified version of `expanded_wall_auto` that doesn't distinguish between wall 
and divertor sections.

# Arguments
 - `vmecSurface::VmecSurface`: the VMEC surface file to use as the expansion base
 - `ζRange::StepRangeLen`: The range of toroidal angles to expand the surface at
 - `θRange::StepRangeLen`: The range of poloidal angles (calculated on the VMEC surface) to expand outwards at
 - `Δ::Float64`: Distance in meters for the boundary (which may be scaled depending on optional arguments)

# Optional Arguments
 - `uniform = false`: If false, boundary will have an adjustable length that scales with the α parameter. If true, it will be uniform distance from the VmecSurface
 - `wf="wall.txt"`: Output file for the resulting surface
 - `α=1`: α value to be used if boundary is non-uniform
 - `showplot = nothing`: For diagnostic purposes, set to an integer it will save a plot to "wallfig.png" of the toroidal _index_ 
"""
function expanded_wall_simple(vmecSurface::VmecSurface, 
        ζRange::StepRangeLen{T, R, S, L}, 
        θRange::StepRangeLen{T, R, S, L},
        Δ::Float64;
        uniform = false, 
        wf = "wall.txt",
        units = "m",
        wrapθ = true,
        α = 1, showplot = nothing) where {T,R,S,L}

    if uniform == true
        expanded_wall_auto(vmecSurface, ζRange, θRange, 0.0, Δ, Δ, α_wall = α, wf=wf, 
                           showplot=showplot, units=units, wrapθ=wrapθ)
    else
        expanded_wall_auto(vmecSurface, ζRange, θRange, 1.0E10, Δ, Δ, α_div = α, wf=wf, 
                           showplot=showplot, units=unirts, wrapθ=wrapθ)

    end
end
