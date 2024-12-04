struct FlareWall{FT, IT}
  nζ::IT
  nθ::IT
  Rp::Array{FT, 2}
  Zp::Array{FT, 2}
  ζs::Array{FT, 1}
  θs::Array{FT, 2}
  nfp::IT
  R::Union{Interpolations.Extrapolation, Nothing}
  Z::Union{Interpolations.Extrapolation, Nothing}
  θparam::String
end

struct StrikePointInfo{FT}
  LC_bwd::Array{FT, 1}
  LC_fwd::Array{FT, 1}
  Lpt_bwd::Array{FT, 1}
  Lpt_fwd::Array{FT, 1}
  R_bwd::Array{FT, 1}
  Z_bwd::Array{FT, 1}
  ϕ_bwd::Array{FT, 1}
  R_fwd::Array{FT, 1}
  Z_fwd::Array{FT, 1}
  ϕ_fwd::Array{FT, 1}
  θflare_bwd::Array{FT, 1}
  ζflare_bwd::Array{FT, 1}
  θflare_fwd::Array{FT, 1}
  ζflare_fwd::Array{FT, 1}
  R_start::Array{FT, 1}
  Z_start::Array{FT, 1}
  ϕ_start::Array{FT, 1}
end

"""
resample_flare_wall(flareWall, θlist, ζlist)

Resamples the flare wall given new lists of poloidal and toroidal values

"""
function resample_flare_wall(fw::FlareWall, θrange, ζrange; write_name = Nothing, units="m")
    nfp = Int64(fw.nfp)
    nζ = length(ζrange)
    nθ = length(θrange)
    ζs = Array{Float64}(undef, nζ)
    θs = Array{Float64}(undef, nθ, nζ)
    Rp = similar(θs)
    Zp = similar(θs)
    for iζ in 1:length(ζrange)
        ζs[iζ] = ζrange[iζ]
        for iθ in 1:length(θrange)
            θs[iθ, iζ] = θrange[iθ]
            println(θs[iθ, iζ], " ",ζs[iζ])
            Rp[iθ, iζ] = fw.R(θs[iθ, iζ], ζs[iζ])
            Zp[iθ, iζ] = fw.Z(θs[iθ, iζ], ζs[iζ])
        end
    end

    (R, Z) = create_wall_splines(ζrange, θrange, Rp, Zp, "none")
    return FlareWall(nζ, nθ, Rp, Zp, ζs, θs, nfp, R, Z, "none")
end

function write_wall(fw::FlareWall, wfname::String; units="m")    
    io = open(wfname,"w")
    write(io,"#Wall\n")
    s = string(fw.nζ)*"\t"*string(fw.nθ)*"\t"*string(fw.nfp)*" 0.0 0.0\n"
    write(io, s)
    for iζ in 1:fw.nζ
        write(io, "  "*string(fw.ζs[iζ]*180/π)*"\n")
        write_segment(fw.Rp[:,iζ], fw.Zp[:,iζ], io, units=units)
    end
    close(io)
end

function write_segment(rs,
                       zs,
                       io::IOStream;
                       units="m")
  mult = 1
  if units == "cm"
    mult = 100
  end

  for (i,r) in enumerate(rs)
    z = zs[i]
    s = string(r*mult)*"\t"*string(z*mult)*"\n"
    write(io, s)
  end
end

"""
This function can read a flare wall and create a wall object, also useful for EMC3 walls

Note: walls must be uniform in poloidal and toroidal directions (for now)

options:
  units: default "m".  either "m" or "cm" if "cm" will convert to m (note EMC3 uses cm, flare uses either)
  angles: default "deg". either "deg" or "rad" if "deg" will convert to "rad" (note EMC3 and flare use degrees)
  make_splines if true, splines will be created.  If not, no splines will be created
  θparam: default "none" none uses a uniform parametrization based on input
"""
function read_flare_wall(filename::String;
                         angles::String = "deg",
                         make_splines::Bool = true,
                         θparam::String = "none",
                         θarray = nothing,
                         units::String = "m")
  c = 1
  data = readlines(filename)

  #skip comment lines
  #and read tpz properly
  skip = false
  while data[c][1] == '#'
    if occursin("NODES", data[c])
      dum = split(data[c])
      nζ = parse(Int32, dum[end-1])
      nθ = parse(Int32, dum[end])
      skip = true
      nfp = nothing
    end
    c += 1
  end
  # this is the master info
  if !skip
    dum = split(data[c])
    nζ = parse(Int32, dum[1])
    nθ = parse(Int32, dum[2])
    nfp = parse(Int32, dum[3])
    c += 1
  end

  if θarray == nothing
      θuniform = range(0, 2π, nθ)
  else
      θuniform = θarray
  end
  #set up the arrays
  ζs = Array{Float64}(undef, nζ)
  θs = Array{Float64}(undef, nθ, nζ)
  Rp = similar(θs)
  Zp = similar(θs)
  for ζi in 1:nζ
    dum = split(data[c])# these should only have one value
    ζ = parse(Float64, dum[1])

    #convert from degrees
    if angles == "deg"
      ζ *= π/180
    end
    
    ζs[ζi] = ζ
    c += 1
    for θi in 1:nθ
      dum = split(data[c])
      Rp[θi, ζi] = parse(Float64, dum[1])
      Zp[θi, ζi] = parse(Float64, dum[2])
      #for now just store the uniform value, can overwrite later
      θs[θi, ζi] = θuniform[θi]
      c += 1
    end
  end

  if units == "cm"
    Rp ./= 100
    Zp ./= 100
  end

  #guess at nfp if it's not there
  if nfp == nothing
    nfp = Int32(div(2π, ζs[end]))
  end

  if make_splines == false
    return FlareWall(nζ, nθ, Rp, Zp, ζs, θs, nfp, nothing, nothing, "")
  end

  ζrange = range(ζs[1], ζs[end], nζ)
  (R, Z) = create_wall_splines(ζrange, θuniform, Rp, Zp, θparam)
  return FlareWall(nζ, nθ, Rp, Zp, ζs, θs, nfp, R, Z, θparam)

end

function create_wall_splines(ζrange, θrange, Rp, Zp, θparam)

  #right now we don't use θparam, save this for later
  knots = (θrange, ζrange)
  #do periodic in theta if the divertor goes from 0 to 2pi
  if θrange[1] == 0 && abs(θrange[end] - 2π) < 1.0E-10
      itp_types = (BSpline(Cubic(Periodic(OnGrid()))),
                   BSpline(Linear(Throw(OnGrid()))))
  else 
      itp_types = (BSpline(Cubic(Throw(OnGrid()))),
                   BSpline(Linear(Throw(OnGrid()))))
  end
  itp = (f) -> scale(interpolate(f, itp_types), knots...)
  extp = (f) -> extrapolate(itp(f), (Periodic(), Throw()))
  R = extp(Rp)
  Z = extp(Zp)
  return R, Z
end

function read_flare_strike(strike_name::String, launch_name::String;
                           units::String = "cm", angles::String = "deg")

  launch_data = read_flare_wall(launch_name, units=units, angles=angles,
                                make_splines = false)

  data = readlines(strike_name)
  c = 1
  npoints = launch_data.nζ * launch_data.nθ
  LC_bwd = Array{Float64}(undef, npoints)
  LC_fwd = similar(LC_bwd)
  Lpt_bwd = similar(LC_bwd)
  Lpt_fwd = similar(LC_bwd)
  R_bwd = similar(LC_bwd)
  Z_bwd = similar(LC_bwd)
  ϕ_bwd = similar(LC_bwd)
  R_fwd = similar(LC_bwd)
  Z_fwd = similar(LC_bwd)
  ϕ_fwd = similar(LC_bwd)
  θflare_bwd = similar(LC_bwd)
  ζflare_bwd = similar(LC_bwd)
  θflare_fwd = similar(LC_bwd)
  ζflare_fwd = similar(LC_bwd)
  R_start = similar(LC_bwd)
  Z_start = similar(LC_bwd)
  ϕ_start = similar(LC_bwd)

  while data[c][1] == '#'
    c += 1
    continue
  end

  for i in 1:npoints
    dum = split(data[c])
    LC_bwd[i] = parse(Float64, dum[1])
    LC_fwd[i] = parse(Float64, dum[2])
    Lpt_bwd[i] = parse(Float64, dum[3])
    Lpt_fwd[i] = parse(Float64, dum[4])
    R_bwd[i] = parse(Float64, dum[7])
    Z_bwd[i] = parse(Float64, dum[8])
    ϕ_bwd[i] = parse(Float64, dum[9])
    R_fwd[i] = parse(Float64, dum[10])
    Z_fwd[i] = parse(Float64, dum[11])
    ϕ_fwd[i] = parse(Float64, dum[12])
    θflare_bwd[i] = parse(Float64, dum[14])
    ζflare_bwd[i] = parse(Float64, dum[15])
    θflare_fwd[i] = parse(Float64, dum[17])
    ζflare_fwd[i] = parse(Float64, dum[18])
    R_start[i] = launch_data.Rp[i]
    Z_start[i] = launch_data.Zp[i]
    ζi = div(i-1, launch_data.nθ) + 1
    ϕ_start[i] = launch_data.ζs[ζi]
    


    c+=1
  end
  return StrikePointInfo(LC_bwd, LC_fwd, Lpt_bwd, Lpt_fwd, R_bwd, Z_bwd, ϕ_bwd,
                         R_fwd, Z_fwd, ϕ_fwd, θflare_bwd, ζflare_bwd,
                         θflare_fwd, ζflare_fwd, R_start, Z_start, ϕ_start)
end


function get_flare_strike_2d(wall::FlareWall{FT,IT}, strike::StrikePointInfo{FT}
                            ) where {FT, IT}
  nζ = wall.nζ
  nθ = wall.nθ
  nfp = wall.nfp
  gi_bwd = findall(!iszero, strike.R_bwd)
  gi_fwd = findall(!iszero, strike.R_fwd)
  θbwd = (strike.θflare_bwd[gi_bwd] ./ wall.nθ) .* 2π
  θfwd = (strike.θflare_fwd[gi_fwd] ./ wall.nθ) .* 2π
  ζbwd = (strike.ζflare_bwd[gi_bwd] ./ wall.nθ) .* (2π/nfp)
  ζfwd = (strike.ζflare_fwd[gi_fwd] ./ wall.nθ) .* (2π/nfp)

  return ζbwd, θbwd, ζfwd, θfwd
end
