module StellaratorGrids

using Requires
using VMEC
using Roots
using Interpolations
using Plots
using PlasmaEquilibriumToolkit
# Write your package code here.

export expanded_wall_auto, expanded_wall_simple

#flare utils
export read_flare_strike, get_flare_strike_2d, read_flare_wall

include("expand_surface.jl")
include("FlareUtils.jl")

end
