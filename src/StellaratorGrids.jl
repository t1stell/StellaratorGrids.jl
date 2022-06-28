module StellaratorGrids

using Requires
using VMEC
using Roots
using Plots
using PlasmaEquilibriumToolkit
# Write your package code here.

export expanded_wall_auto, expanded_wall_simple

include("expand_surface.jl")

end
