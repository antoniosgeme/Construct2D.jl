# Basic usage: mesh an airfoil from a coordinate file and write a Plot3D grid.
#
#   julia --project=. examples/01_basic_naca0012.jl
#
# NACA0012 here has a blunt trailing edge, so Construct2D recommends an O-grid.

using Construct2D

datfile = joinpath(@__DIR__, "..", "test", "fixtures", "naca0012.dat")

res = mesh_airfoil(datfile)

println(res)
println("topology       : ", res.options.topology)
println("grid dimensions: ", res.grid.dims, "  (imax × jmax)")
println("wall spacing y0: ", res.wall_distance)
println("farfield radius: ", round(maximum(hypot.(res.X[:, end] .- 0.5, res.Y[:, end])); digits=2))

# res.X / res.Y are imax×jmax matrices of node coordinates; res.X[:, 1] is the
# airfoil surface, res.X[:, end] the farfield boundary.
outfile = joinpath(@__DIR__, "naca0012.p3d")
write_plot3d(outfile, res.grid)
println("wrote Plot3D grid to ", outfile)
