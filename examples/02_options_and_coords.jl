# Custom options and meshing from in-memory coordinates.
#
#   julia --project=. examples/02_options_and_coords.jl
#
# Builds a sharp-trailing-edge NACA0012 directly in Julia (no file), which gets a
# C-grid, and dials in the wall spacing via a target y+ and Reynolds number.

using Construct2D

# Closed-TE NACA0012 in Selig order (TE → upper → LE → lower → TE).
function naca0012_sharp(n = 100)
    th = range(0, pi; length = n)
    xc = (1 .- cos.(th)) ./ 2                       # cosine clustering
    yt(x) = 0.6 * (0.2969sqrt(x) - 0.1260x - 0.3516x^2 + 0.2843x^3 - 0.1036x^4)
    xu = reverse(xc); yu = yt.(xu)
    xl = xc[2:end];   yl = -yt.(xl)
    x = vcat(xu, xl); y = vcat(yu, yl)
    x[1] = x[end] = 1.0; y[1] = y[end] = 0.0        # exact closure ⇒ sharp TE
    return x, y
end

x, y = naca0012_sharp()

opts = GridOptions(
    jmax = 120,        # points in the wall-normal direction
    radi = 20.0,       # farfield radius (chords)
    ypls = 0.8,        # target y+ for the first cell ...
    recd = 3.0e6,      # ... at this chord Reynolds number
    nwke = 40,         # wake points (C-grid)
)

res = mesh_airfoil((x, y); name = "naca0012_sharp", options = opts)

println(res)
println("topology       : ", res.options.topology, "  (sharp TE ⇒ C-grid)")
println("grid dimensions: ", res.grid.dims)
println("imax = nsrf + 2·nwke = ", length(x), " + 2·", res.options.nwake,
        " = ", res.options.imax)
println("wall spacing y0: ", res.wall_distance)

write_plot3d(joinpath(@__DIR__, "naca0012_sharp.p3d"), res.grid)
println("wrote ", joinpath(@__DIR__, "naca0012_sharp.p3d"))
