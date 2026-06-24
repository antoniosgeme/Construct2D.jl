# Render a generated grid to a standalone SVG (no plotting dependencies), so you
# can eyeball mesh quality near the airfoil.
#
#   julia --project=. examples/03_visualize_svg.jl
#
# Draws every grid line (ξ = const and η = const) into an SVG, zoomed to a window
# around the airfoil.

using Construct2D

"Write grid `g` to `path` as an SVG, clipped to [-win, 1+win] × [-win, win]."
function grid_to_svg(path, g; win = 0.6, W = 1000)
    X, Y = g.X, g.Y
    imax, jmax = g.dims
    xlo, xhi = -win, 1.0 + win
    ylo, yhi = -win, win
    H = round(Int, W * (yhi - ylo) / (xhi - xlo))
    sx(x) = (x - xlo) / (xhi - xlo) * W
    sy(y) = H - (y - ylo) / (yhi - ylo) * H        # flip: SVG y grows downward

    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" """ *
                """viewBox="0 0 $W $H"><rect width="$W" height="$H" fill="white"/>""")
    polyline(pts, color, w) = begin
        s = join(("$(round(sx(x);digits=2)),$(round(sy(y);digits=2))" for (x, y) in pts), " ")
        println(io, """<polyline points="$s" fill="none" stroke="$color" stroke-width="$w"/>""")
    end

    # η = const lines (wrap around the airfoil)
    for j in 1:jmax
        polyline(((X[i, j], Y[i, j]) for i in 1:imax), "#3b7", 0.4)
    end
    # ξ = const lines (shoot outward); subsample for clarity
    step = max(1, imax ÷ 200)
    for i in 1:step:imax
        polyline(((X[i, j], Y[i, j]) for j in 1:jmax), "#37b", 0.4)
    end
    # airfoil surface, emphasized
    polyline(((X[i, 1], Y[i, 1]) for i in 1:imax), "#000", 1.5)

    println(io, "</svg>")
    write(path, take!(io))
    return path
end

res = mesh_airfoil(joinpath(@__DIR__, "..", "test", "fixtures", "naca0012.dat");
                   options = GridOptions(jmax = 80))
out = grid_to_svg(joinpath(@__DIR__, "naca0012_grid.svg"), res.grid)
println(res)
println("wrote ", out, " — open it in a browser to view the mesh.")
