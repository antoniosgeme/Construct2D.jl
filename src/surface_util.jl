# Surface-grid utilities.  <- src/surface_util.f90
#
# Only the routines on the hyperbolic path are ported here: surface normals, the
# y+ wall-distance estimate, and the normal-direction respacing. The elliptic
# solver's `compute_inverse_metrics` / `grid_residual` are deferred with the
# elliptic solver itself.

"""
    surface_normals(grid, options, j) -> normals  (imax×2)

Outward unit normals along grid row `j`, from the surface tangent (central
difference, with topology-specific handling at the ends). Returns an `imax×2`
matrix of `(nx, ny)`.  <- surface_util.f90 :: surface_normals
"""
function surface_normals(grid::SrfGrid, options::MeshOptions, j::Integer)
    imax = grid.imax
    normals = zeros(Float64, imax, 2)
    for i in 1:imax
        if i == 1
            if options.topology == "OGRD"
                vx = grid.x[2, j] - grid.x[imax-1, j]
                vy = grid.y[2, j] - grid.y[imax-1, j]
            else
                vx = grid.x[2, j] - grid.x[1, j]
                vy = grid.y[2, j] - grid.y[1, j]
            end
        elseif i == imax
            if options.topology == "OGRD"
                vx = grid.x[2, j] - grid.x[imax-1, j]
                vy = grid.y[2, j] - grid.y[imax-1, j]
            else
                vx = grid.x[imax, j] - grid.x[imax-1, j]
                vy = grid.y[imax, j] - grid.y[imax-1, j]
            end
        else
            vx = grid.x[i+1, j] - grid.x[i-1, j]
            vy = grid.y[i+1, j] - grid.y[i-1, j]
        end
        len = sqrt(vx^2 + vy^2)
        # Negative reciprocal: outward normal for counter-clockwise points.
        normals[i, 1] = vy / len
        normals[i, 2] = -vx / len
    end
    return normals
end

"""
    wall_distance(yplus, Re, Lref) -> y0

First-layer wall spacing for a target `yplus` at chord Reynolds number `Re`,
using a turbulent flat-plate skin-friction estimate.
<- surface_util.f90 :: wall_distance
"""
function wall_distance(yplus, Re, Lref)
    Cf = (2.0 * log10(Re) - 0.65)^(-2.3)
    return yplus * Lref / (Re * sqrt(0.5 * Cf))
end

"""
    apply_normal_spacing!(grid, options) -> y0

Redistribute interior points along each ξ = const line so the first off-wall
spacing matches the y+-based `y0` and grows geometrically out to the farfield.
Returns `y0` (the first-layer wall distance).
<- surface_util.f90 :: apply_normal_spacing
"""
function apply_normal_spacing!(grid::SrfGrid, options::MeshOptions)
    jmax = grid.jmax
    y0 = wall_distance(options.yplus, options.Re, options.cfrac)
    for i in 1:grid.imax
        x = grid.x[i, :]
        y = grid.y[i, :]
        cdfs = polyline_dist(x, y)                       # arc length along η
        ga = get_growth(cdfs[jmax], y0, jmax - 1)        # growth to fit jmax pts
        cdfj = zeros(Float64, jmax)
        for j in 1:jmax-2
            cdfj[j+1] = xi_length(y0, ga, j)
        end
        cdfj[jmax] = cdfs[jmax]
        pt1 = 1
        for j in 2:jmax-1
            gx, gy, pt1 = interp_spline(cdfj[j], x, y, cdfs, pt1)
            grid.x[i, j] = gx
            grid.y[i, j] = gy
        end
    end
    return y0
end
