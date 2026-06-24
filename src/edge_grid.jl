# Edge / boundary point generation.  <- src/edge_grid.f90
#
# `get_growth`, `xi_length` and `interp_spline` live in math_deps.jl (they are
# shared with the solvers). This file holds the boundary generators: the C-grid
# wake cut, the blunt trailing-edge fillet, and the TE-fillet point spacing.

"""
    normal_spacing(l, i, N, d0) -> width

Width of the `i`-th of `N` intervals spanning length `l`, distributed by a
Gaussian "bump" (denser in the middle) with floor `d0`.
<- edge_grid.f90 :: normal_spacing
"""
function normal_spacing(l, i::Integer, N::Integer, d0)
    mu = 0.5
    sig = sqrt(0.025)
    Gmid = normal_dist(0.5, sig, mu)
    G0 = normal_dist(0.0, sig, mu)
    sumg = 0.0
    for j in 1:N
        sumg += _g_func_norm(j, N, sig, mu, G0, Gmid)
    end
    d1 = (l - N * d0) / sumg + d0
    g = _g_func_norm(i, N, sig, mu, G0, Gmid)
    return d0 + g * (d1 - d0)
end

"""
    fillet_trailing_edge!(foil, nte, surfbounds)

Round a blunt trailing edge by inserting `nte + 1` points along a degree-3
B-spline filleted across the TE gap, rebuilding `foil` in place and recording the
top/bottom corner indices. <- edge_grid.f90 :: fillet_trailing_edge
"""
function fillet_trailing_edge!(foil::AirfoilSurface, nte::Integer, surfbounds::AbstractVector{<:Integer})
    npoints = foil.npoints
    fx = foil.x; fy = foil.y

    # 7 control points filleting the TE and removing the sharp corners.
    CP = zeros(Float64, 2, 7)
    CP[1, 1] = 0.75 * fx[npoints-1] + 0.25 * fx[npoints]
    CP[2, 1] = 0.75 * fy[npoints-1] + 0.25 * fy[npoints]
    CP[1, 2] = 0.5 * (CP[1, 1] + fx[npoints])
    CP[2, 2] = 0.5 * (CP[2, 1] + fy[npoints])
    CP[1, 4] = 0.5 * (fx[npoints] + fx[1])
    CP[2, 4] = 0.5 * (fy[npoints] + fy[1])
    CP[1, 3] = 0.5 * (fx[npoints] + CP[1, 4])
    CP[2, 3] = 0.5 * (fy[npoints] + CP[2, 4])
    CP[1, 5] = 0.5 * (CP[1, 4] + fx[1])
    CP[2, 5] = 0.5 * (CP[2, 4] + fy[1])
    CP[1, 7] = 0.25 * fx[1] + 0.75 * fx[2]
    CP[2, 7] = 0.25 * fy[1] + 0.75 * fy[2]
    CP[1, 6] = 0.5 * (fx[1] + CP[1, 7])
    CP[2, 6] = 0.5 * (fy[1] + CP[2, 7])

    sx, sy, _ = bspline(CP, 3, 100)
    cdfs = polyline_dist(sx, sy)

    # Spacing along the fillet (slightly tighter than uniform).
    d0 = cdfs[100] / (nte + 3) / 1.05
    cdfte = zeros(Float64, nte + 4)
    for i in 1:nte+2
        space = normal_spacing(cdfs[100], i, nte + 3, d0)
        cdfte[i+1] = cdfte[i] + space
    end
    cdfte[nte+4] = cdfs[100]

    tex = zeros(Float64, nte + 4)
    tey = zeros(Float64, nte + 4)
    tex[1] = sx[1]; tey[1] = sy[1]
    pt1 = 1
    for i in 2:nte+4
        tex[i], tey[i], pt1 = interp_spline(cdfte[i], sx, sy, cdfs, pt1)
    end

    # Point that becomes i = 1 (mid-TE), and counts above/below it.
    firstpt = floor(Int, (nte + 4) / 2) + 1
    ntetop = nte + 2 - firstpt + 1
    ntebot = firstpt - 2

    foil.topcorner = ntetop + surfbounds[1]
    foil.botcorner = surfbounds[2] - ntebot

    # Rebuild the surface: TE-top fillet, original interior, TE-bottom fillet.
    newn = npoints + nte + 1
    nx = zeros(Float64, newn)
    ny = zeros(Float64, newn)
    nx[1:ntetop+1] = tex[firstpt:nte+3]
    ny[1:ntetop+1] = tey[firstpt:nte+3]
    nx[ntetop+2:newn-ntebot-1] = fx[2:npoints-1]
    ny[ntetop+2:newn-ntebot-1] = fy[2:npoints-1]
    nx[newn-ntebot:newn] = tex[2:firstpt]
    ny[newn-ntebot:newn] = tey[2:firstpt]

    foil.npoints = newn
    foil.x = nx
    foil.y = ny
    return foil
end

"""
    add_wake_points!(grid, options)

Lay out the C-grid wake cut: a downstream boundary point at `x = radi + 0.5`,
TE-matched initial spacing `d0`, and geometric growth along the TE→back vector,
mirrored onto both wake-cut edges (`srf1-i` and `srf2+i`).
<- edge_grid.f90 :: add_wake_points
"""
function add_wake_points!(grid::SrfGrid, options::MeshOptions)
    grid.x[1, 1] = options.radi + 0.5
    grid.y[1, 1] = 0.0
    grid.x[grid.imax, 1] = grid.x[1, 1]
    grid.y[grid.imax, 1] = grid.y[1, 1]

    srf1 = grid.surfbounds[1]
    srf2 = grid.surfbounds[2]
    d0 = sqrt((grid.x[srf1+1, 1] - grid.x[srf1, 1])^2 +
              (grid.y[srf1+1, 1] - grid.y[srf1, 1])^2)

    wx = grid.x[1, 1] - grid.x[srf1, 1]
    wy = grid.y[1, 1] - grid.y[srf1, 1]
    wlen = sqrt(wx^2 + wy^2)
    wx /= wlen; wy /= wlen

    g1 = get_growth(wlen, d0, options.nwake)

    space = d0
    for i in 1:options.nwake-1
        space = i == 1 ? d0 : space * g1
        grid.x[srf1-i, 1] = grid.x[srf1-i+1, 1] + wx * space
        grid.y[srf1-i, 1] = grid.y[srf1-i+1, 1] + wy * space
        grid.x[srf2+i, 1] = grid.x[srf1-i, 1]
        grid.y[srf2+i, 1] = grid.y[srf1-i, 1]
    end
    return grid
end
