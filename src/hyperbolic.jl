# Hyperbolic surface grid generation.  <- src/hyperbolic_surface_grid.f90
#
# Grows the volume grid outward from the airfoil surface by implicit hyperbolic
# marching: at each η-level it solves a (nearly block-tridiagonal) linear system
# enforcing cell area and orthogonality, with optional area smoothing/blending.
#
# The Fortran kept `cell_sizes` as a module variable allocated at j == 1; here it
# is computed once in `hyperbolic_grid!` and threaded through as an argument
# (same values, no hidden state). The dense solve `lmult` becomes `\`.

"""
    hyperbolic_grid!(grid, options) -> y0

Generate the hyperbolic volume grid in place (surface row `j = 1` must already be
set). Marches all `jmax-1` levels, then applies y+-based normal spacing. Returns
the first-layer wall distance `y0`. <- hyperbolic_surface_grid.f90 :: hyperbolic_grid
"""
function hyperbolic_grid!(grid::SrfGrid, options::MeshOptions)
    imax = grid.imax
    jmax = grid.jmax

    # Cell sizes in the normal direction (geometric growth to the farfield),
    # computed once. <- the j == 1 branch of specify_hyperbolic_area
    lenscale = 1.2
    boundlen = options.radi * lenscale
    y0 = wall_distance(options.yplus, options.Re, options.cfrac)
    ga = get_growth(boundlen, y0, jmax - 1)
    cell_sizes = zeros(Float64, jmax - 1)
    cell_sizes[1] = y0
    for jj in 2:jmax-1
        cell_sizes[jj] = ga * cell_sizes[jj-1]
    end

    area = zeros(Float64, imax)
    for j in 1:jmax-1
        specify_hyperbolic_area!(area, grid, options, j, cell_sizes)
        known_inverse_metrics!(grid, area, imax, j, options.topology)
        LHS, RHS = hyperbolic_system(grid, options, area, j)
        solve_hyperbolic_system!(grid, LHS, RHS, j, imax, options.topology)
    end

    apply_normal_spacing!(grid, options)
    return y0
end

# Specify the target cell area at the next level.
# <- hyperbolic_surface_grid.f90 :: specify_hyperbolic_area
function specify_hyperbolic_area!(area, grid::SrfGrid, options::MeshOptions, j::Integer, cell_sizes)
    imax = grid.imax
    jmax = grid.jmax
    smth_fact = 0.16

    # Surface normals of the current level; remember the wall normals.
    normals = surface_normals(grid, options, j)
    if j == 1
        grid.surfnorm .= normals
    end

    # Offset surface one cell out along the normals.
    xoff = zeros(Float64, imax); yoff = zeros(Float64, imax)
    for i in 1:imax
        xoff[i] = grid.x[i, j] + normals[i, 1] * cell_sizes[j]
        yoff[i] = grid.y[i, j] + normals[i, 2] * cell_sizes[j]
    end

    # Arc lengths of the offset and current surfaces.
    cdfoff = zeros(Float64, imax); cdfj = zeros(Float64, imax)
    for i in 2:imax
        cdfoff[i] = cdfoff[i-1] + sqrt((xoff[i] - xoff[i-1])^2 + (yoff[i] - yoff[i-1])^2)
        cdfj[i] = cdfj[i-1] + sqrt((grid.x[i, j] - grid.x[i-1, j])^2 +
                                   (grid.y[i, j] - grid.y[i-1, j])^2)
    end

    # Interpolate the current-surface spacing onto the offset surface.
    cdfj1 = cdfoff[imax] / cdfj[imax] .* cdfj
    xj1 = zeros(Float64, imax); yj1 = zeros(Float64, imax)
    xj1[1] = xoff[1]; yj1[1] = yoff[1]
    xj1[imax] = xoff[imax]; yj1[imax] = yoff[imax]
    pt1 = 1
    for i in 2:imax-1
        xj1[i], yj1[i], pt1 = interp_spline(cdfj1[i], xoff, yoff, cdfoff, pt1)
    end

    # Cell areas from the Jacobian of the offset map.
    for i in 1:imax
        if i == 1
            if options.topology == "OGRD"
                xz = derv1c(grid.x[i+1, j], grid.x[imax-1, j])
                yz = derv1c(grid.y[i+1, j], grid.y[imax-1, j])
            else
                xz = derv1f(grid.x[i+1, j], grid.x[i, j])
                yz = derv1f(grid.y[i+1, j], grid.y[i, j])
            end
        elseif i == imax
            if options.topology == "OGRD"
                xz = derv1c(grid.x[2, j], grid.x[i-1, j])
                yz = derv1c(grid.y[2, j], grid.y[i-1, j])
            else
                xz = derv1b(grid.x[i-1, j], grid.x[i, j])
                yz = derv1b(grid.y[i-1, j], grid.y[i, j])
            end
        else
            xz = derv1c(grid.x[i+1, j], grid.x[i-1, j])
            yz = derv1c(grid.y[i+1, j], grid.y[i-1, j])
        end
        if j == 1
            xn = derv1f(xj1[i], grid.x[i, j])
            yn = derv1f(yj1[i], grid.y[i, j])
        else
            xn = derv1c(xj1[i], grid.x[i, j-1])
            yn = derv1c(yj1[i], grid.y[i, j-1])
        end
        area[i] = xz * yn - xn * yz
    end

    # Local area smoothing.
    area_smth = similar(area)
    for _ in 1:options.asmt
        for i in 1:imax
            if options.topology == "OGRD"
                if i == 1
                    area_smth[i] = (1.0 - smth_fact) * area[i] +
                                   0.5 * smth_fact * (area[i+1] + area[imax-1])
                elseif i == imax
                    area_smth[i] = (1.0 - smth_fact) * area[i] +
                                   0.5 * smth_fact * (area[2] + area[i-1])
                else
                    area_smth[i] = (1.0 - smth_fact) * area[i] +
                                   0.5 * smth_fact * (area[i+1] + area[i-1])
                end
            else # CGRD
                if i == 1 || i == imax
                    area_smth[i] = area[i]
                else
                    area_smth[i] = (1.0 - smth_fact) * area[i] +
                                   0.5 * smth_fact * (area[i+1] + area[i-1])
                end
            end
        end
        area .= area_smth
    end

    # Blend clustered (near-wall) and uniform (near-farfield) areas.
    area_scale = 0.5 * (sin((j - 2) / (jmax - 2) * pi - 0.5 * pi) + 1.0)
    blend = options.funi * area_scale
    uniform_area = sum(area) / imax
    @. area = blend * uniform_area + (1.0 - blend) * area
    return area
end

# Inverse metrics (xz, yz, xn, yn) at the known level.
# <- hyperbolic_surface_grid.f90 :: known_inverse_metrics
function known_inverse_metrics!(grid::SrfGrid, area, imax::Integer, j::Integer, topology::AbstractString)
    for i in 1:imax-1
        if i == 1
            if topology == "OGRD"
                grid.xz[i, j] = derv1c(grid.x[i+1, j], grid.x[imax-1, j])
                grid.yz[i, j] = derv1c(grid.y[i+1, j], grid.y[imax-1, j])
            else
                grid.xz[i, j] = derv1f(grid.x[i+1, j], grid.x[i, j])
                grid.yz[i, j] = derv1f(grid.y[i+1, j], grid.y[i, j])
            end
        else
            grid.xz[i, j] = derv1c(grid.x[i+1, j], grid.x[i-1, j])
            grid.yz[i, j] = derv1c(grid.y[i+1, j], grid.y[i-1, j])
        end
        denom = grid.xz[i, j]^2 + grid.yz[i, j]^2
        grid.xn[i, j] = -grid.yz[i, j] * area[i] / denom
        grid.yn[i, j] =  grid.xz[i, j] * area[i] / denom
    end
    return grid
end

# Assemble the (nearly block-tridiagonal) system for the next level.
# <- hyperbolic_surface_grid.f90 :: hyperbolic_system
function hyperbolic_system(grid::SrfGrid, options::MeshOptions, area, j::Integer)
    imax = grid.imax
    iimax = 2 * imax - 2
    eps_scale = (j - 2) / (grid.jmax - 2)
    epsi = eps_scale * options.epsi
    epse = eps_scale * options.epse
    alfa = options.alfa

    LHS = zeros(Float64, iimax, iimax)
    RHS = zeros(Float64, iimax)

    for i in 1:imax-1
        iindex = 2 * (i - 1) + 1

        A = [grid.xn[i, j]  grid.yn[i, j];
             grid.yn[i, j] -grid.xn[i, j]]
        B = [grid.xz[i, j]  grid.yz[i, j];
            -grid.yz[i, j]  grid.xz[i, j]]
        Binv = inv(B)
        C = 0.5 * alfa * (Binv * A)

        Bl1 = [-epsi - C[1, 1]   -C[1, 2];
                -C[2, 1]         -epsi - C[2, 2]]
        Bl2 = [1.0 + 2.0 * epsi  0.0;
               0.0               1.0 + 2.0 * epsi]
        Bl3 = [-epsi + C[1, 1]    C[1, 2];
                C[2, 1]          -epsi + C[2, 2]]

        fvec = [0.0, area[i]]
        rhsvec1 = [grid.x[i, j], grid.y[i, j]]
        rhsvec2 = [0.0, 0.0]
        rhsvec3 = [0.0, 0.0]

        if options.topology == "OGRD"
            if i == 1
                LHS[iindex:iindex+1, iimax-1:iimax]   = Bl1
                LHS[iindex:iindex+1, iindex+2:iindex+3] = Bl3
                rhsvec2 = [grid.x[i+1, j] - 2grid.x[i, j] + grid.x[imax-1, j],
                           grid.y[i+1, j] - 2grid.y[i, j] + grid.y[imax-1, j]]
                rhsvec3 = [grid.x[i+1, j] - grid.x[imax-1, j],
                           grid.y[i+1, j] - grid.y[imax-1, j]]
            elseif i == imax - 1
                LHS[iindex:iindex+1, iindex-2:iindex-1] = Bl1
                LHS[iindex:iindex+1, 1:2]               = Bl3
                rhsvec2 = [grid.x[1, j] - 2grid.x[i, j] + grid.x[i-1, j],
                           grid.y[1, j] - 2grid.y[i, j] + grid.y[i-1, j]]
                rhsvec3 = [grid.x[1, j] - grid.x[i-1, j],
                           grid.y[1, j] - grid.y[i-1, j]]
            else
                LHS[iindex:iindex+1, iindex-2:iindex-1] = Bl1
                LHS[iindex:iindex+1, iindex+2:iindex+3] = Bl3
                rhsvec2 = [grid.x[i+1, j] - 2grid.x[i, j] + grid.x[i-1, j],
                           grid.y[i+1, j] - 2grid.y[i, j] + grid.y[i-1, j]]
                rhsvec3 = [grid.x[i+1, j] - grid.x[i-1, j],
                           grid.y[i+1, j] - grid.y[i-1, j]]
            end
        else # CGRD
            if i == 1
                # Constant-plane boundary: Δx = 0, Δy(1) = Δy(2).
                Bl2 = [1.0 0.0; 0.0 1.0]
                Bl3 = [-1.0 0.0; 0.0 -1.0]
                LHS[iindex:iindex+1, iindex+2:iindex+3] = Bl3
                rhsvec1 = [grid.x[i, j] - grid.x[2, j], grid.y[i, j] - grid.y[2, j]]
                rhsvec2 = [0.0, 0.0]; rhsvec3 = [0.0, 0.0]; fvec = [0.0, 0.0]
            elseif i == imax - 1
                LHS[iindex:iindex+1, iindex-2:iindex-1] = Bl1
                # y(imax) = -y(1): flip the y-column of Bl3.
                Bl3[1, 2] = -Bl3[1, 2]
                Bl3[2, 2] = -Bl3[2, 2]
                LHS[iindex:iindex+1, 1:2] = Bl3
                rhsvec2 = [ grid.x[1, j] - 2grid.x[i, j] + grid.x[i-1, j],
                           -grid.y[1, j] - 2grid.y[i, j] + grid.y[i-1, j]]
                rhsvec3 = [ grid.x[1, j] - grid.x[i-1, j],
                           -grid.y[1, j] - grid.y[i-1, j]]
            else
                LHS[iindex:iindex+1, iindex-2:iindex-1] = Bl1
                LHS[iindex:iindex+1, iindex+2:iindex+3] = Bl3
                rhsvec2 = [grid.x[i+1, j] - 2grid.x[i, j] + grid.x[i-1, j],
                           grid.y[i+1, j] - 2grid.y[i, j] + grid.y[i-1, j]]
                rhsvec3 = [grid.x[i+1, j] - grid.x[i-1, j],
                           grid.y[i+1, j] - grid.y[i-1, j]]
            end
        end

        LHS[iindex:iindex+1, iindex:iindex+1] = Bl2
        RHS[iindex:iindex+1] = rhsvec1 - (epsi + epse) * rhsvec2 + C * rhsvec3 + Binv * fvec
    end
    return LHS, RHS
end

# Solve the system and write the next level back into the grid.
# <- hyperbolic_surface_grid.f90 :: solve_hyperbolic_system
function solve_hyperbolic_system!(grid::SrfGrid, LHS, RHS, j::Integer, imax::Integer, topology::AbstractString)
    slnvec = LHS \ RHS
    for i in 1:imax-1
        iindex = 2 * (i - 1) + 1
        grid.x[i, j+1] = slnvec[iindex]
        grid.y[i, j+1] = slnvec[iindex+1]
    end
    grid.x[imax, j+1] = grid.x[1, j+1]
    if topology == "OGRD"
        grid.y[imax, j+1] = grid.y[1, j+1]
    else
        grid.y[imax, j+1] = -grid.y[1, j+1]
    end
    return grid
end
