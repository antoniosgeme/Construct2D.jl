# Surface-grid driver.  <- src/surface_grid.f90

# --- not-yet-ported paths ---------------------------------------------------
# The XFOIL-based surface repaneling (SMTH mode) and the elliptic solver are not
# translated yet; they fail loudly rather than silently producing a wrong grid.

"`apply_foil_spacing!` (SMTH/XFOIL repaneling) is not ported yet — use BUFF mode."
apply_foil_spacing!(::AirfoilSurface, lesp, tesp) =
    error("Surface repaneling (SMTH mode) is not yet ported (needs XFOIL spline " *
          "routines from src/xfoil_deps.f90). Use surface=:buffer for now.")

"`elliptic_grid!` (ELLP solver) is not ported yet — use the hyperbolic solver."
elliptic_grid!(::SrfGrid, ::MeshOptions) =
    error("The elliptic solver (slvr=\"ELLP\") is not yet ported. Use slvr=\"HYPR\".")

"""
    create_grid!(foil, options; smooth=false) -> (grid, y0)

Top-level surface-grid driver. Allocates the grid, sets topology-dependent
surface bounds, lays in the airfoil (filleting a blunt TE and, for `CGRD`, adding
the wake cut), defines the cut boundaries, runs the solver, and stitches the
topology edges. Returns the [`SrfGrid`](@ref) and the first-layer wall distance.
<- surface_grid.f90 :: create_grid
"""
function create_grid!(foil::AirfoilSurface, options::MeshOptions; smooth::Bool=false)
    grid = SrfGrid(options.imax, options.jmax)

    # Surface bounds by topology.
    if options.topology == "OGRD"
        grid.surfbounds[1] = 1
        grid.surfbounds[2] = grid.imax
    else # CGRD: surface sits between the two wake cuts
        grid.surfbounds[1] = options.nwake + 1
        grid.surfbounds[2] = grid.imax - options.nwake
    end
    srf1 = grid.surfbounds[1]
    srf2 = grid.surfbounds[2]

    # Round a blunt trailing edge.
    if foil.tegap
        fillet_trailing_edge!(foil, options.nte, grid.surfbounds)
    end

    # Optional surface point clustering (SMTH; not yet ported).
    if smooth
        apply_foil_spacing!(foil, options.lesp, options.tesp)
    end

    # Place the airfoil into the wall row (j = 1).
    length(foil.x) == srf2 - srf1 + 1 ||
        error("airfoil has $(length(foil.x)) points but the surface needs " *
              "$(srf2 - srf1 + 1) (imax=$(grid.imax), nwake=$(options.nwake), nte=$(options.nte))")
    @views grid.x[srf1:srf2, 1] .= foil.x
    @views grid.y[srf1:srf2, 1] .= foil.y
    if !foil.tegap
        foil.topcorner = srf1 + 1
        foil.botcorner = srf2 - 1
    end

    # Wake cut for the C-grid.
    if options.topology == "CGRD"
        add_wake_points!(grid, options)
    end

    # Cut boundaries.
    grid.xicut .= (1, 0)
    grid.etacut .= (1, 0)
    if options.topology == "OGRD"
        grid.xicut[2] = grid.jmax
    else
        grid.etacut[2] = srf1 - 1
    end

    # Solver.
    y0 = if options.slvr == "HYPR"
        hyperbolic_grid!(grid, options)
    elseif options.slvr == "ELLP"
        elliptic_grid!(grid, options)
    else
        error("solver must be \"HYPR\" or \"ELLP\", got $(repr(options.slvr))")
    end

    copy_edges!(grid, options.topology)
    return grid, y0
end

"""
    copy_edges!(grid, topology) -> grid

Stitch coincident edges. `OGRD`: wrap the last ξ column onto the first. `CGRD`:
mirror the first `srf1` wall points across the branch cut.
<- surface_grid.f90 :: copy_edges
"""
function copy_edges!(grid::SrfGrid, topology::AbstractString)
    imax = grid.imax
    srf1 = grid.surfbounds[1]
    if topology == "OGRD"
        @views grid.x[imax, :] .= grid.x[1, :]
        @views grid.y[imax, :] .= grid.y[1, :]
    elseif topology == "CGRD"
        for i in 1:srf1
            grid.x[imax-i+1, 1] = grid.x[i, 1]
            grid.y[imax-i+1, 1] = grid.y[i, 1]
        end
    end
    return grid
end
