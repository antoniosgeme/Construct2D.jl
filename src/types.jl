# Data structures.  <- src/vardef.f90 :: module vardef
#
# Fortran is 1-based and column-major and so is Julia, so the `grid%x(i,j)`
# indexing carries over unchanged to `grid.x[i,j]`. Fortran pointers/allocatables
# become Julia arrays; component names are preserved.

"""
    AirfoilSurface  <- airfoil_surface_type

Loaded airfoil surface points and trailing-edge metadata. `tegap` is `true` for a
blunt (open) trailing edge.
"""
mutable struct AirfoilSurface
    npoints::Int
    x::Vector{Float64}
    y::Vector{Float64}
    tegap::Bool
    topcorner::Int
    botcorner::Int
end

AirfoilSurface(x::AbstractVector, y::AbstractVector, tegap::Bool) =
    AirfoilSurface(length(x), collect(float.(x)), collect(float.(y)), tegap, 0, 0)

"""
    MeshOptions  <- options_type

Fully-resolved grid-generation parameters. End users normally do not build this
directly — they pass a [`GridOptions`](@ref) to [`mesh_airfoil`](@ref), which is
merged with geometry-aware defaults (`resolve_options`, a port of
`menu.f90::set_defaults`) to produce a `MeshOptions`.
"""
Base.@kwdef mutable struct MeshOptions
    project_name::String = "airfoil"
    imax::Int = 0
    jmax::Int = 100
    nsrf::Int = 0
    nsrfdefault::Int = 250
    lesp::Float64 = 0.0
    tesp::Float64 = 0.0
    nte::Int = 0
    ntedefault::Int = 13
    yplus::Float64 = 0.9
    Re::Float64 = 1.0e6
    cfrac::Float64 = 0.5
    maxsteps::Int = 1000
    fsteps::Int = 20
    radi::Float64 = 15.0
    slvr::String = "HYPR"          # "HYPR" or "ELLP"
    topology::String = "CGRD"      # "OGRD" or "CGRD"
    nwake::Int = 50
    nrmt::Int = 1
    nrmb::Int = 1
    fwkl::Float64 = 1.0
    fwki::Float64 = 10.0
    fdst::Float64 = 1.0
    griddim::Int = 2
    nplanes::Int = 2
    plane_delta::Float64 = 1.0
    alfa::Float64 = 1.0
    epsi::Float64 = 15.0
    epse::Float64 = 0.0
    funi::Float64 = 0.20
    asmt::Int = 20
    f3d_compat::Bool = false
end

"""
    SrfGrid  <- srf_grid_type

The structured grid plus the metric-derivative scratch arrays used by the
solvers. Build an empty grid with `SrfGrid(imax, jmax)` (this replaces
`memory.f90::grid_allocation`). `x[i,j]`/`y[i,j]` are the node coordinates with
`j = 1` the airfoil-surface row.
"""
mutable struct SrfGrid
    imax::Int
    jmax::Int
    x::Matrix{Float64}
    y::Matrix{Float64}
    xz::Matrix{Float64}; yz::Matrix{Float64}        # ∂/∂ξ
    xn::Matrix{Float64}; yn::Matrix{Float64}        # ∂/∂η
    xzz::Matrix{Float64}; yzz::Matrix{Float64}
    xnn::Matrix{Float64}; ynn::Matrix{Float64}
    jac::Matrix{Float64}
    surfnorm::Matrix{Float64}    # (imax, 2): unit surface normal per surface point
    surfbounds::Vector{Int}      # (2,) first/last i index of the airfoil surface
    xicut::Vector{Int}           # (2,) ξ cut boundary
    etacut::Vector{Int}          # (2,) η cut boundary
end

function SrfGrid(imax::Integer, jmax::Integer)
    z() = zeros(Float64, imax, jmax)
    return SrfGrid(imax, jmax,
                   z(), z(),
                   z(), z(), z(), z(),
                   z(), z(), z(), z(), z(),
                   zeros(Float64, imax, 2),
                   zeros(Int, 2), zeros(Int, 2), zeros(Int, 2))
end
