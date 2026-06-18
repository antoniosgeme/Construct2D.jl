# Grid generation options and serialization to Construct2D's `grid_options.in`
# namelist file.
#
# The variable names, namelist groups, and types below mirror exactly what
# Construct2D v2.1.5 reads in `src/menu.f90` (set_defaults):
#
#   namelist /SOPT/ nsrf, lesp, tesp, radi, nwke, fdst, fwkl, fwki
#   namelist /VOPT/ name, jmax, slvr, topo, ypls, recd, stp1, stp2, nrmt, nrmb,
#                   alfa, epsi, epse, funi, asmt, cfrc
#   namelist /OOPT/ gdim, npln, dpln, f3dm
#
# IMPORTANT: the `*_settings.nml` files shipped in the upstream `sample_airfoils`
# directory are STALE (they still use the pre-2.1.0 keys tele/bunc/tept). Do not
# model the file format on them; this serializer follows menu.f90 instead.
#
# Every field is `Union{Nothing, T}` and defaults to `nothing`. A Fortran
# namelist read only assigns variables that appear in the file and leaves the
# rest at the defaults Construct2D itself computed (which depend on nsrf and the
# trailing-edge geometry). So we deliberately emit ONLY the fields the user set,
# and let Construct2D fill in everything else.

"""
    GridOptions(; kwargs...)

Options passed to Construct2D via its `grid_options.in` namelist. Any field left
as `nothing` is omitted from the file, so Construct2D uses its own default.

# Surface grid (`&SOPT`)
- `nsrf::Int`    : number of points distributed over the airfoil surface
- `lesp::Float64`: leading-edge point spacing
- `tesp::Float64`: trailing-edge point spacing
- `radi::Float64`: farfield radius, in chords
- `nwke::Int`    : number of points in the wake (C-grid)
- `fdst`, `fwkl`, `fwki::Float64` : farfield / wake clustering factors

# Volume grid (`&VOPT`)
- `jmax::Int`    : number of points in the normal direction
- `slvr::String` : grid solver, `"HYPR"` (hyperbolic) or `"ELLP"` (elliptic)
- `topo::String` : topology, `"OGRD"` or `"CGRD"`. Leave `nothing` to let
  Construct2D pick the recommended one (O-grid for blunt TE, C-grid for sharp TE).
- `ypls::Float64`: target y+ used to set the first-layer wall spacing
- `recd::Float64`: Reynolds number used together with `ypls`
- `cfrc::Float64`: reference length as a fraction of chord for the y+ calc
- `stp1`, `stp2::Int` : elliptic-solver step counts
- `nrmt`, `nrmb::Int` : normals smoothing iterations (top/bottom)
- `alfa`, `epsi`, `epse`, `funi::Float64` : solver smoothing/uniformity controls
- `asmt::Int`    : number of airfoil-surface smoothing iterations

# Output (`&OOPT`)
- `gdim::Int`    : grid dimension, `2` or `3`
- `npln::Int`    : number of planes when `gdim == 3`
- `dpln::Float64`: spacing between planes when `gdim == 3`
- `f3dm::Bool`   : enable FUN3D compatibility output mode
"""
Base.@kwdef struct GridOptions
    # &SOPT
    nsrf::Union{Nothing,Int}     = nothing
    lesp::Union{Nothing,Float64} = nothing
    tesp::Union{Nothing,Float64} = nothing
    radi::Union{Nothing,Float64} = nothing
    nwke::Union{Nothing,Int}     = nothing
    fdst::Union{Nothing,Float64} = nothing
    fwkl::Union{Nothing,Float64} = nothing
    fwki::Union{Nothing,Float64} = nothing
    # &VOPT
    jmax::Union{Nothing,Int}     = nothing
    slvr::Union{Nothing,String}  = nothing
    topo::Union{Nothing,String}  = nothing
    ypls::Union{Nothing,Float64} = nothing
    recd::Union{Nothing,Float64} = nothing
    cfrc::Union{Nothing,Float64} = nothing
    stp1::Union{Nothing,Int}     = nothing
    stp2::Union{Nothing,Int}     = nothing
    nrmt::Union{Nothing,Int}     = nothing
    nrmb::Union{Nothing,Int}     = nothing
    alfa::Union{Nothing,Float64} = nothing
    epsi::Union{Nothing,Float64} = nothing
    epse::Union{Nothing,Float64} = nothing
    funi::Union{Nothing,Float64} = nothing
    asmt::Union{Nothing,Int}     = nothing
    # &OOPT
    gdim::Union{Nothing,Int}     = nothing
    npln::Union{Nothing,Int}     = nothing
    dpln::Union{Nothing,Float64} = nothing
    f3dm::Union{Nothing,Bool}    = nothing
end

# --- validation -------------------------------------------------------------

function _validate(o::GridOptions)
    if o.slvr !== nothing && uppercase(o.slvr) ∉ ("HYPR", "ELLP")
        throw(ArgumentError("slvr must be \"HYPR\" or \"ELLP\", got $(repr(o.slvr))"))
    end
    if o.topo !== nothing && uppercase(o.topo) ∉ ("OGRD", "CGRD")
        throw(ArgumentError("topo must be \"OGRD\" or \"CGRD\", got $(repr(o.topo))"))
    end
    if o.gdim !== nothing && o.gdim ∉ (2, 3)
        throw(ArgumentError("gdim must be 2 or 3, got $(o.gdim)"))
    end
    return o
end

# --- namelist formatting ----------------------------------------------------

_fmt(v::Bool)            = v ? ".true." : ".false."   # must precede Integer
_fmt(v::Integer)         = string(v)
_fmt(v::AbstractString)  = "'" * String(v) * "'"
function _fmt(v::AbstractFloat)
    s = @sprintf("%.12g", v)
    # Ensure it reads back as a real (Fortran namelist needs a '.' or exponent).
    return occursin(r"[.eEdD]", s) ? s : s * ".0"
end

_emit(io::IO, key, ::Nothing) = nothing
_emit(io::IO, key, val)       = println(io, "  ", key, " = ", _fmt(val))

"""
    write_grid_options(io_or_path, opts::GridOptions; name::AbstractString)

Write `opts` as a Construct2D `grid_options.in` namelist. `name` becomes the
`&VOPT name` entry, which controls the names of the output files
(`<name>.p3d`, `<name>_stats.p3d`, `<name>.nmf`).
"""
function write_grid_options(io::IO, o::GridOptions; name::AbstractString)
    _validate(o)

    println(io, "&SOPT")
    _emit(io, "nsrf", o.nsrf)
    _emit(io, "lesp", o.lesp)
    _emit(io, "tesp", o.tesp)
    _emit(io, "radi", o.radi)
    _emit(io, "nwke", o.nwke)
    _emit(io, "fdst", o.fdst)
    _emit(io, "fwkl", o.fwkl)
    _emit(io, "fwki", o.fwki)
    println(io, "/")

    println(io, "&VOPT")
    _emit(io, "name", name)
    _emit(io, "jmax", o.jmax)
    o.slvr !== nothing && _emit(io, "slvr", uppercase(o.slvr))
    o.topo !== nothing && _emit(io, "topo", uppercase(o.topo))
    _emit(io, "ypls", o.ypls)
    _emit(io, "recd", o.recd)
    _emit(io, "cfrc", o.cfrc)
    _emit(io, "stp1", o.stp1)
    _emit(io, "stp2", o.stp2)
    _emit(io, "nrmt", o.nrmt)
    _emit(io, "nrmb", o.nrmb)
    _emit(io, "alfa", o.alfa)
    _emit(io, "epsi", o.epsi)
    _emit(io, "epse", o.epse)
    _emit(io, "funi", o.funi)
    _emit(io, "asmt", o.asmt)
    println(io, "/")

    println(io, "&OOPT")
    _emit(io, "gdim", o.gdim)
    _emit(io, "npln", o.npln)
    _emit(io, "dpln", o.dpln)
    _emit(io, "f3dm", o.f3dm)
    println(io, "/")
    return nothing
end

function write_grid_options(path::AbstractString, o::GridOptions; name::AbstractString)
    open(io -> write_grid_options(io, o; name=name), path, "w")
    return path
end
