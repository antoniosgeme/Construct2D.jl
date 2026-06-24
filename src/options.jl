# Public options and default resolution.  <- src/menu.f90 :: set_defaults / setup_airfoil_data
#
# `GridOptions` is the user-facing knob set: every field defaults to `nothing`,
# meaning "use Construct2D's geometry-aware default". `resolve_options` ports
# set_defaults: it detects the trailing-edge type, fills the defaults, applies the
# user's overrides, and computes the derived sizes (nsrf, nte, imax).

"""
    GridOptions(; kwargs...)

User-facing grid options. Any field left `nothing` falls back to Construct2D's
default (some defaults depend on the airfoil geometry). Common knobs:

- `jmax`  : points in the wall-normal direction (default 100)
- `radi`  : farfield radius in chords (default 15)
- `slvr`  : `"HYPR"` (hyperbolic, default) — `"ELLP"` not yet ported
- `topo`  : `"OGRD"` or `"CGRD"`; default is geometry-recommended
  (C-grid for a sharp TE, O-grid for a blunt TE)
- `ypls`  : target y+ for the first layer (default 0.9)
- `recd`  : chord Reynolds number used with `ypls` (default 1e6)
- `cfrc`  : reference length as a fraction of chord for the y+ calc (default 0.5)
- `nwke`  : wake points for the C-grid (default 50)
- `nte`   : points inserted to fillet a blunt TE (default 13)
- `alfa`, `epsi`, `epse`, `funi`, `asmt` : hyperbolic solver controls
"""
Base.@kwdef struct GridOptions
    jmax::Union{Nothing,Int}     = nothing
    radi::Union{Nothing,Float64} = nothing
    slvr::Union{Nothing,String}  = nothing
    topo::Union{Nothing,String}  = nothing
    ypls::Union{Nothing,Float64} = nothing
    recd::Union{Nothing,Float64} = nothing
    cfrc::Union{Nothing,Float64} = nothing
    nwke::Union{Nothing,Int}     = nothing
    nte::Union{Nothing,Int}      = nothing
    alfa::Union{Nothing,Float64} = nothing
    epsi::Union{Nothing,Float64} = nothing
    epse::Union{Nothing,Float64} = nothing
    funi::Union{Nothing,Float64} = nothing
    asmt::Union{Nothing,Int}     = nothing
end

function _validate(o::GridOptions)
    o.slvr === nothing || uppercase(o.slvr) in ("HYPR", "ELLP") ||
        throw(ArgumentError("slvr must be \"HYPR\" or \"ELLP\", got $(repr(o.slvr))"))
    o.topo === nothing || uppercase(o.topo) in ("OGRD", "CGRD") ||
        throw(ArgumentError("topo must be \"OGRD\" or \"CGRD\", got $(repr(o.topo))"))
    return o
end

# Trailing-edge classification + recommended topology.
# <- menu.f90 :: setup_airfoil_data. A sharp (closed) TE recommends C-grid; a
# blunt TE recommends O-grid. Closure is detected to a tight tolerance.
function _classify_te(x::Vector{Float64}, y::Vector{Float64})
    tol = 1.0e-12
    n = length(x)
    if abs(x[1] - x[n]) <= tol && abs(y[1] - y[n]) <= tol
        return false, "CGRD"                 # sharp
    elseif abs(y[1] - y[n]) <= tol && abs(x[1] - x[n]) > tol
        x[n] = x[1]                          # "funky geometry" fix
        return false, "CGRD"
    else
        return true, "OGRD"                  # blunt
    end
end

"""
    resolve_options(x, y; name, opts::GridOptions) -> (MeshOptions, AirfoilSurface)

Build the fully-resolved [`MeshOptions`](@ref) and [`AirfoilSurface`](@ref) from
raw coordinates and user overrides — the port of `set_defaults`. Computes the TE
type, the default solver controls, and the derived sizes (`nsrf`, `nte`, `imax`).
"""
function resolve_options(x::AbstractVector, y::AbstractVector;
                         name::AbstractString, opts::GridOptions=GridOptions())
    _validate(opts)
    x = collect(float.(x)); y = collect(float.(y))
    tegap, rec_topo = _classify_te(x, y)

    o = MeshOptions()
    o.project_name = name
    o.nsrfdefault = 250
    o.ntedefault = 13
    o.jmax = 100
    o.radi = 15.0
    o.nwake = 50
    o.yplus = 0.9
    o.Re = 1.0e6
    o.cfrac = 0.5
    o.alfa = 1.0
    o.epsi = 15.0
    o.epse = 0.0
    o.funi = 0.20
    o.asmt = 20
    o.slvr = "HYPR"
    o.topology = rec_topo

    # lesp/tesp (used only by the SMTH path) for fidelity.
    uni = 2.0 / o.nsrfdefault
    o.lesp = uni / 2.0
    o.tesp = tegap ? (y[1] - y[end]) / 10.0 : uni / 1.5

    # Apply user overrides.
    opts.jmax === nothing || (o.jmax = opts.jmax)
    opts.radi === nothing || (o.radi = opts.radi)
    opts.slvr === nothing || (o.slvr = uppercase(opts.slvr))
    opts.ypls === nothing || (o.yplus = opts.ypls)
    opts.recd === nothing || (o.Re = opts.recd)
    opts.cfrc === nothing || (o.cfrac = opts.cfrc)
    opts.nwke === nothing || (o.nwake = opts.nwke)
    opts.alfa === nothing || (o.alfa = opts.alfa)
    opts.epsi === nothing || (o.epsi = opts.epsi)
    opts.epse === nothing || (o.epse = opts.epse)
    opts.funi === nothing || (o.funi = opts.funi)
    opts.asmt === nothing || (o.asmt = opts.asmt)
    if opts.topo !== nothing
        o.topology = uppercase(opts.topo)
    end
    # C-grid uses a much lower uniformity-blend default (unless the user set it).
    if o.topology == "CGRD" && opts.funi === nothing
        o.funi = 0.01
    end

    # Derived sizes.
    o.ntedefault = something(opts.nte, 13)
    if tegap
        o.nte = o.ntedefault
        o.nsrf = length(x) + o.nte + 1
    else
        o.nte = 0
        o.nsrf = length(x)
    end
    o.imax = o.topology == "OGRD" ? o.nsrf : o.nsrf + 2 * o.nwake

    foil = AirfoilSurface(x, y, tegap)
    return o, foil
end
