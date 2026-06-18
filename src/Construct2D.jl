"""
    Construct2D

A cross-platform Julia interface to [Construct2D](https://github.com/furstj/Construct2D),
a structured grid generator for 2D airfoils.

The high-level entry point is [`mesh_airfoil`](@ref), which writes an airfoil and
an options file into a working directory, drives Construct2D headlessly, and
returns the generated grid as Julia arrays. [`run_construct2d`](@ref) is the
lower-level escape hatch that just feeds menu commands to the executable.

# Locating the executable
The binary is resolved, in order, from:
1. a path set via [`set_construct2d_path!`](@ref);
2. the `CONSTRUCT2D_EXE` environment variable;
3. `Construct2D_jll` (once that artifact is installed — see the README).
"""
module Construct2D

using Printf

export GridOptions, MeshResult, Plot3DGrid,
       mesh_airfoil, run_construct2d,
       write_airfoil, read_plot3d, write_grid_options,
       set_construct2d_path!, construct2d_exe

include("options.jl")
include("io.jl")

# --- locating the construct2d executable ------------------------------------

const _EXE_OVERRIDE = Ref{String}("")

"""
    set_construct2d_path!(path)

Set an explicit path to a `construct2d` executable, taking precedence over the
`CONSTRUCT2D_EXE` environment variable and `Construct2D_jll`. Useful during
development against a locally built binary.
"""
function set_construct2d_path!(path::AbstractString)
    isfile(path) || @warn "construct2d executable not found at this path" path
    _EXE_OVERRIDE[] = String(path)
    return path
end

# Set when the optional `Construct2D_jll` provider has been registered (see the
# package extension / README "Wiring the JLL"). Holds a 0-arg function returning
# the executable path, matching the JLL's API.
const _JLL_PROVIDER = Ref{Union{Nothing,Function}}(nothing)

"""
    construct2d_exe() -> String

Return the path to the `construct2d` executable, or throw if none is configured.
"""
function construct2d_exe()
    isempty(_EXE_OVERRIDE[]) || return _EXE_OVERRIDE[]
    haskey(ENV, "CONSTRUCT2D_EXE") && return ENV["CONSTRUCT2D_EXE"]
    _JLL_PROVIDER[] === nothing || return _JLL_PROVIDER[]()
    error("""
          No `construct2d` executable is available. Provide one by either:
            • installing Construct2D_jll (see the README), or
            • setting ENV["CONSTRUCT2D_EXE"] = "/path/to/construct2d", or
            • calling Construct2D.set_construct2d_path!("/path/to/construct2d").
          """)
end

# --- Stage 1: drive the executable directly ---------------------------------

"""
    run_construct2d(workdir, airfoil_file; commands, exe=construct2d_exe())

Run `construct2d <airfoil_file>` with the working directory set to `workdir`,
feeding `commands` (a vector of menu strings) to its stdin, one per line.
Returns `(success::Bool, log::String)` where `log` is the combined stdout/stderr.

`airfoil_file` is resolved relative to `workdir`. The output files Construct2D
writes (`<name>.p3d`, etc.) also land in `workdir`.
"""
function run_construct2d(workdir::AbstractString, airfoil_file::AbstractString;
                         commands::AbstractVector{<:AbstractString},
                         exe::AbstractString=construct2d_exe())
    isdir(workdir) || throw(ArgumentError("workdir does not exist: $workdir"))
    input = isempty(commands) ? "" : join(commands, "\n") * "\n"
    cmd = Cmd(`$exe $airfoil_file`; dir=workdir)
    buf = IOBuffer()
    proc = run(pipeline(cmd; stdin=IOBuffer(input), stdout=buf, stderr=buf); wait=false)
    wait(proc)
    return (success(proc), String(take!(buf)))
end

# --- Stage 2: high-level meshing --------------------------------------------

"""
    mesh_airfoil(airfoil; options=GridOptions(), name=nothing,
                 surface=:smoothed, workdir=mktempdir(), cleanup=false) -> MeshResult

Generate a structured grid around `airfoil` and return a [`MeshResult`](@ref)
holding the parsed grid plus paths to the output files.

`airfoil` may be:
- a path to an existing coordinate `.dat` file, or
- coordinates as `(x, y)`, an `N×2` matrix, or a vector of `(x, y)` points.

Keywords:
- `options`  : a [`GridOptions`](@ref); unset fields use Construct2D's defaults.
- `name`     : project name controlling the output file names. Defaults to the
  airfoil file's base name, or `"airfoil"` for raw coordinates.
- `surface`  : `:smoothed` (XFoil-paneled, default) or `:buffer` (use the loaded
  geometry directly).
- `workdir`  : directory to run in (created fresh by default).
- `cleanup`  : if `true`, delete `workdir` after reading the grid into memory.
"""
function mesh_airfoil(airfoil;
                      options::GridOptions=GridOptions(),
                      name::Union{Nothing,AbstractString}=nothing,
                      surface::Symbol=:smoothed,
                      workdir::AbstractString=mktempdir(),
                      cleanup::Bool=false)
    surface in (:smoothed, :buffer) ||
        throw(ArgumentError("surface must be :smoothed or :buffer, got $(repr(surface))"))
    _validate(options)
    isdir(workdir) || mkpath(workdir)

    # Resolve coordinates + project name and stage the airfoil file in workdir.
    if airfoil isa AbstractString
        isfile(airfoil) || throw(ArgumentError("airfoil file not found: $airfoil"))
        base = something(name, first(splitext(basename(airfoil))))
        datfile = base * ".dat"
        cp(airfoil, joinpath(workdir, datfile); force=true)
        x, y = _read_airfoil_coords(joinpath(workdir, datfile))
    else
        x, y = _coords(airfoil)
        base = something(name, "airfoil")
        datfile = base * ".dat"
        write_airfoil(joinpath(workdir, datfile), x, y; name=base)
    end

    write_grid_options(joinpath(workdir, "grid_options.in"), options; name=base)

    # If the user forces a topology that conflicts with the one Construct2D
    # recommends for this trailing edge, it asks a y/n question at startup
    # (menu.f90 set_defaults). Answer "y" first so the run stays headless.
    commands = String[]
    if options.topo !== nothing
        recommended = _te_is_sharp(x, y) ? "CGRD" : "OGRD"
        uppercase(options.topo) == recommended || push!(commands, "y")
    end
    push!(commands, "GRID")
    push!(commands, surface === :smoothed ? "SMTH" : "BUFF")
    push!(commands, "QUIT")

    ok, log = run_construct2d(workdir, datfile; commands=commands)

    p3d   = joinpath(workdir, base * ".p3d")
    stats = joinpath(workdir, base * "_stats.p3d")
    nmf   = joinpath(workdir, base * ".nmf")
    if !isfile(p3d)
        error("Construct2D did not produce a grid file ($p3d).\n--- log ---\n$log")
    end
    ok || @warn "construct2d exited with a non-zero status; grid file was still produced"

    grid = read_plot3d(p3d)
    result = MeshResult(base, workdir, grid, p3d, stats, nmf, _parse_wall_distance(log), log)

    cleanup && rm(workdir; recursive=true, force=true)
    return result
end

end # module
