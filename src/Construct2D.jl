"""
    Construct2D

A pure-Julia structured grid generator for 2D airfoils — a port of Daniel
Prosser's [Construct2D](https://sourceforge.net/projects/construct2d/). No Fortran
toolchain, no external binary: `]add Construct2D` and mesh an airfoil from Julia.

The high-level entry point is [`mesh_airfoil`](@ref). It reads an airfoil, builds
a hyperbolic C-grid (sharp TE) or O-grid (blunt TE) around it, and returns the
node coordinates as Julia matrices.

```julia
using Construct2D
res = mesh_airfoil("naca0012.dat")     # or mesh_airfoil((x, y))
res.X            # imax×jmax node x-coordinates
res.Y            # imax×jmax node y-coordinates
write_plot3d("naca0012.p3d", res.grid) # optional Plot3D export
```

This package is GPL-3.0 (it is a derivative work of the GPL Construct2D); see
the LICENSE file.
"""
module Construct2D

using LinearAlgebra
using Printf

export GridOptions, MeshResult, Plot3DGrid,
       mesh_airfoil,
       read_airfoil, write_airfoil, read_plot3d, write_plot3d

include("types.jl")
include("math_deps.jl")
include("io.jl")
include("options.jl")
include("surface_util.jl")
include("edge_grid.jl")
include("hyperbolic.jl")
include("surface_grid.jl")

"""
    MeshResult

Returned by [`mesh_airfoil`](@ref).

- `name`          : project name
- `grid`          : the [`Plot3DGrid`](@ref) (node-coordinate matrices)
- `X`, `Y`        : convenience aliases for `grid.X` / `grid.Y` (`imax × jmax`)
- `wall_distance` : first off-wall spacing implied by the y+ target
- `options`       : the fully-resolved [`MeshOptions`](@ref) used
"""
struct MeshResult
    name::String
    grid::Plot3DGrid
    wall_distance::Float64
    options::MeshOptions
end

Base.getproperty(r::MeshResult, s::Symbol) =
    s === :X ? getfield(r, :grid).X :
    s === :Y ? getfield(r, :grid).Y :
    getfield(r, s)

function Base.show(io::IO, r::MeshResult)
    print(io, "MeshResult(\"", r.name, "\", ", r.options.topology,
          " grid ", r.grid.dims, ", wall Δ ≈ ", @sprintf("%.3e", r.wall_distance), ")")
end

"""
    mesh_airfoil(airfoil; options=GridOptions(), name=nothing, surface=:buffer) -> MeshResult

Generate a structured grid around `airfoil` and return a [`MeshResult`](@ref).

`airfoil` may be a path to a labeled `.dat` file, a tuple `(x, y)`, an `N×2`
matrix, or a vector of `(x, y)` points.

Keywords:
- `options` : a [`GridOptions`](@ref); unset fields use geometry-aware defaults.
- `name`    : project name (defaults to the file's base name, else `"airfoil"`).
- `surface` : `:buffer` (use the loaded geometry directly — the default and only
  supported mode) or `:smoothed` (XFOIL repaneling — not yet ported).
"""
function mesh_airfoil(airfoil;
                      options::GridOptions=GridOptions(),
                      name::Union{Nothing,AbstractString}=nothing,
                      surface::Symbol=:buffer)
    surface in (:buffer, :smoothed) ||
        throw(ArgumentError("surface must be :buffer or :smoothed, got $(repr(surface))"))

    if airfoil isa AbstractString
        isfile(airfoil) || throw(ArgumentError("airfoil file not found: $airfoil"))
        fname, x, y = read_airfoil(airfoil)
        base = something(name, first(splitext(basename(airfoil))))
    else
        x, y = _coords(airfoil)
        base = something(name, "airfoil")
    end

    opts, foil = resolve_options(x, y; name=base, opts=options)
    grid, y0 = create_grid!(foil, opts; smooth=(surface === :smoothed))
    return MeshResult(base, Plot3DGrid(grid), y0, opts)
end

end # module
