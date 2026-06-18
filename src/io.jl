# File I/O: writing airfoil coordinate files Construct2D can read, and reading
# back the Plot3D grids it produces.

# --- airfoil coordinate file (XFoil "Selig" labeled format) -----------------
#
# Construct2D's reader (src/util.f90 read_airfoil) skips one header line, then
# reads `x y` pairs (whitespace separated). Points should be ordered like XFoil
# / airfoiltools.com Selig format: from the trailing edge forward over the upper
# surface to the leading edge, then back along the lower surface to the trailing
# edge. Construct2D re-orders to counter-clockwise itself if needed.

"""
    write_airfoil(path, x, y; name="airfoil")
    write_airfoil(path, coords::AbstractMatrix; name="airfoil")   # N×2: column 1 = x

Write airfoil coordinates to `path` in the labeled (Selig) format Construct2D
reads: a header line `name`, followed by `x  y` columns.
"""
function write_airfoil(path::AbstractString, x::AbstractVector{<:Real},
                       y::AbstractVector{<:Real}; name::AbstractString="airfoil")
    length(x) == length(y) ||
        throw(DimensionMismatch("x and y must have equal length, got $(length(x)) and $(length(y))"))
    open(path, "w") do io
        println(io, name)
        for (xi, yi) in zip(x, y)
            @printf(io, "%18.12f  %18.12f\n", xi, yi)
        end
    end
    return path
end

function write_airfoil(path::AbstractString, coords::AbstractMatrix; name::AbstractString="airfoil")
    size(coords, 2) == 2 ||
        throw(DimensionMismatch("coords must be an N×2 matrix (column 1 = x, column 2 = y)"))
    return write_airfoil(path, @view(coords[:, 1]), @view(coords[:, 2]); name=name)
end

# Normalize the various ways a caller may pass coordinates into (x, y) vectors.
_coords(c::Tuple{<:AbstractVector,<:AbstractVector}) = (collect(float.(c[1])), collect(float.(c[2])))
function _coords(c::AbstractMatrix)
    size(c, 2) == 2 || throw(DimensionMismatch("coordinate matrix must be N×2"))
    return (collect(float.(c[:, 1])), collect(float.(c[:, 2])))
end
function _coords(c::AbstractVector)
    # Vector of (x,y) points / 2-element vectors.
    x = [float(p[1]) for p in c]
    y = [float(p[2]) for p in c]
    return (x, y)
end

# Read an airfoil .dat back into (x, y) (header line skipped). Used internally
# for trailing-edge detection when the caller supplied a file path.
function _read_airfoil_coords(path::AbstractString)
    x = Float64[]; y = Float64[]
    open(path, "r") do io
        readline(io)  # header
        for ln in eachline(io)
            t = split(strip(ln))
            length(t) >= 2 || continue
            push!(x, parse(Float64, t[1]))
            push!(y, parse(Float64, t[2]))
        end
    end
    return (x, y)
end

# Sharp trailing edge if the gap between the first and last surface points is
# negligible relative to the chord (mirrors Construct2D's surf%tegap logic).
function _te_is_sharp(x::AbstractVector, y::AbstractVector)
    length(x) >= 2 || return true
    chord = maximum(x) - minimum(x)
    chord = chord == 0 ? 1.0 : chord
    gap = hypot(x[1] - x[end], y[1] - y[end])
    return gap / chord < 1e-5
end

# --- Plot3D grid reader -----------------------------------------------------
#
# Matches src/util.f90 write_srf_grid. 2D (gdim=2): a single line `imax jmax`,
# then imax*jmax x-values followed by imax*jmax y-values, each on its own line
# (es17.8), written in the order ((v(i,j), i=imax,1,-1), j=1,jmax) — i runs
# backwards. 3D (gdim=3, FUN3D mode): a line `1` (block count), then
# `imax kmax jmax`, then x, z, y blocks in the same i-reversed ordering.

"""
    Plot3DGrid

A structured grid read from a Construct2D `.p3d` file.

- `X`, `Y` : node-coordinate arrays. For a 2D grid these are `imax × jmax`
  matrices indexed so that `X[1,j]`/`Y[1,j]` is the first airfoil-surface point
  (the file's i-reversed ordering is undone on read). For a 3D grid they are
  `imax × kmax × jmax` arrays.
- `Z`      : `nothing` for 2D grids, otherwise the spanwise coordinate array.
- `dims`   : `(imax, jmax)` or `(imax, kmax, jmax)`.
"""
struct Plot3DGrid
    X::Array{Float64}
    Y::Array{Float64}
    Z::Union{Nothing,Array{Float64}}
    dims::Tuple{Vararg{Int}}
end

Base.size(g::Plot3DGrid) = g.dims

function _stream_floats(lines, start::Int)
    vals = Float64[]
    for k in start:length(lines)
        for t in split(strip(lines[k]))
            push!(vals, parse(Float64, t))
        end
    end
    return vals
end

"""
    read_plot3d(path) -> Plot3DGrid

Read a Construct2D-produced Plot3D grid file (2D or 3D single block).
"""
function read_plot3d(path::AbstractString)
    lines = readlines(path)
    h = findfirst(l -> !isempty(strip(l)), lines)
    h === nothing && error("empty Plot3D file: $path")
    hdr = split(strip(lines[h]))

    if length(hdr) == 1
        # 3D single-block: header is block count (1), dims on the next line.
        d = findnext(l -> !isempty(strip(l)), lines, h + 1)
        dt = split(strip(lines[d]))
        imax, kmax, jmax = parse(Int, dt[1]), parse(Int, dt[2]), parse(Int, dt[3])
        n = imax * kmax * jmax
        v = _stream_floats(lines, d + 1)
        length(v) >= 3n || error("Plot3D file $path: expected $(3n) values, found $(length(v))")
        X = Array{Float64}(undef, imax, kmax, jmax)
        Z = Array{Float64}(undef, imax, kmax, jmax)
        Y = Array{Float64}(undef, imax, kmax, jmax)
        idx = 0
        for arr in (X, Z, Y)
            for j in 1:jmax, k in 1:kmax, i in imax:-1:1
                idx += 1
                arr[i, k, j] = v[idx]
            end
        end
        return Plot3DGrid(X, Y, Z, (imax, kmax, jmax))
    else
        imax, jmax = parse(Int, hdr[1]), parse(Int, hdr[2])
        n = imax * jmax
        v = _stream_floats(lines, h + 1)
        length(v) >= 2n || error("Plot3D file $path: expected $(2n) values, found $(length(v))")
        X = Matrix{Float64}(undef, imax, jmax)
        Y = Matrix{Float64}(undef, imax, jmax)
        idx = 0
        for arr in (X, Y)
            for j in 1:jmax, i in imax:-1:1
                idx += 1
                arr[i, j] = v[idx]
            end
        end
        return Plot3DGrid(X, Y, nothing, (imax, jmax))
    end
end

# --- result type & stdout parsing -------------------------------------------

"""
    MeshResult

Returned by [`mesh_airfoil`](@ref).

- `name`       : project name (basis of the output file names)
- `workdir`    : directory containing the generated files
- `grid`       : the parsed [`Plot3DGrid`](@ref)
- `p3d`, `stats`, `nmf` : paths to the `.p3d` grid, `_stats.p3d`, and `.nmf` files
- `wall_distance` : first-layer wall spacing parsed from stdout (`nothing` if not found)
- `log`        : full combined stdout/stderr from the Construct2D run
"""
struct MeshResult
    name::String
    workdir::String
    grid::Plot3DGrid
    p3d::String
    stats::String
    nmf::String
    wall_distance::Union{Nothing,Float64}
    log::String
end

function Base.show(io::IO, r::MeshResult)
    print(io, "MeshResult(\"", r.name, "\", grid dims ", r.grid.dims)
    r.wall_distance !== nothing && print(io, ", wall Δ ≈ ", r.wall_distance)
    print(io, ")")
end

# Construct2D prints the first-layer wall distance to stdout when generating a
# mesh (see README release notes for v2.1.4). The exact wording can vary between
# versions, so parse leniently: take the first float on a line mentioning the
# wall distance / first layer spacing.
function _parse_wall_distance(log::AbstractString)
    for ln in split(log, '\n')
        if occursin(r"(?i)(wall\s*(distance|spacing)|first[- ]?layer)", ln)
            m = match(r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?", ln)
            if m !== nothing
                return parse(Float64, replace(m.match, r"[dD]" => "e"))
            end
        end
    end
    return nothing
end
