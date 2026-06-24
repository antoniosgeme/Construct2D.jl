# File I/O: airfoil coordinate files and Plot3D grids.  <- src/util.f90

# --- airfoil coordinate file (XFOIL "Selig" labeled format) -----------------
# One header line, then `x y` pairs, ordered TE -> upper -> LE -> lower -> TE.

"""
    read_airfoil(path) -> (name, x, y)

Read a labeled airfoil `.dat` file: the first line is the name, the rest are
`x  y` columns. <- util.f90 :: read_airfoil
"""
function read_airfoil(path::AbstractString)
    x = Float64[]; y = Float64[]
    name = "airfoil"
    open(path, "r") do io
        name = strip(readline(io))
        for ln in eachline(io)
            t = split(strip(ln))
            length(t) >= 2 || continue
            push!(x, parse(Float64, t[1]))
            push!(y, parse(Float64, t[2]))
        end
    end
    return String(name), x, y
end

"""
    write_airfoil(path, x, y; name="airfoil")

Write airfoil coordinates in the labeled (Selig) format: a header line `name`,
then `x  y` columns.
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

# Normalize the various coordinate inputs into (x, y) vectors.
_coords(c::Tuple{<:AbstractVector,<:AbstractVector}) = (collect(float.(c[1])), collect(float.(c[2])))
function _coords(c::AbstractMatrix)
    size(c, 2) == 2 || throw(DimensionMismatch("coordinate matrix must be N×2"))
    return (collect(float.(c[:, 1])), collect(float.(c[:, 2])))
end
function _coords(c::AbstractVector)
    x = [float(p[1]) for p in c]
    y = [float(p[2]) for p in c]
    return (x, y)
end

# --- Plot3D grid ------------------------------------------------------------

"""
    Plot3DGrid

A structured grid. `X`, `Y` are `imax × jmax` node-coordinate matrices with
`X[i,1]`/`Y[i,1]` on the airfoil surface; `dims = (imax, jmax)`.
"""
struct Plot3DGrid
    X::Matrix{Float64}
    Y::Matrix{Float64}
    dims::Tuple{Int,Int}
end

Plot3DGrid(g::SrfGrid) = Plot3DGrid(copy(g.x), copy(g.y), (g.imax, g.jmax))
Base.size(g::Plot3DGrid) = g.dims

"""
    write_plot3d(path, grid)

Write a 2D Plot3D `.p3d` file: header `imax jmax`, then all x then all y values,
each on its own line in the i-reversed order `((v(i,j), i=imax..1), j=1..jmax)`
that Construct2D uses to keep positive cell volumes.
<- util.f90 :: write_srf_grid (2D path)
"""
function write_plot3d(path::AbstractString, grid::Union{SrfGrid,Plot3DGrid})
    X = grid isa SrfGrid ? grid.x : grid.X
    Y = grid isa SrfGrid ? grid.y : grid.Y
    imax, jmax = grid isa SrfGrid ? (grid.imax, grid.jmax) : grid.dims
    open(path, "w") do io
        println(io, "  ", imax, "  ", jmax)
        for arr in (X, Y)
            for j in 1:jmax, i in imax:-1:1
                @printf(io, "%17.8e\n", arr[i, j])
            end
        end
    end
    return path
end

"""
    read_plot3d(path) -> Plot3DGrid

Read a 2D single-block Plot3D grid written by [`write_plot3d`](@ref) (undoes the
i-reversed ordering).
"""
function read_plot3d(path::AbstractString)
    lines = readlines(path)
    h = findfirst(l -> !isempty(strip(l)), lines)
    h === nothing && error("empty Plot3D file: $path")
    hdr = split(strip(lines[h]))
    imax, jmax = parse(Int, hdr[1]), parse(Int, hdr[2])
    n = imax * jmax
    vals = Float64[]
    for k in h+1:length(lines)
        for t in split(strip(lines[k]))
            push!(vals, parse(Float64, t))
        end
    end
    length(vals) >= 2n || error("Plot3D file $path: expected $(2n) values, found $(length(vals))")
    X = Matrix{Float64}(undef, imax, jmax)
    Y = Matrix{Float64}(undef, imax, jmax)
    idx = 0
    for arr in (X, Y)
        for j in 1:jmax, i in imax:-1:1
            idx += 1
            arr[i, j] = vals[idx]
        end
    end
    return Plot3DGrid(X, Y, (imax, jmax))
end
