# Construct2D.jl

A cross-platform Julia interface to [**Construct2D**](https://github.com/furstj/Construct2D),
a structured grid generator for 2D airfoils (GPL-3.0, originally by Daniel
Prosser; the actively maintained fork is by Jiří Fürst).

The goal: `]add Construct2D` and mesh an airfoil from Julia — no Fortran
compiler, no manual setup — on Linux, macOS, and Windows. The native binary is
delivered as a precompiled artifact via `Construct2D_jll` (see
[Distribution](#distribution)); this package is the idiomatic Julia layer on top.

```julia
using Construct2D

# From a coordinate file (XFoil/Selig labeled format) ...
result = mesh_airfoil("naca0012.dat")

# ... or from coordinates already in Julia, with custom options:
opts = GridOptions(jmax = 120, radi = 20.0, ypls = 0.8, recd = 3.0e6, slvr = "HYPR")
result = mesh_airfoil((x, y); name = "myfoil", options = opts)

result.grid.X      # imax×jmax matrix of node x-coordinates
result.grid.Y      # imax×jmax matrix of node y-coordinates
result.p3d         # path to the Plot3D grid file
result.nmf         # path to the boundary-condition (.nmf) file
result.wall_distance   # first-layer wall spacing parsed from the run log
```

## How it works

Construct2D is a menu-driven Fortran executable, not a library. `Construct2D.jl`
drives it robustly and headlessly via process I/O:

1. writes your airfoil to a `.dat` file and your [`GridOptions`](#options) to a
   `grid_options.in` namelist in a temporary working directory;
2. runs `construct2d <airfoil>.dat`, feeding the menu commands
   (`GRID` → `SMTH`/`BUFF` → `QUIT`) on stdin;
3. reads the resulting `<name>.p3d` Plot3D grid back into Julia arrays.

`mesh_airfoil` is the high-level entry point. `run_construct2d(workdir, datfile;
commands=[...])` is a thin escape hatch if you want to drive the menu yourself.

## Options

[`GridOptions`](src/options.jl) mirrors the `&SOPT` / `&VOPT` / `&OOPT` namelist
groups Construct2D reads. Every field defaults to `nothing` and is **omitted**
from `grid_options.in` when unset, so Construct2D applies its own (geometry-aware)
defaults. Common knobs:

| field  | meaning                                             |
|--------|-----------------------------------------------------|
| `nsrf` | points distributed over the airfoil surface         |
| `jmax` | points in the wall-normal direction                 |
| `radi` | farfield radius (chords)                            |
| `topo` | `"OGRD"` or `"CGRD"` (omit to use the recommended)  |
| `slvr` | `"HYPR"` (hyperbolic) or `"ELLP"` (elliptic)        |
| `ypls` | target y+ (with `recd`, the Reynolds number)        |
| `gdim` | `2` or `3`; `f3dm=true` enables FUN3D output mode   |

See the docstring (`?GridOptions`) for the full list.

## Locating the executable

The binary is resolved in this order:

1. an explicit path from `set_construct2d_path!("/path/to/construct2d")`;
2. the `CONSTRUCT2D_EXE` environment variable;
3. `Construct2D_jll`, once that artifact is installed (see below).

## Distribution

The native binary is built once for every platform with
[BinaryBuilder.jl](https://docs.binarybuilder.org) and published as
`Construct2D_jll` through [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil).
The build recipe lives at [`.ci/build_tarballs.jl`](.ci/build_tarballs.jl).

### Pre-flight CI

The [`JLL preflight`](.github/workflows/jll-preflight.yml) GitHub Action builds the
recipe on Linux and runs the end-to-end test against the freshly built binary (the
test that is skipped when no binary is present). Use it to validate the recipe and
the wrapper before publishing. It does **not** publish the JLL.

### Publishing via Yggdrasil

Copy `.ci/build_tarballs.jl` to `C/Construct2D/build_tarballs.jl` in a fork of
Yggdrasil and open a PR. Yggdrasil's CI builds every platform and registers
`Construct2D_jll` in General.

### Building / testing the JLL locally (BinaryBuilder needs Linux)

On Windows, use WSL2 (or any Linux box). With a recent Julia + `]add BinaryBuilder`:

```bash
julia --color=yes build_tarballs.jl --verbose --debug x86_64-linux-gnu
```

This produces a tarball under `products/` and a local `Construct2D_jll`. Point
this package at the freshly built binary to run the end-to-end tests before the
JLL is registered:

```julia
using Construct2D
Construct2D.set_construct2d_path!("/path/to/products/.../bin/construct2d")
```

### Wiring the JLL (once registered)

After the Yggdrasil PR merges and `Construct2D_jll` is in General, make it the
default provider by adding it as a dependency and registering its `construct2d`
function in `Construct2D`'s `__init__` (one branch in `construct2d_exe`):

```julia
# Project.toml: add Construct2D_jll to [deps] (UUID from the generated JLL)
using Construct2D_jll
function __init__()
    Construct2D._JLL_PROVIDER[] = Construct2D_jll.construct2d
end
```

(Or move this into a package extension so the JLL stays an optional dependency.)

## Status

`v0.1` — staged build. The wrapper API and the BinaryBuilder recipe are complete;
the JLL has not yet been submitted to Yggdrasil. Until then, supply a binary via
`set_construct2d_path!` / `CONSTRUCT2D_EXE`.

## License

`Construct2D.jl` (this wrapper) is MIT-licensed. **Construct2D itself is GPL-3.0**;
this package merely invokes it as a separate executable.
