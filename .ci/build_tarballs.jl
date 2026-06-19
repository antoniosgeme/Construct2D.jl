# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
#
# Build recipe for Construct2D — a structured grid generator for 2D airfoils.
# Upstream (actively maintained fork): https://github.com/furstj/Construct2D  (GPL-3.0)
#
# Submit by opening a PR that adds this file at  C/Construct2D/build_tarballs.jl
# in https://github.com/JuliaPackaging/Yggdrasil . On merge, Yggdrasil's CI
# builds every platform and registers `Construct2D_jll` in the General registry.

using BinaryBuilder

name    = "Construct2D"
version = v"2.1.5"

# Pinned to the v2.1.5 release tag commit of the furstj fork.
sources = [
    GitSource("https://github.com/furstj/Construct2D.git",
              "ac989b5eadb791130db13b8541c85e5a93545fd2"),
]

# Build steps.
#
# Gotcha discovered while writing this recipe: the repo ships THREE makefiles and
# the *system-specific* ones (Makefile_Linux_MacOSX / Makefile_Windows) are STALE
# — they still reference `src/xfoil_deps.f`, but v2.1.5 renamed that file to
# `src/xfoil_deps.f90`, so those makefiles fail. The plain root `Makefile` is the
# up-to-date one (uses `xfoil_deps.f90`) and already honors `$(FC)`, so we just
# point it at the cross-compiler. There is no `install` target, so we copy the
# resulting executable ourselves. Build serially: the root Makefile orders the
# `xfoil_deps` module before its consumers, which `-j` could violate.
#
# Second gotcha: the root Makefile hardcodes `FCFLAGS=... -std=f2023`, but
# BinaryBuilder's gfortran is GCC 8.1, which predates that flag (`-std=f2023`
# arrived in GCC 13). The sources are really Fortran 2018, so downgrade the
# standard rather than forcing a newer, less-portable compiler
# (`preferred_gcc_version`), which would raise the JLL's glibc baseline.
script = raw"""
cd ${WORKSPACE}/srcdir/Construct2D
sed -i 's/-std=f2023/-std=f2018/' Makefile
make FC=${FC}

# Some toolchains (e.g. mingw) append .exe to the linker output automatically.
if [[ -f construct2d.exe ]]; then
    install -Dvm755 construct2d.exe "${bindir}/construct2d.exe"
else
    install -Dvm755 construct2d "${bindir}/construct2d${exeext}"
fi
"""

# Pure-Fortran program with no external library deps. Expand the platform set
# across libgfortran ABI versions so the binary matches whatever libgfortran the
# user's Julia was built against.
platforms = expand_gfortran_versions(supported_platforms())

products = [
    ExecutableProduct("construct2d", :construct2d),
]

# Supplies libgfortran / libquadmath at runtime.
dependencies = [
    Dependency("CompilerSupportLibraries_jll"),
]

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat="1.6")
