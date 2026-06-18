using Construct2D
using Test

const FIX = joinpath(@__DIR__, "fixtures")

# Is a real construct2d executable available for the end-to-end test?
function _exe_available()
    !isempty(Construct2D._EXE_OVERRIDE[]) && return true
    haskey(ENV, "CONSTRUCT2D_EXE") && return true
    return false
end

@testset "Construct2D.jl" begin

    @testset "write_grid_options namelist" begin
        opts = GridOptions(nsrf=200, radi=20.0, slvr="hypr", topo="ogrd",
                           recd=1.0e6, gdim=2, f3dm=false)
        io = IOBuffer()
        write_grid_options(io, opts; name="myfoil")
        s = String(take!(io))

        # Group structure and only-set-fields behavior.
        @test occursin("&SOPT", s)
        @test occursin("&VOPT", s)
        @test occursin("&OOPT", s)
        @test count("/\n", s) == 3
        @test occursin("name = 'myfoil'", s)
        @test occursin("nsrf = 200", s)
        @test occursin("slvr = 'HYPR'", s)   # normalized to upper case
        @test occursin("topo = 'OGRD'", s)
        @test occursin("f3dm = .false.", s)
        # radi is a float -> must be written with a decimal point
        @test occursin(r"radi = 20(\.0|\.)", s)
        @test occursin(r"recd = 1(e\+?0?6|000000\.?0?)", s)
        # unset fields are omitted entirely
        @test !occursin("jmax", s)
        @test !occursin("lesp", s)
    end

    @testset "options validation" begin
        @test_throws ArgumentError write_grid_options(IOBuffer(), GridOptions(slvr="NOPE"); name="x")
        @test_throws ArgumentError write_grid_options(IOBuffer(), GridOptions(topo="XGRD"); name="x")
        @test_throws ArgumentError write_grid_options(IOBuffer(), GridOptions(gdim=4); name="x")
    end

    @testset "airfoil write/read round trip" begin
        x = [1.0, 0.5, 0.0, 0.5, 1.0]
        y = [0.0, 0.08, 0.0, -0.08, 0.0]
        path = tempname() * ".dat"
        Construct2D.write_airfoil(path, x, y; name="diamond")
        @test readlines(path)[1] == "diamond"
        x2, y2 = Construct2D._read_airfoil_coords(path)
        @test x2 ≈ x
        @test y2 ≈ y

        # matrix input writes the same coordinates
        path2 = tempname() * ".dat"
        Construct2D.write_airfoil(path2, hcat(x, y); name="diamond")
        x3, y3 = Construct2D._read_airfoil_coords(path2)
        @test x3 ≈ x && y3 ≈ y
    end

    @testset "trailing-edge sharpness" begin
        # closed (sharp) TE: first point == last point
        @test Construct2D._te_is_sharp([1.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        # open (blunt) TE: finite gap at x=1
        @test !Construct2D._te_is_sharp([1.0, 0.0, 1.0], [0.01, 0.0, -0.01])
    end

    @testset "read_plot3d (2D) matches Construct2D write order" begin
        # Build known node arrays and write them exactly as src/util.f90
        # write_srf_grid does: header "imax jmax", then x then y, each value on
        # its own line, in the order ((v(i,j), i=imax,1,-1), j=1,jmax).
        imax, jmax = 3, 2
        X = Float64[10i + j for i in 1:imax, j in 1:jmax]
        Y = Float64[100i + j for i in 1:imax, j in 1:jmax]
        path = tempname() * ".p3d"
        open(path, "w") do io
            println(io, "  $imax  $jmax")
            for arr in (X, Y), j in 1:jmax, i in imax:-1:1
                println(io, arr[i, j])
            end
        end

        g = read_plot3d(path)
        @test g.dims == (imax, jmax)
        @test g.Z === nothing
        @test g.X == X
        @test g.Y == Y
    end

    @testset "end-to-end mesh (requires construct2d)" begin
        if _exe_available()
            res = mesh_airfoil(joinpath(FIX, "naca0012.dat"); name="naca0012")
            @test isfile(res.p3d)
            @test isfile(res.nmf)
            @test res.grid.dims[1] > 1 && res.grid.dims[2] > 1
            @test size(res.grid.X) == size(res.grid.Y)
        else
            @info "Skipping end-to-end test: no construct2d executable " *
                  "(set CONSTRUCT2D_EXE or call set_construct2d_path!)."
            @test_skip false
        end
    end

end
