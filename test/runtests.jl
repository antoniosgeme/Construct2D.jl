using Construct2D
using Construct2D: angle, growth, get_growth, xi_length, interp1, polyline_dist,
                  bspline, _classify_te, resolve_options, _coords
using Test

const FIX = joinpath(@__DIR__, "fixtures")

# Count genuinely folded cells: quads whose signed area runs against the grid's
# dominant orientation by more than a tiny fraction of the largest cell. (Near-TE
# cells can be ~0 and flip sign at roundoff; those are not folds.)
function count_folds(X, Y)
    imax, jmax = size(X)
    areas = Float64[]
    for j in 1:jmax-1, i in 1:imax-1
        x1, y1 = X[i, j], Y[i, j];         x2, y2 = X[i+1, j], Y[i+1, j]
        x3, y3 = X[i+1, j+1], Y[i+1, j+1]; x4, y4 = X[i, j+1], Y[i, j+1]
        push!(areas, 0.5 * ((x1*y2 - x2*y1) + (x2*y3 - x3*y2) +
                            (x3*y4 - x4*y3) + (x4*y1 - x1*y4)))
    end
    amax = maximum(abs, areas)
    thresh = 1e-8 * amax
    dominant = sum(sign, areas) >= 0 ? 1.0 : -1.0     # most cells' sign
    return count(a -> sign(a) != dominant && abs(a) > thresh, areas)
end

@testset "Construct2D" begin

    @testset "math: angle / growth / interp" begin
        @test angle(1.0, 0.0, 0.0, 1.0, 0.0, 0.0) ≈ 90.0
        @test angle(1.0, 0.0, -1.0, 0.0, 0.0, 0.0) ≈ 180.0
        @test growth(2.0, 0.0, 1.0, 0.0, 0.0, 0.0) ≈ 0.0
        @test growth(3.0, 0.0, 1.0, 0.0, 0.0, 0.0) ≈ 1.0
        @test interp1(0.0, 2.0, 1.0, 10.0, 20.0) ≈ 15.0
        @test polyline_dist([0.0, 3.0, 3.0], [0.0, 0.0, 4.0]) ≈ [0.0, 3.0, 7.0]
    end

    @testset "get_growth inverts xi_length" begin
        # g recovered from a target length reproduces that length.
        for (L, d0, N) in ((10.0, 0.1, 30), (5.0, 0.001, 50), (20.0, 0.05, 40))
            g = get_growth(L, d0, N)
            @test xi_length(d0, g, N) ≈ L rtol = 1e-6
            @test g > 1.0                      # clustered (growing) spacing
        end
    end

    @testset "bspline passes through its end control points" begin
        CP = [0.0 1.0 2.0 3.0;
              0.0 2.0 -1.0 0.5]
        x, y, z = bspline(CP, 3, 50)
        @test length(x) == 50
        @test x[1] == CP[1, 1] && y[1] == CP[2, 1]      # clamped at start
        @test x[end] ≈ CP[1, end] && y[end] ≈ CP[2, end]  # and end
        @test all(==(0.0), z)                            # 2D input ⇒ z = 0
    end

    @testset "airfoil read/write round trip" begin
        x = [1.0, 0.5, 0.0, 0.5, 1.0]
        y = [0.0, 0.08, 0.0, -0.08, 0.0]
        p = tempname() * ".dat"
        write_airfoil(p, x, y; name="diamond")
        nm, x2, y2 = read_airfoil(p)
        @test nm == "diamond"
        @test x2 ≈ x && y2 ≈ y
    end

    @testset "plot3d write/read round trip" begin
        imax, jmax = 4, 3
        X = Float64[10i + j for i in 1:imax, j in 1:jmax]
        Y = Float64[100i + j for i in 1:imax, j in 1:jmax]
        g = Plot3DGrid(X, Y, (imax, jmax))
        p = tempname() * ".p3d"
        write_plot3d(p, g)
        g2 = read_plot3d(p)
        @test g2.dims == (imax, jmax)
        @test g2.X ≈ X && g2.Y ≈ Y
    end

    @testset "TE classification & derived sizes" begin
        # sharp (closed) ⇒ CGRD; blunt ⇒ OGRD
        tegap, topo = _classify_te([1.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        @test tegap == false && topo == "CGRD"
        tegap, topo = _classify_te([1.0, 0.0, 1.0], [0.01, 0.0, -0.01])
        @test tegap == true && topo == "OGRD"

        # imax bookkeeping: OGRD ⇒ nsrf; CGRD ⇒ nsrf + 2·nwke
        _, x, y = read_airfoil(joinpath(FIX, "naca0012.dat"))         # blunt
        o, _ = resolve_options(x, y; name="b")
        @test o.topology == "OGRD"
        @test o.imax == length(x) + o.nte + 1

        _, xs, ys = read_airfoil(joinpath(FIX, "naca0012_sharp.dat")) # sharp
        o2, _ = resolve_options(xs, ys; name="s")
        @test o2.topology == "CGRD"
        @test o2.imax == length(xs) + 2 * o2.nwake
    end

    @testset "end-to-end: blunt TE ⇒ O-grid" begin
        res = mesh_airfoil(joinpath(FIX, "naca0012.dat"); options=GridOptions(jmax=60))
        @test res.options.topology == "OGRD"
        @test res.grid.dims == (res.options.imax, 60)
        @test size(res.X) == size(res.Y) == res.grid.dims
        @test !any(isnan, res.X) && !any(isnan, res.Y)
        @test res.wall_distance > 0
        @test count_folds(res.X, res.Y) == 0           # no folded cells
        # Farfield row roughly at the requested radius (within the solver's scaling).
        rfar = maximum(hypot.(res.X[:, end] .- 0.5, res.Y[:, end]))
        @test 8.0 < rfar < 25.0
    end

    @testset "end-to-end: sharp TE ⇒ C-grid" begin
        # Use jmax=80: a C-grid with ~50 wake points needs enough normal levels
        # to resolve the wake, and folds if starved (e.g. jmax=60). Clean ≥ 80.
        res = mesh_airfoil(joinpath(FIX, "naca0012_sharp.dat"); options=GridOptions(jmax=80))
        @test res.options.topology == "CGRD"
        @test res.grid.dims[1] == res.options.imax
        @test !any(isnan, res.X) && !any(isnan, res.Y)
        @test count_folds(res.X, res.Y) == 0
        # Surface row reaches the leading edge (x ≈ 0) and trailing edge (x ≈ 1).
        @test minimum(res.X[:, 1]) < 0.02
        @test maximum(res.X[:, 1]) > 0.98
    end

    @testset "coordinate input forms agree" begin
        _, x, y = read_airfoil(joinpath(FIX, "naca0012.dat"))
        a = mesh_airfoil((x, y); name="t", options=GridOptions(jmax=40))
        b = mesh_airfoil(hcat(x, y); name="t", options=GridOptions(jmax=40))
        @test a.X ≈ b.X && a.Y ≈ b.Y
    end

    @testset "unsupported paths fail loudly" begin
        @test_throws ErrorException mesh_airfoil(joinpath(FIX, "naca0012.dat"); surface=:smoothed)
        @test_throws ErrorException mesh_airfoil(joinpath(FIX, "naca0012.dat");
                                                 options=GridOptions(slvr="ELLP"))
        @test_throws ArgumentError GridOptions(slvr="NOPE") |> x ->
            mesh_airfoil(joinpath(FIX, "naca0012.dat"); options=x)
    end

end
