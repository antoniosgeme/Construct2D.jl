# Math utilities.  <- src/math_deps.f90
#
# Faithful Julia translation of the math helpers Construct2D's grid generators
# rely on. Finite-difference stencils are kept as small functions (the Fortran
# calls them with step h = 1 throughout). The dense linear solve `lmult` becomes
# Julia's backslash.

# --- finite-difference stencils  <- derv1*/derv2* ---------------------------
# All are called with h = 1 in the solver, but keep h for fidelity.

derv1f(u_plus, u, h=1.0)            = (u_plus - u) / h               # forward
derv1b(u_minus, u, h=1.0)          = (u - u_minus) / h              # backward
derv1c(u_plus, u_minus, h=1.0)     = (u_plus - u_minus) / (2.0 * h) # central
derv2f(u_p2, u_p, u, h=1.0)        = (u - 2.0 * u_p + u_p2) / h^2
derv2b(u_m2, u_m, u, h=1.0)        = (u - 2.0 * u_m + u_m2) / h^2
derv2c(u_plus, u, u_minus, h=1.0)  = (u_plus - 2.0 * u + u_minus) / h^2

# --- small geometry helpers -------------------------------------------------

"""
    angle(x1, y1, x2, y2, x0, y0) -> degrees

Angle at vertex `(x0,y0)` between rays to `(x1,y1)` and `(x2,y2)`.
<- math_deps.f90 :: angle
"""
function angle(x1, y1, x2, y2, x0, y0)
    v1x = x1 - x0; v1y = y1 - y0
    v2x = x2 - x0; v2y = y2 - y0
    mag1 = sqrt(v1x^2 + v1y^2)
    mag2 = sqrt(v2x^2 + v2y^2)
    return acos((v1x * v2x + v1y * v2y) / (mag1 * mag2)) * 180.0 / pi
end

"""
    growth(x1, y1, x0, y0, xm1, ym1) -> ratio

Signed growth ratio of consecutive edges `(xm1,ym1)->(x0,y0)->(x1,y1)`.
<- math_deps.f90 :: growth
"""
function growth(x1, y1, x0, y0, xm1, ym1)
    len1 = sqrt((x1 - x0)^2 + (y1 - y0)^2)
    len0 = sqrt((x0 - xm1)^2 + (y0 - ym1)^2)
    lennorm = len0 <= len1 ? len0 : len1
    return (len1 - len0) / lennorm
end

"`between(a, b, c)` — is `b` in `[a, c]`?  <- math_deps.f90 :: between"
between(a, b, c) = (b >= a) && (b <= c)

"`interp1` — linear interpolation.  <- math_deps.f90 :: interp1"
interp1(x1, x2, x, y1, y2) = y1 + (y2 - y1) * (x - x1) / (x2 - x1)

"""
    polyline_dist(x, y) -> cdf

Cumulative arc length along the polyline `(x,y)` (cdf[1] = 0).
<- math_deps.f90 :: polyline_dist
"""
function polyline_dist(x::AbstractVector, y::AbstractVector)
    n = length(x)
    cdf = zeros(Float64, n)
    for j in 2:n
        cdf[j] = cdf[j-1] + sqrt((x[j] - x[j-1])^2 + (y[j] - y[j-1])^2)
    end
    return cdf
end

# --- geometric-series point spacing  <- get_growth / xi_length / golden_search

"""
    xi_length(d0, ga, N) -> L

Total length of `N` intervals starting at `d0` with geometric ratio `ga`:
`d0 * sum_{i=0}^{N-1} ga^i`.  <- edge_grid.f90 :: xi_length / math_deps.f90 :: lfunc
"""
function xi_length(d0, ga, N::Integer)
    l = 0.0
    for i in 1:N
        l += ga^(i - 1)
    end
    return d0 * l
end

# Objective for the golden search: |lact - xi_length(d0, ga, N)|.
# <- math_deps.f90 :: lfunc
_lfunc(lact, d0, ga, N) = abs(lact - xi_length(d0, ga, N))

"""
    golden_search(bounds, lact, d0, N) -> (gmin, fmin)

Golden-section minimisation of `_lfunc` over `ga ∈ bounds`.
<- math_deps.f90 :: golden_search
"""
function golden_search(bounds::NTuple{2,Float64}, lact, d0, N::Integer)
    tol = 1.0e-9
    imax = 100
    x1, x4 = bounds[1], bounds[2]
    T = (3.0 - sqrt(5.0)) / 2.0               # golden section
    x2 = (1.0 - T) * x1 + T * x4
    x3 = T * x1 + (1.0 - T) * x4
    f1 = _lfunc(lact, d0, x1, N); f2 = _lfunc(lact, d0, x2, N)
    f3 = _lfunc(lact, d0, x3, N); f4 = _lfunc(lact, d0, x4, N)
    if f1 > f4
        xmin, fmin = x4, f4
    else
        xmin, fmin = x1, f1
    end
    i = 2
    dist = x4 - x1
    while i < imax && dist > tol
        if f2 > f3
            xmin, fmin = x3, f3
            x1, x2, f1, f2 = x2, x3, f2, f3
            x3 = T * x1 + (1.0 - T) * x4
            f3 = _lfunc(lact, d0, x3, N)
        else
            xmin, fmin = x2, f2
            x4, x3, f4, f3 = x3, x2, f3, f2
            x2 = (1.0 - T) * x1 + T * x4
            f2 = _lfunc(lact, d0, x2, N)
        end
        dist = x4 - x1
        i += 1
    end
    return xmin, fmin
end

"""
    get_growth(lact, d0, N) -> g

Geometric growth ratio `g` so that `N` intervals starting at `d0` span length
`lact`. Brackets `g` upward from 0 in steps of 0.05, then golden-searches.
<- edge_grid.f90 :: get_growth
"""
function get_growth(lact, d0, N::Integer)
    lo = 0.0
    hi = 0.0
    bracketed = false
    while !bracketed
        hi = lo + 0.05
        ll = xi_length(d0, lo, N)
        lr = xi_length(d0, hi, N)
        bracketed = between(ll, lact, lr)
        if !bracketed
            lo = hi
        end
    end
    g, _ = golden_search((lo, hi), lact, d0, N)
    return g
end

"""
    interp_spline(s, xs, ys, cdfs, pt1) -> (x, y, pt1)

Linearly interpolate `(x,y)` at arc length `s` along the polyline `(xs,ys)` whose
cumulative distance is `cdfs`. `pt1` is a search hint advanced in place (returned).
<- edge_grid.f90 :: interp_spline
"""
function interp_spline(s, xs::AbstractVector, ys::AbstractVector, cdfs::AbstractVector, pt1::Integer)
    npt = length(xs)
    pt1store = pt1
    isbtwn = false
    while !isbtwn && pt1 < npt
        isbtwn = between(cdfs[pt1], s, cdfs[pt1+1])
        if !isbtwn
            pt1 += 1
            if pt1 == npt
                pt1 = pt1store
                error("interp_spline: could not find interpolants (s=$s, cdfmax=$(cdfs[npt]))")
            end
        end
    end
    x = interp1(cdfs[pt1], cdfs[pt1+1], s, xs[pt1], xs[pt1+1])
    y = interp1(cdfs[pt1], cdfs[pt1+1], s, ys[pt1], ys[pt1+1])
    return x, y, pt1
end

# --- normal distribution used by the blunt-TE fillet  <- normal_dist/g_func_norm

"`normal_dist(x, sig, mu)` — Gaussian pdf.  <- math_deps.f90 :: normal_dist"
normal_dist(x, sig, mu) = 1.0 / (sig * sqrt(2.0 * pi)) * exp(-(x - mu)^2 / (2.0 * sig^2))

# <- edge_grid.f90 :: g_func_norm
_g_func_norm(i, N, sig, mu, G0, Gmid) =
    (normal_dist((i - 1) / (N - 1), sig, mu) - G0) / (Gmid - G0)

# --- B-spline  <- bspline / basis_func --------------------------------------

"""
    basis_func(j, i, x, t) -> value

Cox-de Boor B-spline basis function (recursive). `i` is 1-based.
<- math_deps.f90 :: basis_func (recursive)
"""
function basis_func(j::Integer, i::Integer, x::AbstractVector, t)
    m = length(x)
    if j == 0
        if (x[i] <= t) && (t < x[i+1])
            return 1.0
        elseif (x[i] <= t) && (t == x[i+1]) && (x[i+1] == 1.0)
            return 1.0
        else
            return 0.0
        end
    else
        val = 0.0
        if x[i] < x[i+j]
            val = (t - x[i]) / (x[i+j] - x[i]) * basis_func(j - 1, i, x, t)
        end
        if i < m
            if x[i+1] < x[i+j+1]
                val += (x[i+j+1] - t) / (x[i+j+1] - x[i+1]) * basis_func(j - 1, i + 1, x, t)
            end
        end
        return val
    end
end

"""
    bspline(Pin, degree, npoints) -> (x, y, z)

Evaluate a B-spline of the given `degree` through control points `Pin` (a
`d×ncp` matrix, `d` = 2 or 3) at `npoints` uniformly spaced parameter values.
<- math_deps.f90 :: bspline
"""
function bspline(Pin::AbstractMatrix, degree::Integer, npoints::Integer)
    d = size(Pin, 1)
    ncp = size(Pin, 2)
    n = ncp - 1                              # number of control points - 1
    d in (2, 3) || error("bspline: control points must have 2 or 3 spatial dimensions")
    degree <= n || error("bspline: degree must be less than the number of control points")
    degree >= 1 || error("bspline: degree must be >= 1")

    P = zeros(Float64, 3, ncp)               # promote to 3D (z = 0 if 2D input)
    P[1:d, :] .= Pin

    # Knot vector
    T = zeros(Float64, degree + n + 2)
    for i in 1:degree+n+2
        if i < degree + 2
            T[i] = 0.0
        elseif i < n + 2
            T[i] = (i - degree - 1) / (n + 1 - degree)
        else
            T[i] = 1.0
        end
    end

    Q = zeros(Float64, 3, npoints)
    for l in 1:npoints
        tval = (l - 1) / (npoints - 1)
        for i in 0:n
            Nval = basis_func(degree, i + 1, T, tval)
            @views Q[:, l] .+= Nval .* P[:, i+1]
        end
    end
    return Q[1, :], Q[2, :], Q[3, :]
end
