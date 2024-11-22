# This file contains an interface to GeodesicLib through Proj.  It provides some basic interfaces for:
# - `geod_geodesic(eq_radius, flattening)`
# - `geod_direct(g, lat1, lon1, azi1, s12)`
# - `geod_inverse(g, lat1, lon1, lat2, lon2)`
# - `geod_directline(g, lat1, lon1, azi1, s12)`
# - `geod_inverseline(g, lat1, lon1, lat2, lon2)`
# - `geod_position(l, s12)`
# - Many of the `geod_polygon` methods are also wrapped here.

# In addition, we have wrapped some functions for convenience and to minimize allocations:
# - `geod_position(l, s12s::AbstractArray{<:Real})`

# Finally, in order to sample along a geodetic path, we have a convenience function `geod_path(geodesic::geod_geodesic, lat1, lon1, lat2, lon2, npoints)`.

# Simple examples are available in the docstrings for each of these functions, but
# do note that quite a few of them display the C documentation directly.  Look at
# the output of e.g. `methods(geod_geodesic)` for a list of all available methods.

# ## Basic wrappers
# These are some basic wrappers for `geod_direct` and `geod_inverse`.

function geod_direct(g::geod_geodesic, lat::Real, lon::Real, azi::Real, s12::Real)
    lat_out = Ref{Cdouble}(NaN)
    lon_out = Ref{Cdouble}(NaN)
    azi_out = Ref{Cdouble}(NaN)

    geod_direct(Ref(g), lat, lon, azi, s12, lat_out, lon_out, azi_out)

    return lat_out[], lon_out[], azi_out[]
end

function geod_inverse(g::geod_geodesic, lat1::Real, lon1::Real, lat2::Real, lon2::Real)
    s12 = Ref{Cdouble}(NaN)
    azi1 = Ref{Cdouble}(NaN)
    azi2 = Ref{Cdouble}(NaN)

    geod_inverse(Ref(g), lat1, lon1, lat2, lon2, s12, azi1, azi2)

    return s12[], azi1[], azi2[]
end

# ## Constructors and null Constructors

# Since these objects are stack allocated, we create the `_null` function,
# which simply creates an object with its values initialized to `NaN` or 0
# as appropriate.  This is used in the constructors below.

# ### `geod_geodesic`

"""
    Proj._null(geod_*)

A null initializer which returns an object of the given type,
with all floats set to NaN, and all integers set to 0.

Available types are `geod_geodesic`, `geod_geodesicline, `geod_polygon`.
"""
function _null(::Type{geod_geodesic})
    return geod_geodesic(
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        ntuple(_ -> NaN, Val(6)),
        ntuple(_ -> NaN, Val(15)),
        ntuple(_ -> NaN, Val(21)),
    )
end

function geod_geodesic(equatorial_radius::Real, flattening::Real)
    init_obj = _null(geod_geodesic)
    new_objref = Ref(init_obj)
    geod_init(new_objref, Cdouble(equatorial_radius), Cdouble(flattening))
    return new_objref[]
end

# ### `geod_geodesicline`

function _null(::Type{geod_geodesicline})
    return geod_geodesicline(
        (NaN for _ = 1:30)...,
        ntuple(_ -> NaN, 7),
        ntuple(_ -> NaN, 7),
        ntuple(_ -> NaN, 7),
        ntuple(_ -> NaN, 6),
        ntuple(_ -> NaN, 6),
        Cuint(0),
    )
end

function geod_directline(g::geod_geodesic, lat1, lon1, azi1, s12, caps::Cuint = UInt32(0))
    init_obj = _null(geod_geodesicline)
    new_objref = Ref(init_obj)
    geod_directline(new_objref, Ref(g), lat1, lon1, azi1, s12, caps)
    return new_objref[]
end

function geod_inverseline(g::geod_geodesic, lat1, lon1, lat2, lon2, caps::Cuint = UInt32(0))
    init_obj = _null(geod_geodesicline)
    new_objref = Ref(init_obj)

    geod_inverseline(new_objref, Ref(g), lat1, lon1, lat2, lon2, caps)
    return new_objref[]
end

# ## Wrappers for path functions

# These functions wrap path or position calculating functions in GeographicLib.
# in general, `geod_position` accepts real numbers as well as arrays.


# !!! note
#     This returns according to the C order (y, x, az).
#     Do we want this to return by the Julian order (x, y, azi)?
#     If so, should this be a new function?

function geod_position(line::geod_geodesicline, s12::Real)
    lat = Ref{Cdouble}(NaN)
    lon = Ref{Cdouble}(NaN)
    azi = Ref{Cdouble}(NaN)

    geod_position(Ref(line), s12, lat, lon, azi)

    return lat[], lon[], azi[]
end


"""
    geod_position(line::geod_geodesicline, s12::Real)::(lat, lon, azi)
    geod_position(line::geod_geodesicline, s12s::AbstractArray{<: Real})::(lats, lons, azis)

Returns `(lat, lon, azimuth)` at the `s12` distance along the line.
If provided an array, will return three `similar` arrays, which also have

"""
function geod_position(line::geod_geodesicline, s12s::AbstractArray{<:Real})
    result_lon = similar(s12s)
    result_lat = similar(s12s)
    result_azi = similar(s12s)
    lat = Ref{Cdouble}(NaN)
    lon = Ref{Cdouble}(NaN)
    azi = Ref{Cdouble}(NaN)
    for ind in eachindex(s12s)
        geod_position(Ref(line), s12s[ind], lat, lon, azi)
        result_lat[ind] = lat[]
        result_lon[ind] = lon[]
        result_azi[ind] = azi[]
    end
    return result_lat, result_lon, result_azi
end

"""
    geod_position_relative(line::geod_geodesicline, relative_arclength::Real)

Returns `(lat, lon, azimuth)` at the `relative_arclength` between the line's start and end point.
`relative_arclength` can be any real value, but values along the line should be between 0 and 1.

If `relative_arclength` is an Array, then a Tuple of arrays are returned.
"""
function geod_position_relative(
    line::geod_geodesicline,
    relative_arclength::Union{<:Real,<:AbstractArray{<:Real}},
)
    return geod_position(line, line.s13 * relative_arclength)
end

function geod_setdistance(l::geod_geodesicline, s13::Real)
    geod_setdistance(Ref(l), Cdouble(l13))
end

"""
    geod_genposition(l::geod_geodesicline, flags::Union{Cuint, geod_flags}, s12_a12)

Calls the C function `geod_genposition` and returns a tuple of the results, namely:
`(plat2[], plon2[], pazi2[], ps12[], pm12[], pM12[], pM21[], pS12[])`
"""
function geod_genposition(l::geod_geodesicline, flags::Union{Cuint,geod_flags}, s12_a12)
    plat2 = Ref{Cdouble}(NaN)
    plon2 = Ref{Cdouble}(NaN)
    pazi2 = Ref{Cdouble}(NaN)
    ps12 = Ref{Cdouble}(NaN)
    pm12 = Ref{Cdouble}(NaN)
    pM12 = Ref{Cdouble}(NaN)
    pM21 = Ref{Cdouble}(NaN)
    pS12 = Ref{Cdouble}(NaN)

    geod_genposition(
        Ref(l),
        flags,
        s12_a12,
        plat2,
        plon2,
        pazi2,
        ps12,
        pm12,
        pM12,
        pM21,
        pS12,
    )

    return (plat2[], plon2[], pazi2[], ps12[], pm12[], pM12[], pM21[], pS12[])

end


"""
    geod_path(geodesic::geod_geodesic, lat1, lon1, lat2, lon2, npoints = 1000)

Returns a tuple of vectors representing longitude and latitude.

## Example

```julia
geod = Proj.geod_geodesic(6378137, 1/298.257223563)
lats, lons = Proj.geod_path(geod, 40.64, -73.78, 1.36, 103.99)
```
"""
function geod_path(
    geodesic::geod_geodesic,
    lat1,
    lon1,
    lat2,
    lon2,
    npoints = 1000;
    caps = Cuint(0),
)
    @assert npoints > 1

    inverse_line = geod_inverseline(geodesic, lat1, lon1, lat2, lon2, caps)

    lats, lons, azis = geod_position_relative(inverse_line, LinRange(0, 1, npoints))

    return lats, lons
end


# ## Polygon interface

# ### `geod_polygon`

function _null(::Type{geod_polygon})
    return geod_polygon(
        NaN,
        NaN,
        NaN,
        NaN,
        (NaN, NaN),
        (NaN, NaN),
        Cint(0),
        Cint(0),
        Cint(0),
    )
end

function geod_polygon_init(p::geod_polygon, polylinep::Int = 1)
    polylinep_Cint = Cint(polylinep)
    @assert polylinep_Cint == polylinep "Your integer is too large and there will be an under/overflow error.  Please pass an Cint."
    geod_polygon_init(Ref(p), polylinep_Cint)
    return p
end

geod_polygon_clear(p::geod_polygon) = geod_polygon_clear(Ref(p))

function geod_polygon_addpoint(g::geod_geodesic, p::geod_polygon, lat::Real, lon::Real)
    geod_polygon_addpoint(Ref(g), Ref(p), Cdouble(lat), Cdouble(lon))
    return p
end

function geod_polygon_addedge(g::geod_geodesic, p::geod_polygon, azi::Real, s::Real)
    geod_polygon_addedge(Ref(g), Ref(p), Cdouble(azi), Cdouble(s))
    return p
end

"""
    geod_polygon_compute(g::geod_geodesic, p::geod_polygon, reverse::Bool = false, sign::Bool = true)

Returns a tuple of `(area, perimeter)` in m² and m respectively.  Area is only returned if `polyline` is nonzero in the call to `geod_polygon_init`.
"""
function geod_polygon_compute(
    g::geod_geodesic,
    p::geod_polygon,
    reverse::Bool = false,
    sign::Bool = true,
)
    # initializing to zero makes the return not happen for that value, so we need to initialize to 1
    pA = Ref{Cdouble}(1e0)
    pP = Ref{Cdouble}(1e0)

    n = geod_polygon_compute(Ref(g), Ref(p), Cint(reverse), Cint(sign), pA, pP)

    return (pA[], pP[])
end

"""
    geod_polygonarea(g::geod_geodesic, lats::AbstractVector{<: Real}, lons::AbstractVector{<: Real})

    Simple interface to compute the geodesic area of a polygon.
Returns a tuple of `(area, perimeter)` in m² and m respectively.
"""
function geod_polygonarea(
    g::geod_geodesic,
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
)
    @assert length(lats) == length(lons)
    c_lats = convert(Vector{Cdouble}, lats)
    c_lons = convert(Vector{Cdouble}, lons)

    pA = Ref{Cdouble}(1e0)
    pP = Ref{Cdouble}(1e0)

    geod_polygonarea(Ref(g), c_lats, c_lons, length(lats), pA, pP)

    return (pA[], pP[])

end

# TODO: add geod_polygon_addpoint, geod_polygon_addedge, geod_polygon_testedge, geod_polygon_testpoint

# ## GeoInterface wrappers
# What follows are GeoInterface-based wrappers so it's easy to pass through points.  
# ### geod_direct
geod_direct(point, azi::Real, s12::Real; geodesic = geod_geodesic(6378137, 1/298/257223563)) = geod_direct(GeoInterface.geomtrait(point), point, azi, s12; geodesic)
function geod_direct(::GeoInterface.PointTrait, point, azi, s12; geodesic = geod_geodesic(6378137, 1/298/257223563))
    return geod_direct(geodesic, GeoInterface.y(p1), GeoInterface.x(p1), azi, s2)
end
# ### geod_inverse
geod_inverse(p1, p2; geodesic = geod_geodesic(6378137, 1/298/257223563)) = geod_inverse(GeoInterface.geomtrait(p1), GeoInterface.geomtrait(p2), p1, p2; geodesic)
function geod_inverse(::GeoInterface.PointTrait, ::GeoInterface.PointTrait, p1, p2; geodesic = geod_geodesic(6378137, 1/298/257223563))
    return geod_inverse(geodesic, GeoInterface.y(p1), GeoInterface.x(p1), GeoInterface.y(p2), GeoInterface.x(p2))
end
# ### geod_directline
geod_directline(point, azi::Real, s12::Real, caps::Cuint = UInt32(0); geodesic = geod_geodesic(6378137, 1/298/257223563)) = geod_directline(GeoInterface.geomtrait(point), point, azi, s12, caps; geodesic)
function geod_directline(::GeoInterface.PointTrait, point, azi, s12, caps::Cuint = UInt32(0); geodesic = geod_geodesic(6378137, 1/298/257223563))
    return geod_directline(geodesic, GeoInterface.y(p1), GeoInterface.x(p1), azi, s2, caps)
end
# ### geod_inverseline
geod_inverseline(p1, p2, caps::Cuint = UInt32(0); geodesic = geod_geodesic(6378137, 1/298/257223563)) = geod_inverse(GeoInterface.geomtrait(p1), GeoInterface.geomtrait(p2), p1, p2, caps; geodesic)
function geod_inverse(::GeoInterface.PointTrait, ::GeoInterface.PointTrait, p1, p2, caps::Cuint = UInt32(0); geodesic = geod_geodesic(6378137, 1/298/257223563))
    return geod_inverse(geodesic, GeoInterface.y(p1), GeoInterface.x(p1), GeoInterface.y(p2), GeoInterface.x(p2), caps)
end
# ### geod_path
function geod_path(p1, p2, npoints = 1000; geodesic = geod_geodesic(6378137, 1/298/257223563), caps = UInt32(0))
    geod_path(GeoInterface.geomtrait(p1), GeoInterface.geomtrait(p2), p1, p2, npoints; geodesic, caps)
end
function geod_path(::GeoInterface.PointTrait, ::GeoInterface.PointTrait, p1, p2, npoints = 1000; geodesic = geod_geodesic(6378137, 1/298/257223563), caps = UInt32(0))
    geod_path(geodesic, GeoInterface.y(p1), GeoInterface.x(p1), GeoInterface.y(p2), GeoInterface.x(p2), npoints; caps)
end
