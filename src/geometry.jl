# ============================================================================
# Layer geometry: any type with `thickness(g)` and `area_density(g, M)` is one
# ============================================================================

"""
    LayerGeometry(thickness, leaf_area_density)

Finite-volume layer with finite positive `thickness` [m] and finite nonnegative
one-sided `leaf_area_density` [m² leaf m⁻³ bulk canopy volume]. Bulk volume is the
ground footprint times `thickness`, with no implicit porosity factor. Throws an
`ArgumentError` for a nonpositive or nonfinite thickness, a negative or nonfinite
density, or a nonfinite derivative component.

A layer geometry is any type that defines [`thickness`](@ref), [`area_density`](@ref)
and `Base.eltype`; this one has leaves only, so a stem flux on it is a `MethodError`.

# Examples
```jldoctest
julia> g = LayerGeometry(1.5, 0.8);

julia> (thickness(g), area_density(g, LeafArea))
(1.5, 0.8)
```
"""
struct LayerGeometry{FT<:Real} <: ScalarRecord
    thickness::FT
    leaf_area_density::FT
    global function _layer_geometry(dz::FT, lad::FT) where {FT}
        _checked(dz, _positive, "layer thickness must be finite and positive")
        _checked(lad, _nonnegative, "leaf area density must be finite and nonnegative")
        return new{FT}(dz, lad)
    end
end
@inline function LayerGeometry(dz::Real, lad::Real)
    same_precision(dz, lad)
    return _layer_geometry(promote(dz, lad)...)
end

"""
    PlantLayerGeometry(thickness, leaf_area_density, stem_area_density = 0)

As [`LayerGeometry`](@ref), with finite nonnegative one-sided `stem_area_density`
[m² stem m⁻³ bulk canopy volume] (`PAI = LAI + SAI` as in ClimaLand). The two-argument
form sets the stem area density to zero of the leaf area density's type.

# Examples
```jldoctest
julia> g = PlantLayerGeometry(2.0, 1.2, 0.3);

julia> (area_density(g, StemArea), plant_area_density(g))
(0.3, 1.5)
```
"""
struct PlantLayerGeometry{FT<:Real} <: ScalarRecord
    thickness::FT
    leaf_area_density::FT
    stem_area_density::FT
    global function _plant_layer_geometry(dz::FT, lad::FT, sad::FT) where {FT}
        _checked(dz, _positive, "layer thickness must be finite and positive")
        _checked(lad, _nonnegative, "leaf area density must be finite and nonnegative")
        _checked(sad, _nonnegative, "stem area density must be finite and nonnegative")
        return new{FT}(dz, lad, sad)
    end
end
@inline function PlantLayerGeometry(dz::Real, lad::Real, sad::Real)
    same_precision(dz, lad, sad)
    return _plant_layer_geometry(promote(dz, lad, sad)...)
end
@inline PlantLayerGeometry(dz::Real, lad::Real) = PlantLayerGeometry(dz, lad, zero(lad))

const Geometry = Union{LayerGeometry,PlantLayerGeometry}
Base.eltype(::Type{<:LayerGeometry{FT}}) where {FT} = FT
Base.eltype(::Type{<:PlantLayerGeometry{FT}}) where {FT} = FT

# Container promotion between geometries that differ only by derivative information.
for G in (:LayerGeometry, :PlantLayerGeometry)
    @eval begin
        Base.promote_rule(::Type{$G{A}}, ::Type{$G{B}}) where {A,B} =
            _promotable(A, B) ? $G{promote_type(A, B)} : Union{}
        function Base.convert(::Type{$G{FT}}, g::$G) where {FT}
            same_precision(zero(FT), g)
            return $G(map(Base.Fix1(convert, FT), _fields(g))...)
        end
    end
end

"""
    thickness(geometry)

Layer thickness [m] of a layer geometry. Part of the geometry interface; a new
geometry type either has a `thickness` field or adds a method.
"""
@inline thickness(g) = g.thickness

"""
    area_density(geometry, ::Type{M}) where {M<:Measure}

One-sided area density [m² element area m⁻³ bulk canopy volume] of the element class
with measure `M` in `geometry`. Part of the geometry interface: a geometry provides a
method for every element class it contains, and a flux of a class it does not contain
is a `MethodError` in [`volume_source`](@ref).
"""
function area_density end
@inline area_density(g::Geometry, ::Type{LeafArea}) = g.leaf_area_density
@inline area_density(g::PlantLayerGeometry, ::Type{StemArea}) = g.stem_area_density

"""
    plant_area_density(geometry)

Total one-sided plant area density [m² m⁻³]: `LAD` for a [`LayerGeometry`](@ref),
`LAD + SAD` for a [`PlantLayerGeometry`](@ref) (`PAI_k / Δz_k`).

# Examples
```jldoctest
julia> (plant_area_density(LayerGeometry(1.5, 0.8)), plant_area_density(PlantLayerGeometry(1.5, 0.8, 0.2)))
(0.8, 1.0)
```
"""
function plant_area_density end
@inline plant_area_density(g::LayerGeometry) = g.leaf_area_density
@inline plant_area_density(g::PlantLayerGeometry) = g.leaf_area_density + g.stem_area_density

"""
    layer_geometry(thickness, leaf_area_density) -> LayerGeometry
    layer_geometry(thickness, leaf_area_density, stem_area_density) -> PlantLayerGeometry

Function form of the geometry constructors for broadcasts over host fields, where a
type in the function position is rejected. Units and checks as in
[`LayerGeometry`](@ref) and [`PlantLayerGeometry`](@ref).
"""
@inline layer_geometry(dz, lad) = LayerGeometry(dz, lad)
@inline layer_geometry(dz, lad, sad) = PlantLayerGeometry(dz, lad, sad)
