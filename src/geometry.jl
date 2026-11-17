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

It implements the geometry interface ([`CanopyExchange.area_density`](@ref)) for
leaves only, so a stem flux on it is a `MethodError`.

# Examples
```jldoctest
julia> g = LayerGeometry(1.5, 0.8);

julia> (CanopyExchange.thickness(g), CanopyExchange.area_density(g, LeafArea))
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
    PlantLayerGeometry(thickness, leaf_area_density, stem_area_density)
    PlantLayerGeometry(thickness, leaf_area_density)

As [`LayerGeometry`](@ref), with finite nonnegative one-sided `stem_area_density`
[m² stem m⁻³ bulk canopy volume] (`PAI = LAI + SAI` as in ClimaLand). The
two-argument form sets the stem area density to `zero(leaf_area_density)`.

# Examples
```jldoctest
julia> g = PlantLayerGeometry(2.0, 1.2, 0.3);

julia> (CanopyExchange.area_density(g, StemArea), plant_area_density(g))
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
    CanopyExchange.thickness(geometry)

Layer thickness [m] of a layer geometry; part of the geometry interface (see
[`CanopyExchange.area_density`](@ref)). The default reads a `thickness` field.
"""
@inline thickness(g) = g.thickness

"""
    CanopyExchange.area_density(geometry, ::Type{M}) where {M<:Measure}

One-sided area density [m² element area m⁻³ bulk canopy volume] of the element class
with measure `M` in `geometry`.

**Geometry interface.** A layer geometry is any type `G` that provides

- [`CanopyExchange.thickness`](@ref)`(g::G)` [m] (or a field `thickness`), and
- `CanopyExchange.area_density(g::G, ::Type{M})` for every element class `M` it holds.

Both are `public`, not exported: extend them with qualified names. A flux of a class
the geometry does not hold is a `MethodError` in [`volume_source`](@ref).
`Base.eltype` is not required. To use a geometry as a scalar in broadcasts
(`volume_source.(fluxes, g)`), define `Base.Broadcast.broadcastable(g::G) = Ref(g)`
or pass `Ref(g)`. The package cannot validate an external geometry: it must itself
guarantee a finite positive thickness and finite nonnegative area densities, of the
same floating-point type as the fluxes (the kernels check the type, not the values).

# Examples
```jldoctest
julia> struct Orchard{FT}; thickness::FT; leaf_area_density::FT; end

julia> CanopyExchange.area_density(g::Orchard, ::Type{LeafArea}) = g.leaf_area_density;

julia> volume_source(leaf_flux(100.0, 2.0e-5), Orchard(2.0, 0.5)).rates.energy
50.0
```
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
