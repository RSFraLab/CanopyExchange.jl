# ============================================================================
# Element area to bulk volume: one multiplication by the area density
# ============================================================================

"""
    area_source(rate, area_density)

Convert a scalar transfer per m² one-sided element area and second into a source per
m³ bulk canopy volume and second by multiplying by `area_density` [m² m⁻³] once.

Empty-layer rule: when `area_density` is an exact zero (value and every derivative
component zero, see [`CanopyExchange.exact_zero`](@ref)) the result is an exact zero,
even when the unused element `rate` is `NaN` or `Inf`. A zero value with a nonzero
derivative component (`Dual(0.0, 1.0)`) evaluates `rate * area_density`, which keeps
the right derivative `∂(rate a)/∂a = rate` and requires a finite `rate`.

# Examples
```jldoctest
julia> area_source(150.0, 0.8)
120.0

julia> area_source(NaN, 0.0)
0.0
```
"""
@inline function area_source(rate, density)
    same_precision(rate, density)
    source = rate * density
    return ifelse(exact_zero(density), zero(source), source)
end

"""
    volume_source(flux::Exchange{M}, geometry) -> VolumeSource

Volume source of one element class: every component of `flux` times the area density
`area_density(geometry, M)` of its own class, through [`area_source`](@ref) (so the
empty-layer rule applies). Positive into canopy air. A flux of a class the geometry
does not contain is a `MethodError`; a [`FaceFlux`](@ref) or a
[`VolumeSource`](@ref) is not an element flux.

# Examples
```jldoctest
julia> s = volume_source(leaf_flux(100.0, 2.0e-5, (; co2 = -1.0e-5)), LayerGeometry(2.0, 0.5));

julia> (s.rates.energy, s.rates.tracers.co2)
(50.0, -5.0e-6)
```
"""
@inline function volume_source(flux::Exchange{M}, g) where {M}
    same_precision(flux, thickness(g))
    return VolumeSource(map(Base.Fix2(area_source, area_density(g, M)), flux.rates))
end

"""
    leaf_area_source(leaf_rate, geometry)

[`area_source`](@ref) of one scalar leaf rate [per m² leaf s⁻¹] with the leaf area
density of `geometry`: a bulk-volume source [per m³ s⁻¹], positive into canopy air.
Negate it for the leaf donor ledger.

# Examples
```jldoctest
julia> leaf_area_source(2.0e-5, LayerGeometry(1.5, 0.5))
1.0e-5
```
"""
@inline leaf_area_source(rate, g) = area_source(rate, area_density(g, LeafArea))

"""
    stem_area_source(stem_rate, geometry)

[`area_source`](@ref) of one scalar stem rate [per m² stem s⁻¹] with the stem area
density of `geometry`. A geometry without stems is a `MethodError`.

# Examples
```jldoctest
julia> stem_area_source(40.0, PlantLayerGeometry(2.0, 1.2, 0.25))
10.0
```
"""
@inline stem_area_source(rate, g) = area_source(rate, area_density(g, StemArea))

"""
    leaf_air_sources(leaf::LeafFlux, geometry) -> (; air, leaf)

Paired [`VolumeSource`](@ref) ledgers of a leaf layer: `air` is the leaf flux times the
leaf area density ([`volume_source`](@ref)), and the donor ledger `leaf` is its exact
IEEE negation, so `air + leaf` is exactly zero componentwise. The donor ledger records
transfer out of the leaf surface; it is not leaf gas storage or a leaf temperature
tendency. An empty layer returns exact zeros even for nonfinite dummy leaf rates.

# Examples
```jldoctest
julia> p = leaf_air_sources(leaf_flux(100.0, 2.0e-5, (; co2 = -1.0e-5)), LayerGeometry(2.0, 0.5));

julia> (p.air.rates.energy, p.leaf.rates.energy)
(50.0, -50.0)
```
"""
@inline function leaf_air_sources(leaf::LeafFlux, g)
    air = volume_source(leaf, g)
    return (; air, leaf = -air)
end

"""
    plant_air_sources(leaf::LeafFlux, stem::StemFlux, geometry) -> (; air, leaf, stem)

Paired ledgers of a leaf-and-stem layer. The leaf flux is scaled by the leaf area
density and the stem flux by the stem area density, each with its own empty-layer
rule. `air` is the sum of both contributions; `leaf` and `stem` are the exact
negations of their own contributions, so `air + (leaf + stem)` is zero to round-off
(one rounding in the sum `air`). Leaf and stem fluxes must carry the same tracers. A
geometry without stems is a `MethodError`.

# Examples
```jldoctest
julia> lf = leaf_flux(100.0, 2.0e-5, (; co2 = -1.0e-5));

julia> sf = sensible_stem_flux(40.0, (; co2 = 0.0));

julia> s = plant_air_sources(lf, sf, PlantLayerGeometry(2.0, 1.0, 0.25));

julia> (s.air.rates.energy, s.leaf.rates.energy, s.stem.rates.energy)
(110.0, -100.0, -10.0)
```
"""
@inline function plant_air_sources(leaf::LeafFlux, stem::StemFlux, g)
    from_leaf, from_stem = volume_source(leaf, g), volume_source(stem, g)
    return (; air = from_leaf + from_stem, leaf = -from_leaf, stem = -from_stem)
end
