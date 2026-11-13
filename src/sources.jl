"""
    area_source(rate::Real, area_density::Real)

Convert a scalar transfer per unit one-sided plant surface area and second to a
bulk canopy volume source by multiplying by `area_density` once.

When `exact_zero_scalar(area_density)` holds (both primal and derivative partials
are zero), return exact `zero(FT)` without multiplying `rate`, guarding empty
layers against nonfinite dummy element rates on CPU, GPU, and AD paths. When the
primal of `area_density` is zero with a nonzero AD tangent (`Dual(0.0, 1.0)`),
evaluate `rate * area_density` directly (requiring finite `rate`) so the right
directional derivative `∂(rate * a)/∂a = rate` is preserved.

# Examples
```jldoctest
julia> area_source(150.0, 0.8)
120.0

julia> area_source(NaN, 0.0)
0.0
```
"""
@inline function area_source(rate::FT, area_density::FT) where {FT<:Real}
    require_continuous_scalar(FT, area_source, (rate, area_density))
    return exact_zero_scalar(area_density) ? zero(FT) : rate * area_density
end

@inline function area_source(rate::Real, area_density::Real)
    FT = promote_compatible_scalar(area_source, rate, area_density)
    return area_source(convert(FT, rate), convert(FT, area_density))
end

"""
    leaf_area_source(leaf_rate::Real, geometry::AbstractLayerGeometry)

Convert a scalar transfer per m² one-sided leaf area and second to a bulk canopy
volume source [m⁻³ s⁻¹] by multiplying `geometry.leaf_area_density` once.

Positive values denote transfer into canopy air; negate the result for the
paired leaf ledger. For zero leaf-area density (with zero AD partials), return
exact zero even when `leaf_rate` is nonfinite.

# Examples
```jldoctest
julia> leaf_area_source(2.0e-5, LayerGeometry(1.5, 0.5))
1.0e-5
```
"""
@inline function leaf_area_source(rate::Real, geometry::AbstractLayerGeometry)
    return area_source(rate, geometry.leaf_area_density)
end

"""
    stem_area_source(stem_rate::Real, geometry::PlantLayerGeometry)

Convert a scalar transfer per m² one-sided stem area and second to a bulk canopy
volume source [m⁻³ s⁻¹] by multiplying `geometry.stem_area_density` once.

For zero stem-area density (with zero AD partials), return exact zero even when
`stem_rate` is nonfinite.

# Examples
```jldoctest
julia> stem_area_source(40.0, PlantLayerGeometry(2.0, 1.2, 0.25))
10.0
```
"""
@inline function stem_area_source(rate::Real, geometry::PlantLayerGeometry)
    return area_source(rate, geometry.stem_area_density)
end

"""
    plant_area_density(geometry::AbstractLayerGeometry)

Return total one-sided plant area density (`LAD` for [`LayerGeometry`](@ref), or
`LAD + SAD` for [`PlantLayerGeometry`](@ref)) in m² plant area per m³ bulk canopy.

# Examples
```jldoctest
julia> (plant_area_density(LayerGeometry(1.5, 0.8)), plant_area_density(PlantLayerGeometry(1.5, 0.8, 0.2)))
(0.8, 1.0)
```

See also [`LayerGeometry`](@ref), [`PlantLayerGeometry`](@ref), [`plant_air_sources`](@ref).
"""
function plant_area_density end
@inline plant_area_density(geometry::LayerGeometry) = geometry.leaf_area_density
@inline plant_area_density(geometry::PlantLayerGeometry) =
    geometry.leaf_area_density + geometry.stem_area_density

@inline _oxygen_area_source(::NoOxygen{FO}, area_density::FA) where {FO,FA} =
    NoOxygen{promote_compatible_scalar(_oxygen_area_source, zero(FO), area_density)}()
@inline function _oxygen_area_source(
    oxygen::OxygenTransfer{FO},
    area_density::FA,
) where {FO,FA}
    return OxygenTransfer(area_source(oxygen.value, area_density))
end

"""
    sensible_water_flux(sensible_heat, water_flux, carbon_flux, oxygen, vapor_enthalpy)
    sensible_water_flux(sensible_heat, water_flux, carbon_flux, vapor_enthalpy)

Construct a [`LeafFlux`](@ref) from sensible heat [W m⁻² leaf], signed water flux
[kg H₂O m⁻² leaf s⁻¹], signed `CO₂` flux [mol C m⁻² leaf s⁻¹], optional `oxygen`,
and transported water-vapor specific enthalpy `vapor_enthalpy` [J kg⁻¹].

Total leaf-to-air energy transfer is `sensible_heat + water_flux * vapor_enthalpy`.
Positive fluxes enter canopy air (photosynthetic carbon uptake is negative
`carbon_flux`; condensation is negative `water_flux`). Do not add a separate
latent-heat source when coupling this record to a conserved air energy budget.

# Examples
```jldoctest
julia> f = sensible_water_flux(80.0, 4.0e-5, -1.2e-5, 2.5e6);

julia> f.rates.energy
180.0
```
"""
@inline function sensible_water_flux(
    sensible_heat::FT,
    water_flux::FT,
    carbon_flux::FT,
    oxygen::OxygenContent{FT},
    vapor_enthalpy::FT,
) where {FT<:Real}
    require_continuous_scalar(
        FT,
        sensible_water_flux,
        (sensible_heat, water_flux, carbon_flux, oxygen, vapor_enthalpy),
    )
    return LeafFlux(
        ExchangeRates(
            sensible_heat + water_flux * vapor_enthalpy,
            water_flux,
            carbon_flux,
            oxygen,
        ),
    )
end

@inline function sensible_water_flux(
    sensible_heat::Real,
    water_flux::Real,
    carbon_flux::Real,
    oxygen::O,
    vapor_enthalpy::Real,
) where {FO<:Real,O<:OxygenContent{FO}}
    FT = promote_compatible_scalar(
        sensible_water_flux,
        sensible_heat,
        water_flux,
        carbon_flux,
        zero(FO),
        vapor_enthalpy,
    )
    OT = _promote_oxygen_type(O, FT)
    return sensible_water_flux(
        convert(FT, sensible_heat),
        convert(FT, water_flux),
        convert(FT, carbon_flux),
        convert(OT, oxygen),
        convert(FT, vapor_enthalpy),
    )
end

@inline function sensible_water_flux(
    sensible_heat::Real,
    water_flux::Real,
    carbon_flux::Real,
    vapor_enthalpy::Real,
)
    FT = promote_compatible_scalar(
        sensible_water_flux,
        sensible_heat,
        water_flux,
        carbon_flux,
        vapor_enthalpy,
    )
    return sensible_water_flux(
        convert(FT, sensible_heat),
        convert(FT, water_flux),
        convert(FT, carbon_flux),
        NoOxygen{FT}(),
        convert(FT, vapor_enthalpy),
    )
end

"""
    sensible_stem_flux(sensible_heat, oxygen=NoOxygen(FT))
    sensible_stem_flux(sensible_heat, water_flux, carbon_flux, oxygen, vapor_enthalpy)
    sensible_stem_flux(sensible_heat, water_flux, carbon_flux, vapor_enthalpy)

Construct a [`StemFlux`](@ref) per m² one-sided stem area and second.

The 1- and 2-argument forms construct a dry, non-assimilating woody sensible-heat
flux (`water = 0`, `carbon = 0`). When nonzero stem evaporation or respiration is
supplied via the 4- or 5-argument form, `vapor_enthalpy` [J kg⁻¹] is mandatory so
transported water enthalpy `water_flux * vapor_enthalpy` is never silently omitted.

# Examples
```jldoctest
julia> s = sensible_stem_flux(35.0)
StemFlux{Float64, NoOxygen{Float64}}(ExchangeRates{Float64, NoOxygen{Float64}}(35.0, 0.0, 0.0, NoOxygen{Float64}()))
```
"""
@inline function sensible_stem_flux(sensible_heat::FT) where {FT<:Real}
    require_continuous_scalar(FT, sensible_stem_flux, (sensible_heat,))
    z = zero(FT)
    return StemFlux(ExchangeRates(sensible_heat, z, z, NoOxygen{FT}()))
end

@inline function sensible_stem_flux(
    sensible_heat::Real,
    oxygen::O,
) where {FO<:Real,O<:OxygenContent{FO}}
    FT = promote_compatible_scalar(sensible_stem_flux, sensible_heat, zero(FO))
    OT = _promote_oxygen_type(O, FT)
    z = zero(FT)
    return StemFlux(ExchangeRates(convert(FT, sensible_heat), z, z, convert(OT, oxygen)))
end

@inline function sensible_stem_flux(
    sensible_heat::Real,
    water_flux::Real,
    carbon_flux::Real,
    oxygen::O,
    vapor_enthalpy::Real,
) where {FO<:Real,O<:OxygenContent{FO}}
    FT = promote_compatible_scalar(
        sensible_stem_flux,
        sensible_heat,
        water_flux,
        carbon_flux,
        zero(FO),
        vapor_enthalpy,
    )
    OT = _promote_oxygen_type(O, FT)
    sh = convert(FT, sensible_heat)
    wf = convert(FT, water_flux)
    cf = convert(FT, carbon_flux)
    hv = convert(FT, vapor_enthalpy)
    return StemFlux(ExchangeRates(sh + wf * hv, wf, cf, convert(OT, oxygen)))
end

@inline function sensible_stem_flux(
    sensible_heat::Real,
    water_flux::Real,
    carbon_flux::Real,
    vapor_enthalpy::Real,
)
    FT = promote_compatible_scalar(
        sensible_stem_flux,
        sensible_heat,
        water_flux,
        carbon_flux,
        vapor_enthalpy,
    )
    return sensible_stem_flux(
        convert(FT, sensible_heat),
        convert(FT, water_flux),
        convert(FT, carbon_flux),
        NoOxygen{FT}(),
        convert(FT, vapor_enthalpy),
    )
end

"""
    leaf_air_sources(flux::LeafFlux, geometry::AbstractLayerGeometry) -> (air, leaf)

Return paired `(air, leaf)` [`VolumeSource`](@ref) ledgers by multiplying `flux` by `LAD` once.

`leaf` equals `-air` with IEEE sign-bit negation, representing an equal-and-opposite
exchange ledger rather than independent leaf gas storage or leaf temperature
rate of change. When `exact_zero_scalar(geometry.leaf_area_density)` holds, return
exact finite zeros even if `flux` holds nonfinite dummy values.

# Examples
```jldoctest
julia> f = LeafFlux(ExchangeRates(100.0, 2.0e-5, -1.0e-5));

julia> p = leaf_air_sources(f, LayerGeometry(2.0, 0.5));

julia> (p.air.rates.energy, p.leaf.rates.energy)
(50.0, -50.0)
```
"""
@inline function leaf_air_sources(
    flux::LeafFlux{FF,O},
    geometry::AbstractLayerGeometry{FG},
) where {FF<:Real,O<:OxygenContent{FF},FG<:Real}
    FT = promote_compatible_scalar(leaf_air_sources, zero(FF), zero(FG))
    lad = convert(FT, geometry.leaf_area_density)
    r = convert(ExchangeRates{FT,_promote_oxygen_type(O, FT)}, flux.rates)
    rates = ExchangeRates(
        area_source(r.energy, lad),
        area_source(r.water, lad),
        area_source(r.carbon, lad),
        _oxygen_area_source(r.oxygen, lad),
    )
    return (air = VolumeSource(rates), leaf = VolumeSource(-rates))
end

"""
    plant_air_sources(leaf_flux::LeafFlux, stem_flux::StemFlux, geometry::PlantLayerGeometry)

Return conservative `(air, leaf, stem)` [`VolumeSource`](@ref) ledgers for a leaf-and-stem layer.

`leaf_flux` scales by `geometry.leaf_area_density`, `stem_flux` scales by
`geometry.stem_area_density` with independent zero-density guards, `air` equals
`leaf_air + stem_air`, and `leaf` and `stem` hold the opposite donor ledgers so
`air + leaf + stem == zero(air)`.

# Examples
```jldoctest
julia> lf = LeafFlux(ExchangeRates(100.0, 2.0e-5, -1.0e-5));

julia> sf = sensible_stem_flux(40.0);

julia> s = plant_air_sources(lf, sf, PlantLayerGeometry(2.0, 1.0, 0.25));

julia> (s.air.rates.energy, s.leaf.rates.energy, s.stem.rates.energy)
(110.0, -100.0, -10.0)
```
"""
@inline function plant_air_sources(
    leaf_flux::LeafFlux{FL,OL},
    stem_flux::StemFlux{FS,OS},
    geometry::PlantLayerGeometry{FG},
) where {FL<:Real,OL<:OxygenContent{FL},FS<:Real,OS<:OxygenContent{FS},FG<:Real}
    compatible_oxygen(OL, OS) ||
        throw(MethodError(plant_air_sources, (leaf_flux, stem_flux, geometry)))
    FT = promote_compatible_scalar(plant_air_sources, zero(FL), zero(FS), zero(FG))
    OT = _promote_oxygen_type(OL, FT)
    lad = convert(FT, geometry.leaf_area_density)
    sad = convert(FT, geometry.stem_area_density)
    lr = convert(ExchangeRates{FT,OT}, leaf_flux.rates)
    sr = convert(ExchangeRates{FT,OT}, stem_flux.rates)
    leaf_air = ExchangeRates(
        area_source(lr.energy, lad),
        area_source(lr.water, lad),
        area_source(lr.carbon, lad),
        _oxygen_area_source(lr.oxygen, lad),
    )
    stem_air = ExchangeRates(
        area_source(sr.energy, sad),
        area_source(sr.water, sad),
        area_source(sr.carbon, sad),
        _oxygen_area_source(sr.oxygen, sad),
    )
    return (
        air = VolumeSource(leaf_air + stem_air),
        leaf = VolumeSource(-leaf_air),
        stem = VolumeSource(-stem_air),
    )
end

"""
    conserved_tendency(bottom_rate, top_rate, source, geometry::AbstractLayerGeometry)

Compute scalar finite-volume conserved tendency `(bottom_rate - top_rate) / thickness + source`.

`bottom_rate` and `top_rate` are upward-positive fluxes per common horizontal
ground area, and `source` is per bulk canopy volume. Summing thickness-weighted
tendencies across stacked layers cancels interior face fluxes identically.

# Examples
```jldoctest
julia> conserved_tendency(15.0, 12.0, 4.0, LayerGeometry(1.5, 0.8))
6.0
```
"""
@inline function conserved_tendency(
    bottom::FT,
    top::FT,
    source::FT,
    geometry::AbstractLayerGeometry{FT},
) where {FT<:Real}
    require_continuous_scalar(FT, conserved_tendency, (bottom, top, source, geometry))
    return (bottom - top) / geometry.thickness + source
end

@inline function conserved_tendency(
    bottom::Real,
    top::Real,
    source::Real,
    geometry::AbstractLayerGeometry{FG},
) where {FG<:Real}
    FT = promote_compatible_scalar(conserved_tendency, bottom, top, source, zero(FG))
    dz = convert(FT, geometry.thickness)
    return (convert(FT, bottom) - convert(FT, top)) / dz + convert(FT, source)
end

@inline _oxygen_tendency(
    ::NoOxygen{FT},
    ::NoOxygen{FT},
    ::NoOxygen{FT},
    dz::FT,
) where {FT} = NoOxygen{FT}()

@inline function _oxygen_tendency(
    bottom::OxygenTransfer{FT},
    top::OxygenTransfer{FT},
    source::OxygenTransfer{FT},
    dz::FT,
) where {FT}
    return OxygenTransfer{FT}((bottom.value - top.value) / dz + source.value)
end

"""
    layer_tendency(bottom::FaceFlux, top::FaceFlux, source::VolumeSource, geometry::AbstractLayerGeometry)

Compute conservative bulk-volume [`VolumeSource`](@ref) tendency `(bottom - top) / thickness + source`.

All face fluxes are upward-positive across shared horizontal ground area, and all
inputs must share compatible continuous scalar precision and matching oxygen mode.
Summing `thickness * tendency` across adjacent layers telescopes interior faces.

# Examples
```jldoctest
julia> b = FaceFlux(ExchangeRates(30.0, 0.0, 0.0));

julia> t = FaceFlux(ExchangeRates(24.0, 0.0, 0.0));

julia> s = VolumeSource(ExchangeRates(10.0, 0.0, 0.0));

julia> layer_tendency(b, t, s, LayerGeometry(2.0, 0.5)).rates.energy
13.0
```
"""
@inline function layer_tendency(
    bottom::FaceFlux{FT,O},
    top::FaceFlux{FT,O},
    source::VolumeSource{FT,O},
    geometry::AbstractLayerGeometry{FT},
) where {FT<:Real,O<:OxygenContent{FT}}
    dz = geometry.thickness
    return VolumeSource(
        ExchangeRates(
            (bottom.rates.energy - top.rates.energy) / dz + source.rates.energy,
            (bottom.rates.water - top.rates.water) / dz + source.rates.water,
            (bottom.rates.carbon - top.rates.carbon) / dz + source.rates.carbon,
            _oxygen_tendency(
                bottom.rates.oxygen,
                top.rates.oxygen,
                source.rates.oxygen,
                dz,
            ),
        ),
    )
end

@inline function layer_tendency(
    bottom::FaceFlux{FB,OB},
    top::FaceFlux{FTp,OTp},
    source::VolumeSource{FS,OS},
    geometry::AbstractLayerGeometry{FG},
) where {
    FB<:Real,
    OB<:OxygenContent{FB},
    FTp<:Real,
    OTp<:OxygenContent{FTp},
    FS<:Real,
    OS<:OxygenContent{FS},
    FG<:Real,
}
    (compatible_oxygen(OB, OTp) && compatible_oxygen(OB, OS)) ||
        throw(MethodError(layer_tendency, (bottom, top, source, geometry)))
    FT = promote_compatible_scalar(
        layer_tendency,
        zero(FB),
        zero(FTp),
        zero(FS),
        zero(FG),
    )
    OT = _promote_oxygen_type(OB, FT)
    bb = convert(FaceFlux{FT,OT}, bottom)
    tt = convert(FaceFlux{FT,OT}, top)
    ss = convert(VolumeSource{FT,OT}, source)
    dz = convert(FT, geometry.thickness)
    return VolumeSource(
        ExchangeRates(
            (bb.rates.energy - tt.rates.energy) / dz + ss.rates.energy,
            (bb.rates.water - tt.rates.water) / dz + ss.rates.water,
            (bb.rates.carbon - tt.rates.carbon) / dz + ss.rates.carbon,
            _oxygen_tendency(bb.rates.oxygen, tt.rates.oxygen, ss.rates.oxygen, dz),
        ),
    )
end
