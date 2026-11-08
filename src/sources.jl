"""
    leaf_area_source(leaf_rate::FT, geometry::LayerGeometry{FT})

Convert a scalar transfer per one-sided leaf area per second to a source per
bulk canopy volume per second, multiplying LAD once. Positive means into air.
For zero LAD, return exact zero without using an invalid absent-leaf rate.

This extension point supports additional species and host state layouts without
changing `ExchangeRates`. The host declares the species and consistent mass or
molar units; scalar calls do not provide the wrapper types' measure protection.
Use the negative result for the counterparty exchange ledger.
"""
@inline function leaf_area_source(rate::FT, geometry::LayerGeometry{FT}) where {FT}
    return iszero(geometry.leaf_area_density) ? zero(FT) : rate * geometry.leaf_area_density
end

"""
    conserved_tendency(bottom_rate::FT, top_rate::FT, source::FT, geometry)

Scalar conserved-density tendency `(bottom_rate - top_rate)/thickness + source`.
Face rates use the same species/units per common horizontal ground area and are
upward-positive; source uses those units per bulk volume. This is an extension
point for additional conserved tracers in host-owned fields. It does not directly
advance temperature, humidity or mixing ratio. Host transformations, metrics and
thermodynamics remain explicit. Applying the scalar algebra to momentum does not
supply its kinetic-energy/TKE contract. All inputs and geometry use one FT.
"""
@inline function conserved_tendency(
    bottom::FT,
    top::FT,
    source::FT,
    geometry::LayerGeometry{FT},
) where {FT}
    return (bottom - top) / geometry.thickness + source
end

@inline oxygen_source(::NoOxygen{FT}, ::LayerGeometry{FT}) where {FT} = NoOxygen{FT}()
@inline oxygen_source(oxygen::OxygenTransfer{FT}, geometry::LayerGeometry{FT}) where {FT} =
    OxygenTransfer(leaf_area_source(oxygen.value, geometry))
@inline oxygen_tendency(
    ::NoOxygen{FT},
    ::NoOxygen{FT},
    ::NoOxygen{FT},
    ::LayerGeometry{FT},
) where {FT} = NoOxygen{FT}()
@inline oxygen_tendency(
    bottom::OxygenTransfer{FT},
    top::OxygenTransfer{FT},
    source::OxygenTransfer{FT},
    geometry::LayerGeometry{FT},
) where {FT} =
    OxygenTransfer(conserved_tendency(bottom.value, top.value, source.value, geometry))

"""
    sensible_water_flux(sensible_heat, water_flux, carbon_flux, oxygen, vapor_enthalpy)

Construct a `LeafFlux` from sensible heat [W m⁻² leaf], signed water flux
[kg H₂O m⁻² leaf s⁻¹], signed CO₂ flux [mol C m⁻² leaf s⁻¹], optional oxygen,
and caller-supplied transported vapour enthalpy [J kg⁻¹]. Total energy transfer
is `sensible_heat + water_flux * vapor_enthalpy`.

Positive flux is into air: assimilation must be supplied as negative carbon.
Condensation is negative water. The caller chooses the transfer temperature and
enthalpy reference consistently for both reservoirs, including reverse transfer.
This helper neglects trace-gas enthalpy, mass-carried kinetic/potential and carbon chemical energy;
hosts needing those terms construct `LeafFlux(ExchangeRates(...))` directly.
Do not add a second latent-heat source to the air energy equation.
"""
@inline function sensible_water_flux(
    sensible_heat::FT,
    water_flux::FT,
    carbon_flux::FT,
    oxygen::OxygenContent{FT},
    vapor_enthalpy::FT,
) where {FT<:AbstractFloat}
    return LeafFlux(
        ExchangeRates(
            sensible_heat + water_flux * vapor_enthalpy,
            water_flux,
            carbon_flux,
            oxygen,
        ),
    )
end

"""
    leaf_air_sources(flux::LeafFlux, geometry::LayerGeometry)

Return `(air, leaf)` paired `VolumeSource`s by multiplying a per-leaf flux by
leaf-area density exactly once. Leaf is an exchange ledger for energy, water,
carbon and optional O₂; it does not assert independent leaf CO₂/O₂ storage.
Its values are not leaf-temperature tendencies. Biochemical stoichiometry,
total-mass closure, thermodynamics and liquid replenishment belong to the caller.

For zero leaf-area density, return exact finite zeros for computed components,
even if unused dummy leaf fluxes are nonfinite. Optional oxygen remains absent.
For nonempty layers, nonfinite input propagates and `valid_exchange` detects it.
Soil exchange uses a `FaceFlux`; it cannot be passed to this function.
"""
@inline function leaf_air_sources(
    flux::LeafFlux{FT},
    geometry::LayerGeometry{FT},
) where {FT}
    rates = ExchangeRates(
        leaf_area_source(flux.rates.energy, geometry),
        leaf_area_source(flux.rates.water, geometry),
        leaf_area_source(flux.rates.carbon, geometry),
        oxygen_source(flux.rates.oxygen, geometry),
    )
    return (air = VolumeSource(rates), leaf = VolumeSource(scale_rates(rates, -one(FT))))
end

"""
    layer_tendency(bottom::FaceFlux, top::FaceFlux, source::VolumeSource, geometry)

Conservative finite-volume tendency `(bottom - top) / thickness + source`.
All faces are upward-positive, including the soil lower boundary and canopy top.
Adjacent layers must share the **same face flux** and horizontal ground area.
Host adapters own metric and air-volume-fraction conversion. State entries are densities
per bulk volume; converting to temperature or mixing ratio is host-owned.
Summing thickness-weighted tendencies cancels interior faces. No positivity
clipping or timestep choice is hidden here. Face, source and geometry precision,
and oxygen availability, must agree.
"""
@inline function layer_tendency(
    bottom::FaceFlux{FT,O},
    top::FaceFlux{FT,O},
    source::VolumeSource{FT,O},
    geometry::LayerGeometry{FT},
) where {FT,O}
    return VolumeSource(
        ExchangeRates(
            conserved_tendency(
                bottom.rates.energy,
                top.rates.energy,
                source.rates.energy,
                geometry,
            ),
            conserved_tendency(
                bottom.rates.water,
                top.rates.water,
                source.rates.water,
                geometry,
            ),
            conserved_tendency(
                bottom.rates.carbon,
                top.rates.carbon,
                source.rates.carbon,
                geometry,
            ),
            oxygen_tendency(
                bottom.rates.oxygen,
                top.rates.oxygen,
                source.rates.oxygen,
                geometry,
            ),
        ),
    )
end
