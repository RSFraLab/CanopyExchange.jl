# ============================================================================
# Element and boundary fluxes: the only operations that treat energy and water
# specially (vapour enthalpy carried by water, and the energy datum of water)
# ============================================================================

@inline function _sensible_water_rates(sensible_heat, water, vapor_enthalpy, tracers)
    same_precision(sensible_heat, water, vapor_enthalpy, values(tracers)...)
    return ExchangeRates(sensible_heat + water * vapor_enthalpy, water, tracers)
end

"""
    sensible_water_flux(sensible_heat, water_flux, vapor_enthalpy, tracers::NamedTuple = (;)) -> LeafFlux

Leaf flux from sensible heat [W m⁻² leaf], signed water flux [kg H₂O m⁻² leaf s⁻¹],
the specific enthalpy of the transported water vapour [J kg⁻¹] and passive tracers
[mol m⁻² leaf s⁻¹], all per m² one-sided leaf area.

The energy component is `sensible_heat + water_flux * vapor_enthalpy`: total energy
including the enthalpy of the transported vapour. Do not add a separate latent-heat
source when coupling the result to a conserved air energy budget. Positive into the
canopy air (photosynthetic CO₂ uptake is a negative `co2`, condensation a negative
`water_flux`). Pass `vapor_enthalpy` in the datum of the host energy budget, or shift
the stream with [`shift_water_energy_datum`](@ref).

# Examples
```jldoctest
julia> f = sensible_water_flux(80.0, 4.0e-5, 2.5e6, (; co2 = -1.2e-5));

julia> (f.rates.energy, f.rates.tracers.co2)
(180.0, -1.2e-5)
```
"""
@inline sensible_water_flux(sensible_heat, water, vapor_enthalpy, tracers::NamedTuple = (;)) =
    LeafFlux(_sensible_water_rates(sensible_heat, water, vapor_enthalpy, tracers))

"""
    sensible_stem_flux(sensible_heat, tracers::NamedTuple = (;)) -> StemFlux
    sensible_stem_flux(sensible_heat, water_flux, vapor_enthalpy, tracers::NamedTuple = (;)) -> StemFlux

Stem flux per m² one-sided stem area and second, positive from stem into air.

The first form is a dry woody surface: sensible heat [W m⁻² stem] and zero water.
Tracers are whatever the caller passes; a zero stem CO₂ flux must be given explicitly
(`(; co2 = 0.0)`) to combine with a leaf flux that carries `co2`. The second form adds
stem evaporation [kg H₂O m⁻² stem s⁻¹] and makes `vapor_enthalpy` [J kg⁻¹] mandatory,
so the transported enthalpy `water_flux * vapor_enthalpy` is never silently omitted.

# Examples
```jldoctest
julia> s = sensible_stem_flux(35.0, (; co2 = 0.0));

julia> (s.rates.energy, s.rates.water, s.rates.tracers.co2)
(35.0, 0.0, 0.0)
```
"""
@inline sensible_stem_flux(sensible_heat, tracers::NamedTuple = (;)) =
    StemFlux(ExchangeRates(sensible_heat, zero(sensible_heat), tracers))
@inline sensible_stem_flux(sensible_heat, water, vapor_enthalpy, tracers::NamedTuple = (;)) =
    StemFlux(_sensible_water_rates(sensible_heat, water, vapor_enthalpy, tracers))

"""
    surface_static_energy_flux(shf, lhf, evaporation, tracers::NamedTuple = (;)) -> FaceFlux

Upward face flux from SurfaceFluxes.jl static-energy outputs: sensible heat flux `shf`
and latent heat flux `lhf` [W m⁻²], kinematic water mass flux `evaporation`
[kg H₂O m⁻² s⁻¹] and tracers [mol m⁻² s⁻¹], per m² ground area.

In the SurfaceFluxes.jl convention `SHF = -ρ_sfc g_h ΔDSE + VSE_sfc E` and
`LHF = LH_v0 E`, with `LH_v0 + VSE_sfc = Thermodynamics.enthalpy_vapor(T_sfc) + Φ_sfc`.
Hence `shf + lhf` is the moist static enthalpy plus surface-geopotential flux
`-ρ_sfc g_h ΔDSE + E (h_v(T_sfc) + Φ_sfc)` relative to the Thermodynamics.jl datum,
and it becomes the energy component. Conventional eddy-covariance
`c_p ΔT + L_v(T_air) E` does not satisfy this datum without an explicit shift
([`shift_water_energy_datum`](@ref)).

# Examples
```jldoctest
julia> f = surface_static_energy_flux(60.0, 140.0, 5.6e-5, (; co2 = -1.5e-5));

julia> (f.rates.energy, f.rates.water, f.rates.tracers.co2)
(200.0, 5.6e-5, -1.5e-5)
```
"""
@inline function surface_static_energy_flux(shf, lhf, evaporation, tracers::NamedTuple = (;))
    same_precision(shf, lhf, evaporation, values(tracers)...)
    return FaceFlux(ExchangeRates(shf + lhf, evaporation, tracers))
end

"""
    shift_water_energy_datum(stream, delta_specific_energy)

Add `delta_specific_energy * water` to the energy component of `stream` (an
[`ExchangeRates`](@ref) or any [`Exchange`](@ref) record); `delta_specific_energy`
is in J kg⁻¹ H₂O, for example a reference-enthalpy offset `h_v,target - h_v,source` or
a geopotential offset `Φ_target - Φ_source`.

Apply the shift to each water-carrying stream (every face flux and every element or
volume source) **before** [`layer_tendency`](@ref): interior faces then stay identical
on both sides and telescope. Shifting a layer's net divergence afterwards drops the
vertical gradient of the datum across the cell. Positional arguments only, so it
broadcasts over host fields and GPU arrays.

# Examples
```jldoctest
julia> shift_water_energy_datum(face_flux(100.0, 2.0e-4), 5.0e4).rates.energy
110.0
```
"""
@inline function shift_water_energy_datum(r::ExchangeRates, delta)
    same_precision(r, delta)
    return ExchangeRates(r.energy + delta * r.water, r.water, r.tracers)
end
@inline shift_water_energy_datum(x::Exchange{M}, delta) where {M} =
    Exchange{M}(shift_water_energy_datum(x.rates, delta))
