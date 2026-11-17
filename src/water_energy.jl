# ============================================================================
# Energy carried by water: vapour enthalpy of element fluxes, the SurfaceFluxes
# static-energy boundary, and the energy datum of every water-carrying stream.
# These are the only operations that treat energy and water specially.
# ============================================================================

"""
    sensible_water_flux(sensible_heat, water_flux, vapor_enthalpy, tracers::NamedTuple = (;)) -> LeafFlux
    sensible_water_flux(measure::Measure, sensible_heat, water_flux, vapor_enthalpy, tracers::NamedTuple = (;)) -> Exchange{typeof(measure)}

Element flux from sensible heat [W m⁻²], signed water flux [kg H₂O m⁻² s⁻¹], the
specific enthalpy of the transported water vapour [J kg⁻¹] and passive tracers
[mol m⁻² s⁻¹], all per m² of the element's own area. Without a measure the flux is a
[`LeafFlux`](@ref) (per m² one-sided leaf area). The measure is passed as an
instance (`StemArea()`, `GroundArea()` or a measure defined outside the package):
an instance is `isbits` and works in broadcasts over GPU fields, a type does not.

The energy component is `sensible_heat + water_flux * vapor_enthalpy`: total energy
including the enthalpy of the transported vapour. Do not add a separate latent-heat
source when coupling the result to a conserved air energy budget. Evaluate
`vapor_enthalpy` at the temperature of the donor surface (the leaf for transpiration;
the air for condensation onto the surface), in the datum of the host energy budget,
or shift the stream afterwards with [`shift_water_energy_datum`](@ref). Kinetic
energy carried by the transferred mass is excluded. Positive into the canopy air
(photosynthetic CO₂ uptake is a negative `co2`, condensation a negative
`water_flux`).

# Examples
```jldoctest
julia> f = sensible_water_flux(80.0, 4.0e-5, 2.5e6, (; co2 = -1.2e-5));

julia> (f.rates.energy, f.rates.tracers.co2)
(180.0, -1.2e-5)

julia> sensible_water_flux(StemArea(), 20.0, 1.0e-6, 2.5e6) isa StemFlux
true
```
"""
@inline function sensible_water_flux(
    measure::Measure,
    sensible_heat::Real,
    water::Real,
    vapor_enthalpy::Real,
    tracers::NamedTuple = (;),
)
    same_precision(sensible_heat, water, vapor_enthalpy)
    rates = ExchangeRates(sensible_heat + water * vapor_enthalpy, water, tracers)
    return Exchange{typeof(measure)}(rates)
end
@inline sensible_water_flux(
    sensible_heat::Real,
    water::Real,
    vapor_enthalpy::Real,
    tracers::NamedTuple = (;),
) = sensible_water_flux(LeafArea(), sensible_heat, water, vapor_enthalpy, tracers)

"""
    sensible_stem_flux(sensible_heat, tracers::NamedTuple = (;)) -> StemFlux
    sensible_stem_flux(sensible_heat, water_flux, vapor_enthalpy, tracers::NamedTuple = (;)) -> StemFlux

Stem flux per m² one-sided stem area and second, positive from stem into air.

The first form is a dry woody surface: sensible heat [W m⁻² stem] and zero water. Its
tracers are whatever the caller passes. To combine with a leaf flux that carries
`co2`, the stem needs a `co2` key too; passing `(; co2 = 0.0)` asserts that the stem
exchanges no CO₂ (no stem respiration), it does not mean "not computed". The second
form is `sensible_water_flux(StemArea(), ...)`: stem evaporation
[kg H₂O m⁻² stem s⁻¹] with mandatory `vapor_enthalpy` [J kg⁻¹], so the transported
enthalpy `water_flux * vapor_enthalpy` is never silently omitted.

# Examples
```jldoctest
julia> s = sensible_stem_flux(35.0, (; co2 = 0.0));

julia> (s.rates.energy, s.rates.water, s.rates.tracers.co2)
(35.0, 0.0, 0.0)
```
"""
@inline sensible_stem_flux(sensible_heat::Real, tracers::NamedTuple = (;)) =
    StemFlux(ExchangeRates(sensible_heat, zero(sensible_heat), tracers))
@inline sensible_stem_flux(
    sensible_heat::Real,
    water::Real,
    vapor_enthalpy::Real,
    tracers::NamedTuple = (;),
) = sensible_water_flux(StemArea(), sensible_heat, water, vapor_enthalpy, tracers)

"""
    surface_static_energy_flux(shf, lhf, evaporation, tracers::NamedTuple = (;)) -> FaceFlux

Upward face flux from SurfaceFluxes.jl static-energy outputs: sensible heat flux `shf`
and latent heat flux `lhf` [W m⁻²], evaporation mass flux `evaporation`
[kg H₂O m⁻² s⁻¹] and tracers [mol m⁻² s⁻¹], per m² ground area.

For fluxes **computed** by SurfaceFluxes.jl ≥ 1.0 (not prescribed ones),
`SHF = -ρ_sfc g_h ΔDSE + VSE_sfc E` and `LHF = LH_v0 E`, with
`LH_v0 + VSE_sfc = Thermodynamics.enthalpy_vapor(T_sfc) + Φ_sfc`. Hence `shf + lhf` is
the moist static enthalpy plus surface-geopotential flux
`-ρ_sfc g_h ΔDSE + E (h_v(T_sfc) + Φ_sfc)` relative to the Thermodynamics.jl datum,
and it becomes the energy component. Conventional eddy-covariance
`c_p ΔT + L_v(T_air) E` does not satisfy this datum without an explicit shift
([`shift_water_energy_datum`](@ref)). Tracer passing follows the rules of
[`ExchangeRates`](@ref): a zero value asserts a zero transfer.

# Examples
```jldoctest
julia> f = surface_static_energy_flux(60.0, 140.0, 5.6e-5, (; co2 = -1.5e-5));

julia> (f.rates.energy, f.rates.water, f.rates.tracers.co2)
(200.0, 5.6e-5, -1.5e-5)
```
"""
@inline function surface_static_energy_flux(
    shf::Real,
    lhf::Real,
    evaporation::Real,
    tracers::NamedTuple = (;),
)
    same_precision(shf, lhf, evaporation)
    return FaceFlux(ExchangeRates(shf + lhf, evaporation, tracers))
end

"""
    shift_water_energy_datum(stream, delta_specific_energy::Real)

Add `delta_specific_energy * water` to the energy component of `stream` (an
[`ExchangeRates`](@ref) or any [`Exchange`](@ref) record: element fluxes, face fluxes
and volume sources alike); `delta_specific_energy` is in J kg⁻¹ H₂O, for example a
reference-enthalpy offset `h_v,target - h_v,source` or a geopotential offset
`Φ_target - Φ_source`.

Apply the shift to each water-carrying stream (every face flux and every element or
volume source) **before** [`layer_tendency`](@ref): interior faces then stay identical
on both sides and telescope, and the column energy budget changes by exactly the
shift times the water budget. Shifting a layer's net divergence afterwards drops the
vertical gradient of the datum across the cell. Positional arguments only, so it
broadcasts over host fields and GPU arrays.

# Examples
```jldoctest
julia> shift_water_energy_datum(face_flux(100.0, 2.0e-4), 5.0e4).rates.energy
110.0
```
"""
@inline function shift_water_energy_datum(r::ExchangeRates, delta::Real)
    same_precision(r, delta)
    return ExchangeRates(r.energy + delta * r.water, r.water, r.tracers)
end
@inline shift_water_energy_datum(x::Exchange{M}, delta::Real) where {M} =
    Exchange{M}(shift_water_energy_datum(x.rates, delta))
