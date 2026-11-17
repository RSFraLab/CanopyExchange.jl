# ============================================================================
# Bulk-volume tendencies to gas-volume moist-air tendencies
# ============================================================================

"""
    moist_air_tendencies(tendency::VolumeSource, gas_fraction = 1) -> MoistAirTendencies

Divide a bulk-volume tendency [per m³ bulk canopy volume and second] by the gas
volume fraction `gas_fraction` [m³ gas m⁻³ bulk] to obtain source tendencies per m³
gas volume ([`MoistAirTendencies`](@ref)).

Contracts:

- `gas_fraction` is constant in time (`∂ₜ ε_gas = 0`) and lies in `(0, 1]`; a
  nonfinite value, a value outside that interval or a nonfinite derivative component
  throws an `ArgumentError`.
- Dry-air mass is conserved (`S_ρd = 0`), so `rho_total == rho_water ==
  water / gas_fraction`; `rho_energy = energy / gas_fraction`, and every tracer is
  divided by `gas_fraction`. Each component is divided (not multiplied by the
  inverse), which keeps exact identity at `gas_fraction == 1` and tracer ratios.
- Apply any reference-enthalpy or geopotential shift to the individual streams
  ([`shift_water_energy_datum`](@ref)) before [`layer_tendency`](@ref), not here.
- Positional arguments only, so it broadcasts over host fields and GPU arrays.

# Examples
```jldoctest
julia> v = VolumeSource(ExchangeRates(180.0, 3.6e-5, (; co2 = -1.8e-5)));

julia> m = moist_air_tendencies(v, 0.9);

julia> (m.rho_total, m.rho_water, m.rho_energy, m.tracers.co2)
(4.0e-5, 4.0e-5, 200.0, -2.0e-5)
```
"""
@inline function moist_air_tendencies(
    tendency::VolumeSource,
    gas_fraction = one(eltype(tendency)),
)
    _checked(gas_fraction, _unit_interval, "gas_fraction must be finite and in (0, 1]")
    r = tendency.rates / gas_fraction
    return MoistAirTendencies(r.water, r.water, r.energy, r.tracers)
end
