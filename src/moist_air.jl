# ============================================================================
# Gas-volume moist-air tendencies for the host
# ============================================================================

"""
    MoistAirTendencies(rho_total, rho_water, rho_energy, tracers::NamedTuple = (;))
    MoistAirTendencies(; rho_total, rho_water, rho_energy, tracers = (;))

Source tendencies per m³ **gas** volume and second for a host that conserves dry-air
mass (`S_ρd = 0`):

- `rho_total` [kg m⁻³ gas s⁻¹]: total moist-air mass source; equals `rho_water`
  (tracer mass is left out: the trace approximation, see [`ExchangeRates`](@ref)).
- `rho_water` [kg H₂O m⁻³ gas s⁻¹]: total water mass source.
- `rho_energy` [J m⁻³ gas s⁻¹]: total energy source in the host's datum.
- `tracers` [mol m⁻³ gas s⁻¹]: passive molar tracer sources, kept nested so a tracer
  name cannot overwrite a host field.

Positive into the canopy air. Momentum, canopy drag and subgrid turbulent kinetic
energy are excluded. Host adapters own the conversion to temperature and humidity
tendencies. Same precision policy and componentwise algebra as
[`ExchangeRates`](@ref).

# Examples
```jldoctest
julia> m = MoistAirTendencies(4.0e-5, 4.0e-5, 150.0, (; co2 = -1.0e-5));

julia> (m.rho_total, m.rho_energy, m.tracers.co2)
(4.0e-5, 150.0, -1.0e-5)
```

See also [`moist_air_tendencies`](@ref).
"""
struct MoistAirTendencies{FT<:Real,T<:NamedTuple} <: RateRecord
    rho_total::FT
    rho_water::FT
    rho_energy::FT
    tracers::T
    global _unchecked_moist(total::FT, water::FT, energy::FT, tracers::T) where {FT,T} =
        new{FT,T}(total, water, energy, tracers)
end
@inline function MoistAirTendencies(total::Real, water::Real, energy::Real, tracers::NamedTuple)
    (m, w, e), t = _promote_components((total, water, energy), tracers)
    return _unchecked_moist(m, w, e, t)
end
@inline MoistAirTendencies(total::Real, water::Real, energy::Real) =
    MoistAirTendencies(total, water, energy, (;))
MoistAirTendencies(; rho_total, rho_water, rho_energy, tracers = (;)) =
    MoistAirTendencies(rho_total, rho_water, rho_energy, tracers)

# Opt into the componentwise algebra of rate records (rates.jl).
@inline _scalars(r::MoistAirTendencies) = (r.rho_total, r.rho_water, r.rho_energy)
@inline _unchecked(::MoistAirTendencies, s, tracers) = _unchecked_moist(s..., tracers)
@inline _checked_record(::MoistAirTendencies, s, tracers) = MoistAirTendencies(s..., tracers)
@inline Base.:+(a::MoistAirTendencies, b::MoistAirTendencies) = _componentwise(+, a, b)
@inline Base.:-(a::MoistAirTendencies, b::MoistAirTendencies) = _componentwise(-, a, b)
Base.eltype(::Type{<:MoistAirTendencies{FT}}) where {FT} = FT
Base.zero(::Type{MoistAirTendencies{FT,T}}) where {FT,T} =
    MoistAirTendencies(zero(FT), zero(FT), zero(FT), _zero_tracers(T))
Base.promote_rule(
    ::Type{MoistAirTendencies{A,TA}},
    ::Type{MoistAirTendencies{B,TB}},
) where {A,TA,B,TB} = _record_promotion(MoistAirTendencies, A, TA, B, TB)
Base.convert(::Type{MoistAirTendencies{FT,T}}, x::MoistAirTendencies) where {FT,T} =
    _record_convert(FT, T, x)

"""
    moist_air_tendencies(tendency::VolumeSource, gas_fraction = one(eltype(tendency))) -> MoistAirTendencies

Divide a bulk-volume tendency [per m³ bulk canopy volume and second] by the gas
volume fraction `gas_fraction` [m³ gas m⁻³ bulk] to obtain source tendencies per m³
gas volume ([`MoistAirTendencies`](@ref)).

Contracts:

- `gas_fraction` is constant in time (`∂ₜ ε_gas = 0`), lies in `(0, 1]` and has the
  floating-point type of `tendency`; a nonfinite value, a value outside that interval
  or a nonfinite derivative component throws an `ArgumentError`.
- Dry-air mass is conserved (`S_ρd = 0`), so `rho_total == rho_water ==
  water / gas_fraction`; `rho_energy = energy / gas_fraction`, and every tracer is
  divided by `gas_fraction`. Leaving tracer mass out of `rho_total` is the trace
  approximation (about 1.3 % of the water mass source for daytime CO₂ uptake, and
  dominant when transpiration is near zero; see [`ExchangeRates`](@ref)).
- Each component is divided (not multiplied by the inverse), which keeps exact
  identity at `gas_fraction == 1`. Ratios between tracers are kept exactly for ±1:1
  stoichiometry (equal magnitudes) and to rounding otherwise.
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
    gas_fraction::Real = one(eltype(tendency)),
)
    _checked(gas_fraction, _unit_interval, "gas_fraction must be finite and in (0, 1]")
    r = tendency.rates / gas_fraction
    return _unchecked_moist(r.water, r.water, r.energy, r.tracers)
end
