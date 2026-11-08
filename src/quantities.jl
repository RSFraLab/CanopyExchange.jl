"""Optional transported oxygen; absence is distinct from a computed zero."""
abstract type OxygenContent{FT} end

"""Oxygen transfer is not computed. `oxygen_value` returns `nothing`."""
struct NoOxygen{FT<:AbstractFloat} <: OxygenContent{FT} end

"""Signed oxygen transfer in mol O₂, with the area/volume/time measure of its container."""
struct OxygenTransfer{FT<:AbstractFloat} <: OxygenContent{FT}
    value::FT
end

@inline oxygen_value(::NoOxygen) = nothing
@inline oxygen_value(oxygen::OxygenTransfer) = oxygen.value
@inline scale_oxygen(::NoOxygen{FT}, factor::FT) where {FT} = NoOxygen{FT}()
@inline scale_oxygen(oxygen::OxygenTransfer{FT}, factor::FT) where {FT} =
    OxygenTransfer(oxygen.value * factor)
@inline valid_oxygen(::NoOxygen) = true
@inline valid_oxygen(oxygen::OxygenTransfer) = isfinite(oxygen.value)

"""
    ExchangeRates(energy, water, carbon, oxygen)

Conserved transfer components: energy [J], water [kg H₂O], carbon [mol C]
and optional oxygen [mol O₂], each per second and the measure of its wrapper.
All numeric fields and optional oxygen use one floating-point type. Carbon
exchange here is CO₂ transfer (one carbon atom per molecule), not a claim about
the chemical energy of assimilated carbon.

`energy` is the complete energy transfer for a **declared caller convention**,
not automatically sensible heat or latent heat. An adapter must match the host
energy variable and include applicable transported-mass energy. These scalar
records do not implicitly convert to temperature or mixing-ratio tendencies.
"""
struct ExchangeRates{FT<:AbstractFloat,O<:OxygenContent{FT}}
    energy::FT
    water::FT
    carbon::FT
    oxygen::O
end

"""Flux per m² one-sided leaf area per second; positive from leaf into air."""
struct LeafFlux{FT,O<:OxygenContent{FT}}
    rates::ExchangeRates{FT,O}
end

"""Flux per m² ground area per second; positive upward across a horizontal face."""
struct FaceFlux{FT,O<:OxygenContent{FT}}
    rates::ExchangeRates{FT,O}
end

"""Source per m³ bulk canopy volume per second; positive into the named reservoir."""
struct VolumeSource{FT,O<:OxygenContent{FT}}
    rates::ExchangeRates{FT,O}
end

"""
    LayerGeometry(thickness, leaf_area_density)

Positive finite layer thickness [m] and nonnegative finite one-sided leaf area
density [m² leaf / m³ bulk canopy]. Construct during host setup. The volume uses
ground footprint times thickness; no separate canopy-air porosity is assumed.
"""
struct LayerGeometry{FT<:AbstractFloat}
    thickness::FT
    leaf_area_density::FT
    function LayerGeometry(thickness::FT, leaf_area_density::FT) where {FT<:AbstractFloat}
        isfinite(thickness) && thickness > zero(FT) ||
            throw(ArgumentError("layer thickness must be finite and positive"))
        isfinite(leaf_area_density) && leaf_area_density >= zero(FT) ||
            throw(ArgumentError("leaf area density must be finite and nonnegative"))
        return new{FT}(thickness, leaf_area_density)
    end
end

@inline function scale_rates(rates::ExchangeRates{FT}, factor::FT) where {FT}
    return ExchangeRates(
        rates.energy * factor,
        rates.water * factor,
        rates.carbon * factor,
        scale_oxygen(rates.oxygen, factor),
    )
end

"""Whether all computed components are finite. Hosts decide retry/rejection policy."""
@inline valid_exchange(rates::ExchangeRates) =
    isfinite(rates.energy) &&
    isfinite(rates.water) &&
    isfinite(rates.carbon) &&
    valid_oxygen(rates.oxygen)
@inline valid_exchange(flux::Union{LeafFlux,FaceFlux,VolumeSource}) =
    valid_exchange(flux.rates)
