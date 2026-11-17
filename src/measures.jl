# ============================================================================
# Measures: what the rates are "per"
# ============================================================================

"""
    Measure

Supertype of the spatial measures that exchange rates are expressed per:
[`LeafArea`](@ref), [`StemArea`](@ref), [`GroundArea`](@ref), [`BulkVolume`](@ref).
A new element class (fruit, flowers, a second leaf population) is a new subtype; a
geometry that provides [`area_density`](@ref) for it converts its fluxes to volume
sources with [`volume_source`](@ref), without changes to this package.
"""
abstract type Measure end

"""
    LeafArea <: Measure

Per m² one-sided leaf area. A [`LeafFlux`](@ref) is positive from leaf into air.
"""
struct LeafArea <: Measure end

"""
    StemArea <: Measure

Per m² one-sided woody (stem and branch) area. A [`StemFlux`](@ref) is positive from
stem into air.
"""
struct StemArea <: Measure end

"""
    GroundArea <: Measure

Per m² horizontal ground area. A [`FaceFlux`](@ref) is positive upward (+z).
"""
struct GroundArea <: Measure end

"""
    BulkVolume <: Measure

Per m³ bulk canopy volume (air, leaves and stems together; ground footprint times
layer thickness, with no porosity factor). A [`VolumeSource`](@ref) is positive into
the named reservoir.
"""
struct BulkVolume <: Measure end

"""
    Exchange{M<:Measure}(rates::ExchangeRates)

[`ExchangeRates`](@ref) tagged with the spatial measure `M` they are expressed per.
Arithmetic (`+`, `-`, unary `-`, scalar `*` and `/`, `zero`) exists only within one
measure, so leaf, stem, face and volume quantities cannot be added or swapped by
accident (`MethodError`). Equality and hashing include the measure.

# Examples
```jldoctest
julia> f = Exchange{LeafArea}(ExchangeRates(100.0, 3.0e-5, (; co2 = -1.2e-5)));

julia> f isa LeafFlux
true

julia> (2.0 * f).rates.energy
200.0
```
"""
struct Exchange{M<:Measure,R<:ExchangeRates} <: ScalarRecord
    rates::R
end
Exchange{M}(rates::R) where {M<:Measure,R<:ExchangeRates} = Exchange{M,R}(rates)

"""
    LeafFlux(rates::ExchangeRates)

`Exchange{LeafArea}`: rates per m² one-sided leaf area and second, positive from leaf
into air. Usually built with [`sensible_water_flux`](@ref) or [`leaf_flux`](@ref).

See also [`leaf_air_sources`](@ref), [`StemFlux`](@ref).
"""
const LeafFlux = Exchange{LeafArea}

"""
    StemFlux(rates::ExchangeRates)

`Exchange{StemArea}`: rates per m² one-sided stem area and second, positive from stem
into air. Keeping it distinct from [`LeafFlux`](@ref) prevents swapping leaf and stem
arguments in [`plant_air_sources`](@ref) or adding fluxes over unequal area
denominators. Usually built with [`sensible_stem_flux`](@ref) or [`stem_flux`](@ref).
"""
const StemFlux = Exchange{StemArea}

"""
    FaceFlux(rates::ExchangeRates)

`Exchange{GroundArea}`: rates per m² horizontal ground area and second across a layer
face, positive upward. Usually built with [`surface_static_energy_flux`](@ref),
[`face_flux`](@ref) or [`steady_top_flux`](@ref).
"""
const FaceFlux = Exchange{GroundArea}

"""
    VolumeSource(rates::ExchangeRates)

`Exchange{BulkVolume}`: sources per m³ bulk canopy volume and second, positive into
the named reservoir (canopy air for `air`, the plant element for a donor ledger).
Produced by [`volume_source`](@ref), [`leaf_air_sources`](@ref),
[`plant_air_sources`](@ref) and [`layer_tendency`](@ref).
"""
const VolumeSource = Exchange{BulkVolume}

Base.eltype(::Type{<:Exchange{<:Any,R}}) where {R} = eltype(R)
_measure(::Exchange{M}) where {M} = M
# Direct method avoids a recursive call cycle through the generic record `==`.
Base.:(==)(a::Exchange, b::Exchange) = _measure(a) === _measure(b) && a.rates == b.rates

for op in (:+, :-)
    @eval @inline Base.$op(a::Exchange{M}, b::Exchange{M}) where {M} =
        Exchange{M}($op(a.rates, b.rates))
end
@inline Base.:-(a::Exchange{M}) where {M} = Exchange{M}(-a.rates)
@inline Base.:*(a::Exchange{M}, s::Real) where {M} = Exchange{M}(a.rates * s)
@inline Base.:*(s::Real, a::Exchange{M}) where {M} = Exchange{M}(s * a.rates)
@inline Base.:/(a::Exchange{M}, s::Real) where {M} = Exchange{M}(a.rates / s)
@inline Base.zero(a::Exchange{M}) where {M} = Exchange{M}(zero(a.rates))
Base.zero(::Type{Exchange{M,R}}) where {M,R} = Exchange{M}(zero(R))
@inline valid_exchange(a::Exchange) = valid_exchange(a.rates)

function Base.promote_rule(::Type{Exchange{M,RA}}, ::Type{Exchange{M,RB}}) where {M,RA,RB}
    R = promote_rule(RA, RB)
    return R === Union{} ? Union{} : Exchange{M,R}
end
Base.convert(::Type{Exchange{M,R}}, x::Exchange{M}) where {M,R} =
    Exchange{M}(convert(R, x.rates))

# Function forms of the constructors. Host-field broadcasts on GPU reject a type in
# the function position, so construction by broadcast goes through these.

"""
    leaf_flux(energy, water, tracers::NamedTuple = (;)) -> LeafFlux

Function form of `LeafFlux(ExchangeRates(energy, water, tracers))` for broadcasts:
energy [W m⁻² leaf] including transported vapour enthalpy, water
[kg H₂O m⁻² leaf s⁻¹], tracers [mol m⁻² leaf s⁻¹]; positive from leaf into air.
"""
@inline leaf_flux(energy, water, tracers::NamedTuple = (;)) =
    LeafFlux(ExchangeRates(energy, water, tracers))

"""
    stem_flux(energy, water, tracers::NamedTuple = (;)) -> StemFlux

Function form of `StemFlux(ExchangeRates(energy, water, tracers))` for broadcasts:
per m² one-sided stem area and second, positive from stem into air.
"""
@inline stem_flux(energy, water, tracers::NamedTuple = (;)) =
    StemFlux(ExchangeRates(energy, water, tracers))

"""
    face_flux(energy, water, tracers::NamedTuple = (;)) -> FaceFlux

Function form of `FaceFlux(ExchangeRates(energy, water, tracers))` for broadcasts:
per m² ground area and second, positive upward.
"""
@inline face_flux(energy, water, tracers::NamedTuple = (;)) =
    FaceFlux(ExchangeRates(energy, water, tracers))
