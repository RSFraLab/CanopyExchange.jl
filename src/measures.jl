# ============================================================================
# Measures: what the rates are "per"
# ============================================================================

"""
    Measure

Supertype of the spatial measures that exchange rates are expressed per:
[`LeafArea`](@ref), [`StemArea`](@ref), [`GroundArea`](@ref), [`BulkVolume`](@ref).

A new kind of element surface (fruit, flowers) is a new subtype; a geometry that
provides [`CanopyExchange.area_density`](@ref) for it converts its fluxes to volume
sources with [`volume_source`](@ref), without changes to this package.

Several populations of the same element (sunlit and shaded leaves, cohorts) share
`LeafArea`; they are not new subtypes. Their fluxes are per m² of their own leaf area,
so the layer source is the sum over classes `c` of
`volume_source(F_c, LayerGeometry(Δz, LAD_c))`. Adding the `LeafFlux`es of different
populations without their area weights is wrong.
"""
abstract type Measure end
# A measure instance is a scalar in broadcasts (e.g. `sensible_water_flux.(StemArea(), ...)`).
Base.Broadcast.broadcastable(m::Measure) = Ref(m)

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
    function Exchange{M,R}(rates::ExchangeRates) where {M<:Measure,R<:ExchangeRates}
        isconcretetype(M) || throw(ArgumentError("the measure must be a concrete Measure type"))
        return new{M,R}(_as(R, rates))
    end
end
@inline _as(::Type{R}, rates::R) where {R} = rates
@inline _as(::Type{R}, rates) where {R} = convert(R, rates)
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
# Direct methods avoid a recursive call cycle through the generic record methods.
@inline _record_equal(compare, a::Exchange, b::Exchange) =
    _measure(a) === _measure(b) && _record_equal(compare, a.rates, b.rates)
Base.hash(a::Exchange, h::UInt) = hash(a.rates, hash(nameof(_measure(a)), hash(:Exchange, h)))

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
`energy` [W m⁻² leaf] is the total transfer including the enthalpy of transported
vapour, water [kg H₂O m⁻² leaf s⁻¹], tracers [mol m⁻² leaf s⁻¹]; positive from leaf
into air. To build `energy` from sensible heat use [`sensible_water_flux`](@ref).
"""
@inline leaf_flux(energy, water, tracers::NamedTuple = (;)) =
    LeafFlux(ExchangeRates(energy, water, tracers))

"""
    stem_flux(energy, water, tracers::NamedTuple = (;)) -> StemFlux

Function form of `StemFlux(ExchangeRates(energy, water, tracers))` for broadcasts:
per m² one-sided stem area and second, positive from stem into air. `energy` is the
total transfer including the enthalpy of transported vapour; to build it from
sensible heat use [`sensible_stem_flux`](@ref) or [`sensible_water_flux`](@ref).
"""
@inline stem_flux(energy, water, tracers::NamedTuple = (;)) =
    StemFlux(ExchangeRates(energy, water, tracers))

"""
    face_flux(energy, water, tracers::NamedTuple = (;)) -> FaceFlux

Function form of `FaceFlux(ExchangeRates(energy, water, tracers))` for broadcasts:
per m² ground area and second, positive upward. `energy` is the total transfer
including the enthalpy of transported vapour; to build it from sensible heat use
[`sensible_water_flux`](@ref) with `GroundArea()`, or
[`surface_static_energy_flux`](@ref) for SurfaceFluxes.jl output.
"""
@inline face_flux(energy, water, tracers::NamedTuple = (;)) =
    FaceFlux(ExchangeRates(energy, water, tracers))
