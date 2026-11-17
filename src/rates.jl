# ============================================================================
# Rate records: energy and water are thermodynamically coupled; tracers are passive
# ============================================================================

# Internal supertype of the two rate records (ExchangeRates, MoistAirTendencies):
# componentwise algebra is defined once for it.
abstract type RateRecord <: ScalarRecord end

"""
    ExchangeRates(energy, water, tracers::NamedTuple = (;))
    ExchangeRates(; energy, water, tracers = (;))

Conserved transfer rates, each per second and per measure of the enclosing
[`Exchange`](@ref) (one-sided leaf or stem area, ground area or bulk volume):

- `energy` [J s⁻¹ per measure]: total energy transfer under one caller-declared
  thermodynamic datum, not sensible or latent heat in isolation. It includes the
  enthalpy carried by transported water vapour (see [`sensible_water_flux`](@ref));
  do not add a second latent-heat term when coupling to a conserved air energy budget.
- `water` [kg H₂O s⁻¹ per measure]: total water mass transfer.
- `tracers` [mol s⁻¹ per measure]: passive molar tracers, named by the keys of a
  `NamedTuple`, e.g. `(; co2 = -1.2e-5)` for CO₂ in mol CO₂ (one mol C per mol CO₂;
  not the biochemical enthalpy of assimilated carbon).

Sign: positive in the direction of the enclosing measure (element into air, upward
across a face, into the named reservoir for a volume source). Host adapters own the
conversions between conserved density tendencies and temperature or humidity.

Tracer rules:

- Tracers are molar and passive. They carry no energy and no mass in these records.
  Leaving their mass out of the total-mass source is the trace approximation: with
  12 µmol CO₂ m⁻² s⁻¹ against 4e-5 kg H₂O m⁻² s⁻¹ of transpiration, `M_CO₂ × S_CO₂` is
  about 1.3 % of the water mass source (about 0.3 % net when the opposite O₂ flux is
  counted), and it dominates when transpiration is near zero (at night).
- Water isotopologues (HDO, H₂¹⁸O) are subsets of `water`. Carried as tracers, they
  never enter total mass, which is `water` alone.
- Bulk O₂ is not a trace constituent. An `o2` key is a transfer ledger; it becomes a
  prognostic oxygen budget only if the host has a separate oxygen store contract.
- A tracer that is not computed is not a key: absence is part of the type, never a
  stored zero. Passing a key with value zero asserts a zero transfer (for example no
  stem respiration); it does not mean "not computed". Records with different tracer
  names cannot be combined (`ArgumentError`).
- Tracer order is part of the type: `(; co2, o2)` and `(; o2, co2)` are different
  record types and cannot be combined either.
- The records assume a zero dry-air mass source: the total-mass source is the water
  source (see [`moist_air_tendencies`](@ref)).

All scalars, tracer values included, must be real and carry one floating-point type.
`Float32` with `Float64` throws an `ArgumentError`, integers are rejected, and
`Float64` with a `Float64`-valued dual number promotes to the dual (with the
ForwardDiff extension loaded).

Records support componentwise `+`, `-`, unary `-`, scalar `*` and `/`, `zero`, value
equality `==` (`-0.0 == 0.0`, `NaN != NaN`; derivative components compared too, see
[`CanopyExchange.scalar_equal`](@ref)), `isequal` and a consistent `hash`. In
broadcasts a record is a scalar. Inference and zero allocations are tested up to 12
tracers; beyond 31 tracers Julia's tuple-inference limit makes the kernels type
unstable.

# Examples
```jldoctest
julia> r = ExchangeRates(120.0, 2.0e-5, (; co2 = -1.0e-5));

julia> (r.energy, r.water, r.tracers.co2)
(120.0, 2.0e-5, -1.0e-5)

julia> keys((r + r).tracers)
(:co2,)
```

See also [`Tracers`](@ref), [`Exchange`](@ref), [`valid_exchange`](@ref).
"""
struct ExchangeRates{FT<:Real,T<:NamedTuple} <: RateRecord
    energy::FT
    water::FT
    tracers::T
    global _unchecked_rates(energy::FT, water::FT, tracers::T) where {FT,T} =
        new{FT,T}(energy, water, tracers)
end
@inline function ExchangeRates(energy::Real, water::Real, tracers::NamedTuple)
    (e, w), t = _promote_components((energy, water), tracers)
    return _unchecked_rates(e, w, t)
end
@inline ExchangeRates(energy::Real, water::Real) = ExchangeRates(energy, water, (;))
ExchangeRates(; energy, water, tracers = (;)) = ExchangeRates(energy, water, tracers)

"""
    Tracers(names::Symbol...)
    Tracers(:co2, :o2)(co2_value, o2_value) -> (; co2 = co2_value, o2 = o2_value)

Name tracer values positionally. `Tracers(:co2, :o2)` is a singleton callable, so it
can be used inside broadcasts over GPU arrays and host fields, where a keyword
argument or a `NamedTuple` type in the function position cannot. Values are in
mol s⁻¹ per measure of the record that receives them (see [`ExchangeRates`](@ref)).
Names must be unique.

# Examples
```jldoctest
julia> Tracers(:co2, :o2)(-1.2e-5, 1.3e-5)
(co2 = -1.2e-5, o2 = 1.3e-5)
```
"""
struct Tracers{names} end
function Tracers(names::Symbol...)
    allunique(names) || throw(ArgumentError("tracer names must be unique"))
    return Tracers{names}()
end
@inline (::Tracers{names})(values::Vararg{Real,N}) where {names,N} = NamedTuple{names}(values)

# ============================================================================
# Checked promotion of record components (constructors)
# ============================================================================

# Check the precision policy and promote every component to one scalar type. A tracer
# tuple of one type (the usual case) is checked and converted without recursion.
@inline function _promote_components(scalars::Tuple, tracers::NamedTuple{names}) where {names}
    same_precision(scalars...)
    _same_tracer_precision(_float_type(first(scalars)), values(tracers))
    FT = _promoted_type(scalars, values(tracers))
    return _convert_all(FT, scalars), NamedTuple{names}(_convert_all(FT, values(tracers)))
end
@inline _promoted_type(scalars, ::Tuple{}) = promote_type(map(typeof, scalars)...)
@inline _promoted_type(scalars, ::Tuple{T,Vararg{T}}) where {T} =
    promote_type(map(typeof, scalars)..., T)
@inline _promoted_type(scalars, t::Tuple) =
    promote_type(map(typeof, scalars)..., map(typeof, t)...)
@inline _convert_all(::Type{FT}, t::Tuple{Vararg{FT}}) where {FT} = t
@inline _convert_all(::Type{FT}, t::Tuple) where {FT} = map(Base.Fix1(convert, FT), t)

@inline _tracer_type(::Type{<:NamedTuple{names}}, ::Type{FT}) where {names,FT} =
    NamedTuple{names,NTuple{length(names),FT}}

# ============================================================================
# Componentwise algebra, defined once for both rate records
# ============================================================================

@inline _scalars(r::ExchangeRates) = (r.energy, r.water)
@inline _unchecked(::ExchangeRates, s, tracers) = _unchecked_rates(s..., tracers)
@inline _checked_record(::ExchangeRates, s, tracers) = ExchangeRates(s..., tracers)
@inline _tracer_values(r) = values(r.tracers)

@inline function _same_tracers(a, bs...)
    all(b -> keys(b.tracers) === keys(a.tracers), bs) || throw(ArgumentError(_TRACER_SETS))
    return nothing
end

# Apply `f` componentwise to records of one kind with the same tracers in the same
# order; their scalar types may differ only by derivative information. Results of one
# scalar type are assembled without re-checking.
@inline function _componentwise(f::F, a::RateRecord, bs::Vararg{RateRecord,N}) where {F,N}
    same_precision(a, bs...)
    _same_tracers(a, bs...)
    scalars = map(f, _scalars(a), map(_scalars, bs)...)
    tracers = NamedTuple{keys(a.tracers)}(map(f, _tracer_values(a), map(_tracer_values, bs)...))
    return _assemble(a, scalars, tracers)
end
@inline _assemble(a, scalars::Tuple{FT,Vararg{FT}}, tracers::NamedTuple{K,<:Tuple{Vararg{FT}}}) where {FT,K} =
    _unchecked(a, scalars, tracers)
@inline _assemble(a, scalars, tracers) = _checked_record(a, scalars, tracers)

@inline Base.:+(a::ExchangeRates, b::ExchangeRates) = _componentwise(+, a, b)
@inline Base.:-(a::ExchangeRates, b::ExchangeRates) = _componentwise(-, a, b)
@inline Base.:-(a::RateRecord) = _componentwise(-, a)
@inline Base.:*(a::RateRecord, s::Real) = (same_precision(a, s); _componentwise(Base.Fix2(*, s), a))
@inline Base.:*(s::Real, a::RateRecord) = (same_precision(a, s); _componentwise(Base.Fix1(*, s), a))
@inline Base.:/(a::RateRecord, s::Real) = (same_precision(a, s); _componentwise(Base.Fix2(/, s), a))
@inline Base.zero(a::RateRecord) = _componentwise(zero, a)

Base.eltype(::Type{<:ExchangeRates{FT}}) where {FT} = FT
@inline _zero_tracers(::Type{NamedTuple{names,V}}) where {names,V} =
    NamedTuple{names}(map(zero, fieldtypes(V)))
Base.zero(::Type{ExchangeRates{FT,T}}) where {FT,T} =
    ExchangeRates(zero(FT), zero(FT), _zero_tracers(T))

# Container promotion: [record, dual_record] gets a concrete element type. A
# conversion returns exactly the requested type or throws.
@inline function _record_promotion(R, ::Type{A}, ::Type{TA}, ::Type{B}, ::Type{TB}) where {A,TA,B,TB}
    (_promotable(A, B) && fieldnames(TA) === fieldnames(TB)) || return Union{}
    FT = promote_type(A, B)
    return R{FT,_tracer_type(TA, FT)}
end
function _record_convert(::Type{FT}, ::Type{T}, x::RateRecord) where {FT,T}
    same_precision(zero(FT), x)
    fieldnames(T) === keys(x.tracers) || throw(ArgumentError(_TRACER_SETS))
    T === _tracer_type(T, FT) ||
        throw(ArgumentError("tracer values must carry the record scalar type"))
    return _componentwise(Base.Fix1(convert, FT), x)
end
Base.promote_rule(::Type{ExchangeRates{A,TA}}, ::Type{ExchangeRates{B,TB}}) where {A,TA,B,TB} =
    _record_promotion(ExchangeRates, A, TA, B, TB)
Base.convert(::Type{ExchangeRates{FT,T}}, x::ExchangeRates) where {FT,T} =
    _record_convert(FT, T, x)

"""
    valid_exchange(x) -> Bool

Whether every computed component of `x` is finite, including derivative components
of dual numbers. Accepts [`ExchangeRates`](@ref), any [`Exchange`](@ref) record and
[`MoistAirTendencies`](@ref). Absent tracers are not keys and are not checked. Hosts
remain responsible for thermodynamic admissibility and retry policy.

# Examples
```jldoctest
julia> valid_exchange(leaf_flux(100.0, 2.0e-5, (; co2 = -1.0e-5)))
true

julia> valid_exchange(leaf_flux(NaN, 2.0e-5, (; co2 = -1.0e-5)))
false
```
"""
@inline valid_exchange(a::RateRecord) =
    all(all_finite, _scalars(a)) && all(all_finite, _tracer_values(a))
