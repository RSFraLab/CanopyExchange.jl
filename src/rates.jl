# ============================================================================
# Rate records: energy and water are thermodynamically coupled; tracers are passive
# ============================================================================

# Internal supertype of every record this package defines. Value equality, hashing
# and broadcasting as a scalar are defined once for it (end of this file).
abstract type ScalarRecord end

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
across a face, into the named reservoir for a volume source).

Tracer rules:

- Tracers are molar and passive. They carry no energy and no mass in these records.
- Water isotopologues (HDO, H₂¹⁸O) are subsets of `water`. Carried as tracers, they
  never enter total mass, which is `water` alone.
- An `o2` key carries the O₂ transfer; it does not create an oxygen inventory.
- A tracer that is not computed is not a key: absence is part of the type, never a
  stored zero. Records with different tracer names cannot be combined
  (`ArgumentError`).
- Tracer order is part of the type: `(; co2, o2)` and `(; o2, co2)` are different
  record types and cannot be combined either.
- The records assume a zero dry-air mass source: the total-mass source is the water
  source (see [`moist_air_tendencies`](@ref)).

All scalars must carry one floating-point type. `Float32` with `Float64` throws an
`ArgumentError`, integers are rejected, and `Float64` with a `Float64`-valued dual
number promotes to the dual (with the ForwardDiff extension loaded).

Records support componentwise `+`, `-`, unary `-`, scalar `*` and `/`, `zero`,
`map(f, rates...)`, value equality `==` (`-0.0 == 0.0`, `NaN != NaN`), `isequal`
and a consistent `hash`. In broadcasts a record is a scalar.

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
struct ExchangeRates{FT<:Real,T<:NamedTuple} <: ScalarRecord
    energy::FT
    water::FT
    tracers::T
    global _unchecked_rates(energy::FT, water::FT, tracers::T) where {FT,T} =
        new{FT,T}(energy, water, tracers)
end
@inline function ExchangeRates(energy::Real, water::Real, tracers::NamedTuple{names}) where {names}
    same_precision(energy, water, values(tracers)...)
    promoted = promote(energy, water, values(tracers)...)
    return _unchecked_rates(promoted[1], promoted[2], NamedTuple{names}(Base.tail(Base.tail(promoted))))
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

# Examples
```jldoctest
julia> Tracers(:co2, :o2)(-1.2e-5, 1.3e-5)
(co2 = -1.2e-5, o2 = 1.3e-5)
```
"""
struct Tracers{names} end
Tracers(names::Symbol...) = Tracers{names}()
@inline (::Tracers{names})(values::Vararg{Real,N}) where {names,N} = NamedTuple{names}(values)

"""
    MoistAirTendencies(rho_total, rho_water, rho_energy, tracers::NamedTuple = (;))
    MoistAirTendencies(; rho_total, rho_water, rho_energy, tracers = (;))

Source tendencies per m³ **gas** volume and second for a host that conserves dry-air
mass (`S_ρd = 0`):

- `rho_total` [kg m⁻³ gas s⁻¹]: total moist-air mass source; equals `rho_water`.
- `rho_water` [kg H₂O m⁻³ gas s⁻¹]: total water mass source.
- `rho_energy` [J m⁻³ gas s⁻¹]: total energy source in the host's datum.
- `tracers` [mol m⁻³ gas s⁻¹]: passive molar tracer sources, kept nested so a tracer
  name cannot overwrite a host field.

Positive into the canopy air. Momentum, canopy drag and subgrid turbulent kinetic
energy are excluded. Same precision policy and componentwise algebra as
[`ExchangeRates`](@ref).

# Examples
```jldoctest
julia> m = MoistAirTendencies(4.0e-5, 4.0e-5, 150.0, (; co2 = -1.0e-5));

julia> (m.rho_total, m.rho_energy, m.tracers.co2)
(4.0e-5, 150.0, -1.0e-5)
```

See also [`moist_air_tendencies`](@ref).
"""
struct MoistAirTendencies{FT<:Real,T<:NamedTuple} <: ScalarRecord
    rho_total::FT
    rho_water::FT
    rho_energy::FT
    tracers::T
    global _unchecked_moist(total::FT, water::FT, energy::FT, tracers::T) where {FT,T} =
        new{FT,T}(total, water, energy, tracers)
end
@inline function MoistAirTendencies(
    total::Real,
    water::Real,
    energy::Real,
    tracers::NamedTuple{names},
) where {names}
    same_precision(total, water, energy, values(tracers)...)
    promoted = promote(total, water, energy, values(tracers)...)
    return _unchecked_moist(
        promoted[1],
        promoted[2],
        promoted[3],
        NamedTuple{names}(Base.tail(Base.tail(Base.tail(promoted)))),
    )
end
@inline MoistAirTendencies(total::Real, water::Real, energy::Real) =
    MoistAirTendencies(total, water, energy, (;))
MoistAirTendencies(; rho_total, rho_water, rho_energy, tracers = (;)) =
    MoistAirTendencies(rho_total, rho_water, rho_energy, tracers)

# ============================================================================
# Componentwise algebra, defined once for both rate records
# ============================================================================

const Rates = Union{ExchangeRates,MoistAirTendencies}
@inline _scalars(r::ExchangeRates) = (r.energy, r.water)
@inline _scalars(r::MoistAirTendencies) = (r.rho_total, r.rho_water, r.rho_energy)
@inline _rebuild(::ExchangeRates, s, tracers) = ExchangeRates(s..., tracers)
@inline _rebuild(::MoistAirTendencies, s, tracers) = MoistAirTendencies(s..., tracers)
@inline _tracers(r) = r.tracers

Base.eltype(::Type{<:ExchangeRates{FT}}) where {FT} = FT
Base.eltype(::Type{<:MoistAirTendencies{FT}}) where {FT} = FT

@inline _tracer_type(::Type{<:NamedTuple{names}}, ::Type{FT}) where {names,FT} =
    NamedTuple{names,NTuple{length(names),FT}}

# Records must be of one kind and hold the same tracers in the same order; their
# scalar types may differ only by derivative information.
for R in (:ExchangeRates, :MoistAirTendencies)
    @eval begin
        @inline function Base.map(f::F, a::$R, bs::Vararg{$R,N}) where {F,N}
            same_precision(a, bs...)
            scalars = map(f, _scalars(a), map(_scalars, bs)...)
            tracers = map(f, a.tracers, map(_tracers, bs)...)
            return _rebuild(a, scalars, tracers)
        end
        @inline Base.:+(a::$R, b::$R) = map(+, a, b)
        @inline Base.:-(a::$R, b::$R) = map(-, a, b)

        # Container promotion: [record, dual_record] gets a concrete element type.
        function Base.promote_rule(::Type{$R{A,TA}}, ::Type{$R{B,TB}}) where {A,TA,B,TB}
            (_promotable(A, B) && fieldnames(TA) === fieldnames(TB)) || return Union{}
            FT = promote_type(A, B)
            return $R{FT,_tracer_type(TA, FT)}
        end
        function Base.convert(::Type{$R{FT,T}}, x::$R) where {FT,T}
            same_precision(zero(FT), x)
            fieldnames(T) === keys(x.tracers) ||
                throw(ArgumentError("tracer sets of canopy exchange records differ"))
            return map(Base.Fix1(convert, FT), x)
        end
    end
end
@inline _zero_tracers(::Type{NamedTuple{names,V}}) where {names,V} =
    NamedTuple{names}(map(zero, fieldtypes(V)))
Base.zero(::Type{ExchangeRates{FT,T}}) where {FT,T} =
    ExchangeRates(zero(FT), zero(FT), _zero_tracers(T))
Base.zero(::Type{MoistAirTendencies{FT,T}}) where {FT,T} =
    MoistAirTendencies(zero(FT), zero(FT), zero(FT), _zero_tracers(T))

@inline Base.:-(a::Rates) = map(-, a)
@inline Base.:*(a::Rates, s::Real) = (same_precision(a, s); map(Base.Fix2(*, s), a))
@inline Base.:*(s::Real, a::Rates) = (same_precision(a, s); map(Base.Fix1(*, s), a))
@inline Base.:/(a::Rates, s::Real) = (same_precision(a, s); map(Base.Fix2(/, s), a))
@inline Base.zero(a::Rates) = map(zero, a)

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
@inline valid_exchange(a::Rates) = all(all_finite, _scalars(a)) && all(all_finite, a.tracers)

# ============================================================================
# Value equality, hashing and broadcasting for every record
# ============================================================================

@inline _fields(x::T) where {T} = ntuple(i -> getfield(x, i), Val(fieldcount(T)))
_tag(x) = (nameof(typeof(x)), _measure(x))
_measure(_) = nothing

Base.Broadcast.broadcastable(x::ScalarRecord) = Ref(x)
_same_kind(a, b) = nameof(typeof(a)) === nameof(typeof(b)) && _measure(a) === _measure(b)
Base.:(==)(a::ScalarRecord, b::ScalarRecord) =
    _same_kind(a, b) && _all_equal(_fields(a), _fields(b))
# Componentwise scalar `==` that returns a `Bool` (tuple `==` may return `missing`).
_all_equal(a::Tuple, b::Tuple) = length(a) == length(b) && all(map(_value_equal, a, b))
_value_equal(x::NamedTuple, y::NamedTuple) =
    keys(x) === keys(y) && _all_equal(values(x), values(y))
_value_equal(x, y) = (x == y)::Bool
Base.isequal(a::ScalarRecord, b::ScalarRecord) =
    _same_kind(a, b) && isequal(_fields(a), _fields(b))
# Hash the primal value with zero sign normalised, so records that compare equal by
# value (`-0.0 == 0.0`, or a dual with zero derivative and a plain float) hash equal.
Base.hash(a::ScalarRecord, h::UInt) = _hash_values(_fields(a), hash(_tag(a), h))
_hash_values(values, h::UInt) = foldl((acc, x) -> _hash_value(x, acc), values; init = h)
_hash_value(x::ScalarRecord, h::UInt) = hash(x, h)
_hash_value(x::NamedTuple, h::UInt) = _hash_values(values(x), hash(keys(x), h))
_hash_value(x, h::UInt) = (p = primal(x); hash(ifelse(iszero(p), zero(p), p), h))
