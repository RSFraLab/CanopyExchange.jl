"""
    OxygenContent{FT<:Real}

Represent optional transported molecular oxygen in scalar exchange records.

Absence ([`NoOxygen`](@ref)) is distinct from a computed zero transfer
([`OxygenTransfer`](@ref)). See [`oxygen_value`](@ref).
"""
abstract type OxygenContent{FT<:Real} end

"""
    NoOxygen{FT}()
    NoOxygen(::Type{FT})

Mark oxygen transfer as uncomputed for continuous scalar type `FT`.

[`oxygen_value`](@ref) returns `nothing` on `NoOxygen` instances.

# Examples
```jldoctest
julia> oxygen_value(NoOxygen(Float64)) === nothing
true
```
"""
struct NoOxygen{FT<:Real} <: OxygenContent{FT}
    function NoOxygen{FT}() where {FT<:Real}
        require_continuous_scalar(FT, NoOxygen{FT}, ())
        return new{FT}()
    end
end
@inline NoOxygen(::Type{FT}) where {FT<:Real} = NoOxygen{FT}()

"""
    OxygenTransfer(value)
    OxygenTransfer{FT}(value)

Wrap signed molecular oxygen transfer `value` in mol O₂ per container measure and second.

# Examples
```jldoctest
julia> oxygen_value(OxygenTransfer(1.25e-5))
1.25e-5
```
"""
struct OxygenTransfer{FT<:Real} <: OxygenContent{FT}
    value::FT
    function OxygenTransfer{FT}(value::Real) where {FT<:Real}
        require_continuous_scalar(FT, OxygenTransfer{FT}, (value,))
        return new{FT}(convert(FT, value))
    end
end
@inline function OxygenTransfer(value::FT) where {FT<:Real}
    require_continuous_scalar(FT, OxygenTransfer, (value,))
    return OxygenTransfer{FT}(value)
end

"""
    oxygen_value(oxygen::OxygenContent)

Return `nothing` for [`NoOxygen`](@ref) or the scalar mol O₂ transfer for [`OxygenTransfer`](@ref).

# Examples
```jldoctest
julia> (oxygen_value(NoOxygen(Float32)), oxygen_value(OxygenTransfer(2.5f-5)))
(nothing, 2.5f-5)
```

See also [`NoOxygen`](@ref), [`OxygenTransfer`](@ref).
"""
function oxygen_value end
@inline oxygen_value(::NoOxygen) = nothing
@inline oxygen_value(oxygen::OxygenTransfer) = oxygen.value

"""
    compatible_oxygen(::Type{O1}, ::Type{O2}) -> Bool

Return `true` if `O1` and `O2` represent the same oxygen mode over compatible continuous scalars.

Both types must be [`NoOxygen`](@ref) or both must be [`OxygenTransfer`](@ref),
and their scalar parameters must satisfy [`CanopyExchange.compatible_scalar`](@ref).

# Examples
```jldoctest
julia> CanopyExchange.compatible_oxygen(NoOxygen{Float64}, NoOxygen{Float64})
true

julia> CanopyExchange.compatible_oxygen(NoOxygen{Float64}, OxygenTransfer{Float64})
false
```

See also [`CanopyExchange.compatible_scalar`](@ref), [`OxygenContent`](@ref).
"""
function compatible_oxygen end
@inline compatible_oxygen(::Type{NoOxygen{A}}, ::Type{NoOxygen{B}}) where {A,B} =
    compatible_scalar(A, B)
@inline compatible_oxygen(::Type{OxygenTransfer{A}}, ::Type{OxygenTransfer{B}}) where {A,B} =
    compatible_scalar(A, B)
@inline compatible_oxygen(::Type, ::Type) = false

@inline _promote_oxygen_type(::Type{NoOxygen{A}}, ::Type{FT}) where {A,FT} = NoOxygen{FT}
@inline _promote_oxygen_type(::Type{OxygenTransfer{A}}, ::Type{FT}) where {A,FT} =
    OxygenTransfer{FT}

@inline _negate_oxygen(::NoOxygen{FT}) where {FT} = NoOxygen{FT}()
@inline _negate_oxygen(oxygen::OxygenTransfer{FT}) where {FT} = OxygenTransfer{FT}(-oxygen.value)

@inline _scale_oxygen_left(factor::FT, ::NoOxygen{FT}) where {FT} = NoOxygen{FT}()
@inline _scale_oxygen_left(factor::FT, oxygen::OxygenTransfer{FT}) where {FT} =
    OxygenTransfer{FT}(factor * oxygen.value)

@inline _scale_oxygen_right(::NoOxygen{FT}, factor::FT) where {FT} = NoOxygen{FT}()
@inline _scale_oxygen_right(oxygen::OxygenTransfer{FT}, factor::FT) where {FT} =
    OxygenTransfer{FT}(oxygen.value * factor)

@inline _div_oxygen(::NoOxygen{FT}, denom::FT) where {FT} = NoOxygen{FT}()
@inline _div_oxygen(oxygen::OxygenTransfer{FT}, denom::FT) where {FT} =
    OxygenTransfer{FT}(oxygen.value / denom)

@inline _add_oxygen(::NoOxygen{FT}, ::NoOxygen{FT}) where {FT} = NoOxygen{FT}()
@inline _add_oxygen(a::OxygenTransfer{FT}, b::OxygenTransfer{FT}) where {FT} =
    OxygenTransfer{FT}(a.value + b.value)

@inline _sub_oxygen(::NoOxygen{FT}, ::NoOxygen{FT}) where {FT} = NoOxygen{FT}()
@inline _sub_oxygen(a::OxygenTransfer{FT}, b::OxygenTransfer{FT}) where {FT} =
    OxygenTransfer{FT}(a.value - b.value)

@inline valid_oxygen(::NoOxygen) = true
@inline valid_oxygen(oxygen::OxygenTransfer) = finite_scalar(oxygen.value)

"""
    ExchangeRates(energy, water, carbon, oxygen=NoOxygen(FT))
    ExchangeRates(; energy, water, carbon, oxygen=nothing)

Bundle conserved transfer rates for energy [J], water [kg H₂O], carbon [mol C],
and optional oxygen [mol O₂], each per second and per wrapper spatial measure.

All numeric fields and `oxygen` share one continuous scalar type `FT` (or
promote cleanly across compatible primal and automatic-differentiation scalars).
Passing `oxygen=nothing` in the keyword constructor selects `NoOxygen{FT}()`.
Carbon exchange represents `CO₂` molar transfer (`1` mol C per mol `CO₂`), not
biochemical enthalpy of assimilated carbohydrate.

`energy` represents total energy transfer under a caller-declared thermodynamic
convention rather than sensible or latent heat in isolation. Host adapters own
conversions between conserved density tendencies and temperature or humidity.

# Examples
```jldoctest
julia> r = ExchangeRates(120.0, 2.0e-5, -1.0e-5)
ExchangeRates{Float64, NoOxygen{Float64}}(120.0, 2.0e-5, -1.0e-5, NoOxygen{Float64}())

julia> valid_exchange(r)
true
```

See also [`LeafFlux`](@ref), [`StemFlux`](@ref), [`FaceFlux`](@ref), [`VolumeSource`](@ref).
"""
struct ExchangeRates{FT<:Real,O<:OxygenContent{FT}}
    energy::FT
    water::FT
    carbon::FT
    oxygen::O
    function ExchangeRates{FT,O}(
        energy::Real,
        water::Real,
        carbon::Real,
        oxygen::OxygenContent,
    ) where {FT<:Real,O<:OxygenContent{FT}}
        require_continuous_scalar(FT, ExchangeRates{FT,O}, (energy, water, carbon, oxygen))
        compatible_oxygen(typeof(oxygen), O) ||
            throw(MethodError(ExchangeRates{FT,O}, (energy, water, carbon, oxygen)))
        return new{FT,O}(
            convert(FT, energy),
            convert(FT, water),
            convert(FT, carbon),
            convert(O, oxygen),
        )
    end
end

@inline function ExchangeRates(
    energy::FT,
    water::FT,
    carbon::FT,
    oxygen::O,
) where {FT<:Real,O<:OxygenContent{FT}}
    require_continuous_scalar(FT, ExchangeRates, (energy, water, carbon, oxygen))
    return ExchangeRates{FT,O}(energy, water, carbon, oxygen)
end

@inline function ExchangeRates(
    energy::Real,
    water::Real,
    carbon::Real,
    oxygen::O,
) where {FO<:Real,O<:OxygenContent{FO}}
    FT = promote_compatible_scalar(ExchangeRates, energy, water, carbon, zero(FO))
    OT = _promote_oxygen_type(O, FT)
    return ExchangeRates{FT,OT}(energy, water, carbon, oxygen)
end

@inline function ExchangeRates(energy::Real, water::Real, carbon::Real)
    FT = promote_compatible_scalar(ExchangeRates, energy, water, carbon)
    return ExchangeRates{FT,NoOxygen{FT}}(energy, water, carbon, NoOxygen{FT}())
end

@inline function ExchangeRates(;
    energy::Real,
    water::Real,
    carbon::Real,
    oxygen::Union{Nothing,OxygenContent} = nothing,
)
    if oxygen === nothing
        return ExchangeRates(energy, water, carbon)
    end
    return ExchangeRates(energy, water, carbon, oxygen)
end

"""
    AbstractExchangeFlux{FT<:Real,O<:OxygenContent{FT}}

Supertype for measure-tagged exchange wrappers ([`LeafFlux`](@ref),
[`StemFlux`](@ref), [`FaceFlux`](@ref), and [`VolumeSource`](@ref)).

Distinct concrete subtypes prevent accidental addition or argument substitution
across one-sided leaf area, one-sided stem area, horizontal ground area, and
bulk canopy volume.
"""
abstract type AbstractExchangeFlux{FT<:Real,O<:OxygenContent{FT}} end

"""
    LeafFlux(rates::ExchangeRates)

Wrap `rates` as a flux per m² one-sided leaf area per second, positive from leaf into air.

# Examples
```jldoctest
julia> f = LeafFlux(ExchangeRates(100.0, 3.0e-5, -1.2e-5));

julia> f.rates.carbon
-1.2e-5
```

See also [`sensible_water_flux`](@ref), [`leaf_air_sources`](@ref), [`StemFlux`](@ref).
"""
struct LeafFlux{FT<:Real,O<:OxygenContent{FT}} <: AbstractExchangeFlux{FT,O}
    rates::ExchangeRates{FT,O}
    function LeafFlux{FT,O}(rates::ExchangeRates) where {FT<:Real,O<:OxygenContent{FT}}
        return new{FT,O}(convert(ExchangeRates{FT,O}, rates))
    end
end
@inline LeafFlux(rates::ExchangeRates{FT,O}) where {FT<:Real,O<:OxygenContent{FT}} =
    LeafFlux{FT,O}(rates)

"""
    StemFlux(rates::ExchangeRates)

Wrap `rates` as a flux per m² one-sided woody/stem area per second, positive from stem into air.

Keeping `StemFlux` distinct from [`LeafFlux`](@ref) prevents swapping leaf and
stem arguments in [`plant_air_sources`](@ref) or adding fluxes across unequal
leaf (`LAD`) and stem (`SAD`) area denominators.

# Examples
```jldoctest
julia> s = StemFlux(ExchangeRates(25.0, 0.0, 0.0));

julia> s.rates.energy
25.0
```

See also [`sensible_stem_flux`](@ref), [`plant_air_sources`](@ref), [`LeafFlux`](@ref).
"""
struct StemFlux{FT<:Real,O<:OxygenContent{FT}} <: AbstractExchangeFlux{FT,O}
    rates::ExchangeRates{FT,O}
    function StemFlux{FT,O}(rates::ExchangeRates) where {FT<:Real,O<:OxygenContent{FT}}
        return new{FT,O}(convert(ExchangeRates{FT,O}, rates))
    end
end
@inline StemFlux(rates::ExchangeRates{FT,O}) where {FT<:Real,O<:OxygenContent{FT}} =
    StemFlux{FT,O}(rates)

"""
    FaceFlux(rates::ExchangeRates)

Wrap `rates` as a flux per m² horizontal ground area per second, positive upward.

# Examples
```jldoctest
julia> ff = FaceFlux(ExchangeRates(180.0, 5.0e-5, -1.5e-5));

julia> ff.rates.energy
180.0
```

See also [`surface_static_energy_flux`](@ref), [`layer_tendency`](@ref), [`steady_top_flux`](@ref).
"""
struct FaceFlux{FT<:Real,O<:OxygenContent{FT}} <: AbstractExchangeFlux{FT,O}
    rates::ExchangeRates{FT,O}
    function FaceFlux{FT,O}(rates::ExchangeRates) where {FT<:Real,O<:OxygenContent{FT}}
        return new{FT,O}(convert(ExchangeRates{FT,O}, rates))
    end
end
@inline FaceFlux(rates::ExchangeRates{FT,O}) where {FT<:Real,O<:OxygenContent{FT}} =
    FaceFlux{FT,O}(rates)

"""
    VolumeSource(rates::ExchangeRates)

Wrap `rates` as a source per m³ bulk canopy volume per second, positive into the named reservoir.

# Examples
```jldoctest
julia> vs = VolumeSource(ExchangeRates(50.0, 1.0e-5, -0.5e-5));

julia> vs.rates.water
1.0e-5
```

See also [`leaf_air_sources`](@ref), [`plant_air_sources`](@ref), [`moist_air_tendencies`](@ref).
"""
struct VolumeSource{FT<:Real,O<:OxygenContent{FT}} <: AbstractExchangeFlux{FT,O}
    rates::ExchangeRates{FT,O}
    function VolumeSource{FT,O}(rates::ExchangeRates) where {FT<:Real,O<:OxygenContent{FT}}
        return new{FT,O}(convert(ExchangeRates{FT,O}, rates))
    end
end
@inline VolumeSource(rates::ExchangeRates{FT,O}) where {FT<:Real,O<:OxygenContent{FT}} =
    VolumeSource{FT,O}(rates)

"""
    MoistAirTendencies(rho_total, rho_water, rho_energy, carbon, oxygen=NoOxygen(FT))
    MoistAirTendencies(; rho_total, rho_water, rho_energy, carbon, oxygen=nothing)

Store conserved atmospheric gas-volume source tendencies per m³ gas volume per second.

Fields map directly to conserved moist-air prognostic budgets (`ClimaAtmos` /
porous-canopy hosts with conserved dry-air mass `S_ρd = 0`):
- `rho_total`: total moist-air mass source [kg m⁻³ gas s⁻¹] (`== rho_water` when `S_ρd = 0`),
- `rho_water`: total water mass source [kg H₂O m⁻³ gas s⁻¹],
- `rho_energy`: total energy density source [J m⁻³ gas s⁻¹] in the host's datum,
- `carbon`: passive molar `CO₂` source [mol C m⁻³ gas s⁻¹],
- `oxygen`: optional molar `O₂` source [mol O₂ m⁻³ gas s⁻¹] (`nothing` keyword selects `NoOxygen{FT}()`).

Momentum, canopy drag, and subgrid turbulent kinetic energy are excluded.

# Examples
```jldoctest
julia> m = MoistAirTendencies(; rho_total=4.0e-5, rho_water=4.0e-5, rho_energy=150.0, carbon=-1.0e-5);

julia> (m.rho_total, m.rho_energy)
(4.0e-5, 150.0)
```

See also [`moist_air_tendencies`](@ref), [`VolumeSource`](@ref).
"""
struct MoistAirTendencies{FT<:Real,O<:OxygenContent{FT}}
    rho_total::FT
    rho_water::FT
    rho_energy::FT
    carbon::FT
    oxygen::O
    function MoistAirTendencies{FT,O}(
        rho_total::Real,
        rho_water::Real,
        rho_energy::Real,
        carbon::Real,
        oxygen::OxygenContent,
    ) where {FT<:Real,O<:OxygenContent{FT}}
        require_continuous_scalar(
            FT,
            MoistAirTendencies{FT,O},
            (rho_total, rho_water, rho_energy, carbon, oxygen),
        )
        compatible_oxygen(typeof(oxygen), O) || throw(
            MethodError(
                MoistAirTendencies{FT,O},
                (rho_total, rho_water, rho_energy, carbon, oxygen),
            ),
        )
        return new{FT,O}(
            convert(FT, rho_total),
            convert(FT, rho_water),
            convert(FT, rho_energy),
            convert(FT, carbon),
            convert(O, oxygen),
        )
    end
end

@inline function MoistAirTendencies(
    rho_total::FT,
    rho_water::FT,
    rho_energy::FT,
    carbon::FT,
    oxygen::O,
) where {FT<:Real,O<:OxygenContent{FT}}
    require_continuous_scalar(
        FT,
        MoistAirTendencies,
        (rho_total, rho_water, rho_energy, carbon, oxygen),
    )
    return MoistAirTendencies{FT,O}(rho_total, rho_water, rho_energy, carbon, oxygen)
end

@inline function MoistAirTendencies(
    rho_total::Real,
    rho_water::Real,
    rho_energy::Real,
    carbon::Real,
    oxygen::O,
) where {FO<:Real,O<:OxygenContent{FO}}
    FT = promote_compatible_scalar(
        MoistAirTendencies,
        rho_total,
        rho_water,
        rho_energy,
        carbon,
        zero(FO),
    )
    OT = _promote_oxygen_type(O, FT)
    return MoistAirTendencies{FT,OT}(rho_total, rho_water, rho_energy, carbon, oxygen)
end

@inline function MoistAirTendencies(
    rho_total::Real,
    rho_water::Real,
    rho_energy::Real,
    carbon::Real,
)
    FT = promote_compatible_scalar(
        MoistAirTendencies,
        rho_total,
        rho_water,
        rho_energy,
        carbon,
    )
    return MoistAirTendencies{FT,NoOxygen{FT}}(
        rho_total,
        rho_water,
        rho_energy,
        carbon,
        NoOxygen{FT}(),
    )
end

@inline function MoistAirTendencies(;
    rho_total::Real,
    rho_water::Real,
    rho_energy::Real,
    carbon::Real,
    oxygen::Union{Nothing,OxygenContent} = nothing,
)
    if oxygen === nothing
        return MoistAirTendencies(rho_total, rho_water, rho_energy, carbon)
    end
    return MoistAirTendencies(rho_total, rho_water, rho_energy, carbon, oxygen)
end

"""
    AbstractLayerGeometry{FT<:Real}

Supertype for canopy finite-volume layer geometry descriptors ([`LayerGeometry`](@ref)
and [`PlantLayerGeometry`](@ref)).
"""
abstract type AbstractLayerGeometry{FT<:Real} end

"""
    LayerGeometry(thickness, leaf_area_density)

Store finite positive layer `thickness` [m] and finite nonnegative one-sided
`leaf_area_density` [m² leaf m⁻³ bulk canopy].

Bulk volume equals horizontal ground footprint multiplied by `thickness` with no
implicit air-porosity factor. Throw `ArgumentError` if `thickness ≤ 0`,
`leaf_area_density < 0`, or any scalar component (including dual partials) is
nonfinite.

# Examples
```jldoctest
julia> g = LayerGeometry(1.5, 0.8)
LayerGeometry{Float64}(1.5, 0.8)
```
"""
struct LayerGeometry{FT<:Real} <: AbstractLayerGeometry{FT}
    thickness::FT
    leaf_area_density::FT
    function LayerGeometry{FT}(thickness::Real, leaf_area_density::Real) where {FT<:Real}
        require_continuous_scalar(FT, LayerGeometry{FT}, (thickness, leaf_area_density))
        dz = convert(FT, thickness)
        lad = convert(FT, leaf_area_density)
        P = primal_float_type(FT)
        finite_scalar(dz) && primal_scalar(dz) > zero(P) ||
            throw(ArgumentError("layer thickness must be finite and positive"))
        finite_scalar(lad) && primal_scalar(lad) >= zero(P) ||
            throw(ArgumentError("leaf area density must be finite and nonnegative"))
        return new{FT}(dz, lad)
    end
end

@inline function LayerGeometry(thickness::FT, leaf_area_density::FT) where {FT<:Real}
    return LayerGeometry{FT}(thickness, leaf_area_density)
end

@inline function LayerGeometry(thickness::Real, leaf_area_density::Real)
    FT = promote_compatible_scalar(LayerGeometry, thickness, leaf_area_density)
    return LayerGeometry{FT}(thickness, leaf_area_density)
end

"""
    PlantLayerGeometry(thickness, leaf_area_density, stem_area_density=zero(FT))

Store finite positive `thickness` [m] with nonnegative one-sided `leaf_area_density`
[m² leaf m⁻³] and `stem_area_density` [m² stem m⁻³] (`PAI = LAI + SAI` as in `ClimaLand`).

Throw `ArgumentError` if `thickness ≤ 0`, either area density is negative, or any
scalar component (including dual partials) is nonfinite.

# Examples
```jldoctest
julia> g = PlantLayerGeometry(2.0, 1.2, 0.3)
PlantLayerGeometry{Float64}(2.0, 1.2, 0.3)

julia> plant_area_density(g)
1.5
```
"""
struct PlantLayerGeometry{FT<:Real} <: AbstractLayerGeometry{FT}
    thickness::FT
    leaf_area_density::FT
    stem_area_density::FT
    function PlantLayerGeometry{FT}(
        thickness::Real,
        leaf_area_density::Real,
        stem_area_density::Real,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            PlantLayerGeometry{FT},
            (thickness, leaf_area_density, stem_area_density),
        )
        dz = convert(FT, thickness)
        lad = convert(FT, leaf_area_density)
        sad = convert(FT, stem_area_density)
        P = primal_float_type(FT)
        finite_scalar(dz) && primal_scalar(dz) > zero(P) ||
            throw(ArgumentError("layer thickness must be finite and positive"))
        finite_scalar(lad) && primal_scalar(lad) >= zero(P) ||
            throw(ArgumentError("leaf area density must be finite and nonnegative"))
        finite_scalar(sad) && primal_scalar(sad) >= zero(P) ||
            throw(ArgumentError("stem area density must be finite and nonnegative"))
        return new{FT}(dz, lad, sad)
    end
end

@inline function PlantLayerGeometry(
    thickness::FT,
    leaf_area_density::FT,
    stem_area_density::FT,
) where {FT<:Real}
    return PlantLayerGeometry{FT}(thickness, leaf_area_density, stem_area_density)
end

@inline function PlantLayerGeometry(
    thickness::Real,
    leaf_area_density::Real,
    stem_area_density::Real,
)
    FT = promote_compatible_scalar(
        PlantLayerGeometry,
        thickness,
        leaf_area_density,
        stem_area_density,
    )
    return PlantLayerGeometry{FT}(thickness, leaf_area_density, stem_area_density)
end

@inline function PlantLayerGeometry(thickness::Real, leaf_area_density::Real)
    FT = promote_compatible_scalar(PlantLayerGeometry, thickness, leaf_area_density)
    return PlantLayerGeometry{FT}(thickness, leaf_area_density, zero(FT))
end

# Invariant Base.eltype definitions (zero unbound parameters on Union{})
Base.eltype(::Type{NoOxygen{FT}}) where {FT<:Real} = FT
Base.eltype(::Type{OxygenTransfer{FT}}) where {FT<:Real} = FT
Base.eltype(::Type{OxygenContent{FT}}) where {FT<:Real} = FT
Base.eltype(::Type{ExchangeRates{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}} = FT
Base.eltype(::Type{LeafFlux{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}} = FT
Base.eltype(::Type{StemFlux{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}} = FT
Base.eltype(::Type{FaceFlux{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}} = FT
Base.eltype(::Type{VolumeSource{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}} = FT
Base.eltype(::Type{AbstractExchangeFlux{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}} = FT
Base.eltype(::Type{MoistAirTendencies{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}} = FT
Base.eltype(::Type{LayerGeometry{FT}}) where {FT<:Real} = FT
Base.eltype(::Type{PlantLayerGeometry{FT}}) where {FT<:Real} = FT
Base.eltype(::Type{AbstractLayerGeometry{FT}}) where {FT<:Real} = FT

const _ExchangeScalarRecord = Union{
    OxygenContent,
    ExchangeRates,
    AbstractExchangeFlux,
    MoistAirTendencies,
    AbstractLayerGeometry,
}

Base.eltype(x::_ExchangeScalarRecord) = eltype(typeof(x))
broadcastable(x::_ExchangeScalarRecord) = Ref(x)

# Compatible-scalar Base.convert methods
Base.convert(::Type{NoOxygen{FT}}, x::NoOxygen{FT}) where {FT<:Real} = x
function Base.convert(::Type{NoOxygen{FT}}, x::NoOxygen{FT2}) where {FT<:Real,FT2<:Real}
    compatible_scalar(FT, FT2) || throw(MethodError(convert, (NoOxygen{FT}, x)))
    return NoOxygen{FT}()
end

Base.convert(::Type{OxygenTransfer{FT}}, x::OxygenTransfer{FT}) where {FT<:Real} = x
function Base.convert(
    ::Type{OxygenTransfer{FT}},
    x::OxygenTransfer{FT2},
) where {FT<:Real,FT2<:Real}
    compatible_scalar(FT, FT2) || throw(MethodError(convert, (OxygenTransfer{FT}, x)))
    return OxygenTransfer{FT}(convert(FT, x.value))
end

Base.convert(
    ::Type{ExchangeRates{FT,O}},
    x::ExchangeRates{FT,O},
) where {FT<:Real,O<:OxygenContent{FT}} = x
function Base.convert(
    ::Type{ExchangeRates{FT,O}},
    x::ExchangeRates{FT2,O2},
) where {FT<:Real,O<:OxygenContent{FT},FT2<:Real,O2<:OxygenContent{FT2}}
    compatible_oxygen(O, O2) || throw(MethodError(convert, (ExchangeRates{FT,O}, x)))
    return ExchangeRates{FT,O}(
        convert(FT, x.energy),
        convert(FT, x.water),
        convert(FT, x.carbon),
        convert(O, x.oxygen),
    )
end

for FluxT in (:LeafFlux, :StemFlux, :FaceFlux, :VolumeSource)
    @eval begin
        Base.convert(
            ::Type{$FluxT{FT,O}},
            x::$FluxT{FT,O},
        ) where {FT<:Real,O<:OxygenContent{FT}} = x
        function Base.convert(
            ::Type{$FluxT{FT,O}},
            x::$FluxT{FT2,O2},
        ) where {FT<:Real,O<:OxygenContent{FT},FT2<:Real,O2<:OxygenContent{FT2}}
            compatible_oxygen(O, O2) || throw(MethodError(convert, ($FluxT{FT,O}, x)))
            return $FluxT{FT,O}(convert(ExchangeRates{FT,O}, x.rates))
        end
    end
end

Base.convert(
    ::Type{MoistAirTendencies{FT,O}},
    x::MoistAirTendencies{FT,O},
) where {FT<:Real,O<:OxygenContent{FT}} = x
function Base.convert(
    ::Type{MoistAirTendencies{FT,O}},
    x::MoistAirTendencies{FT2,O2},
) where {FT<:Real,O<:OxygenContent{FT},FT2<:Real,O2<:OxygenContent{FT2}}
    compatible_oxygen(O, O2) || throw(MethodError(convert, (MoistAirTendencies{FT,O}, x)))
    return MoistAirTendencies{FT,O}(
        convert(FT, x.rho_total),
        convert(FT, x.rho_water),
        convert(FT, x.rho_energy),
        convert(FT, x.carbon),
        convert(O, x.oxygen),
    )
end

Base.convert(::Type{LayerGeometry{FT}}, x::LayerGeometry{FT}) where {FT<:Real} = x
function Base.convert(
    ::Type{LayerGeometry{FT}},
    x::LayerGeometry{FT2},
) where {FT<:Real,FT2<:Real}
    compatible_scalar(FT, FT2) || throw(MethodError(convert, (LayerGeometry{FT}, x)))
    return LayerGeometry{FT}(convert(FT, x.thickness), convert(FT, x.leaf_area_density))
end

Base.convert(::Type{PlantLayerGeometry{FT}}, x::PlantLayerGeometry{FT}) where {FT<:Real} = x
function Base.convert(
    ::Type{PlantLayerGeometry{FT}},
    x::PlantLayerGeometry{FT2},
) where {FT<:Real,FT2<:Real}
    compatible_scalar(FT, FT2) || throw(MethodError(convert, (PlantLayerGeometry{FT}, x)))
    return PlantLayerGeometry{FT}(
        convert(FT, x.thickness),
        convert(FT, x.leaf_area_density),
        convert(FT, x.stem_area_density),
    )
end

# Compatible-scalar Base.promote_rule methods (concrete AD container promotion; Union{} on mismatch)
function promote_rule(::Type{NoOxygen{A}}, ::Type{NoOxygen{B}}) where {A<:Real,B<:Real}
    compatible_scalar(A, B) || return Union{}
    FT = promote_type(A, B)
    is_continuous_scalar(FT) || return Union{}
    return NoOxygen{FT}
end

function promote_rule(
    ::Type{OxygenTransfer{A}},
    ::Type{OxygenTransfer{B}},
) where {A<:Real,B<:Real}
    compatible_scalar(A, B) || return Union{}
    FT = promote_type(A, B)
    is_continuous_scalar(FT) || return Union{}
    return OxygenTransfer{FT}
end

for RecT in (
    :ExchangeRates,
    :LeafFlux,
    :StemFlux,
    :FaceFlux,
    :VolumeSource,
    :MoistAirTendencies,
)
    @eval begin
        function promote_rule(
            ::Type{$RecT{FA,OA}},
            ::Type{$RecT{FB,OB}},
        ) where {FA<:Real,OA<:OxygenContent{FA},FB<:Real,OB<:OxygenContent{FB}}
            compatible_oxygen(OA, OB) || return Union{}
            FT = promote_type(FA, FB)
            is_continuous_scalar(FT) || return Union{}
            return $RecT{FT,_promote_oxygen_type(OA, FT)}
        end
    end
end

for GeomT in (:LayerGeometry, :PlantLayerGeometry)
    @eval begin
        function promote_rule(
            ::Type{$GeomT{A}},
            ::Type{$GeomT{B}},
        ) where {A<:Real,B<:Real}
            compatible_scalar(A, B) || return Union{}
            FT = promote_type(A, B)
            is_continuous_scalar(FT) || return Union{}
            return $GeomT{FT}
        end
    end
end

# Unified ==, isequal, and hash (both isequal(a, b) => hash(a) == hash(b) and a == b => hash(a) == hash(b))
Base.:(==)(::NoOxygen, ::NoOxygen) = true
Base.isequal(::NoOxygen, ::NoOxygen) = true
Base.hash(::NoOxygen, h::UInt) = hash(:NoOxygen, h)

Base.:(==)(a::OxygenTransfer, b::OxygenTransfer) = a.value == b.value
Base.isequal(a::OxygenTransfer, b::OxygenTransfer) = isequal(a.value, b.value)
Base.hash(a::OxygenTransfer, h::UInt) = _hash_scalar(a.value, hash(:OxygenTransfer, h))

Base.:(==)(::OxygenContent, ::OxygenContent) = false
Base.isequal(::OxygenContent, ::OxygenContent) = false

function Base.:(==)(a::ExchangeRates, b::ExchangeRates)
    return (a.oxygen == b.oxygen) &&
           (a.energy == b.energy) &&
           (a.water == b.water) &&
           (a.carbon == b.carbon)
end
function Base.isequal(a::ExchangeRates, b::ExchangeRates)
    return isequal(a.oxygen, b.oxygen) &&
           isequal(a.energy, b.energy) &&
           isequal(a.water, b.water) &&
           isequal(a.carbon, b.carbon)
end
function Base.hash(a::ExchangeRates, h::UInt)
    h = hash(:ExchangeRates, h)
    h = _hash_scalar(a.energy, h)
    h = _hash_scalar(a.water, h)
    h = _hash_scalar(a.carbon, h)
    return hash(a.oxygen, h)
end

for (FluxT, Sym) in (
    (:LeafFlux, :(:LeafFlux)),
    (:StemFlux, :(:StemFlux)),
    (:FaceFlux, :(:FaceFlux)),
    (:VolumeSource, :(:VolumeSource)),
)
    @eval begin
        Base.:(==)(a::$FluxT, b::$FluxT) = a.rates == b.rates
        Base.isequal(a::$FluxT, b::$FluxT) = isequal(a.rates, b.rates)
        Base.hash(a::$FluxT, h::UInt) = hash(a.rates, hash($Sym, h))
    end
end
Base.:(==)(::AbstractExchangeFlux, ::AbstractExchangeFlux) = false
Base.isequal(::AbstractExchangeFlux, ::AbstractExchangeFlux) = false

function Base.:(==)(a::MoistAirTendencies, b::MoistAirTendencies)
    return (a.oxygen == b.oxygen) &&
           (a.rho_total == b.rho_total) &&
           (a.rho_water == b.rho_water) &&
           (a.rho_energy == b.rho_energy) &&
           (a.carbon == b.carbon)
end
function Base.isequal(a::MoistAirTendencies, b::MoistAirTendencies)
    return isequal(a.oxygen, b.oxygen) &&
           isequal(a.rho_total, b.rho_total) &&
           isequal(a.rho_water, b.rho_water) &&
           isequal(a.rho_energy, b.rho_energy) &&
           isequal(a.carbon, b.carbon)
end
function Base.hash(a::MoistAirTendencies, h::UInt)
    h = hash(:MoistAirTendencies, h)
    h = _hash_scalar(a.rho_total, h)
    h = _hash_scalar(a.rho_water, h)
    h = _hash_scalar(a.rho_energy, h)
    h = _hash_scalar(a.carbon, h)
    return hash(a.oxygen, h)
end

function Base.:(==)(a::LayerGeometry, b::LayerGeometry)
    return (a.thickness == b.thickness) && (a.leaf_area_density == b.leaf_area_density)
end
function Base.isequal(a::LayerGeometry, b::LayerGeometry)
    return isequal(a.thickness, b.thickness) &&
           isequal(a.leaf_area_density, b.leaf_area_density)
end
function Base.hash(a::LayerGeometry, h::UInt)
    h = hash(:LayerGeometry, h)
    h = _hash_scalar(a.thickness, h)
    return _hash_scalar(a.leaf_area_density, h)
end

function Base.:(==)(a::PlantLayerGeometry, b::PlantLayerGeometry)
    return (a.thickness == b.thickness) &&
           (a.leaf_area_density == b.leaf_area_density) &&
           (a.stem_area_density == b.stem_area_density)
end
function Base.isequal(a::PlantLayerGeometry, b::PlantLayerGeometry)
    return isequal(a.thickness, b.thickness) &&
           isequal(a.leaf_area_density, b.leaf_area_density) &&
           isequal(a.stem_area_density, b.stem_area_density)
end
function Base.hash(a::PlantLayerGeometry, h::UInt)
    h = hash(:PlantLayerGeometry, h)
    h = _hash_scalar(a.thickness, h)
    h = _hash_scalar(a.leaf_area_density, h)
    return _hash_scalar(a.stem_area_density, h)
end
Base.:(==)(::AbstractLayerGeometry, ::AbstractLayerGeometry) = false
Base.isequal(::AbstractLayerGeometry, ::AbstractLayerGeometry) = false

# Base.zero on flux and tendency records (omitted on positive-thickness geometries)
Base.zero(::Type{NoOxygen{FT}}) where {FT<:Real} = NoOxygen{FT}()
Base.zero(x::NoOxygen) = zero(typeof(x))

Base.zero(::Type{OxygenTransfer{FT}}) where {FT<:Real} = OxygenTransfer{FT}(zero(FT))
Base.zero(x::OxygenTransfer) = zero(typeof(x))

function Base.zero(::Type{ExchangeRates{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}}
    isconcretetype(O) || throw(MethodError(zero, (ExchangeRates{FT,O},)))
    return ExchangeRates{FT,O}(zero(FT), zero(FT), zero(FT), zero(O))
end
Base.zero(x::ExchangeRates) = zero(typeof(x))

for FluxT in (:LeafFlux, :StemFlux, :FaceFlux, :VolumeSource)
    @eval begin
        function Base.zero(::Type{$FluxT{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}}
            return $FluxT{FT,O}(zero(ExchangeRates{FT,O}))
        end
        Base.zero(x::$FluxT) = zero(typeof(x))
    end
end

function Base.zero(::Type{MoistAirTendencies{FT,O}}) where {FT<:Real,O<:OxygenContent{FT}}
    isconcretetype(O) || throw(MethodError(zero, (MoistAirTendencies{FT,O},)))
    return MoistAirTendencies{FT,O}(zero(FT), zero(FT), zero(FT), zero(FT), zero(O))
end
Base.zero(x::MoistAirTendencies) = zero(typeof(x))

# Measure-safe linear algebra (preserving unary sign-bit negation and fieldwise division)
@inline function Base.:-(r::ExchangeRates{FT,O}) where {FT,O}
    return ExchangeRates{FT,O}(-r.energy, -r.water, -r.carbon, _negate_oxygen(r.oxygen))
end

@inline function Base.:+(
    a::ExchangeRates{FA,OA},
    b::ExchangeRates{FB,OB},
) where {FA,OA,FB,OB}
    compatible_oxygen(OA, OB) || throw(MethodError(+, (a, b)))
    FT = promote_compatible_scalar(+, zero(FA), zero(FB))
    OT = _promote_oxygen_type(OA, FT)
    aa = convert(ExchangeRates{FT,OT}, a)
    bb = convert(ExchangeRates{FT,OT}, b)
    return ExchangeRates{FT,OT}(
        aa.energy + bb.energy,
        aa.water + bb.water,
        aa.carbon + bb.carbon,
        _add_oxygen(aa.oxygen, bb.oxygen),
    )
end

@inline function Base.:-(
    a::ExchangeRates{FA,OA},
    b::ExchangeRates{FB,OB},
) where {FA,OA,FB,OB}
    compatible_oxygen(OA, OB) || throw(MethodError(-, (a, b)))
    FT = promote_compatible_scalar(-, zero(FA), zero(FB))
    OT = _promote_oxygen_type(OA, FT)
    aa = convert(ExchangeRates{FT,OT}, a)
    bb = convert(ExchangeRates{FT,OT}, b)
    return ExchangeRates{FT,OT}(
        aa.energy - bb.energy,
        aa.water - bb.water,
        aa.carbon - bb.carbon,
        _sub_oxygen(aa.oxygen, bb.oxygen),
    )
end

@inline function Base.:*(α::Real, r::ExchangeRates{FR,OR}) where {FR,OR}
    FT = promote_compatible_scalar(*, α, zero(FR))
    OT = _promote_oxygen_type(OR, FT)
    s = convert(FT, α)
    rr = convert(ExchangeRates{FT,OT}, r)
    return ExchangeRates{FT,OT}(
        s * rr.energy,
        s * rr.water,
        s * rr.carbon,
        _scale_oxygen_left(s, rr.oxygen),
    )
end

@inline function Base.:*(r::ExchangeRates{FR,OR}, α::Real) where {FR,OR}
    FT = promote_compatible_scalar(*, zero(FR), α)
    OT = _promote_oxygen_type(OR, FT)
    s = convert(FT, α)
    rr = convert(ExchangeRates{FT,OT}, r)
    return ExchangeRates{FT,OT}(
        rr.energy * s,
        rr.water * s,
        rr.carbon * s,
        _scale_oxygen_right(rr.oxygen, s),
    )
end

@inline function Base.:/(r::ExchangeRates{FR,OR}, α::Real) where {FR,OR}
    FT = promote_compatible_scalar(/, zero(FR), α)
    OT = _promote_oxygen_type(OR, FT)
    s = convert(FT, α)
    rr = convert(ExchangeRates{FT,OT}, r)
    return ExchangeRates{FT,OT}(
        rr.energy / s,
        rr.water / s,
        rr.carbon / s,
        _div_oxygen(rr.oxygen, s),
    )
end

for FluxT in (:LeafFlux, :StemFlux, :FaceFlux, :VolumeSource)
    @eval begin
        @inline Base.:-(x::$FluxT) = $FluxT(-x.rates)
        @inline Base.:+(a::$FluxT, b::$FluxT) = $FluxT(a.rates + b.rates)
        @inline Base.:-(a::$FluxT, b::$FluxT) = $FluxT(a.rates - b.rates)
        @inline Base.:*(α::Real, x::$FluxT) = $FluxT(α * x.rates)
        @inline Base.:*(x::$FluxT, α::Real) = $FluxT(x.rates * α)
        @inline Base.:/(x::$FluxT, α::Real) = $FluxT(x.rates / α)
    end
end

@inline function Base.:-(m::MoistAirTendencies{FT,O}) where {FT,O}
    return MoistAirTendencies{FT,O}(
        -m.rho_total,
        -m.rho_water,
        -m.rho_energy,
        -m.carbon,
        _negate_oxygen(m.oxygen),
    )
end

@inline function Base.:+(
    a::MoistAirTendencies{FA,OA},
    b::MoistAirTendencies{FB,OB},
) where {FA,OA,FB,OB}
    compatible_oxygen(OA, OB) || throw(MethodError(+, (a, b)))
    FT = promote_compatible_scalar(+, zero(FA), zero(FB))
    OT = _promote_oxygen_type(OA, FT)
    aa = convert(MoistAirTendencies{FT,OT}, a)
    bb = convert(MoistAirTendencies{FT,OT}, b)
    return MoistAirTendencies{FT,OT}(
        aa.rho_total + bb.rho_total,
        aa.rho_water + bb.rho_water,
        aa.rho_energy + bb.rho_energy,
        aa.carbon + bb.carbon,
        _add_oxygen(aa.oxygen, bb.oxygen),
    )
end

@inline function Base.:-(
    a::MoistAirTendencies{FA,OA},
    b::MoistAirTendencies{FB,OB},
) where {FA,OA,FB,OB}
    compatible_oxygen(OA, OB) || throw(MethodError(-, (a, b)))
    FT = promote_compatible_scalar(-, zero(FA), zero(FB))
    OT = _promote_oxygen_type(OA, FT)
    aa = convert(MoistAirTendencies{FT,OT}, a)
    bb = convert(MoistAirTendencies{FT,OT}, b)
    return MoistAirTendencies{FT,OT}(
        aa.rho_total - bb.rho_total,
        aa.rho_water - bb.rho_water,
        aa.rho_energy - bb.rho_energy,
        aa.carbon - bb.carbon,
        _sub_oxygen(aa.oxygen, bb.oxygen),
    )
end

@inline function Base.:*(α::Real, m::MoistAirTendencies{FM,OM}) where {FM,OM}
    FT = promote_compatible_scalar(*, α, zero(FM))
    OT = _promote_oxygen_type(OM, FT)
    s = convert(FT, α)
    mm = convert(MoistAirTendencies{FT,OT}, m)
    return MoistAirTendencies{FT,OT}(
        s * mm.rho_total,
        s * mm.rho_water,
        s * mm.rho_energy,
        s * mm.carbon,
        _scale_oxygen_left(s, mm.oxygen),
    )
end

@inline function Base.:*(m::MoistAirTendencies{FM,OM}, α::Real) where {FM,OM}
    FT = promote_compatible_scalar(*, zero(FM), α)
    OT = _promote_oxygen_type(OM, FT)
    s = convert(FT, α)
    mm = convert(MoistAirTendencies{FT,OT}, m)
    return MoistAirTendencies{FT,OT}(
        mm.rho_total * s,
        mm.rho_water * s,
        mm.rho_energy * s,
        mm.carbon * s,
        _scale_oxygen_right(mm.oxygen, s),
    )
end

@inline function Base.:/(m::MoistAirTendencies{FM,OM}, α::Real) where {FM,OM}
    FT = promote_compatible_scalar(/, zero(FM), α)
    OT = _promote_oxygen_type(OM, FT)
    s = convert(FT, α)
    mm = convert(MoistAirTendencies{FT,OT}, m)
    return MoistAirTendencies{FT,OT}(
        mm.rho_total / s,
        mm.rho_water / s,
        mm.rho_energy / s,
        mm.carbon / s,
        _div_oxygen(mm.oxygen, s),
    )
end

"""
    valid_exchange(x) -> Bool

Return `true` if every computed scalar component of `x` (including dual partials) is finite.

Accept [`ExchangeRates`](@ref), any [`AbstractExchangeFlux`](@ref), and
[`MoistAirTendencies`](@ref). [`NoOxygen`](@ref) is treated as valid; hosts
remain responsible for thermodynamic admissibility and retry policy.

# Examples
```jldoctest
julia> valid_exchange(LeafFlux(ExchangeRates(100.0, 2.0e-5, -1.0e-5)))
true

julia> valid_exchange(LeafFlux(ExchangeRates(NaN, 2.0e-5, -1.0e-5)))
false
```

See also [`CanopyExchange.finite_scalar`](@ref), [`ExchangeRates`](@ref), [`MoistAirTendencies`](@ref).
"""
function valid_exchange end

@inline valid_exchange(rates::ExchangeRates) =
    finite_scalar(rates.energy) &&
    finite_scalar(rates.water) &&
    finite_scalar(rates.carbon) &&
    valid_oxygen(rates.oxygen)

@inline valid_exchange(flux::AbstractExchangeFlux) = valid_exchange(flux.rates)

@inline valid_exchange(m::MoistAirTendencies) =
    finite_scalar(m.rho_total) &&
    finite_scalar(m.rho_water) &&
    finite_scalar(m.rho_energy) &&
    finite_scalar(m.carbon) &&
    valid_oxygen(m.oxygen)
