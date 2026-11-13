"""
    surface_static_energy_flux(shf, lhf, evaporation, carbon_flux=zero(FT), oxygen=NoOxygen(FT))

Construct an upward-positive [`FaceFlux`](@ref) from `SurfaceFluxes.jl` static-energy outputs.

In CliMA's `SurfaceFluxes.jl` convention (`SHF = -ρ_sfc g_h ΔDSE + VSE_sfc E` and
`LHF = LH_v0 E` with `LH_v0 + VSE_sfc == Thermodynamics.enthalpy_vapor(T_sfc) + Φ_sfc`),
`shf + lhf` equals the total moist static enthalpy plus surface-geopotential flux
`-ρ_sfc g_h ΔDSE + E * (h_v(T_sfc) + Φ_sfc)` relative to `Thermodynamics.jl`'s
zero-energy datum, where `evaporation` [kg H₂O m⁻² s⁻¹] is the kinematic water mass
flux. Conventional eddy-covariance `c_p ΔT + L_v(T_air) E` does **not** satisfy this
static-energy datum without an explicit shift via [`shift_water_energy_datum`](@ref).

# Examples
```jldoctest
julia> f = surface_static_energy_flux(60.0, 140.0, 5.6e-5, -1.5e-5);

julia> (f.rates.energy, f.rates.water, f.rates.carbon)
(200.0, 5.6e-5, -1.5e-5)
```

See also [`shift_water_energy_datum`](@ref), [`sensible_water_flux`](@ref), [`FaceFlux`](@ref).
"""
@inline function surface_static_energy_flux(
    shf::Real,
    lhf::Real,
    evaporation::Real,
    carbon_flux::Real,
    oxygen::O,
) where {FO<:Real,O<:OxygenContent{FO}}
    FT = promote_compatible_scalar(
        surface_static_energy_flux,
        shf,
        lhf,
        evaporation,
        carbon_flux,
        zero(FO),
    )
    OT = _promote_oxygen_type(O, FT)
    s = convert(FT, shf)
    l = convert(FT, lhf)
    e = convert(FT, evaporation)
    c = convert(FT, carbon_flux)
    return FaceFlux(ExchangeRates(s + l, e, c, convert(OT, oxygen)))
end

@inline function surface_static_energy_flux(
    shf::Real,
    lhf::Real,
    evaporation::Real,
    carbon_flux::Real,
)
    FT = promote_compatible_scalar(
        surface_static_energy_flux,
        shf,
        lhf,
        evaporation,
        carbon_flux,
    )
    return surface_static_energy_flux(shf, lhf, evaporation, carbon_flux, NoOxygen{FT}())
end

@inline function surface_static_energy_flux(shf::Real, lhf::Real, evaporation::Real)
    FT = promote_compatible_scalar(surface_static_energy_flux, shf, lhf, evaporation)
    return surface_static_energy_flux(shf, lhf, evaporation, zero(FT), NoOxygen{FT}())
end

"""
    shift_water_energy_datum(stream, delta_specific_energy::Real)

Shift the energy component of `stream` by `delta_specific_energy * water` [J kg⁻¹ H₂O].

Apply reference-enthalpy offsets (`h_v,target - h_v,source`) or geopotential
offsets (`Φ_target - Φ_source`) directly to individual [`ExchangeRates`](@ref),
[`LeafFlux`](@ref), [`StemFlux`](@ref), [`FaceFlux`](@ref), or
[`VolumeSource`](@ref) streams **before** [`layer_tendency`](@ref) so interior
face fluxes remain identical on shared interfaces and telescope conservatively
across multi-layer columns. Purely positional arguments support `@.` broadcasting
over `ClimaCore.Fields.Field` and `CuArray`.

# Examples
```jldoctest
julia> f = FaceFlux(ExchangeRates(100.0, 2.0e-4, 0.0));

julia> shift_water_energy_datum(f, 5.0e4).rates.energy
110.0
```
"""
@inline function shift_water_energy_datum(
    rates::ExchangeRates{FR,OR},
    delta_specific_energy::Real,
) where {FR<:Real,OR<:OxygenContent{FR}}
    FT = promote_compatible_scalar(shift_water_energy_datum, zero(FR), delta_specific_energy)
    OT = _promote_oxygen_type(OR, FT)
    rr = convert(ExchangeRates{FT,OT}, rates)
    de = convert(FT, delta_specific_energy)
    return ExchangeRates{FT,OT}(
        rr.energy + de * rr.water,
        rr.water,
        rr.carbon,
        rr.oxygen,
    )
end

for FluxT in (:LeafFlux, :StemFlux, :FaceFlux, :VolumeSource)
    @eval begin
        @inline function shift_water_energy_datum(
            stream::$FluxT,
            delta_specific_energy::Real,
        )
            return $FluxT(shift_water_energy_datum(stream.rates, delta_specific_energy))
        end
    end
end

"""
    layer_flux_increment(source::VolumeSource, geometry::AbstractLayerGeometry) -> FaceFlux

Integrate bulk-volume `source` across `geometry.thickness` to a [`FaceFlux`](@ref) increment.

# Examples
```jldoctest
julia> s = VolumeSource(ExchangeRates(40.0, 1.0e-5, -0.5e-5));

julia> layer_flux_increment(s, LayerGeometry(2.5, 0.8)).rates.energy
100.0
```
"""
@inline function layer_flux_increment(
    source::VolumeSource{FS,O},
    geometry::AbstractLayerGeometry{FG},
) where {FS<:Real,O<:OxygenContent{FS},FG<:Real}
    FT = promote_compatible_scalar(layer_flux_increment, zero(FS), zero(FG))
    OT = _promote_oxygen_type(O, FT)
    dz = convert(FT, geometry.thickness)
    sr = convert(ExchangeRates{FT,OT}, source.rates)
    return FaceFlux(dz * sr)
end

"""
    steady_top_flux(bottom::FaceFlux, source::VolumeSource, geometry::AbstractLayerGeometry) -> FaceFlux

Return the upward [`FaceFlux`](@ref) at the layer top under zero air-storage tendency.

Satisfies `layer_tendency(bottom, steady_top_flux(bottom, source, geometry), source, geometry) ≈ zero(source)`,
bridging diagnostic zero-storage canopy models (`ClimaLand` big-leaf / diagnostic
multi-layer turbulence) with prognostic finite-volume columns.

# Examples
```jldoctest
julia> b = FaceFlux(ExchangeRates(20.0, 0.0, 0.0));

julia> s = VolumeSource(ExchangeRates(15.0, 0.0, 0.0));

julia> g = LayerGeometry(2.0, 0.6);

julia> t = steady_top_flux(b, s, g);

julia> (t.rates.energy, layer_tendency(b, t, s, g).rates.energy)
(50.0, 0.0)
```
"""
@inline function steady_top_flux(
    bottom::FaceFlux{FB,OB},
    source::VolumeSource{FS,OS},
    geometry::AbstractLayerGeometry{FG},
) where {FB<:Real,OB<:OxygenContent{FB},FS<:Real,OS<:OxygenContent{FS},FG<:Real}
    compatible_oxygen(OB, OS) ||
        throw(MethodError(steady_top_flux, (bottom, source, geometry)))
    FT = promote_compatible_scalar(steady_top_flux, zero(FB), zero(FS), zero(FG))
    OT = _promote_oxygen_type(OB, FT)
    bb = convert(FaceFlux{FT,OT}, bottom)
    inc = layer_flux_increment(convert(VolumeSource{FT,OT}, source), geometry)
    return bb + inc
end

"""
    accumulate_steady_fluxes!(faces, bottom, sources, geometries)

Overwrite `faces` (`N + 1` elements) in place from `bottom` upward across `N` layers and return `faces`.

`sources` and `geometries` must share identical `axes(sources, 1) == axes(geometries, 1)`
(supporting `OffsetArray` and `@view`), and `faces` must satisfy
`length(faces) == length(sources) + 1`, throwing `DimensionMismatch` otherwise.
At completion, `faces[end]` holds the emergent canopy-top [`FaceFlux`](@ref)
equal to `bottom + sum_k thickness_k * sources_k`.

# Examples
```jldoctest
julia> b = FaceFlux(ExchangeRates(10.0, 0.0, 0.0));

julia> s = [VolumeSource(ExchangeRates(6.0, 0.0, 0.0))];

julia> g = [LayerGeometry(2.0, 0.5)];

julia> faces = Vector{typeof(b)}(undef, 2);

julia> accumulate_steady_fluxes!(faces, b, s, g);

julia> [f.rates.energy for f in faces]
2-element Vector{Float64}:
 10.0
 22.0
```

See also [`accumulate_steady_fluxes`](@ref), [`steady_top_flux`](@ref), [`layer_flux_increment`](@ref).
"""
function accumulate_steady_fluxes!(
    faces::AbstractVector{<:FaceFlux},
    bottom::FaceFlux,
    sources::AbstractVector{<:VolumeSource},
    geometries::AbstractVector{<:AbstractLayerGeometry},
)
    axes(sources, 1) == axes(geometries, 1) || throw(
        DimensionMismatch("sources and geometries must share identical first-dimension axes"),
    )
    length(faces) == length(sources) + 1 || throw(
        DimensionMismatch("faces length must equal length(sources) + 1"),
    )
    FaceOut = eltype(faces)
    fi = firstindex(faces)
    curr = convert(FaceOut, bottom)
    @inbounds faces[fi] = curr
    for idx in eachindex(sources, geometries)
        fi += 1
        curr = convert(FaceOut, steady_top_flux(curr, sources[idx], geometries[idx]))
        @inbounds faces[fi] = curr
    end
    return faces
end

"""
    accumulate_steady_fluxes(bottom::FaceFlux, sources, geometries) -> Vector{<:FaceFlux}

Allocate and return a 1-based `N + 1` steady-state [`FaceFlux`](@ref) profile from `bottom` to canopy top.

Pass a preallocated staggered vector (including `OffsetVector`) to
[`accumulate_steady_fluxes!`](@ref) when custom face index axes are required.

# Examples
```jldoctest
julia> b = FaceFlux(ExchangeRates(10.0, 0.0, 0.0));

julia> s = [VolumeSource(ExchangeRates(5.0, 0.0, 0.0)), VolumeSource(ExchangeRates(8.0, 0.0, 0.0))];

julia> g = [LayerGeometry(2.0, 0.5), LayerGeometry(1.5, 0.4)];

julia> faces = accumulate_steady_fluxes(b, s, g);

julia> [f.rates.energy for f in faces]
3-element Vector{Float64}:
 10.0
 20.0
 32.0
```

See also [`accumulate_steady_fluxes!`](@ref), [`steady_top_flux`](@ref).
"""
function accumulate_steady_fluxes(
    bottom::FaceFlux{FB,OB},
    sources::AbstractVector{<:VolumeSource{FS,OS}},
    geometries::AbstractVector{<:AbstractLayerGeometry{FG}},
) where {FB<:Real,OB<:OxygenContent{FB},FS<:Real,OS<:OxygenContent{FS},FG<:Real}
    axes(sources, 1) == axes(geometries, 1) || throw(
        DimensionMismatch("sources and geometries must share identical first-dimension axes"),
    )
    compatible_oxygen(OB, OS) ||
        throw(MethodError(accumulate_steady_fluxes, (bottom, sources, geometries)))
    FT = promote_compatible_scalar(
        accumulate_steady_fluxes,
        zero(FB),
        zero(FS),
        zero(FG),
    )
    OT = _promote_oxygen_type(OB, FT)
    faces = Vector{FaceFlux{FT,OT}}(undef, length(sources) + 1)
    return accumulate_steady_fluxes!(faces, bottom, sources, geometries)
end

"""
    moist_air_tendencies(tendency::VolumeSource, gas_fraction::Real=one(eltype(tendency))) -> MoistAirTendencies

Convert a bulk-volume [`VolumeSource`](@ref) into gas-volume [`MoistAirTendencies`](@ref) by dividing by `gas_fraction`.

Contracts:
- `gas_fraction` is the time-invariant air porosity `ε_gas ∈ (0, 1]` (`∂_t ε_gas = 0`);
  throw `ArgumentError` if `gas_fraction ≤ 0`, `gas_fraction > 1`, or nonfinite.
- Dry-air mass is conserved (`S_ρd = 0`), so `rho_total == rho_water == water / gas_fraction`.
- Fieldwise division `/ gas_fraction` preserves exact identity at `gas_fraction == 1`
  and zero-ULP molar `O₂ / CO₂` stoichiometry.
- Apply any reference-enthalpy or geopotential shift to individual streams via
  [`shift_water_energy_datum`](@ref) **prior** to [`layer_tendency`](@ref).
- Purely positional arguments support `@.` broadcasting over `ClimaCore.Fields.Field`
  and `CuArray`.

# Examples
```jldoctest
julia> v = VolumeSource(ExchangeRates(180.0, 3.6e-5, -1.8e-5));

julia> m = moist_air_tendencies(v, 0.9);

julia> (m.rho_total, m.rho_water, m.rho_energy, m.carbon)
(4.0e-5, 4.0e-5, 200.0, -2.0e-5)
```
"""
@inline function moist_air_tendencies(
    tendency::VolumeSource{FV,OV},
    gas_fraction::Real = one(FV),
) where {FV<:Real,OV<:OxygenContent{FV}}
    FT = promote_compatible_scalar(moist_air_tendencies, zero(FV), gas_fraction)
    OT = _promote_oxygen_type(OV, FT)
    eps_gas = convert(FT, gas_fraction)
    P = primal_float_type(FT)
    finite_scalar(eps_gas) &&
        (zero(P) < primal_scalar(eps_gas)) &&
        (primal_scalar(eps_gas) <= one(P)) ||
        throw(ArgumentError("gas_fraction must be finite and in (0, 1]"))
    r = convert(ExchangeRates{FT,OT}, tendency.rates)
    w = r.water / eps_gas
    return MoistAirTendencies{FT,OT}(
        w,
        w,
        r.energy / eps_gas,
        r.carbon / eps_gas,
        _div_oxygen(r.oxygen, eps_gas),
    )
end
