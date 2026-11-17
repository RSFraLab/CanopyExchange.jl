# ============================================================================
# Finite-volume flux divergence and steady (zero-storage) face fluxes
# ============================================================================

"""
    conserved_tendency(bottom, top, source, geometry)

Finite-volume tendency `(bottom - top) / thickness + source` of one conserved scalar.
`bottom` and `top` are upward-positive face fluxes per m² ground area and second;
`source` and the result are per m³ bulk canopy volume and second. Summing
`thickness * tendency` over stacked layers cancels every interior face flux.

# Examples
```jldoctest
julia> conserved_tendency(15.0, 12.0, 4.0, LayerGeometry(1.5, 0.8))
6.0
```
"""
@inline function conserved_tendency(bottom, top, source, g)
    same_precision(bottom, top, source, thickness(g))
    return (bottom - top) / thickness(g) + source
end

# Componentwise `conserved_tendency` of one layer (a functor, not a closure, so it
# stays a concrete isbits callable in device kernels).
struct LayerDivergence{G}
    geometry::G
end
@inline (d::LayerDivergence)(bottom, top, source) =
    conserved_tendency(bottom, top, source, d.geometry)

"""
    layer_tendency(bottom::FaceFlux, top::FaceFlux, source::VolumeSource, geometry) -> VolumeSource

Bulk-volume tendency `(bottom - top) / thickness + source` of every component
([`conserved_tendency`](@ref) applied componentwise). Face fluxes are upward positive
per m² ground area; all records must carry the same tracers and precision. Summing
`thickness * tendency` over adjacent layers telescopes the interior faces. Apply any
datum shift to the individual streams before this call
([`shift_water_energy_datum`](@ref)).

# Examples
```jldoctest
julia> b, t = face_flux(30.0, 0.0), face_flux(24.0, 0.0);

julia> layer_tendency(b, t, VolumeSource(ExchangeRates(10.0, 0.0)), LayerGeometry(2.0, 0.5)).rates.energy
13.0
```
"""
@inline function layer_tendency(bottom::FaceFlux, top::FaceFlux, source::VolumeSource, g)
    rates = map(LayerDivergence(g), bottom.rates, top.rates, source.rates)
    return VolumeSource(rates)
end

"""
    layer_flux_increment(source::VolumeSource, geometry) -> FaceFlux

Integrate a bulk-volume source over the layer thickness: `thickness * source`, a face
flux increment per m² ground area and second.

# Examples
```jldoctest
julia> s = VolumeSource(ExchangeRates(40.0, 1.0e-5));

julia> layer_flux_increment(s, LayerGeometry(2.5, 0.8)).rates.energy
100.0
```
"""
@inline layer_flux_increment(source::VolumeSource, g) = FaceFlux(thickness(g) * source.rates)

"""
    steady_top_flux(bottom::FaceFlux, source::VolumeSource, geometry) -> FaceFlux

Upward face flux at the layer top under zero air storage:
`bottom + thickness * source`, so that
`layer_tendency(bottom, steady_top_flux(bottom, source, g), source, g)` is zero to
round-off. Bridges diagnostic zero-storage canopy models (ClimaLand big-leaf or
diagnostic multilayer turbulence) and prognostic finite-volume columns.

# Examples
```jldoctest
julia> b, s, g = face_flux(20.0, 0.0), VolumeSource(ExchangeRates(15.0, 0.0)), LayerGeometry(2.0, 0.6);

julia> t = steady_top_flux(b, s, g);

julia> (t.rates.energy, layer_tendency(b, t, s, g).rates.energy)
(50.0, 0.0)
```
"""
@inline steady_top_flux(bottom::FaceFlux, source::VolumeSource, g) =
    bottom + layer_flux_increment(source, g)

"""
    accumulate_steady_fluxes!(faces, bottom::FaceFlux, sources, geometries) -> faces

Overwrite `faces` (`N + 1` elements) from `bottom` upward across `N` layers with
[`steady_top_flux`](@ref) and return `faces`. `faces[end]` is the emergent canopy-top
flux `bottom + Σ_k thickness_k * sources_k`. `sources` and `geometries` must have
identical first axes (`OffsetArray` and views are supported) and
`length(faces) == length(sources) + 1`; otherwise a `DimensionMismatch` is thrown.
Values are converted to `eltype(faces)`.

# Examples
```jldoctest
julia> b = face_flux(10.0, 0.0);

julia> faces = accumulate_steady_fluxes!(Vector{typeof(b)}(undef, 2), b, [VolumeSource(ExchangeRates(6.0, 0.0))], [LayerGeometry(2.0, 0.5)]);

julia> [f.rates.energy for f in faces]
2-element Vector{Float64}:
 10.0
 22.0
```

See also [`accumulate_steady_fluxes`](@ref).
"""
function accumulate_steady_fluxes!(
    faces::AbstractVector{<:FaceFlux},
    bottom::FaceFlux,
    sources::AbstractVector{<:VolumeSource},
    geometries::AbstractVector,
)
    axes(sources, 1) == axes(geometries, 1) || throw(
        DimensionMismatch("sources and geometries must share identical first-dimension axes"),
    )
    length(faces) == length(sources) + 1 ||
        throw(DimensionMismatch("faces length must equal length(sources) + 1"))
    FaceOut = eltype(faces)
    index = firstindex(faces)
    current = convert(FaceOut, bottom)
    @inbounds faces[index] = current
    for k in eachindex(sources, geometries)
        index += 1
        current = convert(FaceOut, steady_top_flux(current, sources[k], geometries[k]))
        @inbounds faces[index] = current
    end
    return faces
end

"""
    accumulate_steady_fluxes(bottom::FaceFlux, sources, geometries) -> Vector{<:FaceFlux}

Allocate and fill a 1-based `N + 1` steady face-flux profile from `bottom` to the
canopy top (see [`accumulate_steady_fluxes!`](@ref) for the contract). The element
type promotes `bottom`, the sources and the geometries, so a dual-valued source or
geometry gives dual-valued faces. Pass a preallocated vector (for example an
`OffsetVector`) to [`accumulate_steady_fluxes!`](@ref) for other face index axes.

# Examples
```jldoctest
julia> s = [VolumeSource(ExchangeRates(5.0, 0.0)), VolumeSource(ExchangeRates(8.0, 0.0))];

julia> faces = accumulate_steady_fluxes(face_flux(10.0, 0.0), s, [LayerGeometry(2.0, 0.5), LayerGeometry(1.5, 0.4)]);

julia> [f.rates.energy for f in faces]
3-element Vector{Float64}:
 10.0
 20.0
 32.0
```
"""
function accumulate_steady_fluxes(
    bottom::FaceFlux,
    sources::AbstractVector{<:VolumeSource},
    geometries::AbstractVector,
)
    axes(sources, 1) == axes(geometries, 1) || throw(
        DimensionMismatch("sources and geometries must share identical first-dimension axes"),
    )
    FB, FS, FG = eltype(bottom), eltype(eltype(sources)), eltype(eltype(geometries))
    same_precision(zero(FB), zero(FS), zero(FG))
    FT = promote_type(FB, FS, FG)
    R = ExchangeRates{FT,_tracer_type(typeof(bottom.rates.tracers), FT)}
    faces = Vector{FaceFlux{R}}(undef, length(sources) + 1)
    return accumulate_steady_fluxes!(faces, bottom, sources, geometries)
end
