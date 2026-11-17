# CanopyExchange.jl

*Measure-safe, conservative scalar exchange contracts for layered plant canopies and atmospheric or land host models, with no runtime dependencies.*

`CanopyExchange` tags rates of energy, water and an open set of passive molar tracers
([`ExchangeRates`](@ref)) with the measure they are expressed per ([`Exchange`](@ref)),
and provides the conservative arithmetic that connects leaf and stem physiology to
canopy-air and atmospheric column budgets. Kernels are allocation-free in `Float32`
and `Float64` (tested with up to 12 tracers), broadcast over `CuArray` and `ClimaCore.Fields.Field`, and propagate
`ForwardDiff.Dual` derivatives (through a package extension).

## Spatial measures and sign conventions

| Record | Measure | Positive direction | Constructors and producers |
| --- | --- | --- | --- |
| [`LeafFlux`](@ref) | m² one-sided leaf area (`LAD`) | leaf → air | [`sensible_water_flux`](@ref), [`leaf_flux`](@ref) |
| [`StemFlux`](@ref) | m² one-sided stem area (`SAD`) | stem → air | [`sensible_stem_flux`](@ref), [`stem_flux`](@ref) |
| [`FaceFlux`](@ref) | m² horizontal ground area | upward (+z) | [`surface_static_energy_flux`](@ref), [`face_flux`](@ref), [`steady_top_flux`](@ref) |
| [`VolumeSource`](@ref) | m³ bulk canopy volume | into the named ledger | [`volume_source`](@ref), [`leaf_air_sources`](@ref), [`plant_air_sources`](@ref), [`layer_tendency`](@ref) |
| [`MoistAirTendencies`](@ref) | m³ gas volume (ε_gas ∈ (0, 1]) | into canopy air | [`moist_air_tendencies`](@ref) |

Records of different measures cannot be added or swapped (`MethodError`); records with
different tracer sets or floating-point precisions cannot be combined
(`ArgumentError`).

## Quickstart

```jldoctest
julia> using CanopyExchange

julia> geom = PlantLayerGeometry(2.0, 1.2, 0.25); # thickness [m], LAD [m² m⁻³], SAD [m² m⁻³]

julia> T = Tracers(:co2); # passive molar tracers carried by every record

julia> lf = sensible_water_flux(95.0, 4.0e-5, 2.547e6, T(-1.4e-5)); # H, E, h_v(T_leaf), CO₂

julia> sf = sensible_stem_flux(20.0, T(0.0)); # dry woody sensible heat [W m⁻² stem]

julia> srcs = plant_air_sources(lf, sf, geom); # conservative (; air, leaf, stem) ledgers

julia> srcs.air + (srcs.leaf + srcs.stem) == zero(srcs.air)
true

julia> soil = surface_static_energy_flux(45.0, 90.0, 3.6e-5, T(2.0e-6)); # SurfaceFluxes boundary

julia> top = steady_top_flux(soil, srcs.air, geom); # zero-storage diagnostic top FaceFlux

julia> tend = layer_tendency(soil, top, srcs.air, geom); # (bottom - top)/dz + air

julia> moist = moist_air_tendencies(tend, 0.96); # per m³ gas volume (ε_gas = 0.96)

julia> abs(moist.rho_energy) <= 1.0e-12
true
```

## Documentation & API Reference

| Topic | Conceptual Guides | Structured API Reference |
| --- | --- | --- |
| **Measures, Tracers, Donor Ledgers & Finite-Volume Conservation** | [Measures & Conservation](pages/measures_and_conservation.md) | [Records, Measures & Layer Geometries](api/quantities_and_geometries.md), [Area Sources & Divergence](api/sources_and_tendencies.md) |
| **CliMA (`SurfaceFluxes`, `Thermodynamics`, `ClimaLand`, `ClimaAtmos`, `ClimaCore`)** | [CliMA Interoperability](pages/clima_interoperability.md) | [Boundary Fluxes, Steady Columns, Moist Air & Scalar Hooks](api/clima_and_scalars.md) |
| **Frozen Bit Records, AD & GPU Verification** | [Validation & Bit Parity](pages/validation.md) | [Full Module & Symbol Index](api/index.md) |
