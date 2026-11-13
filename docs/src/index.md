# CanopyExchange.jl

*Zero-dependency, measure-safe scalar exchange contracts for layered plant canopies and atmospheric or land host models.*

`CanopyExchange` owns spatial measure wrappers (`LeafFlux`, `StemFlux`, `FaceFlux`,
`VolumeSource`, `MoistAirTendencies`) and conservative transfer arithmetic connecting
leaf/stem physiology to canopy-air and atmospheric column budgets. It carries zero
runtime dependencies, runs allocation-free across `Float32`, `Float64`, `CuArray`,
and `ClimaCore.Fields.Field` broadcasting, and propagates `ForwardDiff.Dual`
derivatives cleanly.

## Spatial Measures & Sign Conventions

| Type | Spatial Denominator | Positive Direction | Primary Constructors & Producers |
| --- | --- | --- | --- |
| [`LeafFlux`](@ref) | $\mathrm{m^2}$ one-sided leaf area (`LAD`) | Leaf $\to$ air | [`LeafFlux`](@ref), [`sensible_water_flux`](@ref) |
| [`StemFlux`](@ref) | $\mathrm{m^2}$ one-sided stem area (`SAD`) | Stem $\to$ air | [`StemFlux`](@ref), [`sensible_stem_flux`](@ref) |
| [`FaceFlux`](@ref) | $\mathrm{m^2}$ horizontal ground area | Upward ($+z$) | [`FaceFlux`](@ref), [`surface_static_energy_flux`](@ref), [`steady_top_flux`](@ref) |
| [`VolumeSource`](@ref) | $\mathrm{m^3}$ bulk canopy volume | Into named ledger | [`leaf_air_sources`](@ref), [`plant_air_sources`](@ref), [`layer_tendency`](@ref) |
| [`MoistAirTendencies`](@ref) | $\mathrm{m^3}$ gas volume ($\varepsilon_{\mathrm{gas}} \in (0, 1]$) | Into canopy air | [`moist_air_tendencies`](@ref) |

Distinct concrete types forbid accidental cross-measure addition (`LeafFlux + StemFlux`
or `LeafFlux + FaceFlux` throws `MethodError`) and prevent swapping leaf and stem
arguments in [`plant_air_sources`](@ref).

## Quickstart

```jldoctest
julia> using CanopyExchange

julia> geom = PlantLayerGeometry(2.0, 1.2, 0.25); # thickness [m], LAD [m² m⁻³], SAD [m² m⁻³]

julia> lf = sensible_water_flux(95.0, 4.0e-5, -1.4e-5, 2.547e6); # H, E, CO₂, h_v(T_leaf)

julia> sf = sensible_stem_flux(20.0); # dry woody sensible heat [W m⁻² stem]

julia> srcs = plant_air_sources(lf, sf, geom); # conservative (; air, leaf, stem) ledgers

julia> srcs.air + (srcs.leaf + srcs.stem) == zero(srcs.air)
true

julia> soil = surface_static_energy_flux(45.0, 90.0, 3.6e-5, 2.0e-6); # SurfaceFluxes boundary

julia> top = steady_top_flux(soil, srcs.air, geom); # zero-storage diagnostic top FaceFlux

julia> tend = layer_tendency(soil, top, srcs.air, geom); # (bottom - top)/dz + air

julia> moist = moist_air_tendencies(tend, 0.96); # per m³ gas volume (ε_gas = 0.96)

julia> abs(moist.rho_energy) <= 1.0e-12
true
```

## Documentation & API Reference

| Topic | Conceptual Guides | Structured API Reference |
| --- | --- | --- |
| **Measures, Donor Ledgers & Finite-Volume Conservation** | [Measures & Conservation](pages/measures_and_conservation.md) | [Quantities & Layer Geometries](api/quantities_and_geometries.md), [Area Sources & Divergence](api/sources_and_tendencies.md) |
| **CliMA (`SurfaceFluxes`, `Thermodynamics`, `ClimaLand`, `ClimaAtmos`, `ClimaCore`)** | [CliMA Interoperability](pages/clima_interoperability.md) | [CliMA Contracts & Scalar Traits](api/clima_and_scalars.md) |
| **Extraction Bit Parity, AD & GPU Verification** | [Validation & Bit Parity](pages/validation.md) | [Full Module & Symbol Index](api/index.md) |
