# API Reference

This section documents all exported and `public` types and functions in
`CanopyExchange`, grouped by physical domain rather than alphabetical order.

```@docs
CanopyExchange
```

## Subsystem Map

| Physical Domain | API Reference Page | Conceptual Guides |
| --- | --- | --- |
| Optional oxygen modes, conserved `ExchangeRates`, measure-tagged wrappers (`LeafFlux`, `StemFlux`, `FaceFlux`, `VolumeSource`), `MoistAirTendencies`, layer geometries (`LayerGeometry`, `PlantLayerGeometry`), and finite-value predicates | [Quantities & Layer Geometries](quantities_and_geometries.md) | [Measures & Conservation](../pages/measures_and_conservation.md) |
| Area-to-volume source conversions (`area_source`, `leaf_area_source`, `stem_area_source`, `plant_area_density`), enthalpy-coupled flux constructors (`sensible_water_flux`, `sensible_stem_flux`), paired donor ledgers (`leaf_air_sources`, `plant_air_sources`), and finite-volume divergences (`conserved_tendency`, `layer_tendency`) | [Area Sources & Divergence](sources_and_tendencies.md) | [Measures & Conservation](../pages/measures_and_conservation.md) |
| `SurfaceFluxes` boundary assembly (`surface_static_energy_flux`), water-energy datum shifts (`shift_water_energy_datum`), diagnostic steady-air flux profiles (`layer_flux_increment`, `steady_top_flux`, `accumulate_steady_fluxes!`, `accumulate_steady_fluxes`), porous moist-air density tendencies (`moist_air_tendencies`), and `public` continuous-scalar traits | [CliMA Contracts & Scalar Traits](clima_and_scalars.md) | [CliMA Interoperability](../pages/clima_interoperability.md), [Validation & Bit Parity](../pages/validation.md) |

## Conventions

- **Exported vs. `public` bindings**: Symbols shown unqualified (`LeafFlux`,
  `plant_air_sources`, `moist_air_tendencies`) are exported by `using CanopyExchange`.
  Continuous-scalar trait and leaf-float predicate functions declared `public` on
  Julia 1.11+ (`CanopyExchange.primal_float_type`, `CanopyExchange.finite_scalar`,
  `CanopyExchange.exact_zero_scalar`, etc.) are part of the supported API when
  qualified with `CanopyExchange.` or explicitly imported.
- **Precision, device & AD contract**: All structs and functions preserve caller
  floating-point precision (`Float32`, `Float64`, `BigFloat`), allocate zero heap
  memory on scalar paths, broadcast via `@.` over `CuArray` and `ClimaCore.Fields.Field`,
  and propagate `ForwardDiff.Dual` tangents while rejecting mixed `Float32`/`Float64`
  precisions and discrete integer/rational types.

## Section Contents

```@contents
Pages = [
    "quantities_and_geometries.md",
    "sources_and_tendencies.md",
    "clima_and_scalars.md",
]
Depth = 2
```

## Alphabetical Symbol Index

```@index
Pages = [
    "index.md",
    "quantities_and_geometries.md",
    "sources_and_tendencies.md",
    "clima_and_scalars.md",
]
```
