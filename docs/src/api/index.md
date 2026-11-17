# API Reference

This section documents all exported and `public` names of `CanopyExchange`, grouped by
physical process.

```@docs
CanopyExchange
```

## Subsystem map

| Process | API page | Concept guide |
| --- | --- | --- |
| Rates and passive tracers (`ExchangeRates`, `Tracers`), measures and measure-tagged records (`Exchange`, `LeafFlux`, `StemFlux`, `FaceFlux`, `VolumeSource`), layer geometry interface (`LayerGeometry`, `PlantLayerGeometry`, `thickness`, `area_density`) | [Records, Measures & Layer Geometries](quantities_and_geometries.md) | [Measures & Conservation](../pages/measures_and_conservation.md) |
| Area-to-volume conversion (`area_source`, `volume_source`), paired ledgers (`leaf_air_sources`, `plant_air_sources`), finite-volume divergence (`conserved_tendency`, `layer_tendency`) | [Area Sources & Divergence](sources_and_tendencies.md) | [Measures & Conservation](../pages/measures_and_conservation.md) |
| Element and boundary fluxes (`sensible_water_flux`, `sensible_stem_flux`, `surface_static_energy_flux`), energy datum (`shift_water_energy_datum`), steady columns, `MoistAirTendencies`, scalar hooks | [Boundary Fluxes, Steady Columns, Moist Air & Scalar Hooks](clima_and_scalars.md) | [CliMA Interoperability](../pages/clima_interoperability.md), [Validation](../pages/validation.md) |

## Conventions

- **Exported and `public` names**: names shown unqualified are exported by
  `using CanopyExchange`. The scalar hooks (`CanopyExchange.primal_type`,
  `CanopyExchange.primal`, `CanopyExchange.all_finite`, `CanopyExchange.exact_zero`)
  are `public` on Julia 1.11+ and not exported.
- **Precision**: every scalar that meets in a record or kernel carries one
  floating-point type. `Float32` with `Float64`, integers, rationals and irrationals
  throw an `ArgumentError`; `Float64` with a `Float64`-valued dual promotes.
- **Kernels**: positional arguments only, no keyword arguments, static error messages,
  zero heap allocations; they broadcast over `CuArray` and `ClimaCore.Fields.Field`.
  Construction inside such broadcasts uses the function forms (`leaf_flux`,
  `face_flux`, `layer_geometry`, `Tracers(:co2)`), not type constructors.

## Section contents

```@contents
Pages = [
    "quantities_and_geometries.md",
    "sources_and_tendencies.md",
    "clima_and_scalars.md",
]
Depth = 2
```

## Alphabetical index

```@index
Pages = [
    "index.md",
    "quantities_and_geometries.md",
    "sources_and_tendencies.md",
    "clima_and_scalars.md",
]
```
