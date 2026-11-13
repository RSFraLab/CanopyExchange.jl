# Changelog

## v0.1.1

### Added
- Continuous-scalar `{FT<:Real}` support across all exchange records, wrappers, and
  layer geometries (`Float32`, `Float64`, `BigFloat`, and 1st/2nd-order
  `ForwardDiff.Dual`), backed by zero-allocation `@generated` structural predicates
  (`CanopyExchange.finite_scalar`, `CanopyExchange.exact_zero_scalar`) that reject
  nonfinite derivative partials while preserving empty-layer `NaN` guards and
  right-directional AD tangents.
- Strict woody/stem measure separation via `StemFlux`, `PlantLayerGeometry`
  (`PAI = LAI + SAI`), `area_source`, `stem_area_source`, `plant_area_density`,
  `sensible_stem_flux`, and conservative `(air, leaf, stem)` donor ledgers in
  `plant_air_sources`.
- CliMA boundary, datum-shift, steady-column, and porous moist-air contracts
  (`surface_static_energy_flux`, `shift_water_energy_datum`, `layer_flux_increment`,
  `steady_top_flux`, `accumulate_steady_fluxes!`, `accumulate_steady_fluxes`,
  `MoistAirTendencies`, and `moist_air_tendencies`), supporting `@.` broadcasting
  over `ClimaCore.Fields.Field` and `CuArray` as well as `OffsetArray` and `@view`
  slices.
- Complete `Base` interfaces (`eltype`, `convert`, `promote_rule`, `==`, `isequal`,
  `hash`, `zero`, `broadcastable`, and measure-safe linear algebra `+`, `-`, unary `-`,
  scalar `*`, scalar `/`).
- Structured `Documenter.jl` site and domain-organized API reference (`docs/make.jl`,
  `docs/src/pages/`, `docs/src/api/`) with `checkdocs = :all` and doctests, plus
  `Aqua.jl` and `ExplicitImports.jl` quality-assurance gates.

### Preserved
- Zero runtime `[deps]` in `Project.toml` (`[compat] julia = "1.10"`) and exact
  IEEE-754 primal `Float32`/`Float64` bit parity (`SHA-256 19bcfb73...`) with `v0.1.0`.

## v0.1.0

### Added
- Initial extraction of conservative scalar canopy exchange contracts (`OxygenContent`,
  `NoOxygen`, `OxygenTransfer`, `oxygen_value`, `ExchangeRates`, `LeafFlux`,
  `FaceFlux`, `VolumeSource`, `LayerGeometry`, `leaf_area_source`,
  `sensible_water_flux`, `leaf_air_sources`, `conserved_tendency`, `layer_tendency`,
  `valid_exchange`) from `CanopyColumn.ConservativeExchange` (`05723a9`) with CPU and
  CUDA broadcast verification.
