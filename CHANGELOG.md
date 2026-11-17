# Changelog

## Unreleased (0.2.0)

The version in `Project.toml` is not yet bumped.

### Breaking
- `ExchangeRates(energy, water, tracers::NamedTuple)` replaces the fixed `carbon` field
  and the optional oxygen: CO₂, O₂ and any further passive molar tracer are keys of
  `tracers`. An uncomputed tracer is not a key; records with different tracer sets
  (or tracer order) cannot be combined.
- Removed: `OxygenContent`, `NoOxygen`, `OxygenTransfer`, `oxygen_value`,
  `compatible_oxygen`, the `carbon` and `oxygen` fields, `AbstractExchangeFlux`,
  `AbstractLayerGeometry`, and the scalar reflection helpers `primal_float_type`,
  `primal_scalar`, `is_continuous_scalar`, `compatible_scalar`, `finite_scalar`,
  `exact_zero_scalar`.
- `LeafFlux`, `StemFlux`, `FaceFlux`, `VolumeSource` are aliases of
  `Exchange{LeafArea}`, `Exchange{StemArea}`, `Exchange{GroundArea}`,
  `Exchange{BulkVolume}`; their type parameters are now `{M,R}` instead of `{FT,O}`.
- New argument order: `sensible_water_flux(H, E, h_v, tracers)` and
  `sensible_stem_flux(H, E, h_v, tracers)` (vapour enthalpy before the tracers);
  `surface_static_energy_flux(shf, lhf, E, tracers)`.
- Mixed floating-point precision, integer, rational and irrational inputs, and
  mismatched tracer sets throw `ArgumentError` (previously `MethodError`).
- Keyword constructors take `tracers` instead of `carbon`/`oxygen`.

### Added
- Open passive tracer set (`Tracers(:co2, :o2)` callable for broadcasts), the
  `Measure` types and `Exchange{M}`, `volume_source`, duck-typed geometry interface
  (`thickness`, `area_density`), function-form constructors (`leaf_flux`, `stem_flux`,
  `face_flux`, `layer_geometry`) for GPU field broadcasts.
- `public` scalar hooks `CanopyExchange.primal_type`, `primal`, `all_finite`,
  `exact_zero`, and the package extension `CanopyExchangeForwardDiffExt` (ForwardDiff
  0.10 and 1.x as a weak dependency; runtime `[deps]` stay empty).
- Frozen v0.1.1 bit record (`test/reference/exchange_v0_1_1.toml`) of the leaf, stem,
  datum-shift, steady-column, surface and moist-air outputs; the new core reproduces
  it bit for bit. Tests for the precision policy, a second tracer set, a third element
  class defined outside the package, `@inferred` and zero allocations of every
  exported kernel; opt-in `test/run_fields.jl` (ClimaCore fields on CPU and CUDA).

### Migration

| v0.1.x | 0.2.0 |
| --- | --- |
| `ExchangeRates(e, w, c)` | `ExchangeRates(e, w, (; co2 = c))` |
| `ExchangeRates(e, w, c, OxygenTransfer(o))` | `ExchangeRates(e, w, (; co2 = c, o2 = o))` |
| `r.carbon`, `oxygen_value(r.oxygen)` | `r.tracers.co2`, `r.tracers.o2` (absent: `!haskey(r.tracers, :o2)`) |
| `NoOxygen{FT}()` | omit the `o2` key |
| `sensible_water_flux(H, E, c, h_v)` | `sensible_water_flux(H, E, h_v, (; co2 = c))` |
| `sensible_water_flux(H, E, c, OxygenTransfer(o), h_v)` | `sensible_water_flux(H, E, h_v, (; co2 = c, o2 = o))` |
| `sensible_stem_flux(H)` | `sensible_stem_flux(H, (; co2 = zero(H)))` (or no tracers) |
| `sensible_stem_flux(H, E, c, h_v)` | `sensible_stem_flux(H, E, h_v, (; co2 = c))` |
| `surface_static_energy_flux(shf, lhf, E, c)` | `surface_static_energy_flux(shf, lhf, E, (; co2 = c))` |
| `surface_static_energy_flux(shf, lhf, E)` (carbon 0) | `surface_static_energy_flux(shf, lhf, E, (; co2 = zero(E)))` |
| `MoistAirTendencies(...).carbon` | `.tracers.co2` |
| `f isa AbstractExchangeFlux` | `f isa Exchange` |
| `g isa AbstractLayerGeometry` | any type with `thickness` and `area_density` |
| `CanopyExchange.finite_scalar(x)` / `exact_zero_scalar(x)` | `CanopyExchange.all_finite(x)` / `CanopyExchange.exact_zero(x)` |
| `CanopyExchange.primal_float_type(T)` / `primal_scalar(x)` | `CanopyExchange.primal_type(T)` / `CanopyExchange.primal(x)` |
| `catch e::MethodError` for mixed precision or oxygen mode | `ArgumentError` |

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
