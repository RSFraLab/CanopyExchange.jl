# CanopyExchange.jl

Zero-dependency, measure-safe scalar exchange contracts for layered plant canopies
and atmospheric or land host models (`ClimaLand`, `SurfaceFluxes`, `ClimaAtmos`,
and standalone finite-volume columns). Originally extracted from
`CanopyColumn.ConservativeExchange` (`05723a9`) with Apache-2.0 attribution and
exact IEEE-754 primal bit parity.

## Source layout

- `src/scalars.jl`: continuous-scalar `{FT<:Real}` traits and zero-allocation
  `@generated` leaf-float predicates (`finite_scalar`, `exact_zero_scalar`) for
  `Float32`, `Float64`, `BigFloat`, and `ForwardDiff.Dual`.
- `src/quantities.jl`: measure-tagged wrappers (`LeafFlux`, `StemFlux`, `FaceFlux`,
  `VolumeSource`), conserved gas-volume `MoistAirTendencies`, layer geometries
  (`LayerGeometry`, `PlantLayerGeometry`), optional oxygen (`NoOxygen`,
  `OxygenTransfer`), and `Base` interfaces (`eltype`, `broadcastable`, `convert`,
  `==`, `isequal`, `hash`, `zero`, linear algebra).
- `src/sources.jl`: one-sided leaf (`LAD`) and stem (`SAD`) area-to-volume source
  conversions, paired donor ledgers (`leaf_air_sources`, `plant_air_sources`), and
  conservative finite-volume divergence (`conserved_tendency`, `layer_tendency`).
- `src/clima_contracts.jl`: CliMA boundary and host contracts
  (`surface_static_energy_flux`, `shift_water_energy_datum`, `layer_flux_increment`,
  `steady_top_flux`, `accumulate_steady_fluxes!`, `accumulate_steady_fluxes`,
  `moist_air_tendencies`).
- `test/`: conservation, migration bit-parity, `Aqua`, `ExplicitImports`,
  `Documenter` doctests, 1st/2nd-order `ForwardDiff` AD, `OffsetArrays`, and
  optional `CUDA` broadcast tests (`test/run_cuda.jl`).

## Physical measures and ownership

- `LeafFlux`: per m² **one-sided leaf area** per second, positive from leaf into air.
- `StemFlux`: per m² **one-sided woody/stem area** per second, positive from stem into air.
- `FaceFlux`: per m² **horizontal ground area** per second, positive upward.
- `VolumeSource`: per m³ **bulk canopy volume** per second, positive into the named reservoir.
- `MoistAirTendencies`: per m³ **gas volume** per second (`ε_gas ∈ (0, 1]`, conserved
  dry air `S_ρd = 0`, passive molar `CO₂`/`O₂`, excluding momentum/drag/TKE).

`leaf_air_sources` and `plant_air_sources` multiply element area densities (`LAD`,
`SAD`) once and return equal-and-opposite donor ledgers without asserting leaf gas
storage. `shift_water_energy_datum` shifts water reference enthalpy or geopotential
on individual streams **prior** to `layer_tendency` so shared interior face fluxes
telescope conservatively across multi-layer columns.

No leaf biochemistry, radiative transfer, turbulence closure, equation of state, or
time integrator belongs in `CanopyExchange`.

## Verification

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
```

For GPU verification on an environment providing `CUDA.jl`, run `test/run_cuda.jl`
under the shared GPU lock (`CUDA.allowscalar(false)`).
