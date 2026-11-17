# CanopyExchange.jl

Measure-safe, conservative scalar exchange contracts for layered plant canopies and
atmospheric or land host models (`ClimaLand`, `SurfaceFluxes`, `Thermodynamics`,
`ClimaAtmos`, `ClimaCore`, and standalone finite-volume columns), with no runtime
dependencies. Extracted from `CanopyColumn.ConservativeExchange` (`05723a9`) with
Apache-2.0 attribution; the arithmetic is fixed by two frozen IEEE-754 bit records
(leaf record `SHA-256 19bcfb73...`, and the extended v0.1.1 record in
`test/reference/`).

`CanopyExchange` owns spatial area/volume measures and conservative transfer
arithmetic for energy, water and an open set of passive molar tracers, across
`Float32`, `Float64`, `BigFloat`, `CuArray`, `ClimaCore.Fields.Field` broadcasting,
and `ForwardDiff` automatic differentiation (through a package extension). Leaf
biochemistry, radiative transfer, equations of state, turbulence closures, and time
integrators compose outside this package.

## Quickstart

```julia
using CanopyExchange

# 1. Layer geometry (thickness [m], LAD [m² leaf m⁻³], SAD [m² stem m⁻³])
geom = PlantLayerGeometry(2.0, 1.2, 0.25)

# 2. Passive molar tracers: an uncomputed tracer is simply not a key
T = Tracers(:co2, :o2)

# 3. One-sided leaf & stem fluxes -> paired (; air, leaf, stem) bulk-volume sources
lf   = sensible_water_flux(95.0, 4.0e-5, 2.547e6, T(-1.4e-5, 1.5e-5)) # H, E, h_v(T_leaf), tracers
sf   = sensible_stem_flux(20.0, T(0.0, 0.0))                          # dry stem sensible heat
srcs = plant_air_sources(lf, sf, geom)                                # air + leaf + stem == 0

# 4. SurfaceFluxes static-energy boundary & diagnostic zero-storage top flux
soil = surface_static_energy_flux(45.0, 90.0, 3.6e-5, T(2.0e-6, -2.0e-6))
top  = steady_top_flux(soil, srcs.air, geom)

# 5. Finite-volume bulk divergence & porous moist-air density tendencies (ε_gas = 0.96)
tend  = layer_tendency(soil, top, srcs.air, geom)
moist = moist_air_tendencies(tend, 0.96)
moist.tracers.o2
```

## Records, measures and sign conventions

`ExchangeRates(energy, water, tracers)` holds energy [J s⁻¹] including transported
vapour enthalpy (never add a second latent-heat term), total water [kg s⁻¹] and a
`NamedTuple` of passive molar tracers [mol s⁻¹], per measure. `Exchange{M}(rates)`
tags them with the measure `M`:

| Record | Measure | Positive sign | Typical producer |
| --- | --- | --- | --- |
| `LeafFlux` = `Exchange{LeafArea}` | m² one-sided leaf area (`LAD`) | leaf → air | `sensible_water_flux`, `leaf_flux` |
| `StemFlux` = `Exchange{StemArea}` | m² one-sided stem area (`SAD`) | stem → air | `sensible_stem_flux`, `stem_flux` |
| `FaceFlux` = `Exchange{GroundArea}` | m² horizontal ground area | upward (+z) | `surface_static_energy_flux`, `face_flux`, `steady_top_flux` |
| `VolumeSource` = `Exchange{BulkVolume}` | m³ bulk canopy volume | into named ledger | `volume_source`, `leaf_air_sources`, `plant_air_sources`, `layer_tendency` |
| `MoistAirTendencies` | m³ gas volume (ε_gas ∈ (0, 1]) | into canopy air | `moist_air_tendencies` |

Cross-measure addition (`LeafFlux + StemFlux`) is a `MethodError`. Mixed
floating-point precision (`Float32` with `Float64`), integer inputs, and records with
different tracer sets (`(; co2)` with `(; co2, o2)`, or another order) are
`ArgumentError`s. Tracers are molar and passive; water isotopologues are subsets of
`water` and never enter total mass; an `o2` key does not create an oxygen inventory.
Empty layers with exact zero `LAD` or `SAD` (and zero AD partials) turn nonfinite
dummy element fluxes into exact zeros, while perturbed dual densities
(`Dual(0.0, 1.0)`) keep the right derivative.

A layer geometry is any type with `thickness(g)`, `area_density(g, M)` and
`Base.eltype`; a new element class is a new `Measure` subtype. Neither needs changes
to this package (see `test/test_extension_points.jl`).

## Where to find things

| Question | Location |
| --- | --- |
| API reference for all exported and `public` names | [docs/src/api/index.md](docs/src/api/index.md) |
| Scalar hooks and the precision policy | [src/precision.jl](src/precision.jl), [ext/CanopyExchangeForwardDiffExt.jl](ext/CanopyExchangeForwardDiffExt.jl) |
| `ExchangeRates`, `Tracers`, `MoistAirTendencies`, record algebra, equality | [src/rates.jl](src/rates.jl) |
| Measures and `Exchange` records, function-form constructors | [src/measures.jl](src/measures.jl) |
| Layer geometry interface | [src/geometry.jl](src/geometry.jl) |
| Area-to-volume sources and paired ledgers | [src/area_to_volume.jl](src/area_to_volume.jl), [docs/src/pages/measures_and_conservation.md](docs/src/pages/measures_and_conservation.md) |
| Flux divergence and steady columns | [src/flux_divergence.jl](src/flux_divergence.jl) |
| Vapour enthalpy, SurfaceFluxes boundary, energy datum | [src/boundary_fluxes.jl](src/boundary_fluxes.jl), [docs/src/pages/clima_interoperability.md](docs/src/pages/clima_interoperability.md) |
| Gas-volume moist-air tendencies | [src/moist_air.jl](src/moist_air.jl) |
| Tests: invariants, contracts, extension points, inference, QA and AD | [test/runtests.jl](test/runtests.jl) |
| Frozen bit records and device scripts | [docs/validation.md](docs/validation.md), [test/reference/](test/reference/), [test/run_cuda.jl](test/run_cuda.jl), [test/run_fields.jl](test/run_fields.jl) |

## Verification and documentation build

```sh
# Full CPU test suite (Aqua, ExplicitImports, doctests, frozen bit records, AD, OffsetArrays)
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'

# Opt-in device checks (each builds its own temporary environment)
julia --startup-file=no test/run_cuda.jl     # CuArray broadcasts
julia --startup-file=no test/run_fields.jl   # ClimaCore fields on CPU and CUDA

# Documenter.jl site (checkdocs = :all, doctest = true)
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
```
