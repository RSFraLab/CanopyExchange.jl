# CanopyExchange.jl

Zero-dependency, measure-safe scalar exchange contracts for layered plant canopies
and atmospheric or land host models (`ClimaLand`, `SurfaceFluxes`, `Thermodynamics`,
`ClimaAtmos`, `ClimaCore`, and standalone finite-volume columns). Extracted from
`CanopyColumn.ConservativeExchange` (`05723a9`) with Apache-2.0 attribution and exact
IEEE-754 primal bit parity (`SHA-256 19bcfb73...`).

`CanopyExchange` owns spatial area/volume measures and conservative transfer
arithmetic across `Float32`, `Float64`, `BigFloat`, `CuArray`, `ClimaCore.Fields.Field`
broadcasting, and `ForwardDiff` automatic differentiation. Leaf biochemistry,
radiative transfer, equations of state, turbulence closures, and time integrators
compose outside this package.

## Quickstart

```julia
using CanopyExchange

# 1. Layer geometry (thickness [m], LAD [m² leaf m⁻³], SAD [m² stem m⁻³])
geom = PlantLayerGeometry(2.0, 1.2, 0.25)

# 2. One-sided leaf & stem fluxes -> paired (; air, leaf, stem) bulk-volume sources
lf   = sensible_water_flux(95.0, 4.0e-5, -1.4e-5, 2.547e6) # H, E, CO₂, h_v(T_leaf)
sf   = sensible_stem_flux(20.0)                            # dry stem sensible heat
srcs = plant_air_sources(lf, sf, geom)                     # air + leaf + stem == 0

# 3. SurfaceFluxes static-energy boundary & diagnostic zero-storage top flux
soil = surface_static_energy_flux(45.0, 90.0, 3.6e-5, 2.0e-6)
top  = steady_top_flux(soil, srcs.air, geom)

# 4. Finite-volume bulk divergence & porous moist-air density tendencies (ε_gas = 0.96)
tend  = layer_tendency(soil, top, srcs.air, geom)
moist = moist_air_tendencies(tend, 0.96)
```

## Spatial measures and sign conventions

| Wrapper | Denominator | Positive sign | Typical producer |
| --- | --- | --- | --- |
| `LeafFlux` | $\mathrm{m^2}$ one-sided leaf area (`LAD`) | Leaf $\to$ air | `sensible_water_flux` |
| `StemFlux` | $\mathrm{m^2}$ one-sided woody/stem area (`SAD`) | Stem $\to$ air | `sensible_stem_flux` |
| `FaceFlux` | $\mathrm{m^2}$ horizontal ground area | Upward ($+z$) | `surface_static_energy_flux`, `steady_top_flux` |
| `VolumeSource` | $\mathrm{m^3}$ bulk canopy volume | Into named ledger | `leaf_air_sources`, `plant_air_sources`, `layer_tendency` |
| `MoistAirTendencies` | $\mathrm{m^3}$ gas volume ($\varepsilon_{\mathrm{gas}} \in (0, 1]$) | Into canopy air | `moist_air_tendencies` |

Cross-measure addition (`LeafFlux + StemFlux`, `LeafFlux + FaceFlux`) and mixed
floating-point precisions (`Float32 + Float64`) throw `MethodError`. Empty layers
with exact zero `LAD` or `SAD` (and zero AD partials) short-circuit nonfinite dummy
element fluxes to exact zero, while perturbed dual densities (`Dual(0.0, 1.0)`)
preserve exact right-directional sensitivities.

## Where to find things

| Question | Location |
| --- | --- |
| Structured API Reference for all exported and `public` symbols | [docs/src/api/index.md](docs/src/api/index.md) |
| Continuous-scalar `{FT<:Real}` traits & `@generated` leaf-float predicates | [src/scalars.jl](src/scalars.jl), [docs/src/api/clima_and_scalars.md](docs/src/api/clima_and_scalars.md) |
| Measure-tagged wrappers, `MoistAirTendencies`, layer geometries & `Base` interfaces | [src/quantities.jl](src/quantities.jl), [docs/src/api/quantities_and_geometries.md](docs/src/api/quantities_and_geometries.md) |
| `LAD`/`SAD` area-to-volume ledgers & finite-volume flux divergence | [src/sources.jl](src/sources.jl), [docs/src/pages/measures_and_conservation.md](docs/src/pages/measures_and_conservation.md) |
| CliMA (`SurfaceFluxes`, `Thermodynamics`, `ClimaLand`, `ClimaAtmos`, `ClimaCore`) coupling | [src/clima_contracts.jl](src/clima_contracts.jl), [docs/src/pages/clima_interoperability.md](docs/src/pages/clima_interoperability.md) |
| QA (`Aqua`, `ExplicitImports`, doctests), 1st/2nd-order `ForwardDiff` AD & `OffsetArray` tests | [test/qa_and_ad.jl](test/qa_and_ad.jl), [test/test_conservative_exchange.jl](test/test_conservative_exchange.jl) |
| Extraction bit-parity validation & NVIDIA CUDA broadcast suite | [docs/validation.md](docs/validation.md), [test/run_cuda.jl](test/run_cuda.jl) |

## Verification and documentation build

```sh
# Run full CPU test suite (Aqua, ExplicitImports, doctests, bit parity, AD, OffsetArrays)
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'

# Build local Documenter.jl HTML site (checkdocs = :all, doctest = true)
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
```

For GPU verification on an environment providing `CUDA.jl`, run `test/run_cuda.jl`
(`CUDA.allowscalar(false)`).
