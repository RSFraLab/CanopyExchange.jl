# CanopyExchange.jl

Small scalar contracts for conservative leaf–air and vertical face transfers.
This local, private-first package has no runtime dependencies. It was extracted
from `CanopyColumn.ConservativeExchange` at CanopyColumn commit `05723a9`,
preserving the arithmetic and Apache-2.0 attribution.

## Finding the code

- `src/quantities.jl`: typed leaf-area fluxes, ground-area face fluxes, bulk-volume
  sources, layer geometry and optional oxygen.
- `src/sources.jl`: area conversion, paired exchange and finite-volume divergence.
- `test/parameters/`: manufactured inputs and numerical tolerances in TOML.
- `test/runtests.jl`: independent conservation, precision and invalid-input checks.
- `test/run_cuda.jl`: optional CUDA broadcast verification; CUDA is supplied by the
  caller's test environment, not a dependency of the package.

All scientific values use one floating-point type; Float32 and Float64 are tested.
The immutable scalar records can be broadcast on CPU or GPU. Device arrays, grids,
solvers and physics laws remain the caller's choice.

## Measures and ownership

`LeafFlux` is per square metre of **one-sided leaf area**. `FaceFlux` is per square
metre of common horizontal ground area and positive upward. `VolumeSource` is per
cubic metre of **bulk canopy**, not physical gas volume. `leaf_air_sources` applies
LAD exactly once and returns opposite air/leaf ledgers. `layer_tendency` adds
`(bottom - top)/thickness`. Hosts own metrics, porosity and conversion to their
conserved variables. A leaf ledger does not imply a leaf CO₂/O₂ storage model.

`NoOxygen` means not computed, distinctly from a computed zero. Scalar
`leaf_area_source` and `conserved_tendency` support additional tracers without
changing the convenience records. `valid_exchange` tests finiteness only; it does
not establish admissible thermodynamics, positive inventories or stoichiometry.

Energy has a caller-declared reference and definition. `sensible_water_flux`
constructs `Q + E*h_v`; it excludes mass-carried kinetic/potential energy,
trace-gas enthalpy and biochemical energy. Do not add that helper's result to a
separate complete thermal packet. An equal/opposite ledger does not by itself
prove compatibility with a host's energy variable.

No leaf biochemistry, radiation, thermodynamics, turbulence, host fields or time
integration belongs here. The original CanopyColumn namespace is retained as a
compatibility alias in the integration branch, with one canonical type owner.
Existing serialized objects from the old module are not a supported migration
format.

## Local verification

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

For CUDA, develop this package in an environment containing CUDA and run
`test/run_cuda.jl` under the shared GPU allocation policy. Tests exercise kernels;
they are not throughput benchmarks or complete atmosphere-coupling validation.
