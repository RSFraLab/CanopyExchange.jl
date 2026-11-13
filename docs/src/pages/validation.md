# Validation & Bit Parity

## 1. Canonical Extraction & IEEE-754 Bit Parity

`CanopyExchange` was extracted from `CanopyColumn.ConservativeExchange`
(`05723a9eec3baa96da7185dd164a88d622ee72d1`) and subsequently extended in `v0.1.1`
with continuous-scalar `{FT<:Real}` automatic differentiation, `StemFlux` /
`PlantLayerGeometry`, and CliMA boundary/tendency contracts.

Homogeneous primal `Float32` and `Float64` evaluations preserve **exact IEEE-754 bit
parity** with the original donor implementation:

- `test/write_arithmetic_record.jl canonical <out.toml>` verifies that all primal
  `Float32` and `Float64` bitstrings across computed and absent oxygen match the
  locked SHA-256 digest `19bcfb73ba7adb1e88f285bbc2679182a861902729216e4948b9d3d350f0d08d`.
- Unary negation uses `-field` rather than `zero(FT) - field`, preserving IEEE `-0.0`
  donor sign bits on empty layers.
- Scalar division uses fieldwise `/ α` rather than `* inv(α)`, preserving exact
  identity at `gas_fraction == 1` and zero-ULP molar `O₂ / CO₂` stoichiometry.

## 2. Continuous-Scalar AD, QA & Device Verification Gates

`Pkg.test()` (`test/runtests.jl`) executes deterministic quality-assurance and
numerical invariants in a clean sandboxed test environment:

- **`Aqua.test_all(CanopyExchange)`:** zero method ambiguities, zero unbound type
  parameters (including `Base.eltype` and `Base.promote_rule` on `Union{}`), zero
  undefined exports, zero piracies, and clean `Project.toml` `[compat]`.
- **`ExplicitImports.jl` (Julia $\ge 1.11$):** passes all seven checks
  (`check_no_implicit_imports`, `check_no_stale_explicit_imports`,
  `check_all_explicit_imports_via_owners`, `check_all_explicit_imports_are_public`,
  `check_all_qualified_accesses_via_owners`, `check_all_qualified_accesses_are_public`,
  and `check_no_self_qualified_accesses`).
- **Documentation & doctests:** zero `Docs.undocumented_names(CanopyExchange)` across
  exported and `public` bindings, and `Documenter.doctest(CanopyExchange)` verified
  on both `Pkg.test()` and `docs/make.jl` (`checkdocs = :all`).
- **1st- and 2nd-order `ForwardDiff` vs. `BigFloat`:** end-to-end gradients and
  Hessians of coupled `PlantLayerGeometry` + `plant_air_sources` + `steady_top_flux` +
  `layer_tendency` + `moist_air_tendencies` pipelines match `BigFloat` central
  differences (`rtol = 1e-11` gradients, `rtol = 1e-9` Hessians), while rejecting
  `NaN`/`Inf` in dual partials via `@generated` [`CanopyExchange.finite_scalar`](@ref).
- **Generic indexing & zero allocations:** `accumulate_steady_fluxes!` and scalar
  kernels allocate `0` bytes after warmup and support `OffsetArray`, `@view`, and
  0-layer columns.
- **GPU & CliMA co-load verification:** verified on NVIDIA L40S GPUs
  (`test/run_cuda.jl`, `CUDA.allowscalar(false)`) in `Float32` and `Float64`, and
  co-loaded with `ClimaAtmos`, `ClimaCore`, `SurfaceFluxes`, and `Thermodynamics`
  (`ClimaCore.Fields.Field` `@.` broadcasting in `Float32` and `Float64`).
