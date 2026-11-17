# Validation & Bit Parity

## 1. Frozen reference records

Two committed records fix the arithmetic, so a change that moves every code path
together is still detected:

- **Leaf record** (`test/write_arithmetic_record.jl canonical <out.toml>`): IEEE-754
  bit strings of leaf fluxes, paired leaf sources and layer tendencies, `Float32` and
  `Float64`, with and without O₂. Its SHA-256 is locked at
  `19bcfb73ba7adb1e88f285bbc2679182a861902729216e4948b9d3d350f0d08d`, unchanged since
  the extraction from `CanopyColumn.ConservativeExchange` (`05723a9`).
- **Extended record** (`test/reference/exchange_v0_1_1.toml`, written by
  `test/write_reference_record.jl` from v0.1.1, commit `c35ad8d`, and never
  regenerated from later code): additionally the stem path with dry and transpiring
  stems (`plant_air_sources`), the water-energy datum shift of leaf, stem, face and
  volume streams, `steady_top_flux`, `accumulate_steady_fluxes`,
  `surface_static_energy_flux` and `moist_air_tendencies` at gas fraction 1 and below
  1. `test/test_reference_record.jl` recomputes it and compares bit for bit.

The 0.2.0 open-tracer core reproduces both records exactly: CO₂ and O₂ are carried as
the tracers `co2` and `o2` in the order of the old `carbon` and `oxygen` fields. Unary
negation uses `-field` (keeping IEEE `-0.0` donor sign bits on empty layers), and
division by the gas fraction is fieldwise `/ ε` (exact at `ε == 1`).

## 2. Test gates

`Pkg.test()` (`test/runtests.jl`) runs:

- **Physical invariants** (`test_conservative_exchange.jl`,
  `test_migration_invariants.jl`): interior-face cancellation on a nonuniform grid,
  donor ledgers, datum shifts on faces and leaves with an operand-based round-off
  bound, layer subdivision, absent versus computed-zero tracers, empty layers.
- **Contracts** (`test_contracts.jl`): the precision policy in every constructor and
  kernel (`ArgumentError` for `Float32` with `Float64` and for integers), measure
  safety, tracer rules, the counterexamples of the first review round, and AD
  contracts (mixed plain and dual inputs, right derivative at an empty layer,
  rejection of nonfinite derivative components, gradient against central differences).
- **Extension points** (`test_extension_points.jl`): a second tracer set
  `(co2, cos, ch4)` through the tagged records, bit-identical to the scalar path, and a
  third element class with its own geometry defined only in the test.
- **Inference and allocations** (`test_inference_and_allocations.jl`): `@inferred`,
  result element type, `isbits` results and zero allocations for every exported
  kernel and constructor in `Float32` and `Float64`.
- **QA and AD** (`qa_and_ad.jl`): `Aqua.test_all`, `ExplicitImports` (Julia ≥ 1.11),
  no undocumented names, doctests, first and second derivatives of a coupled
  leaf + stem + column + moist-air pipeline against `BigFloat` central differences,
  `BigFloat` duals, container promotion, value equality and hashing, `OffsetArray`
  and view columns, zero-allocation kernels.

Opt-in device scripts build their own temporary environment:
`test/run_cuda.jl` (`CuArray` broadcasts, `CUDA.allowscalar(false)`) and
`test/run_fields.jl` (`ClimaCore.Fields.Field` on CPU and CUDA, `Float32` and
`Float64`). GPU arithmetic may contract `a*b + c`, so sums are compared to rounding
tolerance there and pure products exactly.
