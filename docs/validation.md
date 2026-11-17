# Extraction validation, 2026-11-08 (historical, v0.1.0/v0.1.1)

This section records the extraction checkpoint. The file names and hashes below refer
to the v0.1.x sources (`src/quantities.jl`, `src/sources.jl`), which the 0.2.0 core
replaced; the leaf bit record and its SHA-256 still hold.

This package extracts the conservative scalar contracts from CanopyColumn
`05723a9eec3baa96da7185dd164a88d622ee72d1`. The original checkout is preserved.
The compatibility branch uses `import CanopyExchange` and
`const ConservativeExchange = CanopyExchange`; there is only one owner of the
new runtime types. Historical serialized object identities are outside this claim.

Both `src/quantities.jl` and `src/sources.jl` are byte-for-byte unchanged:

| File | SHA-256 |
|---|---|
| quantities.jl | d83101868892c09adbc33124e9ce10e3fb17811dcee0ef8d8781339d7e5b4710 |
| sources.jl | 99b2263e064324784ea6f49dfe5507a1b47b8edc3675f47e2148a2a69aeb0836 |

Executed using Julia 1.12.6:

- `Pkg.test()` passed 14 additional-tracer, 228 inherited exchange and 50
  independent physical-measure checks (292 total).
- Optional CUDA broadcast verification passed 14,420 checks, Float32/Float64,
  computed and absent oxygen, on physical GPU 1 under the shared GPU lock.
- Actual CanopyColumn compatibility tests passed in separate processes in both
  load orders: 31 identity/mixed-call checks plus the 242 inherited checks each.
- Full CanopyColumn `Pkg.test()` passed 1,415 checks, including the alias identity
  assertion registered in its normal test suite.
- Before/after arithmetic records from the old actual package and the canonical
  package are byte-identical, including IEEE bitstrings of fluxes, paired sources
  and tendencies in both precisions and oxygen modes. Their common SHA-256 is
  `19bcfb73ba7adb1e88f285bbc2679182a861902729216e4948b9d3d350f0d08d`.

The independent additions exercise nonunit ground footprint, unequal layer
thicknesses, one-sided LAD, face-only transfer signs, internal-face cancellation,
layer subdivision, water-energy reference changes on both faces and leaves,
absent versus computed-zero oxygen, and invalid active versus absent-leaf inputs.
The initial new Float32 datum test used a relative-only comparison of a nearly
cancelling net difference and failed. Its corrected assertion uses an operand-
based forward-roundoff bound with a declared operation budget, and checks that
bound is smaller than half the resolved signal. Production arithmetic and inherited
physical/numerical acceptance thresholds were not changed. The failure log is retained.

Both independent code/structure and physics/biophysics reviewers cleared this
package and the alias. These checks establish the scalar extraction and source
compatibility. They do not establish a new canopy transport closure, a ClimaAtmos
lifecycle adapter, a wet-leaf model or a global throughput result. Consumer migration
and native-package co-loading are separate integration gates.

## Reproduction

Run `julia --project=. -e 'using Pkg; Pkg.test()'`. For the optional device check,
develop this package in an environment with CUDA and execute `test/run_cuda.jl`.
The compatibility branch in CanopyColumn supplies `test/test_exchange_alias.jl`:
run it with `canonical-first` and `prototype-first` in fresh processes, then run
that package's full tests.

For numerical comparison, execute `test/write_arithmetic_record.jl legacy before.toml`
in an environment holding the original CanopyColumn and the same script with
`canonical after.toml` in an environment holding this package. Compare the files
byte-for-byte. No Julia serialization or module-name comparison is used as a
substitute for the numerical comparison.

Local logs for this checkpoint are under `/tmp/canopy-exchange-20261108/logs`;
`package-test-initial-cancellation.log` preserves the initial assertion failure.
The CPU environments were cloned from previously accepted environments and
resolved using `Pkg.develop`, without modifying the original manifests. Julia 1.10
compatibility is retained from the source package; this checkpoint executed 1.12.6.

# Open tracer core (Unreleased 0.2.0), 2026-11-17

Before `src/` was changed, `test/write_reference_record.jl` wrote
`test/reference/exchange_v0_1_1.toml` from v0.1.1 (commit `c35ad8d`): IEEE bit strings
of leaf fluxes, leaf air and donor sources, layer tendencies, the stem path (dry and
transpiring stems), datum shifts of leaf, stem, face and volume streams,
`steady_top_flux`, `accumulate_steady_fluxes`, `surface_static_energy_flux` and
`moist_air_tendencies` at gas fraction 1 and 0.96, in `Float32` and `Float64`, with
and without O₂. That file is never regenerated from later code.

Executed on the refactored code (two NVIDIA L40S):

- `Pkg.test()` on Julia 1.12.6 (ForwardDiff 1.5.0) and on Julia 1.10.12 (ForwardDiff
  1.5.0) passed; the same suite with ForwardDiff 0.10.22 forced in a separate
  environment passed. The frozen record (105 checks) and the leaf SHA-256 record
  `19bcfb73...` are reproduced bit for bit.
- `test/run_cuda.jl` (CUDA 6.5.1) passed 14,456 checks; `test/run_fields.jl`
  (CUDA 6.4.1, ClimaCore 1.0.8) passed on `CuArray` and on ClimaCore fields on CPU and
  CUDA, `Float32` and `Float64`.
- The Documenter build (`checkdocs = :all`, doctests) completed without errors.

Reproduction: the commands in the README; for ForwardDiff 0.10, develop this package
into an environment with the test dependencies and `ForwardDiff = "0.10"` pinned and
run `test/runtests.jl` with that project.
