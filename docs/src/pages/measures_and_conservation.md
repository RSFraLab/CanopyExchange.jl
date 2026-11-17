# Measures & Conservation

`CanopyExchange` enforces measure safety and discrete finite-volume conservation with
immutable, `isbits` records whose scalars share one floating-point type.

## 1. Conserved quantities and passive tracers (`ExchangeRates`)

[`ExchangeRates`](@ref) holds, per second and per measure of the enclosing record:

1. `energy` [$\mathrm{J\,s^{-1}}$]: total energy transfer under one caller-declared
   thermodynamic datum, including the enthalpy carried by transported water vapour.
   No second latent-heat term may be added when coupling to a conserved energy budget.
2. `water` [$\mathrm{kg\,H_2O\,s^{-1}}$]: total water mass transfer.
3. `tracers` [$\mathrm{mol\,s^{-1}}$]: passive molar tracers in a `NamedTuple`, for
   example `(; co2, o2)` or `(; co2, cos, ch4)`.

Energy and water are named fields because two operations treat them specially: the
water-energy datum shift and the total-mass source of the host (equal to the water
source with zero dry-air mass source). Everything else applies componentwise to the
tracers present. The tracer rules:

- Tracers are molar and passive; they carry no energy and no mass in these records.
- Water isotopologues are subsets of `water`; as tracers they never enter total mass.
- An `o2` key carries the O₂ transfer; it does not create an oxygen inventory.
- A tracer that is not computed is **not a key**. Records with and without it are
  different types and cannot be combined (`ArgumentError`), so an uncomputed tracer
  can never masquerade as a measured zero.
- Tracer order is part of the type.

## 2. Measures (`Exchange{M}`)

| Record | Measure | Positive sign |
| --- | --- | --- |
| [`LeafFlux`](@ref) = `Exchange{LeafArea}` | m² one-sided leaf area | leaf → air |
| [`StemFlux`](@ref) = `Exchange{StemArea}` | m² one-sided stem area | stem → air |
| [`FaceFlux`](@ref) = `Exchange{GroundArea}` | m² horizontal ground area | upward |
| [`VolumeSource`](@ref) = `Exchange{BulkVolume}` | m³ bulk canopy volume | into the named reservoir |

Arithmetic exists only within one measure; `LeafFlux + StemFlux` and
`LeafFlux + FaceFlux` are `MethodError`s. A new element class is a new
[`Measure`](@ref) subtype; see section 5.

## 3. Paired area-to-volume ledgers

A layer geometry is any type with [`thickness`](@ref) [m], [`area_density`](@ref)
[m² m⁻³] for each element class it holds, and `Base.eltype`. The package provides
[`LayerGeometry`](@ref) (leaves) and [`PlantLayerGeometry`](@ref) (leaves and stems;
[`plant_area_density`](@ref) $= \mathrm{LAD}_k + \mathrm{SAD}_k$).

[`volume_source`](@ref) multiplies an element flux by its own area density once.
[`leaf_air_sources`](@ref) returns `(air = +S_leaf, leaf = -S_leaf)` with IEEE sign-bit
negation; [`plant_air_sources`](@ref) scales leaf and stem fluxes independently and
returns `(air, leaf, stem)` with

```math
S_{\mathrm{air},k} + \left(S_{\mathrm{leaf},k} + S_{\mathrm{stem},k}\right) = 0
```

to round-off. Each donor ledger is the exact negative of its own contribution; it
records transfer out of the plant surface, not plant gas storage.

## 4. Empty layers and derivatives

Layers with zero foliage or zero woody area often hold `NaN` or `Inf` dummy rates.
[`area_source`](@ref) asks [`CanopyExchange.exact_zero`](@ref):

- **Exact zero density** (value and every derivative component zero): an exact zero,
  without contamination by `0 * NaN`, on CPU, `CuArray` and ClimaCore fields.
- **Perturbed zero density** (`ForwardDiff.Dual(0.0, 1.0)`): `rate * area_density`, so
  the right derivative $\partial(\text{rate}\cdot a)/\partial a|_{a=0^+} = \text{rate}$
  is kept (requires a finite rate).

## 5. Extension without changing the package

A third element class (for example fruit) needs a `Measure` subtype, a geometry type
with `thickness`, `area_density(g, ::Type{FruitArea})` and `Base.eltype`, and
`Exchange{FruitArea}(rates)` fluxes; [`volume_source`](@ref), [`layer_tendency`](@ref)
and the steady-column functions then work unchanged. The test suite
(`test/test_extension_points.jl`) does this and also carries a second tracer set
`(co2, cos, ch4)` through the records.

## 6. Telescoping finite-volume divergence

Across $N$ stacked layers with shared upward-positive face fluxes
$F_1, \dots, F_{N+1}$ ([`FaceFlux`](@ref)), [`layer_tendency`](@ref) evaluates

```math
\left(\frac{\partial \psi}{\partial t}\right)_{\mathrm{bulk},k}
= \frac{F_k - F_{k+1}}{\Delta z_k} + S_{\mathrm{air},k},
```

and the column integral cancels every interior face:

```math
\sum_{k=1}^N \Delta z_k \left(\frac{\partial \psi}{\partial t}\right)_{\mathrm{bulk},k}
= F_1 - F_{N+1} + \sum_{k=1}^N \Delta z_k \, S_{\mathrm{air},k}.
```
