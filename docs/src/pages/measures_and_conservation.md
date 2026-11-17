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
   example `(; co2, o2)` or `(; co2, cos, ch4)`. CO₂ is in mol CO₂ (one mol C per mol
   CO₂); photosynthetic uptake from air into the leaf is a negative `co2` in a
   [`LeafFlux`](@ref), and condensation onto a surface a negative `water`.

Host adapters own the conversions between conserved density tendencies and
temperature or humidity.

Energy and water are named fields because two operations treat them specially: the
water-energy datum shift and the total-mass source of the host (equal to the water
source with zero dry-air mass source). Everything else applies componentwise to the
tracers present. The tracer rules:

- Tracers are molar and passive; they carry no energy and no mass in these records.
  Leaving their mass out of the total-mass source is the **trace approximation**. With
  the test fixture's leaf values (12 µmol CO₂ m⁻² s⁻¹ against 4·10⁻⁵ kg H₂O m⁻² s⁻¹),
  $M_{\mathrm{CO_2}} S_{\mathrm{CO_2}}$ is about 1.3 % of the water mass source, about
  0.3 % net when the opposite O₂ flux is counted, and it dominates when transpiration
  is near zero (at night).
- Water isotopologues are subsets of `water`; as tracers they never enter total mass.
- Bulk O₂ is not a trace constituent. `tracers.o2` is a transfer ledger; it becomes a
  prognostic oxygen budget only if the host has a separate oxygen store contract.
- A tracer that is not computed is **not a key**. Records with and without it are
  different types and cannot be combined (`ArgumentError`), so an uncomputed tracer
  can never masquerade as a measured zero. Conversely, a key with value zero (for
  example `(; co2 = 0.0)` on a dry stem) asserts a zero transfer; it does not mean
  "not computed".
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

A layer geometry is any type with [`CanopyExchange.thickness`](@ref) [m] and
[`CanopyExchange.area_density`](@ref) [m² m⁻³] for each element class it holds (the
interface is documented in one place, the `area_density` docstring). The package
provides
[`LayerGeometry`](@ref) (leaves) and [`PlantLayerGeometry`](@ref) (leaves and stems;
[`plant_area_density`](@ref) $= \mathrm{LAD}_k + \mathrm{SAD}_k$).

[`volume_source`](@ref) multiplies an element flux by its own area density once.
[`leaf_air_sources`](@ref) returns `(air = +S_leaf, leaf = -S_leaf)` with IEEE sign-bit
negation; [`plant_air_sources`](@ref) scales leaf and stem fluxes independently and
returns `(air, leaf, stem)` with

```math
S_{\mathrm{air},k} + \left(S_{\mathrm{leaf},k} + S_{\mathrm{stem},k}\right) = 0
```

exactly, with the parentheses as written; other groupings are zero to rounding.

Several populations of one element (sunlit and shaded leaves, cohorts) share
[`LeafArea`](@ref): each flux is per m² of its own leaf area, so the layer source is
$\sum_c$ `volume_source(F_c, LayerGeometry(Δz, LAD_c))`. Adding the `LeafFlux`es of
different populations without their area weights is wrong. Each donor ledger is the exact negative of its own contribution; it
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

A third element class (for example fruit) needs a `Measure` subtype `FruitArea`, a
geometry type with a `thickness` (field or method) and
`CanopyExchange.area_density(g, ::Type{FruitArea})`, and fluxes built with
`sensible_water_flux(FruitArea(), H, E, h_v, tracers)` (which enforces the vapour
enthalpy) or `Exchange{FruitArea}(rates)`. [`volume_source`](@ref),
[`layer_tendency`](@ref) and both steady-column functions then work unchanged;
`Base.eltype` is not needed. Define `Base.Broadcast.broadcastable(g) = Ref(g)` to use
the geometry as a scalar in broadcasts. The package cannot validate such a geometry:
it must itself guarantee a finite positive thickness and finite nonnegative area
densities. The test suite (`test/test_extension_points.jl`) does this and also
carries a second tracer set `(co2, cos, ch4)` through the records.

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
