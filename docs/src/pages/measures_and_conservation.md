# Measures & Conservation

`CanopyExchange` enforces dimensional measure safety and discrete finite-volume
conservation through immutable, `isbitstype` wrappers parameterized on a continuous
scalar type `FT <: Real`.

## 1. Conserved Budget Quantities (`ExchangeRates`)

[`ExchangeRates`](@ref) bundles four extensively conserved transfer rates per second
and per wrapper spatial denominator:

1. `energy` [$\mathrm{J\,s^{-1}}$ per measure]: total energy transfer under a
   caller-declared thermodynamic datum (not sensible or latent heat in isolation).
2. `water` [$\mathrm{kg\,H_2O\,s^{-1}}$ per measure]: signed water mass transfer.
3. `carbon` [$\mathrm{mol\,C\,s^{-1}}$ per measure]: signed molar $\mathrm{CO_2}$
   transfer ($1\,\mathrm{mol\,C}$ per $\mathrm{mol\,CO_2}$; photosynthetic uptake
   from air into leaf is negative in `LeafFlux`).
4. `oxygen` ([`OxygenContent`](@ref)): either [`NoOxygen`](@ref) (`oxygen_value === nothing`,
   zero storage overhead) or [`OxygenTransfer`](@ref) [$\mathrm{mol\,O_2\,s^{-1}}$
   per measure].

Adding or combining `NoOxygen` with `OxygenTransfer` throws `MethodError`, so an
uncomputed oxygen budget can never silently masquerade as zero physical $\mathrm{O_2}$
exchange.

## 2. Paired Area-to-Volume Ledgers (`LAD` and `SAD`)

Canopy layers hold finite positive `thickness` $\Delta z_k > 0$ [$\mathrm{m}$] and
nonnegative one-sided area densities:
- [`LayerGeometry`](@ref): leaf-only layers with `leaf_area_density` ($\mathrm{LAD}_k \ge 0$).
- [`PlantLayerGeometry`](@ref): coupled foliage and woody layers with
  `leaf_area_density` ($\mathrm{LAD}_k \ge 0$) and `stem_area_density`
  ($\mathrm{SAD}_k \ge 0$), where [`plant_area_density`](@ref) returns
  $\mathrm{PAI}_k / \Delta z_k = \mathrm{LAD}_k + \mathrm{SAD}_k$.

[`leaf_air_sources`](@ref) multiplies a [`LeafFlux`](@ref) by $\mathrm{LAD}_k$ once
and returns `(air = +S_leaf, leaf = -S_leaf)` using unary IEEE-754 sign-bit negation.
[`plant_air_sources`](@ref) scales `leaf_flux` by $\mathrm{LAD}_k$ and `stem_flux`
by $\mathrm{SAD}_k$ independently, returning `(air, leaf, stem)` such that:

```math
S_{\mathrm{air},k} + \left(S_{\mathrm{leaf},k} + S_{\mathrm{stem},k}\right) = 0.
```

The negative donor ledgers (`leaf`, `stem`) record instantaneous transfer out of the
plant surface rather than asserting internal leaf gas storage.

## 3. Empty-Layer Nonfinite Guards vs. Directional AD Tangents

In batched CPU and GPU arrays, layers with zero foliage ($\mathrm{LAD}_k = 0$) or
zero woody area ($\mathrm{SAD}_k = 0$) frequently hold `NaN` or `Inf` sentinel rates
in unused slots. [`area_source`](@ref) inspects [`CanopyExchange.exact_zero_scalar`](@ref):

- **Exact zero density (primal and all AD partials equal zero):** returns exact
  `zero(FT)` without multiplying `rate`, preventing `0 * NaN == NaN` contamination on
  CPU, `CuArray`, and `ClimaCore.Fields.Field` kernels.
- **Perturbed zero density (`ForwardDiff.Dual(0.0, 1.0)`):** evaluates `rate * area_density`
  directly (requiring finite `rate`) so the right directional derivative
  $\partial(\text{rate}\cdot a)/\partial a\big|_{a=0^+} = \text{rate}$ is preserved.

## 4. Telescoping Finite-Volume Divergence

Across $N$ stacked layers with shared upward-positive horizontal interface fluxes
$F_{1}, \dots, F_{N+1}$ ([`FaceFlux`](@ref)), [`layer_tendency`](@ref) evaluates:

```math
\left(\frac{\partial \psi}{\partial t}\right)_{\mathrm{bulk},k}
= \frac{F_k - F_{k+1}}{\Delta z_k} + S_{\mathrm{air},k}.
```

Integrating over the column $\sum_{k=1}^N \Delta z_k (\partial_t \psi)_{\mathrm{bulk},k}$
cancels every interior face flux $F_2, \dots, F_N$ identically:

```math
\sum_{k=1}^N \Delta z_k \left(\frac{\partial \psi}{\partial t}\right)_{\mathrm{bulk},k}
= F_1 - F_{N+1} + \sum_{k=1}^N \Delta z_k \, S_{\mathrm{air},k}.
```
