# Area Sources & Finite-Volume Divergence

Source and divergence kernels map one-sided leaf (`LAD`) and woody/stem (`SAD`)
fluxes into paired donor-conserving bulk-volume sources and evaluate finite-volume
vertical flux divergences across shared horizontal interfaces.

```@contents
Pages = ["sources_and_tendencies.md"]
Depth = 2
```

## Scalar Area-to-Volume Primitives

```@docs
area_source
leaf_area_source
stem_area_source
plant_area_density
```

## Enthalpy-Coupled Element Flux Constructors

```@docs
sensible_water_flux
sensible_stem_flux
```

## Paired Leaf, Stem & Canopy-Air Ledgers

```@docs
leaf_air_sources
plant_air_sources
```

## Finite-Volume Vertical Divergence

```@docs
conserved_tendency
layer_tendency
```
