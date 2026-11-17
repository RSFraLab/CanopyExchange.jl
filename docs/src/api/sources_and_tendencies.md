# Area Sources & Finite-Volume Divergence

Element fluxes become bulk-volume sources through one multiplication by their own
area density, with paired donor ledgers; face fluxes and sources become conservative
layer tendencies.

```@contents
Pages = ["sources_and_tendencies.md"]
Depth = 2
```

## Element area to bulk volume

```@docs
area_source
volume_source
leaf_area_source
stem_area_source
```

## Paired leaf, stem and canopy-air ledgers

```@docs
leaf_air_sources
plant_air_sources
```

## Finite-volume vertical divergence

```@docs
conserved_tendency
layer_tendency
```
