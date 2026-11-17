# Records, Measures & Layer Geometries

Exchange rates ([`ExchangeRates`](@ref): energy, water and an open set of passive
molar tracers) are separated from the spatial measure they are expressed per
([`Exchange`](@ref) with a [`Measure`](@ref)) and from the finite-volume layer
descriptors. Gas-volume host tendencies are [`MoistAirTendencies`](@ref).

All records are immutable and `isbits` for plain floats, act as scalars in broadcasts,
compare by value and derivative components (`-0.0 == 0.0`, `NaN != NaN`; `isequal`
and `hash` consistent), and support componentwise `+`, `-`, unary `-`, scalar `*`
and `/` and `zero`. They are not collections (no `map`, iteration or `length`). A vector of a
plain and a dual-valued record promotes to a concrete dual-valued element type.

```@contents
Pages = ["quantities_and_geometries.md"]
Depth = 2
```

## Rates and tracers

```@docs
ExchangeRates
Tracers
valid_exchange
```

## Measures and measure-tagged records

```@docs
Measure
LeafArea
StemArea
GroundArea
BulkVolume
Exchange
LeafFlux
StemFlux
FaceFlux
VolumeSource
leaf_flux
stem_flux
face_flux
```

## Layer geometries

```@docs
LayerGeometry
PlantLayerGeometry
layer_geometry
CanopyExchange.area_density
CanopyExchange.thickness
plant_area_density
```
