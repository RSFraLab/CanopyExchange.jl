# Quantities & Layer Geometries

Scalar exchange records in `CanopyExchange` separate conserved chemical/energetic
payloads ([`ExchangeRates`](@ref)) from their spatial integration denominators
([`LeafFlux`](@ref), [`StemFlux`](@ref), [`FaceFlux`](@ref), [`VolumeSource`](@ref),
and [`MoistAirTendencies`](@ref)) and finite-volume layer descriptors
([`LayerGeometry`](@ref) and [`PlantLayerGeometry`](@ref)).

All wrapper types support `eltype`, `convert` and `promote_rule` across compatible
continuous scalars (`Float64` with `ForwardDiff.Dual{Tag,Float64,N}`), `==`,
`isequal`, `hash`, `broadcastable`, and measure-safe linear algebra (`+`, `-`,
unary `-`, scalar `*`, scalar `/`). Flux and tendency records also define `zero`.

```@contents
Pages = ["quantities_and_geometries.md"]
Depth = 2
```

## Optional Molecular Oxygen

```@docs
OxygenContent
NoOxygen
OxygenTransfer
oxygen_value
```

## Conserved Exchange Rates & Validity

```@docs
ExchangeRates
valid_exchange
```

## Measure-Tagged Flux & Source Wrappers

```@docs
AbstractExchangeFlux
LeafFlux
StemFlux
FaceFlux
VolumeSource
MoistAirTendencies
```

## Finite-Volume Layer Geometries

```@docs
AbstractLayerGeometry
LayerGeometry
PlantLayerGeometry
```
