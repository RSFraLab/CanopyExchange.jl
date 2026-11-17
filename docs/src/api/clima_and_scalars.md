# Boundary Fluxes, Steady Columns, Moist Air & Scalar Hooks

Constructors that couple energy and water (transported vapour enthalpy, the
SurfaceFluxes static-energy convention, the energy datum of water), diagnostic
steady-air flux profiles, the gas-volume host tendencies, and the four `public` scalar
hooks through which number types with derivative information take part.

```@contents
Pages = ["clima_and_scalars.md"]
Depth = 2
```

## Element and boundary fluxes

```@docs
sensible_water_flux
sensible_stem_flux
surface_static_energy_flux
shift_water_energy_datum
```

## Diagnostic steady-air flux profiles

```@docs
layer_flux_increment
steady_top_flux
accumulate_steady_fluxes!
accumulate_steady_fluxes
```

## Gas-volume moist-air tendencies

```@docs
MoistAirTendencies
moist_air_tendencies
```

## Scalar hooks (`public`)

The defaults are exact for plain floating-point numbers. The package extension
`CanopyExchangeForwardDiffExt` (loaded with ForwardDiff 0.10 or 1.x) defines them for
`ForwardDiff.Dual`; another number type with derivative information takes part by
defining the same five methods.

```@docs
CanopyExchange.primal_type
CanopyExchange.primal
CanopyExchange.all_finite
CanopyExchange.exact_zero
CanopyExchange.scalar_equal
```
