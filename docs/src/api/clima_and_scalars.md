# CliMA Contracts & Continuous-Scalar Traits

Boundary, column-accumulation, and porous-medium tendency contracts bridge
`CanopyExchange` with `SurfaceFluxes.jl`, `Thermodynamics.jl`, `ClimaLand.jl`,
`ClimaAtmos.jl`, and `ClimaCore.jl` without introducing runtime package dependencies.
Semi-public continuous-scalar trait functions (`public` on Julia 1.11+) govern
precision safety, `ForwardDiff.Dual` promotion, and structural leaf-float predicates.

```@contents
Pages = ["clima_and_scalars.md"]
Depth = 2
```

## Boundary Static-Energy & Datum Shifts

```@docs
surface_static_energy_flux
shift_water_energy_datum
```

## Diagnostic Steady-Air Flux Profiles

```@docs
layer_flux_increment
steady_top_flux
accumulate_steady_fluxes!
accumulate_steady_fluxes
```

## Conserved Gas-Volume Moist-Air Tendencies

```@docs
moist_air_tendencies
```

## Continuous-Scalar & Oxygen Traits (`public`)

```@docs
CanopyExchange.primal_float_type
CanopyExchange.primal_scalar
CanopyExchange.is_continuous_scalar
CanopyExchange.compatible_scalar
CanopyExchange.compatible_oxygen
CanopyExchange.finite_scalar
CanopyExchange.exact_zero_scalar
```
