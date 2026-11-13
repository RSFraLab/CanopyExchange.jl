"""
    CanopyExchange

Provide zero-dependency, measure-safe scalar exchange contracts for layered plant
canopies and atmospheric / land host models (including `ClimaLand`, `SurfaceFluxes`,
and `ClimaAtmos`).

`CanopyExchange` converts one-sided leaf (`LAD`) and stem (`SAD`) fluxes into
paired bulk-volume sources, evaluates conservative finite-volume vertical flux
divergences, bridges diagnostic steady-air flux profiles with prognostic air
columns, and maps bulk sources to conserved gas-volume moist-air tendencies across
`Float32`, `Float64`, GPU arrays (`CuArray`), `ClimaCore.Fields.Field` broadcasting,
and `ForwardDiff` automatic differentiation. It owns no state, grid, turbulence
closure, equation of state, or time integrator.
"""
module CanopyExchange

import Base: ==, +, -, *, /, convert, eltype, hash, isequal, zero
import Base.Broadcast: broadcastable

export OxygenContent,
    NoOxygen,
    OxygenTransfer,
    oxygen_value,
    ExchangeRates,
    AbstractExchangeFlux,
    LeafFlux,
    StemFlux,
    FaceFlux,
    VolumeSource,
    MoistAirTendencies,
    AbstractLayerGeometry,
    LayerGeometry,
    PlantLayerGeometry,
    area_source,
    leaf_area_source,
    stem_area_source,
    plant_area_density,
    sensible_water_flux,
    sensible_stem_flux,
    leaf_air_sources,
    plant_air_sources,
    conserved_tendency,
    layer_tendency,
    surface_static_energy_flux,
    shift_water_energy_datum,
    layer_flux_increment,
    steady_top_flux,
    accumulate_steady_fluxes!,
    accumulate_steady_fluxes,
    moist_air_tendencies,
    valid_exchange

@static if VERSION >= v"1.11"
    eval(
        Expr(
            :public,
            :primal_float_type,
            :primal_scalar,
            :is_continuous_scalar,
            :compatible_scalar,
            :compatible_oxygen,
            :finite_scalar,
            :exact_zero_scalar,
        ),
    )
end

include("scalars.jl")
include("quantities.jl")
include("sources.jl")
include("clima_contracts.jl")

end
