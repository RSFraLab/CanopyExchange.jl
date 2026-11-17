"""
    CanopyExchange

Measure-safe, conservative scalar exchange contracts between layered plant canopies
and atmospheric or land host models (ClimaLand, SurfaceFluxes, ClimaAtmos and
standalone finite-volume columns), with no runtime dependencies.

Rates of energy (including transported vapour enthalpy), water and an open set of
passive molar tracers ([`ExchangeRates`](@ref)) are tagged with the measure they are
expressed per ([`Exchange`](@ref): leaf area, stem area, ground area, bulk volume).
The package converts element fluxes to paired bulk-volume sources, evaluates
conservative finite-volume flux divergences, bridges diagnostic steady-air profiles
and prognostic columns, and maps bulk sources to gas-volume moist-air tendencies. It
works with `Float32` and `Float64`, inside broadcasts over GPU arrays and ClimaCore
fields, and with ForwardDiff dual numbers (through a package extension). It owns no
state, grid, turbulence closure, equation of state or time integrator.

Source files, one per process: `precision.jl` (scalar hooks and the precision
policy), `records.jl` (equality, hashing and broadcasting of every record),
`rates.jl`, `measures.jl`, `geometry.jl`, `area_to_volume.jl`, `flux_divergence.jl`,
`water_energy.jl` (vapour enthalpy, surface static energy, energy datum) and
`moist_air.jl`.
"""
module CanopyExchange

# Records and measures
export ExchangeRates, Tracers, MoistAirTendencies
export Measure, LeafArea, StemArea, GroundArea, BulkVolume
export Exchange, LeafFlux, StemFlux, FaceFlux, VolumeSource
export leaf_flux, stem_flux, face_flux, valid_exchange

# Layer geometries (the interface functions `thickness` and `area_density` are public)
export LayerGeometry, PlantLayerGeometry, layer_geometry, plant_area_density

# Area to volume and paired ledgers
export area_source, volume_source, leaf_area_source, stem_area_source
export leaf_air_sources, plant_air_sources

# Flux divergence and steady columns
export conserved_tendency, layer_tendency, layer_flux_increment, steady_top_flux
export accumulate_steady_fluxes!, accumulate_steady_fluxes

# Element and boundary fluxes, energy datum, moist air
export sensible_water_flux, sensible_stem_flux, surface_static_energy_flux
export shift_water_energy_datum, moist_air_tendencies

# Interfaces extended by other code: the geometry interface and the scalar hooks for
# number types with derivative information (public, not exported)
@static if VERSION >= v"1.11"
    eval(
        Expr(
            :public,
            :thickness,
            :area_density,
            :primal_type,
            :primal,
            :all_finite,
            :exact_zero,
            :scalar_equal,
        ),
    )
end

include("precision.jl")
include("records.jl")
include("rates.jl")
include("measures.jl")
include("geometry.jl")
include("area_to_volume.jl")
include("flux_divergence.jl")
include("water_energy.jl")
include("moist_air.jl")

end
