"""
    CanopyExchange

Scalar exchange contracts for a layered canopy. This module owns conversion of
one-sided leaf-area fluxes to paired volume sources and conservative face-flux
divergence. It owns no state, grid, thermodynamics, turbulence or time integrator.
Host adapters own state layout, thermodynamics and time integration.
"""
module CanopyExchange

export OxygenContent,
    NoOxygen,
    OxygenTransfer,
    oxygen_value,
    ExchangeRates,
    LeafFlux,
    FaceFlux,
    VolumeSource,
    LayerGeometry,
    leaf_air_sources,
    layer_tendency,
    leaf_area_source,
    conserved_tendency,
    sensible_water_flux,
    valid_exchange

include("quantities.jl")
include("sources.jl")

end
