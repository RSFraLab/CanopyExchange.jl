# Host setup shared by CPU and optional CUDA verification; no physics lives here.
# Record construction and component access for the TOML fixture live here, so the
# record files (write_arithmetic_record.jl, reference_record.jl) only call operations.
using TOML

fixture_config() =
    TOML.parsefile(joinpath(@__DIR__, "parameters", "conservative_exchange.toml"))

# Tracers of the fixture: CO2, plus O2 when `enabled` (an uncomputed tracer is not a key).
fixture_tracers(::Type{FT}, enabled::Bool, carbon, oxygen) where {FT} =
    enabled ? Tracers(:co2, :o2)(FT(carbon), FT(oxygen)) : Tracers(:co2)(FT(carbon))

function exchange_fixture(::Type{FT}, enabled::Bool) where {FT}
    config = fixture_config()
    layers, faces = config["layers"], config["faces"]
    geometry =
        LayerGeometry.(FT.(layers["thickness_m"]), FT.(layers["leaf_area_density_m2_m3"]))
    leaves = [
        sensible_water_flux(
            FT(layers["sensible_heat_W_m2_leaf"][i]),
            FT(layers["water_kg_m2_leaf_s"][i]),
            FT(layers["vapor_enthalpy_J_kg"][i]),
            fixture_tracers(
                FT,
                enabled,
                layers["carbon_mol_m2_leaf_s"][i],
                layers["oxygen_mol_m2_leaf_s"][i],
            ),
        ) for i in eachindex(geometry)
    ]
    boundaries = [
        face_flux(
            FT(faces["energy_W_m2_ground"][i]),
            FT(faces["water_kg_m2_ground_s"][i]),
            fixture_tracers(
                FT,
                enabled,
                faces["carbon_mol_m2_ground_s"][i],
                faces["oxygen_mol_m2_ground_s"][i],
            ),
        ) for i in eachindex(faces["energy_W_m2_ground"])
    ]
    precision = FT === Float32 ? "float32" : "float64"
    return (;
        config,
        geometry,
        leaves,
        boundaries,
        rtol = FT(config["validation"][precision*"_relative_tolerance"]),
        atol = FT(config["validation"][precision*"_scaled_absolute_tolerance"]),
    )
end

# (energy, water, carbon[, oxygen]) in the order of the v0.1.x records.
exchange_components(rates::ExchangeRates) = (rates.energy, rates.water, values(rates.tracers)...)

# Stem, soil-surface and gas-fraction inputs for the frozen reference record.
function plant_fixture(::Type{FT}, enabled::Bool) where {FT}
    config = fixture_config()
    layers, stems, surface = config["layers"], config["stems"], config["surface"]
    geometry = PlantLayerGeometry.(
        FT.(layers["thickness_m"]),
        FT.(layers["leaf_area_density_m2_m3"]),
        FT.(stems["stem_area_density_m2_m3"]),
    )
    # A dry stem has zero water and, as in v0.1.x, an explicitly zero CO2 transfer.
    dry_stems = [
        sensible_stem_flux(
            FT(stems["sensible_heat_W_m2_stem"][i]),
            fixture_tracers(FT, enabled, zero(FT), stems["oxygen_mol_m2_stem_s"][i]),
        ) for i in eachindex(geometry)
    ]
    wet_stems = [
        sensible_stem_flux(
            FT(stems["sensible_heat_W_m2_stem"][i]),
            FT(stems["water_kg_m2_stem_s"][i]),
            FT(layers["vapor_enthalpy_J_kg"][i]),
            fixture_tracers(
                FT,
                enabled,
                stems["carbon_mol_m2_stem_s"][i],
                stems["oxygen_mol_m2_stem_s"][i],
            ),
        ) for i in eachindex(geometry)
    ]
    soil = surface_static_energy_flux(
        FT(surface["sensible_heat_flux_W_m2"]),
        FT(surface["latent_heat_flux_W_m2"]),
        FT(surface["evaporation_kg_m2_s"]),
        fixture_tracers(FT, enabled, surface["carbon_mol_m2_s"], surface["oxygen_mol_m2_s"]),
    )
    return (;
        geometry,
        dry_stems,
        wet_stems,
        soil,
        shift = FT(config["enthalpy_reference"]["shift_J_kg"]),
        gas_fraction = FT(config["air"]["gas_fraction_m3_gas_m3_bulk"]),
    )
end

# (rho_total, rho_water, rho_energy, carbon[, oxygen]) in the order of the v0.1.x records.
moist_components(m::MoistAirTendencies) =
    (m.rho_total, m.rho_water, m.rho_energy, values(m.tracers)...)
