# Host setup shared by CPU and optional CUDA verification; no physics lives here.
using TOML

function exchange_fixture(::Type{FT}, enabled::Bool) where {FT}
    config = TOML.parsefile(joinpath(@__DIR__, "parameters", "conservative_exchange.toml"))
    layers, faces = config["layers"], config["faces"]
    oxygen(value) = enabled ? OxygenTransfer(FT(value)) : NoOxygen{FT}()
    geometry =
        LayerGeometry.(FT.(layers["thickness_m"]), FT.(layers["leaf_area_density_m2_m3"]))
    leaves = [
        sensible_water_flux(
            FT(layers["sensible_heat_W_m2_leaf"][i]),
            FT(layers["water_kg_m2_leaf_s"][i]),
            FT(layers["carbon_mol_m2_leaf_s"][i]),
            oxygen(layers["oxygen_mol_m2_leaf_s"][i]),
            FT(layers["vapor_enthalpy_J_kg"][i]),
        ) for i in eachindex(geometry)
    ]
    boundaries = [
        FaceFlux(
            ExchangeRates(
                FT(faces["energy_W_m2_ground"][i]),
                FT(faces["water_kg_m2_ground_s"][i]),
                FT(faces["carbon_mol_m2_ground_s"][i]),
                oxygen(faces["oxygen_mol_m2_ground_s"][i]),
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

function exchange_components(rates::ExchangeRates)
    oxygen = oxygen_value(rates.oxygen)
    return oxygen === nothing ? (rates.energy, rates.water, rates.carbon) :
           (rates.energy, rates.water, rates.carbon, oxygen)
end

# Stem, soil-surface and gas-fraction inputs for the frozen reference record.
function plant_fixture(::Type{FT}, enabled::Bool) where {FT}
    config = TOML.parsefile(joinpath(@__DIR__, "parameters", "conservative_exchange.toml"))
    layers, stems, surface = config["layers"], config["stems"], config["surface"]
    oxygen(value) = enabled ? OxygenTransfer(FT(value)) : NoOxygen{FT}()
    geometry = PlantLayerGeometry.(
        FT.(layers["thickness_m"]),
        FT.(layers["leaf_area_density_m2_m3"]),
        FT.(stems["stem_area_density_m2_m3"]),
    )
    dry_stems = [
        sensible_stem_flux(
            FT(stems["sensible_heat_W_m2_stem"][i]),
            oxygen(stems["oxygen_mol_m2_stem_s"][i]),
        ) for i in eachindex(geometry)
    ]
    wet_stems = [
        sensible_stem_flux(
            FT(stems["sensible_heat_W_m2_stem"][i]),
            FT(stems["water_kg_m2_stem_s"][i]),
            FT(stems["carbon_mol_m2_stem_s"][i]),
            oxygen(stems["oxygen_mol_m2_stem_s"][i]),
            FT(layers["vapor_enthalpy_J_kg"][i]),
        ) for i in eachindex(geometry)
    ]
    soil = surface_static_energy_flux(
        FT(surface["sensible_heat_flux_W_m2"]),
        FT(surface["latent_heat_flux_W_m2"]),
        FT(surface["evaporation_kg_m2_s"]),
        FT(surface["carbon_mol_m2_s"]),
        oxygen(surface["oxygen_mol_m2_s"]),
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

function moist_components(m::MoistAirTendencies)
    oxygen = oxygen_value(m.oxygen)
    return oxygen === nothing ? (m.rho_total, m.rho_water, m.rho_energy, m.carbon) :
           (m.rho_total, m.rho_water, m.rho_energy, m.carbon, oxygen)
end
