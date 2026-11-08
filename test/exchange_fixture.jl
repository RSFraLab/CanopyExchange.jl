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
