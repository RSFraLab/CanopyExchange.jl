# Operations covered by the frozen reference record in test/reference/.
# Record construction and component access live in exchange_fixture.jl; this file
# only calls the operations, so it must stay valid across API revisions.
bit_record(values) = collect(bitstring.(values))
exchange_record(xs) = [bit_record(exchange_components(x.rates)) for x in xs]
moist_record(xs) = [bit_record(moist_components(x)) for x in xs]

function reference_profiles()
    profiles = Dict{String,Any}[]
    for FT in (Float32, Float64), enabled in (false, true)
        case = exchange_fixture(FT, enabled)
        plant = plant_fixture(FT, enabled)
        bottom, top = case.boundaries[1:(end-1)], case.boundaries[2:end]
        paired = leaf_air_sources.(case.leaves, case.geometry)
        air = getproperty.(paired, :air)
        tendency = layer_tendency.(bottom, top, air, case.geometry)
        dry = plant_air_sources.(case.leaves, plant.dry_stems, plant.geometry)
        wet = plant_air_sources.(case.leaves, plant.wet_stems, plant.geometry)
        wet_air = getproperty.(wet, :air)
        push!(
            profiles,
            Dict(
                "precision" => string(FT),
                "oxygen" => enabled,
                "leaf_flux" => exchange_record(case.leaves),
                "leaf_air_source" => exchange_record(air),
                "leaf_donor_source" => exchange_record(getproperty.(paired, :leaf)),
                "layer_tendency" => exchange_record(tendency),
                "dry_stem_flux" => exchange_record(plant.dry_stems),
                "dry_plant_air_source" => exchange_record(getproperty.(dry, :air)),
                "dry_plant_leaf_source" => exchange_record(getproperty.(dry, :leaf)),
                "dry_plant_stem_source" => exchange_record(getproperty.(dry, :stem)),
                "wet_stem_flux" => exchange_record(plant.wet_stems),
                "wet_plant_air_source" => exchange_record(wet_air),
                "wet_plant_leaf_source" => exchange_record(getproperty.(wet, :leaf)),
                "wet_plant_stem_source" => exchange_record(getproperty.(wet, :stem)),
                "plant_layer_tendency" =>
                    exchange_record(layer_tendency.(bottom, top, wet_air, plant.geometry)),
                "shifted_leaf_flux" =>
                    exchange_record(shift_water_energy_datum.(case.leaves, plant.shift)),
                "shifted_stem_flux" =>
                    exchange_record(shift_water_energy_datum.(plant.wet_stems, plant.shift)),
                "shifted_face_flux" =>
                    exchange_record(shift_water_energy_datum.(case.boundaries, plant.shift)),
                "shifted_volume_source" =>
                    exchange_record(shift_water_energy_datum.(wet_air, plant.shift)),
                "steady_top_flux" =>
                    exchange_record(steady_top_flux.(bottom, air, case.geometry)),
                "accumulated_steady_fluxes" => exchange_record(
                    accumulate_steady_fluxes(plant.soil, wet_air, plant.geometry),
                ),
                "surface_static_energy_flux" => exchange_record([plant.soil]),
                "moist_air_unit_gas_fraction" =>
                    moist_record(moist_air_tendencies.(tendency, one(FT))),
                "moist_air_gas_fraction" =>
                    moist_record(moist_air_tendencies.(tendency, plant.gas_fraction)),
            ),
        )
    end
    return profiles
end
