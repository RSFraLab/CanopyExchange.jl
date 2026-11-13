# Optional CUDA check; CUDA belongs to the caller's verification environment.
using CUDA, Test
using CanopyExchange
include("exchange_fixture.jl")

CUDA.functional() || error("CUDA is required for this verification")
CUDA.allowscalar(false)

@testset "Conservative exchange CUDA Float32/Float64" begin
    for FT in (Float32, Float64), enabled in (false, true)
        case = exchange_fixture(FT, enabled)
        tracer = case.config["additional_tracer"]
        tracer_leaf = FT.(tracer["leaf_flux_mol_CH4_m2_leaf_s"])
        tracer_face = FT.(tracer["face_flux_mol_CH4_m2_ground_s"])
        tracer_leaf[1] = FT(NaN)
        tracer_source = leaf_area_source.(CuArray(tracer_leaf), CuArray(case.geometry))
        tracer_tendency = Array(
            conserved_tendency.(
                CuArray(tracer_face[1:(end-1)]),
                CuArray(tracer_face[2:end]),
                tracer_source,
                CuArray(case.geometry),
            ),
        )
        tracer_expected_source = leaf_area_source.(tracer_leaf, case.geometry)
        tracer_expected = conserved_tendency.(
            tracer_face[1:(end-1)],
            tracer_face[2:end],
            tracer_expected_source,
            case.geometry,
        )
        @test all(isfinite, tracer_tendency)
        @test all(
            isapprox.(tracer_tendency, tracer_expected; rtol = case.rtol, atol = zero(FT)),
        )
        layers = case.config["layers"]
        device_heat = CuArray(FT.(layers["sensible_heat_W_m2_leaf"]))
        device_water = CuArray(FT.(layers["water_kg_m2_leaf_s"]))
        device_carbon = CuArray(FT.(layers["carbon_mol_m2_leaf_s"]))
        device_oxygen = CuArray([leaf.rates.oxygen for leaf in case.leaves])
        device_enthalpy = CuArray(FT.(layers["vapor_enthalpy_J_kg"]))
        constructed = Array(
            sensible_water_flux.(
                device_heat,
                device_water,
                device_carbon,
                device_oxygen,
                device_enthalpy,
            ),
        )
        for i in eachindex(constructed)
            for (received, expected) in zip(
                exchange_components(constructed[i].rates),
                exchange_components(case.leaves[i].rates),
            )
                @test isapprox(
                    received,
                    expected;
                    rtol = case.rtol,
                    atol = case.atol * abs(expected),
                )
            end
        end
        count = case.config["validation"]["gpu_repetitions"]
        leaves = repeat(case.leaves, count)
        geometry = repeat(case.geometry, count)
        bottom = repeat(case.boundaries[1:(end-1)], count)
        top = repeat(case.boundaries[2:end], count)
        dummy_oxygen = enabled ? OxygenTransfer(FT(NaN)) : NoOxygen{FT}()
        dummy = LeafFlux(ExchangeRates(FT(NaN), FT(Inf), FT(NaN), dummy_oxygen))
        for i in eachindex(leaves)
            iszero(geometry[i].leaf_area_density) && (leaves[i] = dummy)
        end
        reference_pairs = leaf_air_sources.(leaves, geometry)
        reference_air = getproperty.(reference_pairs, :air)
        reference = layer_tendency.(bottom, top, reference_air, geometry)
        device_leaves, device_geometry = CuArray(leaves), CuArray(geometry)
        device_bottom, device_top = CuArray(bottom), CuArray(top)
        device_pairs = leaf_air_sources.(device_leaves, device_geometry)
        device_air = (pair -> pair.air).(device_pairs)
        device_leaf = (pair -> pair.leaf).(device_pairs)
        device_tendency =
            layer_tendency.(device_bottom, device_top, device_air, device_geometry)
        @test all(Array(valid_exchange.(device_tendency)))
        @test all(Array(valid_exchange.(device_leaf)))
        actual = Array(device_tendency)
        @test eltype(actual) === eltype(reference)
        @test Array(device_leaf) == getproperty.(reference_pairs, :leaf)
        for i in eachindex(reference)
            for (received, expected) in zip(
                exchange_components(actual[i].rates),
                exchange_components(reference[i].rates),
            )
                @test isapprox(
                    received,
                    expected;
                    rtol = case.rtol,
                    atol = case.atol * abs(expected),
                )
            end
        end
        @test (oxygen_value(actual[1].rates.oxygen) === nothing) == !enabled

        # Plant (LAD + SAD), datum shift, steady flux, and moist-air CUDA broadcast
        device_pg = CuArray(
            PlantLayerGeometry.(
                getproperty.(geometry, :thickness),
                getproperty.(geometry, :leaf_area_density),
                FT(0.2),
            ),
        )
        stem_o2 = enabled ? OxygenTransfer(zero(FT)) : NoOxygen{FT}()
        device_stems = CuArray(fill(sensible_stem_flux(FT(15), stem_o2), length(leaves)))
        device_plant = plant_air_sources.(device_leaves, device_stems, device_pg)
        device_plant_air = (p -> p.air).(device_plant)
        device_shifted = shift_water_energy_datum.(device_bottom, FT(500))
        device_steady_top = steady_top_flux.(device_shifted, device_plant_air, device_pg)
        device_plant_tend =
            layer_tendency.(device_shifted, device_steady_top, device_plant_air, device_pg)
        device_moist = moist_air_tendencies.(device_plant_tend, FT(0.92))
        host_moist = Array(device_moist)
        @test eltype(host_moist) === MoistAirTendencies{FT,typeof(stem_o2)}
        @test all(valid_exchange, host_moist)
        @test maximum(m -> abs(m.rho_energy), host_moist) <= FT(100) * eps(FT) * FT(1e4)
    end
end
