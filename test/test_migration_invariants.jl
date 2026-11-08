using Test, TOML, CanopyExchange

@testset "Independent physical-measure invariants" begin
    config = TOML.parsefile(joinpath(@__DIR__, "parameters", "migration.toml"))
    for FT in (Float32, Float64)
        area = FT(config["geometry"]["footprint_m2"])
        dz = FT.(config["geometry"]["thickness_m"])
        lad = FT(config["geometry"]["leaf_area_density_m2_m3"])
        transfer = config["transfer"]
        pulse = FT(transfer["face_pulse"])
        rtol = FT(config["validation"]["relative_tolerance_eps"]) * eps(FT)
        geometry = LayerGeometry.(dz, lad)
        volume = area .* dz
        z = zero(FT)
        for face in 1:3
            faces = zeros(FT, 3)
            faces[face] = pulse
            tendencies = conserved_tendency.(faces[1:2], faces[2:3], z, geometry)
            inventory = volume .* tendencies
            expected = face == 1 ? (area * pulse, z) :
                       face == 2 ? (-area * pulse, area * pulse) :
                                   (z, -area * pulse)
            @test all(isapprox.(inventory, expected; rtol, atol = z))
            @test isapprox(sum(inventory), area * (faces[1] - faces[end]);
                           rtol, atol = rtol * area * pulse)
        end

        # Refining a uniform leaf population preserves physical transfer totals.
        leaf_rate = FT(transfer["leaf_energy_W_m2_leaf"])
        faces = FT.(transfer["face_energy_W_m2_ground"])
        split_source = leaf_area_source.(leaf_rate, geometry)
        split_tendency = conserved_tendency.(faces[1:2], faces[2:3], split_source, geometry)
        whole = LayerGeometry(sum(dz), lad)
        whole_source = leaf_area_source(leaf_rate, whole)
        whole_tendency = conserved_tendency(first(faces), last(faces), whole_source, whole)
        @test isapprox(sum(volume .* split_source), area * sum(dz) * lad * leaf_rate; rtol)
        @test isapprox(sum(volume .* split_tendency), area * sum(dz) * whole_tendency; rtol)

        # Shift energy on every water-carrying face AND leaf stream, not only leaves.
        water = FT(transfer["leaf_water_kg_m2_leaf_s"])
        water_faces = FT.(transfer["face_water_kg_m2_ground_s"])
        datum = FT(transfer["energy_datum_shift_J_kg"])
        energy_tendency = conserved_tendency.(faces[1:2], faces[2:3], split_source, geometry)
        shifted_faces = faces .+ datum .* water_faces
        shifted_source = leaf_area_source.(leaf_rate + datum * water, geometry)
        shifted_tendency = conserved_tendency.(shifted_faces[1:2], shifted_faces[2:3], shifted_source, geometry)
        water_tendency = (water_faces[1:2] .- water_faces[2:3]) ./ dz .+ lad * water
        # This normal-range manufactured case subtracts large energy terms to
        # recover a small datum shift. Bound rounded arithmetic using INPUT
        # magnitudes, including cancellation inside the divergence and source.
        n = FT(config["validation"]["roundoff_operation_budget"])
        gamma = n * eps(FT) / (one(FT) - n * eps(FT))
        magnitude = abs.(faces[1:2]) .+ abs.(faces[2:3]) .+
                    abs.(shifted_faces[1:2]) .+ abs.(shifted_faces[2:3]) .+
                    dz .* (abs.(split_source) .+ abs.(shifted_source))
        water_magnitude = abs(datum) .* (abs.(water_faces[1:2]) .+
                          abs.(water_faces[2:3]) .+ dz .* lad .* abs(water))
        local_bound = gamma .* (magnitude .+ water_magnitude) ./ dz
        @test all(abs.((shifted_tendency .- energy_tendency) .- datum .* water_tendency) .<= local_bound)
        expected_inventory_shift = datum * area * (first(water_faces) - last(water_faces) + sum(dz) * lad * water)
        independent_magnitude = area * abs(datum) * (abs(first(water_faces)) +
            abs(last(water_faces)) + sum(dz) * lad * abs(water))
        bound = gamma * (area * sum(magnitude) + independent_magnitude)
        @test bound < FT(config["validation"]["maximum_roundoff_fraction_of_signal"]) * abs(expected_inventory_shift)
        @test abs(sum(volume .* (shifted_tendency .- energy_tendency)) - expected_inventory_shift) <= bound

        carbon = FT(transfer["leaf_carbon_mol_m2_leaf_s"])
        absent = ExchangeRates(leaf_rate, water, carbon, NoOxygen{FT}())
        computed = ExchangeRates(leaf_rate, water, carbon, OxygenTransfer(z))
        absent_pair = leaf_air_sources(LeafFlux(absent), whole)
        computed_pair = leaf_air_sources(LeafFlux(computed), whole)
        @test typeof(absent_pair) !== typeof(computed_pair)
        @test oxygen_value(absent_pair.air.rates.oxygen) === nothing
        @test oxygen_value(computed_pair.air.rates.oxygen) === z
        for name in (:energy, :water, :carbon)
            @test getproperty(absent_pair.air.rates, name) === getproperty(computed_pair.air.rates, name)
        end
        @test_throws MethodError layer_tendency(FaceFlux(absent), FaceFlux(computed), absent_pair.air, whole)
        empty = LayerGeometry(sum(dz), z)
        for invalid in (FT(NaN), FT(Inf))
            bad = ExchangeRates(invalid, water, carbon, NoOxygen{FT}())
            @test !valid_exchange(leaf_air_sources(LeafFlux(bad), whole).air)
            empty_source = leaf_air_sources(LeafFlux(bad), empty).air
            @test valid_exchange(empty_source)
            @test !valid_exchange(layer_tendency(FaceFlux(bad), FaceFlux(absent), empty_source, empty))
        end
        large_geometry = LayerGeometry(one(FT), floatmax(FT))
        @test !valid_exchange(leaf_air_sources(LeafFlux(absent), large_geometry).air)
    end
end
