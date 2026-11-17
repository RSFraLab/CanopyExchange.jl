using Test
using CanopyExchange
@isdefined(exchange_fixture) || include("exchange_fixture.jl")

function exchange_allocation_check(leaf, geometry, bottom, top)
    paired = leaf_air_sources(leaf, geometry)
    layer_tendency(bottom, top, paired.air, geometry)
    return @allocated layer_tendency(
        bottom,
        top,
        leaf_air_sources(leaf, geometry).air,
        geometry,
    )
end

@testset "Additional tracer without changing the convenience records" begin
    for FT in (Float32, Float64)
        case = exchange_fixture(FT, false)
        tracer = case.config["additional_tracer"]
        leaves = FT.(tracer["leaf_flux_mol_CH4_m2_leaf_s"])
        faces = FT.(tracer["face_flux_mol_CH4_m2_ground_s"])
        source = leaf_area_source.(leaves, case.geometry)
        donor = -source
        tendency =
            conserved_tendency.(faces[1:(end-1)], faces[2:end], source, case.geometry)
        thickness = getproperty.(case.geometry, :thickness)
        area = thickness .* getproperty.(case.geometry, :leaf_area_density)
        scale = maximum(abs, faces) + sum(abs.(area .* leaves))
        @test isapprox(
            sum(thickness .* tendency),
            faces[1] - faces[end] + sum(area .* leaves);
            rtol = case.rtol,
            atol = case.atol * scale,
        )
        @test isapprox(
            sum(thickness .* (tendency .+ donor)),
            faces[1] - faces[end];
            rtol = case.rtol,
            atol = case.atol * scale,
        )
        @test all(iszero, source .+ donor)
        @test @inferred(leaf_area_source(FT(NaN), case.geometry[1])) === zero(FT)
        @test @inferred(
            conserved_tendency(faces[1], faces[2], source[1], case.geometry[1])
        ) isa FT
        @test source[2] > zero(FT) && source[3] < zero(FT)
        other = FT === Float32 ? Float64 : Float32
        @test_throws ArgumentError conserved_tendency(
            other(faces[1]),
            faces[2],
            source[1],
            case.geometry[1],
        )
    end
end

@testset "Conservative exchange Float32/Float64" begin
    for FT in (Float32, Float64), enabled in (false, true)
        case = exchange_fixture(FT, enabled)
        (; geometry, leaves, boundaries, rtol, atol, config) = case
        layers, faces = config["layers"], config["faces"]
        paired = leaf_air_sources.(leaves, geometry)
        air = getproperty.(paired, :air)
        donor = getproperty.(paired, :leaf)
        tendencies =
            layer_tendency.(boundaries[1:(end-1)], boundaries[2:end], air, geometry)

        for i in eachindex(geometry)
            @test @inferred(leaf_air_sources(leaves[i], geometry[i])) === paired[i]
            @test @inferred(
                layer_tendency(boundaries[i], boundaries[i+1], air[i], geometry[i])
            ) === tendencies[i]
            @test isbitstype(typeof(tendencies[i]))
            @test all(x -> x isa FT, exchange_components(tendencies[i].rates))
            @test valid_exchange(tendencies[i])
            @test exchange_components(air[i].rates) ==
                  map(-, exchange_components(donor[i].rates))
            # Independent source expression in input units, not a second call to core.
            expected = (
                FT(layers["leaf_area_density_m2_m3"][i]) * (
                    FT(layers["sensible_heat_W_m2_leaf"][i]) +
                    FT(layers["water_kg_m2_leaf_s"][i]) *
                    FT(layers["vapor_enthalpy_J_kg"][i])
                ),
                FT(layers["leaf_area_density_m2_m3"][i]) *
                FT(layers["water_kg_m2_leaf_s"][i]),
                FT(layers["leaf_area_density_m2_m3"][i]) *
                FT(layers["carbon_mol_m2_leaf_s"][i]),
            )
            @test exchange_components(air[i].rates)[1:3] == expected
        end

        # Integrate independent densities on a nonuniform grid. Interior faces
        # disappear; the soil boundary enters once even though its layer has no leaves.
        for (j, face_key, leaf_key) in (
            (1, "energy_W_m2_ground", "sensible_heat_W_m2_leaf"),
            (2, "water_kg_m2_ground_s", "water_kg_m2_leaf_s"),
            (3, "carbon_mol_m2_ground_s", "carbon_mol_m2_leaf_s"),
            (4, "oxygen_mol_m2_ground_s", "oxygen_mol_m2_leaf_s"),
        )
            j == 4 && !enabled && continue
            flux = FT.(faces[face_key])
            leaf_input = FT.(layers[leaf_key])
            if j == 1
                leaf_input .+=
                    FT.(layers["water_kg_m2_leaf_s"]) .* FT.(layers["vapor_enthalpy_J_kg"])
            end
            source_integral = sum(
                FT.(layers["thickness_m"]) .* FT.(layers["leaf_area_density_m2_m3"]) .*
                leaf_input,
            )
            integrated = sum(
                geometry[i].thickness * exchange_components(tendencies[i].rates)[j] for
                i in eachindex(geometry)
            )
            expected = flux[1] - flux[end] + source_integral
            scale =
                maximum(abs, flux) + sum(
                    abs.(
                        FT.(layers["thickness_m"]) .*
                        FT.(layers["leaf_area_density_m2_m3"]) .* leaf_input,
                    ),
                )
            @test isapprox(integrated, expected; rtol, atol = atol * scale)
            total = sum(
                geometry[i].thickness * (
                    exchange_components(tendencies[i].rates)[j] +
                    exchange_components(donor[i].rates)[j]
                ) for i in eachindex(geometry)
            )
            @test isapprox(total, flux[1] - flux[end]; rtol, atol = atol * scale)
        end

        # Changing the common water-enthalpy datum changes the paired energy
        # transfer by exactly datum_shift * water; it must not create extra heat.
        shift = FT(config["enthalpy_reference"]["shift_J_kg"])
        for i in eachindex(leaves)
            original = leaves[i].rates
            shifted = sensible_water_flux(
                FT(layers["sensible_heat_W_m2_leaf"][i]),
                original.water,
                FT(layers["vapor_enthalpy_J_kg"][i]) + shift,
                original.tracers,
            )
            difference =
                leaf_air_sources(shifted, geometry[i]).air.rates.energy -
                air[i].rates.energy
            @test isapprox(
                difference,
                shift * air[i].rates.water;
                rtol,
                atol = atol * max(abs(air[i].rates.energy), abs(difference)),
            )
        end

        dummy = LeafFlux(ExchangeRates(FT(NaN), FT(Inf), fixture_tracers(FT, enabled, NaN, NaN)))
        empty = leaf_air_sources(dummy, geometry[1])
        @test all(iszero, exchange_components(empty.air.rates))
        @test all(iszero, exchange_components(empty.leaf.rates))
        @test valid_exchange(empty.air)
        @test !valid_exchange(leaf_air_sources(dummy, geometry[2]).air)
        @test haskey(empty.air.rates.tracers, :o2) == enabled
        @test exchange_allocation_check(
            leaves[2],
            geometry[2],
            boundaries[2],
            boundaries[3],
        ) == config["validation"]["allocation_bytes"]
        @test_throws MethodError leaf_air_sources(boundaries[1], geometry[1])
        @test_throws MethodError layer_tendency(
            leaves[1],
            boundaries[2],
            air[1],
            geometry[1],
        )
        other = FT === Float32 ? Float64 : Float32
        # Mixed precision and mismatched tracer sets are ArgumentErrors (v0.1.x: MethodError).
        @test_throws ArgumentError leaf_air_sources(
            leaves[1],
            LayerGeometry(one(other), zero(other)),
        )
        incompatible = FaceFlux(
            ExchangeRates(zero(FT), zero(FT), fixture_tracers(FT, !enabled, zero(FT), zero(FT))),
        )
        @test_throws ArgumentError layer_tendency(
            incompatible,
            boundaries[2],
            air[1],
            geometry[1],
        )
        @test_throws ArgumentError ExchangeRates(one(FT), one(other), (; co2 = one(FT)))
        for thickness in (zero(FT), -one(FT), FT(NaN), FT(Inf))
            @test_throws ArgumentError LayerGeometry(thickness, zero(FT))
        end
        for density in (-one(FT), FT(NaN), FT(Inf))
            @test_throws ArgumentError LayerGeometry(one(FT), density)
        end
    end
end

@testset "Datum shift, surface flux, wet stem and gas fraction against independent expressions" begin
    for FT in (Float32, Float64), enabled in (false, true)
        case = exchange_fixture(FT, enabled)
        plant = plant_fixture(FT, enabled)
        (; geometry, leaves, boundaries, rtol, atol, config) = case
        layers, faces, stems, surface = config["layers"], config["faces"], config["stems"], config["surface"]
        dz, lad = FT.(layers["thickness_m"]), FT.(layers["leaf_area_density_m2_m3"])
        hv = FT.(layers["vapor_enthalpy_J_kg"])
        leaf_water, face_water = FT.(layers["water_kg_m2_leaf_s"]), FT.(faces["water_kg_m2_ground_s"])
        shift = plant.shift

        # Datum shift on every stream before the divergence changes the column energy
        # budget by exactly shift x column water budget (to rounding).
        bottom, top = boundaries[1:(end-1)], boundaries[2:end]
        air = getproperty.(leaf_air_sources.(leaves, geometry), :air)
        shifted_leaves = shift_water_energy_datum.(leaves, shift)
        shifted_faces = shift_water_energy_datum.(boundaries, shift)
        shifted_air = getproperty.(leaf_air_sources.(shifted_leaves, geometry), :air)
        energy(x) = x.rates.energy
        original = layer_tendency.(bottom, top, air, geometry)
        shifted = layer_tendency.(shifted_faces[1:(end-1)], shifted_faces[2:end], shifted_air, geometry)
        change = sum(dz .* (energy.(shifted) .- energy.(original)))
        water_budget = face_water[1] - face_water[end] + sum(dz .* lad .* leaf_water)
        magnitude = sum(abs.(energy.(boundaries))) + sum(abs.(energy.(shifted_faces))) +
                    sum(dz .* (abs.(energy.(air)) .+ abs.(energy.(shifted_air))))
        @test isapprox(change, shift * water_budget; rtol, atol = atol * magnitude)

        # SurfaceFluxes boundary: energy = shf + lhf, water = evaporation, tracers as given
        soil = plant.soil
        @test soil.rates.energy == FT(surface["sensible_heat_flux_W_m2"]) + FT(surface["latent_heat_flux_W_m2"])
        @test soil.rates.water == FT(surface["evaporation_kg_m2_s"])
        @test soil.rates.tracers.co2 == FT(surface["carbon_mol_m2_s"])
        @test haskey(soil.rates.tracers, :o2) == enabled

        # Transpiring stem: energy = H + E h_v per m² stem, scaled by SAD into the air
        sad = FT.(stems["stem_area_density_m2_m3"])
        stem_heat, stem_water = FT.(stems["sensible_heat_W_m2_stem"]), FT.(stems["water_kg_m2_stem_s"])
        for i in eachindex(plant.wet_stems)
            stem = plant.wet_stems[i]
            @test stem.rates.energy == stem_heat[i] + stem_water[i] * hv[i]
            @test stem.rates.water == stem_water[i]
            pair = plant_air_sources(leaves[i], stem, plant.geometry[i])
            @test pair.stem.rates.energy == -(sad[i] * (stem_heat[i] + stem_water[i] * hv[i]))
            @test pair.stem.rates.water == -(sad[i] * stem_water[i])
        end

        # Gas fraction below one: multiplying back by the gas fraction recovers the bulk tendency
        ε = plant.gas_fraction
        @test ε < one(FT)
        for t in original
            m = moist_air_tendencies(t, ε)
            @test m.rho_total === m.rho_water
            @test isapprox(m.rho_water * ε, t.rates.water; rtol)
            @test isapprox(m.rho_energy * ε, t.rates.energy; rtol)
            @test all(isapprox.(values(m.tracers) .* ε, values(t.rates.tracers); rtol))
        end
    end
end
