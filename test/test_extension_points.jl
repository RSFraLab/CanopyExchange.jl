# Extension without touching the package: a second tracer set through the tagged
# records, and a third element class with its own geometry defined only here.
using Test, CanopyExchange
@isdefined(exchange_fixture) || include("exchange_fixture.jl")

struct FruitArea <: CanopyExchange.Measure end

struct OrchardGeometry{FT}
    thickness::FT
    leaf_area_density::FT
    fruit_area_density::FT
end
Base.eltype(::Type{OrchardGeometry{FT}}) where {FT} = FT
CanopyExchange.area_density(g::OrchardGeometry, ::Type{LeafArea}) = g.leaf_area_density
CanopyExchange.area_density(g::OrchardGeometry, ::Type{FruitArea}) = g.fruit_area_density

@testset "Second tracer set (co2, cos, ch4) through the tagged records" begin
    for FT in (Float32, Float64)
        case = exchange_fixture(FT, false)
        config = case.config
        layers, faces = config["layers"], config["faces"]
        ch4, cos = config["additional_tracer"], config["second_tracer_set"]
        T3 = Tracers(:co2, :cos, :ch4)
        leaf_tracers = T3.(
            FT.(layers["carbon_mol_m2_leaf_s"]),
            FT.(cos["leaf_flux_mol_COS_m2_leaf_s"]),
            FT.(ch4["leaf_flux_mol_CH4_m2_leaf_s"]),
        )
        face_tracers = T3.(
            FT.(faces["carbon_mol_m2_ground_s"]),
            FT.(cos["face_flux_mol_COS_m2_ground_s"]),
            FT.(ch4["face_flux_mol_CH4_m2_ground_s"]),
        )
        leaves = sensible_water_flux.(
            FT.(layers["sensible_heat_W_m2_leaf"]),
            FT.(layers["water_kg_m2_leaf_s"]),
            FT.(layers["vapor_enthalpy_J_kg"]),
            leaf_tracers,
        )
        boundaries = face_flux.(
            FT.(faces["energy_W_m2_ground"]),
            FT.(faces["water_kg_m2_ground_s"]),
            face_tracers,
        )
        geometry = case.geometry
        paired = leaf_air_sources.(leaves, geometry)
        air = getproperty.(paired, :air)
        tendency = layer_tendency.(boundaries[1:(end-1)], boundaries[2:end], air, geometry)
        # The co2-only records of the fixture are unchanged by the extra tracers.
        reference_air = getproperty.(leaf_air_sources.(case.leaves, geometry), :air)
        reference = layer_tendency.(
            case.boundaries[1:(end-1)],
            case.boundaries[2:end],
            reference_air,
            geometry,
        )
        for i in eachindex(geometry)
            @test keys(tendency[i].rates.tracers) == (:co2, :cos, :ch4)
            @test eltype(tendency[i]) === FT && isbitstype(typeof(tendency[i]))
            @test @inferred(layer_tendency(boundaries[i], boundaries[i+1], air[i], geometry[i])) ===
                  tendency[i]
            @test exchange_components(tendency[i].rates)[1:3] ===
                  exchange_components(reference[i].rates)
            # Each tracer equals the scalar path bit for bit.
            for name in (:cos, :ch4)
                source = leaf_area_source(leaves[i].rates.tracers[name], geometry[i])
                @test air[i].rates.tracers[name] === source
                @test paired[i].leaf.rates.tracers[name] === -source
                @test tendency[i].rates.tracers[name] === conserved_tendency(
                    boundaries[i].rates.tracers[name],
                    boundaries[i+1].rates.tracers[name],
                    source,
                    geometry[i],
                )
            end
            moist = moist_air_tendencies(tendency[i], FT(config["air"]["gas_fraction_m3_gas_m3_bulk"]))
            @test keys(moist.tracers) == (:co2, :cos, :ch4)
            @test moist.rho_total === moist.rho_water
        end
        # Tracers are passive: the water-energy datum shift leaves them unchanged.
        shifted = shift_water_energy_datum.(leaves, FT(config["enthalpy_reference"]["shift_J_kg"]))
        @test getproperty.(getproperty.(shifted, :rates), :tracers) == leaf_tracers
        # A record with a different tracer set cannot be combined.
        @test_throws ArgumentError air[2] + reference_air[2]
    end
end

@testset "Third element class and geometry defined outside the package" begin
    for FT in (Float32, Float64)
        config = fixture_config()
        layers, faces, fruit = config["layers"], config["faces"], config["third_element_class"]
        dz = FT.(layers["thickness_m"])
        geometry = OrchardGeometry.(dz, FT.(layers["leaf_area_density_m2_m3"]),
            FT.(fruit["fruit_area_density_m2_m3"]))
        T = Tracers(:co2)
        leaves = sensible_water_flux.(
            FT.(layers["sensible_heat_W_m2_leaf"]),
            FT.(layers["water_kg_m2_leaf_s"]),
            FT.(layers["vapor_enthalpy_J_kg"]),
            T.(FT.(layers["carbon_mol_m2_leaf_s"])),
        )
        fruits = Exchange{FruitArea}.(
            ExchangeRates.(
                FT.(fruit["sensible_heat_W_m2_fruit"]) .+
                FT.(fruit["water_kg_m2_fruit_s"]) .* FT.(layers["vapor_enthalpy_J_kg"]),
                FT.(fruit["water_kg_m2_fruit_s"]),
                T.(FT.(fruit["carbon_mol_m2_fruit_s"])),
            ),
        )
        boundaries = face_flux.(
            FT.(faces["energy_W_m2_ground"]),
            FT.(faces["water_kg_m2_ground_s"]),
            T.(FT.(faces["carbon_mol_m2_ground_s"])),
        )
        from_leaf = volume_source.(leaves, geometry)
        from_fruit = volume_source.(fruits, geometry)
        air = from_leaf .+ from_fruit
        for i in eachindex(geometry)
            @test @inferred(volume_source(fruits[i], geometry[i])) isa VolumeSource
            @test from_fruit[i].rates.energy ===
                  area_source(fruits[i].rates.energy, geometry[i].fruit_area_density)
            @test from_leaf[i] == leaf_air_sources(leaves[i], geometry[i]).air
            ledger = air[i] + (-from_leaf[i] + -from_fruit[i])
            @test abs(ledger.rates.energy) <= eps(FT) * abs(air[i].rates.energy)
        end
        tendency = layer_tendency.(boundaries[1:(end-1)], boundaries[2:end], air, geometry)
        @test eltype(eltype(tendency)) === FT
        # Column budget: interior faces cancel, both element classes enter once.
        integral = sum(dz .* getproperty.(getproperty.(tendency, :rates), :energy))
        expected = boundaries[1].rates.energy - boundaries[end].rates.energy +
                   sum(dz .* getproperty.(getproperty.(air, :rates), :energy))
        rtol = FT(config["validation"][FT === Float32 ? "float32_relative_tolerance" :
                                                         "float64_relative_tolerance"])
        @test isapprox(integral, expected; rtol)
        # Measure safety for the new class
        @test_throws MethodError leaves[2] + fruits[2]
        @test_throws MethodError volume_source(fruits[2], LayerGeometry(dz[2], one(FT)))
        @test_throws MethodError leaf_air_sources(fruits[2], geometry[2])
        # Empty fruit layer: nonfinite dummy fruit flux contributes exact zeros
        dummy = Exchange{FruitArea}(ExchangeRates(FT(NaN), FT(Inf), T(FT(NaN))))
        @test valid_exchange(volume_source(dummy, geometry[1]))
        steady = accumulate_steady_fluxes(boundaries[1], air, geometry)
        @test length(steady) == length(geometry) + 1
    end
end
