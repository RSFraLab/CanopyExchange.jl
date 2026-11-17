# Extension without touching the package: a second tracer set through the tagged
# records, and a third element class with its own geometry defined only here.
using Test, CanopyExchange
@isdefined(exchange_fixture) || include("exchange_fixture.jl")

struct FruitArea <: CanopyExchange.Measure end

# The geometry interface: `thickness` (here a field), `area_density` per element class,
# and `broadcastable` to act as a scalar in broadcasts. No `Base.eltype` is needed.
struct OrchardGeometry{FT}
    thickness::FT
    leaf_area_density::FT
    fruit_area_density::FT
end
Base.Broadcast.broadcastable(g::OrchardGeometry) = Ref(g)
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
        precision = FT === Float32 ? "float32" : "float64"
        rtol = FT(config["validation"][precision*"_relative_tolerance"])
        atol = FT(config["validation"][precision*"_scaled_absolute_tolerance"])
        dz = FT.(layers["thickness_m"])
        lad, fad = FT.(layers["leaf_area_density_m2_m3"]), FT.(fruit["fruit_area_density_m2_m3"])
        hv = FT.(layers["vapor_enthalpy_J_kg"])
        leaf_heat, leaf_water, leaf_co2 = (FT.(layers[k]) for k in
            ("sensible_heat_W_m2_leaf", "water_kg_m2_leaf_s", "carbon_mol_m2_leaf_s"))
        fruit_heat, fruit_water, fruit_co2 = (FT.(fruit[k]) for k in
            ("sensible_heat_W_m2_fruit", "water_kg_m2_fruit_s", "carbon_mol_m2_fruit_s"))
        face_energy, face_water, face_co2 = (FT.(faces[k]) for k in
            ("energy_W_m2_ground", "water_kg_m2_ground_s", "carbon_mol_m2_ground_s"))
        geometry = OrchardGeometry.(dz, lad, fad)
        T = Tracers(:co2)
        leaves = sensible_water_flux.(leaf_heat, leaf_water, hv, T.(leaf_co2))
        # The enthalpy-enforcing constructor takes a measure instance (broadcast-safe).
        fruits = sensible_water_flux.(FruitArea(), fruit_heat, fruit_water, hv, T.(fruit_co2))
        boundaries = face_flux.(face_energy, face_water, T.(face_co2))
        @test eltype(fruits) <: Exchange{FruitArea}
        from_leaf = volume_source.(leaves, geometry)
        from_fruit = volume_source.(fruits, geometry)
        air = from_leaf .+ from_fruit
        # Independent per-layer sources from the TOML inputs: density x (H + E h_v), E, CO2
        components(x) = (x.rates.energy, x.rates.water, x.rates.tracers.co2)
        for i in eachindex(geometry)
            fruit_expected = (
                fad[i] * (fruit_heat[i] + fruit_water[i] * hv[i]),
                fad[i] * fruit_water[i],
                fad[i] * fruit_co2[i],
            )
            leaf_expected = (
                lad[i] * (leaf_heat[i] + leaf_water[i] * hv[i]),
                lad[i] * leaf_water[i],
                lad[i] * leaf_co2[i],
            )
            @test components(from_fruit[i]) == fruit_expected
            @test components(from_leaf[i]) == leaf_expected
            @test @inferred(volume_source(fruits[i], geometry[i])) isa VolumeSource
        end
        # Column budgets of energy, water and CO2: interior faces cancel, both element
        # classes enter once.
        tendency = layer_tendency.(boundaries[1:(end-1)], boundaries[2:end], air, geometry)
        @test eltype(eltype(tendency)) === FT
        steady = accumulate_steady_fluxes(boundaries[1], air, geometry)
        for (j, face_values, leaf_rate, fruit_rate) in (
            (1, face_energy, leaf_heat .+ leaf_water .* hv, fruit_heat .+ fruit_water .* hv),
            (2, face_water, leaf_water, fruit_water),
            (3, face_co2, leaf_co2, fruit_co2),
        )
            plant = sum(dz .* (lad .* leaf_rate .+ fad .* fruit_rate))
            scale = maximum(abs, face_values) + sum(abs.(dz .* (lad .* leaf_rate .+ fad .* fruit_rate)))
            integral = sum(dz[i] * components(tendency[i])[j] for i in eachindex(dz))
            @test isapprox(integral, face_values[1] - face_values[end] + plant; rtol, atol = atol * scale)
            # Steady column: canopy-top flux = bottom flux + all plant sources
            @test isapprox(components(steady[end])[j], face_values[1] + plant; rtol, atol = atol * scale)
        end
        # Measure safety for the new class
        @test_throws MethodError leaves[2] + fruits[2]
        @test_throws MethodError volume_source(fruits[2], LayerGeometry(dz[2], one(FT)))
        @test_throws MethodError leaf_air_sources(fruits[2], geometry[2])
        # Empty fruit layer: nonfinite dummy fruit flux contributes exact zeros
        dummy = Exchange{FruitArea}(ExchangeRates(FT(NaN), FT(Inf), T(FT(NaN))))
        @test components(volume_source(dummy, geometry[1])) == (zero(FT), zero(FT), zero(FT))
        # A geometry with `broadcastable` is a scalar in broadcasts
        @test volume_source.(fruits, geometry[2]) == [volume_source(f, geometry[2]) for f in fruits]
    end
end
