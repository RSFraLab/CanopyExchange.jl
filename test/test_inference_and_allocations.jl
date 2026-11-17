# `@inferred`, result element type, isbits results and zero allocations for every
# exported kernel and constructor, in Float32 and Float64.
using Test, CanopyExchange

# Helper functions keep @allocated out of @testset scope (Julia 1.10).
allocations(f::F, args::Vararg{Any,N}) where {F,N} = (f(args...); @allocated f(args...))
# A type passed through `Vararg{Any}` is not specialized on Julia 1.10; wrap it.
stem_area_density(g) = CanopyExchange.area_density(g, StemArea)

function kernel_cases(::Type{FT}, ntracers::Int) where {FT}
    T2 = Tracers(ntuple(i -> Symbol(:t, i), ntracers)...)
    tr = T2(ntuple(i -> FT(i) * FT(1e-6), ntracers)...)
    rates = ExchangeRates(FT(100), FT(2e-5), tr)
    leaf = sensible_water_flux(FT(80), FT(4e-5), FT(2.5e6), tr)
    stem = sensible_stem_flux(FT(20), FT(1e-6), FT(2.5e6), tr)
    face = face_flux(FT(40), FT(3e-5), tr)
    g = LayerGeometry(FT(2), FT(1.2))
    pg = PlantLayerGeometry(FT(2), FT(1.2), FT(0.3))
    source = leaf_air_sources(leaf, g).air
    moist = moist_air_tendencies(source, FT(0.96))
    return Any[
        (ExchangeRates, (FT(1), FT(2), tr)),
        (ExchangeRates, (FT(1), FT(2))),
        (MoistAirTendencies, (FT(1), FT(1), FT(2), tr)),
        (T2, values(tr)),
        (Exchange{LeafArea}, (rates,)),
        (LeafFlux, (rates,)),
        (leaf_flux, (FT(1), FT(2), tr)),
        (stem_flux, (FT(1), FT(2), tr)),
        (face_flux, (FT(1), FT(2), tr)),
        (LayerGeometry, (FT(2), FT(1))),
        (PlantLayerGeometry, (FT(2), FT(1), FT(0.5))),
        (PlantLayerGeometry, (FT(2), FT(1))),
        (layer_geometry, (FT(2), FT(1))),
        (layer_geometry, (FT(2), FT(1), FT(0.5))),
        (CanopyExchange.thickness, (pg,)),
        (stem_area_density, (pg,)),
        (plant_area_density, (pg,)),
        (area_source, (FT(150), FT(0.8))),
        (area_source, (FT(NaN), zero(FT))),
        (volume_source, (stem, pg)),
        (leaf_area_source, (FT(2e-5), g)),
        (stem_area_source, (FT(40), pg)),
        (leaf_air_sources, (leaf, g)),
        (plant_air_sources, (leaf, stem, pg)),
        (conserved_tendency, (FT(15), FT(12), FT(4), g)),
        (layer_tendency, (face, face, source, g)),
        (layer_flux_increment, (source, g)),
        (steady_top_flux, (face, source, g)),
        (sensible_water_flux, (FT(80), FT(4e-5), FT(2.5e6), tr)),
        (sensible_water_flux, (FT(80), FT(4e-5), FT(2.5e6))),
        (sensible_stem_flux, (FT(20), tr)),
        (sensible_stem_flux, (FT(20),)),
        (sensible_stem_flux, (FT(20), FT(1e-6), FT(2.5e6), tr)),
        (surface_static_energy_flux, (FT(45), FT(95), FT(3.8e-5), tr)),
        (surface_static_energy_flux, (FT(45), FT(95), FT(3.8e-5))),
        (shift_water_energy_datum, (rates, FT(170000))),
        (shift_water_energy_datum, (face, FT(170000))),
        (moist_air_tendencies, (source,)),
        (moist_air_tendencies, (source, FT(0.96))),
        (valid_exchange, (leaf,)),
        (valid_exchange, (moist,)),
        (+, (leaf, leaf)),
        (-, (rates, rates)),
        (-, (source,)),
        (*, (FT(2), leaf)),
        (*, (moist, FT(2))),
        (/, (face, FT(2))),
        (zero, (moist,)),
        (zero, (typeof(leaf),)),
        (==, (leaf, leaf)),
        (==, (g, g)),
        (isequal, (leaf, leaf)),
        (isequal, (moist, moist)),
        (hash, (leaf,)),
        (hash, (pg,)),
        (hash, (moist,)),
        (sensible_water_flux, (StemArea(), FT(20), FT(1e-6), FT(2.5e6), tr)),
        (CanopyExchange.all_finite, (FT(1),)),
        (CanopyExchange.exact_zero, (FT(0),)),
        (CanopyExchange.primal, (FT(1),)),
    ]
end

element_type(x::Bool) = Bool
element_type(x::UInt) = UInt
element_type(x::Type) = Type
element_type(x::NamedTuple) = eltype(first(values(x)))
element_type(x) = eltype(x)

# Two tracers (co2, o2) and twelve tracers; Julia's tuple limit is 31 (see ExchangeRates).
@testset "Inference, result type and allocations: $FT, $n tracers" for FT in (Float32, Float64),
    n in (2, 12)

    for (f, args) in kernel_cases(FT, n)
        result = @inferred f(args...)
        @test element_type(result) in (FT, Bool, UInt)
        @test isbits(result)
        @test allocations(f, args...) == 0
    end
    # Steady accumulation over vectors: inferred result, no allocation in place
    face = face_flux(FT(40), FT(3e-5), (; co2 = FT(1e-6)))
    sources = fill(VolumeSource(ExchangeRates(FT(5), FT(1e-5), (; co2 = FT(-1e-6)))), 3)
    geometries = fill(PlantLayerGeometry(FT(2), FT(1.2), FT(0.3)), 3)
    faces = @inferred accumulate_steady_fluxes(face, sources, geometries)
    @test eltype(eltype(faces)) === FT
    @test allocations(accumulate_steady_fluxes!, faces, face, sources, geometries) == 0
    @test CanopyExchange.primal_type(FT) === FT
end
