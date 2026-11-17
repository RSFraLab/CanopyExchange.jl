# Precision policy, measure safety, tracer rules and the counterexamples with which
# the first review round broke the v1 prototype (canopy-refactor/prototype/check.jl).
using Test, ForwardDiff, CanopyExchange

# An external layer geometry with a thickness only and no Base.eltype.
struct BareGeometry
    thickness::Float64
end

@testset "Precision policy in constructors and every kernel" begin
    T = Tracers(:co2)
    r64 = ExchangeRates(100.0, 2e-5, T(-1e-5))
    r32 = ExchangeRates(100.0f0, 2.0f-5, T(-1.0f-5))
    g64, g32 = LayerGeometry(2.0, 1.0), LayerGeometry(2.0f0, 1.0f0)
    pg32 = PlantLayerGeometry(2.0f0, 1.0f0, 0.25f0)
    @test_throws ArgumentError ExchangeRates(1.0f0, 2.0)
    @test_throws ArgumentError ExchangeRates(1.0f0, 2.0f0, T(3.0))
    @test_throws ArgumentError MoistAirTendencies(1.0f0, 1.0f0, 2.0, T(3.0f0))
    @test_throws ArgumentError LayerGeometry(2.0f0, 1.0)
    @test_throws ArgumentError PlantLayerGeometry(2.0f0, 1.0f0, 0.5)
    @test_throws ArgumentError r64 + r32
    @test_throws ArgumentError r32 * 2.0
    @test_throws ArgumentError 2.0 * r32
    @test_throws ArgumentError r32 / 2.0
    @test_throws ArgumentError LeafFlux(r32) * 2.0
    @test_throws ArgumentError leaf_air_sources(LeafFlux(r32), g64)
    @test_throws ArgumentError leaf_air_sources(LeafFlux(r64), g32)
    @test_throws ArgumentError plant_air_sources(LeafFlux(r32), StemFlux(r32), PlantLayerGeometry(2.0, 1.0, 0.25))
    @test_throws ArgumentError plant_air_sources(LeafFlux(r32), StemFlux(r64), pg32)
    @test_throws ArgumentError area_source(1.0f0, 2.0)
    @test_throws ArgumentError leaf_area_source(1.0f0, g64)
    @test_throws ArgumentError conserved_tendency(1.0f0, 2.0, 3.0f0, g32)
    @test_throws ArgumentError conserved_tendency(1.0f0, 2.0f0, 3.0f0, g64)
    @test_throws ArgumentError shift_water_energy_datum(LeafFlux(r32), 500.0)
    @test_throws ArgumentError moist_air_tendencies(VolumeSource(r32), 0.96)
    @test_throws ArgumentError layer_tendency(FaceFlux(r32), FaceFlux(r32), VolumeSource(r64), g64)
    @test_throws ArgumentError layer_tendency(FaceFlux(r32), FaceFlux(r32), VolumeSource(r32), g64)
    @test_throws ArgumentError steady_top_flux(FaceFlux(r32), VolumeSource(r32), g64)
    @test_throws ArgumentError layer_flux_increment(VolumeSource(r32), g64)
    # constructors that combine components before building the record
    @test_throws ArgumentError sensible_water_flux(80.0f0, 4e-5, 2.5e6)
    @test_throws ArgumentError sensible_water_flux(80.0f0, 4.0f-5, 2.5e6, T(1.0f0))
    @test_throws ArgumentError sensible_water_flux(80, 4e-5, 2.5e6)
    @test_throws ArgumentError sensible_stem_flux(20.0f0, 1e-6, 2.5e6)
    @test_throws ArgumentError sensible_stem_flux(20, (; co2 = 0.0))
    @test_throws ArgumentError surface_static_energy_flux(45.0f0, 95.0, 3.8e-5)
    @test_throws ArgumentError surface_static_energy_flux(45, 95, 0.0)
    @test_throws ArgumentError leaf_flux(1.0f0, 2.0)
    @test_throws ArgumentError layer_geometry(1.0f0, 0.5)
    # integers and booleans are not exchange scalars, in any argument position
    @test_throws ArgumentError ExchangeRates(1, 2)
    @test_throws ArgumentError("canopy exchange scalars must be floating point") ExchangeRates(1.0, 2)
    @test_throws ArgumentError("canopy exchange scalars must be floating point") ExchangeRates(1.0, 2.0, T(1))
    @test_throws ArgumentError("canopy exchange scalars must be floating point") conserved_tendency(1.0, 2, 3.0, g64)
    @test_throws ArgumentError("mixed floating-point precision in canopy exchange") ExchangeRates(1.0, 2.0f0)
    # tracer values must be real scalars
    for bad in ([1.0, 2.0], (1.0, 2.0), nothing, "x", 1.0 + 0im, missing, r64)
        @test_throws ArgumentError("tracer values must be real scalars") ExchangeRates(1.0, 2.0, (; a = bad))
    end
    @test_throws ArgumentError Tracers(:co2, :co2)
    @test_throws ArgumentError ExchangeRates(true, false)
    @test_throws ArgumentError LayerGeometry(1, 0)
    @test_throws ArgumentError moist_air_tendencies(VolumeSource(r64), 1)
end

@testset "Measure safety, tracer rules, empty layer, equality, algebra" begin
    T = Tracers(:co2)
    r64 = ExchangeRates(100.0, 2e-5, T(-1e-5))
    g64 = LayerGeometry(2.0, 1.0)
    pg = PlantLayerGeometry(2.0, 1.0, 0.25)
    leaf, stem = LeafFlux(r64), StemFlux(ExchangeRates(40.0, 0.0, T(0.0)))
    @test_throws MethodError leaf + stem
    @test_throws MethodError leaf + FaceFlux(r64)
    @test_throws MethodError leaf_air_sources(stem, g64)
    @test_throws MethodError leaf_air_sources(r64, g64)
    # a stem flux needs a geometry that has stems
    @test_throws MethodError plant_air_sources(leaf, stem, g64)
    @test_throws MethodError volume_source(stem, g64)
    @test_throws MethodError volume_source(FaceFlux(r64), pg)
    @test_throws MethodError volume_source(VolumeSource(r64), pg)
    # scalar kernels take scalars only, so a record cannot lose its measure tag
    @test_throws MethodError area_source(leaf, 0.5)
    @test_throws MethodError leaf_area_source(leaf, g64)
    @test_throws MethodError stem_area_source(stem, pg)
    @test_throws MethodError conserved_tendency(FaceFlux(r64), FaceFlux(r64), FaceFlux(r64), g64)
    @test_throws MethodError conserved_tendency(VolumeSource(r64), VolumeSource(r64), VolumeSource(r64), g64)
    # a measure-tagged record needs a concrete measure
    @test_throws ArgumentError Exchange{CanopyExchange.Measure}(r64)
    @test_throws ArgumentError Exchange{Union{}}(r64)
    # conversion returns exactly the requested type or throws
    @test convert(typeof(r64), r64) === r64
    @test_throws ArgumentError convert(ExchangeRates{Float64,NamedTuple{(:co2,),Tuple{Float32}}}, r64)
    s = plant_air_sources(leaf, stem, pg)
    @test iszero(s.air.rates.energy + s.leaf.rates.energy + s.stem.rates.energy)
    @test volume_source(leaf, pg) == -s.leaf

    # absent tracer is not a zero tracer; tracer order is part of the type
    with_o2 = LeafFlux(ExchangeRates(100.0, 2e-5, Tracers(:co2, :o2)(-1e-5, 0.0)))
    o2_first = LeafFlux(ExchangeRates(100.0, 2e-5, Tracers(:o2, :co2)(0.0, -1e-5)))
    @test typeof(leaf) !== typeof(with_o2)
    @test typeof(with_o2) !== typeof(o2_first)
    @test_throws ArgumentError("tracer sets of canopy exchange records differ") leaf + with_o2
    @test_throws ArgumentError("tracer sets of canopy exchange records differ") with_o2 + o2_first
    @test with_o2 != o2_first
    @test_throws ArgumentError plant_air_sources(with_o2, stem, pg)
    # tracer names cannot overwrite host fields
    odd = VolumeSource(ExchangeRates(10.0, 1e-5, (; rho_total = 7.0)))
    @test moist_air_tendencies(odd).rho_total === 1e-5
    @test moist_air_tendencies(odd).tracers.rho_total === 7.0
    # no tracers at all is a valid tracer set
    bare = moist_air_tendencies(VolumeSource(ExchangeRates(10.0, 1e-5)), 0.5)
    @test bare.tracers === (;) && bare.rho_energy === 20.0

    # empty layer
    empty = LayerGeometry(1.0, 0.0)
    bad = LeafFlux(ExchangeRates(NaN, Inf, T(NaN)))
    @test valid_exchange(leaf_air_sources(bad, empty).air)
    @test leaf_air_sources(bad, empty).air == zero(leaf_air_sources(bad, empty).air)
    @test !valid_exchange(leaf_air_sources(bad, g64).air)
    @test area_source(NaN, -0.0) === 0.0

    # value equality and hashing
    z, mz = ExchangeRates(0.0, 0.0), ExchangeRates(-0.0, 0.0)
    @test z == mz && !isequal(z, mz) && hash(z) == hash(mz)
    @test ExchangeRates(NaN, 0.0) != ExchangeRates(NaN, 0.0)
    @test isequal(ExchangeRates(NaN, 0.0), ExchangeRates(NaN, 0.0))
    @test hash(r64) == hash(ExchangeRates(100.0, 2e-5, T(-1e-5)))
    @test LeafFlux(r64) != StemFlux(r64) && LeafFlux(r64) == leaf
    @test hash(LeafFlux(r64)) != hash(StemFlux(r64))
    @test LayerGeometry(2.0, 1.0) == g64 && LayerGeometry(2.0, 1.0) != PlantLayerGeometry(2.0, 1.0, 0.0)
    @test length(Set([leaf, LeafFlux(r64), StemFlux(r64)])) == 2

    # algebra on tagged records and reductions
    @test (2.0 * leaf).rates.energy === 200.0 && (leaf / 2.0).rates.energy === 50.0
    @test sum([VolumeSource(r64), VolumeSource(r64)]).rates.energy === 200.0
    @test zero(typeof(leaf)) === zero(leaf)
    @test Exchange{LeafArea,typeof(r64)}(r64) === leaf
    # Records are not collections: no map, iteration or length
    @test_throws MethodError map(x -> 2x, r64)
    @test_throws MethodError collect(r64)
    m = moist_air_tendencies(VolumeSource(r64), 0.5)
    @test (m + m).rho_energy === 400.0 && (2.0 * m).rho_total === 2 * m.rho_total
    @test moist_air_tendencies(VolumeSource(r64)) === moist_air_tendencies(VolumeSource(r64), 1.0)

    # records are scalars in broadcasts
    @test (r -> r.energy).(LeafFlux(r64).rates) === 100.0
    @test getproperty.(leaf_air_sources.([leaf, leaf], g64), :air) == fill(volume_source(leaf, g64), 2)
end

@testset "Steady accumulation accepts what the in-place version accepts" begin
    T = Tracers(:co2)
    face = face_flux(10.0, 1e-5, T(1e-6))
    source = VolumeSource(ExchangeRates(5.0, 2e-5, T(-1e-6)))
    sources = [source, source, source]
    leaf_only, plant = LayerGeometry(2.0, 0.5), PlantLayerGeometry(2.0, 0.5, 0.25)
    expected = [10.0, 20.0, 30.0, 40.0]
    energies(faces) = [f.rates.energy for f in faces]
    for geometries in (
        [leaf_only, plant, leaf_only],
        LayerGeometry[leaf_only, leaf_only, leaf_only],
        Any[leaf_only, plant, leaf_only],
        [(; thickness = 2.0, leaf_area_density = 0.5) for _ in 1:3],
    )
        faces = accumulate_steady_fluxes(face, sources, geometries)
        @test energies(faces) == expected
        @test eltype(faces) === typeof(face)
        inplace = accumulate_steady_fluxes!(similar(faces), face, sources, geometries)
        @test inplace == faces
    end
    @test energies(accumulate_steady_fluxes(face, VolumeSource[source, source, source], [plant, plant, plant])) == expected
    @test accumulate_steady_fluxes(face, VolumeSource[], LayerGeometry[]) == [face]
    @test accumulate_steady_fluxes(face, typeof(source)[], Any[]) == [face]
    # external geometry without Base.eltype
    @test energies(accumulate_steady_fluxes(face, sources, fill(BareGeometry(2.0), 3))) == expected
    # dual thickness promotes the faces
    d(v, p) = ForwardDiff.Dual{Nothing}(v, p)
    dual_faces = accumulate_steady_fluxes(face, sources, [LayerGeometry(d(2.0, 1.0), 0.5) for _ in 1:3])
    @test dual_faces[end].rates.energy == d(40.0, 15.0)
    # mixed precision is still rejected
    @test_throws ArgumentError accumulate_steady_fluxes(face, sources, fill(LayerGeometry(2.0f0, 0.5f0), 3))
    @test_throws ArgumentError accumulate_steady_fluxes(face, sources, [(; thickness = 2) for _ in 1:3])
end

@testset "Automatic differentiation contracts" begin
    @test Base.get_extension(CanopyExchange, :CanopyExchangeForwardDiffExt) isa Module
    D = ForwardDiff.Dual
    d(v, p) = D{Nothing}(v, p)
    # mixed plain and dual inputs promote; Float32 with a Float64 dual does not
    lf = sensible_water_flux(d(80.0, 1.0), 4e-5, 2.5e6, (; co2 = -1.2e-5))
    @test eltype(lf) === typeof(d(1.0, 1.0)) && isbitstype(typeof(lf))
    @test LayerGeometry(2.0, d(0.6, 0.5)) isa LayerGeometry{typeof(d(1.0, 1.0))}
    @test_throws ArgumentError ExchangeRates(d(1.0, 1.0), 2.0f0)
    @test_throws ArgumentError ExchangeRates(D{Nothing}(1, 1), D{Nothing}(1, 1))
    @test_throws ArgumentError LeafFlux(ExchangeRates(1.0f0, 2.0f0)) * d(2.0, 1.0)
    # dual flux on a plain geometry, plain flux on a dual geometry
    f(h) = leaf_air_sources(sensible_water_flux(h, 4e-5, 2.5e6, (; co2 = -1.2e-5)),
        LayerGeometry(2.0, 1.0)).air.rates.energy
    @test ForwardDiff.derivative(f, 80.0) === 1.0
    g(lad) = leaf_air_sources(leaf_flux(180.0, 4e-5, (; co2 = -1.2e-5)),
        LayerGeometry(2.0, lad)).air.rates.energy
    @test ForwardDiff.derivative(g, 0.5) === 180.0
    # empty-layer guard and its right derivative
    @test area_source(d(NaN, NaN), d(0.0, 0.0)) === d(0.0, 0.0)
    @test ForwardDiff.derivative(a -> area_source(25.0, a), 0.0) === 25.0
    @test ForwardDiff.derivative(a -> area_source(BigFloat(25), a), BigFloat(0)) == 25
    # nonfinite derivative components are rejected
    @test !valid_exchange(ExchangeRates(d(1.0, NaN), d(1.0, 0.0)))
    @test !valid_exchange(ExchangeRates(1.0, 0.0, (; co2 = d(1.0, Inf))))
    @test_throws ArgumentError LayerGeometry(d(1.0, Inf), d(0.5, 0.0))
    # checks compare values, not derivatives (identical on ForwardDiff 0.10 and 1.x)
    v = VolumeSource(ExchangeRates(10.0, 1e-5))
    @test moist_air_tendencies(v, d(1.0, 1.0)).rho_energy == d(10.0, -10.0)
    @test PlantLayerGeometry(d(1.0, 0.0), d(0.0, -1.0), d(0.0, 0.0)) isa PlantLayerGeometry
    @test_throws ArgumentError moist_air_tendencies(v, d(1.5, 0.0))
    @test_throws ArgumentError LayerGeometry(d(0.0, 1.0), 1.0)
    # record equality compares values and derivative components on every ForwardDiff
    # version; a plain number has zero derivative components
    plain = leaf_flux(1.0, 2.0, (; co2 = 3.0))
    @test leaf_flux(d(1.0, 5.0), 2.0, (; co2 = 3.0)) != plain
    @test plain != leaf_flux(d(1.0, 5.0), 2.0, (; co2 = 3.0))
    @test leaf_flux(d(1.0, 0.0), 2.0, (; co2 = 3.0)) == plain
    @test plain == leaf_flux(d(1.0, -0.0), 2.0, (; co2 = 3.0))
    @test !isequal(plain, leaf_flux(d(1.0, -0.0), 2.0, (; co2 = 3.0)))
    @test leaf_flux(d(1.0, 5.0), 2.0) != leaf_flux(d(1.0, 4.0), 2.0)
    @test CanopyExchange.scalar_equal(==, d(1.0, 0.0), 1.0)
    @test !CanopyExchange.scalar_equal(==, d(1.0, 1.0), 1.0)
    @test length(Set([plain, leaf_flux(d(1.0, 5.0), 2.0, (; co2 = 3.0))])) == 2
    # gradient against central differences and a finite Hessian through the pipeline
    function pipeline(x)
        sh_leaf, evap, hv, sh_stem, lad, sad, dz, ε = x
        T = Tracers(:co2, :o2)
        leaf = sensible_water_flux(sh_leaf, evap, hv, T(-1.2e-5, 1.3e-5))
        stem = sensible_stem_flux(sh_stem, 0.2 * evap, hv, T(1e-6, -1e-6))
        geometry = PlantLayerGeometry(dz, lad, sad)
        src = plant_air_sources(leaf, stem, geometry)
        bottom = shift_water_energy_datum(face_flux(100.0, 2.8e-5, T(2e-6, -2e-6)), 0.0)
        top = steady_top_flux(bottom, 0.5 * src.air, geometry)
        m = moist_air_tendencies(layer_tendency(bottom, top, src.air, geometry), ε)
        return m.rho_energy + 1e6 * m.rho_water + 1e5 * m.tracers.o2
    end
    x0 = [80.0, 4e-5, 2.5e6, 20.0, 1.2, 0.3, 2.0, 0.96]
    grad = ForwardDiff.gradient(pipeline, x0)
    fd = map(eachindex(x0)) do i
        h = 1e-6 * max(abs(x0[i]), 1e-6)
        e = zeros(8)
        e[i] = h
        (pipeline(x0 .+ e) - pipeline(x0 .- e)) / 2h
    end
    @test isapprox(grad, fd; rtol = 1e-5)
    @test all(isfinite, ForwardDiff.hessian(pipeline, x0))
end
