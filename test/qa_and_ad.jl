using Aqua
using CanopyExchange
using Documenter
using ExplicitImports
using ForwardDiff
using OffsetArrays
using SHA
using Test

@testset "Package QA (Aqua, ExplicitImports, Docs, Doctests)" begin
    Aqua.test_all(CanopyExchange)

    if VERSION >= v"1.11"
        @test check_no_implicit_imports(CanopyExchange) === nothing
        @test check_no_stale_explicit_imports(CanopyExchange) === nothing
        @test check_all_explicit_imports_via_owners(CanopyExchange) === nothing
        # The ForwardDiff extension imports `Dual`, `value` and `partials`, which
        # ForwardDiff documents but does not declare `public`.
        @test check_all_explicit_imports_are_public(
            CanopyExchange;
            ignore = (:Dual, :value, :partials),
        ) === nothing
        @test check_all_qualified_accesses_via_owners(CanopyExchange) === nothing
        @test check_all_qualified_accesses_are_public(CanopyExchange) === nothing
        @test check_no_self_qualified_accesses(CanopyExchange) === nothing
        @test isempty(Docs.undocumented_names(CanopyExchange))
    end

    DocMeta.setdocmeta!(
        CanopyExchange,
        :DocTestSetup,
        :(using CanopyExchange);
        recursive = true,
    )
    doctest(CanopyExchange; manual = false)
end

@testset "Exact IEEE-754 primal bit parity (SHA-256 19bcfb73...)" begin
    mktemp() do path, io
        close(io)
        script = joinpath(@__DIR__, "write_arithmetic_record.jl")
        proj = dirname(@__DIR__)
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$proj $script canonical $path`
        run(pipeline(cmd; stdout = devnull))
        digest = bytes2hex(sha256(read(path)))
        @test digest == "19bcfb73ba7adb1e88f285bbc2679182a861902729216e4948b9d3d350f0d08d"
    end
end

@testset "Floating-point scalar guards, measure separation, Base interfaces" begin
    # Reject discrete, irrational and discrete-dual scalars (v0.1.x: MethodError)
    dd_int = ForwardDiff.Dual{Nothing}(1, 1)
    @test dd_int isa ForwardDiff.Dual{Nothing,Int,1}
    @test_throws ArgumentError ExchangeRates(1, 2)
    @test_throws ArgumentError ExchangeRates(true, false)
    @test_throws ArgumentError ExchangeRates(1.0, 2.0, (; o2 = 1))
    @test_throws ArgumentError ExchangeRates(1.0, 2.0, (; o2 = 1 // 2))
    @test_throws ArgumentError ExchangeRates(1.0, 2.0, (; o2 = pi))
    @test_throws ArgumentError ExchangeRates(pi, 1.0)
    @test_throws ArgumentError ExchangeRates(dd_int, dd_int)
    @test_throws ArgumentError ExchangeRates(1.0, 2.0, (; o2 = dd_int))
    @test_throws ArgumentError ExchangeRates(1, 2, (; co2 = 3))
    @test_throws ArgumentError LayerGeometry(1, 0)
    @test_throws ArgumentError PlantLayerGeometry(1, 0, 0)

    # Reject mixed float precisions and mismatched tracer sets
    r64 = ExchangeRates(100.0, 2.0e-5, (; co2 = -1.0e-5))
    r32 = ExchangeRates(100.0f0, 2.0f-5, (; co2 = -1.0f-5))
    r64_o2 = ExchangeRates(100.0, 2.0e-5, (; co2 = -1.0e-5, o2 = 1.2e-5))
    @test_throws ArgumentError (r64 + r32)
    @test_throws ArgumentError (r64 + r64_o2)
    @test_throws ArgumentError convert(typeof(r32), r64)
    @test_throws ArgumentError convert(typeof(r64_o2), r64)

    # Strict measure separation across LeafFlux, StemFlux, FaceFlux, VolumeSource
    lf = LeafFlux(r64)
    sf = StemFlux(r64)
    ff = FaceFlux(r64)
    vs = VolumeSource(r64)
    pg = PlantLayerGeometry(2.0, 1.2, 0.3)
    @test_throws MethodError (lf + sf)
    @test_throws MethodError (lf + ff)
    @test_throws MethodError (sf + vs)
    @test_throws MethodError (ff - vs)
    @test_throws MethodError plant_air_sources(sf, lf, pg)
    @test_throws MethodError stem_area_source(10.0, LayerGeometry(2.0, 1.2))
    @test lf != sf
    @test !isequal(lf, sf)

    # Base.eltype, zero, keyword constructors, linear algebra
    for FT in (Float32, Float64)
        kw_no = ExchangeRates(; energy = FT(10), water = FT(1), tracers = (; co2 = FT(-2)))
        kw_o2 = ExchangeRates(;
            energy = FT(10),
            water = FT(1),
            tracers = (; co2 = FT(-2), o2 = FT(3)),
        )
        @test eltype(kw_no) === FT
        @test eltype(typeof(kw_o2)) === FT
        @test eltype(Exchange{LeafArea,typeof(kw_o2)}) === FT
        @test eltype(LayerGeometry{FT}) === FT
        @test eltype(PlantLayerGeometry{FT}) === FT
        @test eltype(MoistAirTendencies{FT,typeof(kw_o2.tracers)}) === FT
        @test zero(kw_no) == ExchangeRates(zero(FT), zero(FT), (; co2 = zero(FT)))
        @test zero(kw_o2) ==
              ExchangeRates(zero(FT), zero(FT), (; co2 = zero(FT), o2 = zero(FT)))
        @test zero(typeof(kw_o2)) === zero(kw_o2)
        for W in (LeafFlux, StemFlux, FaceFlux, VolumeSource)
            w = W(kw_o2)
            @test eltype(w) === FT
            @test zero(w) == W(zero(kw_o2))
            @test zero(typeof(w)) === zero(w)
            @test (w + (-w)) == zero(w)
            @test (FT(2) * w) / FT(2) == w
            @test (w * FT(4)) / FT(4) == w
            @test (w - w) == zero(w)
        end
        m = MoistAirTendencies(;
            rho_total = FT(1),
            rho_water = FT(1),
            rho_energy = FT(50),
            tracers = (; co2 = FT(-0.5), o2 = FT(0.6)),
        )
        @test eltype(m) === FT
        @test valid_exchange(m)
        @test (m + (-m)) == zero(m)
        @test zero(typeof(m)) === zero(m)
        @test (FT(2) * m) / FT(2) == m
        @test (m * FT(2) - m) == m
        @test_throws MethodError zero(LayerGeometry{FT})
        @test_throws MethodError zero(PlantLayerGeometry{FT})
    end

    # ==, isequal, and hash invariants across ±0.0, NaN, and ForwardDiff.Dual
    mkdual(v, p) = ForwardDiff.Dual{Nothing}(v, p)
    a_neg = LeafFlux(ExchangeRates(-0.0, 1.0, (; co2 = -2.0, o2 = -0.0)))
    a_pos = LeafFlux(ExchangeRates(0.0, 1.0, (; co2 = -2.0, o2 = 0.0)))
    @test a_neg == a_pos
    @test !isequal(a_neg, a_pos)
    @test hash(a_neg) == hash(a_pos)

    d_neg = LeafFlux(ExchangeRates(mkdual(-0.0, 1.0), mkdual(1.0, 0.0), (; co2 = mkdual(-2.0, 0.0))))
    d_pos = LeafFlux(ExchangeRates(mkdual(0.0, 1.0), mkdual(1.0, 0.0), (; co2 = mkdual(-2.0, 0.0))))
    @test d_neg == d_pos
    @test hash(d_neg) == hash(d_pos)

    # Record equality is componentwise scalar `==`. ForwardDiff 1.x compares partials
    # (so a nonzero partial differs from a plain zero); ForwardDiff 0.10 compares values.
    d_f64 = LeafFlux(ExchangeRates(0.0, 1.0, (; co2 = -2.0)))
    @test (d_pos != d_f64) == (mkdual(0.0, 1.0) != 0.0)
    @test d_pos != d_f64 || hash(d_pos) == hash(d_f64)
    d_zero_partial = LeafFlux(
        ExchangeRates(mkdual(0.0, 0.0), mkdual(1.0, 0.0), (; co2 = mkdual(-2.0, 0.0))),
    )
    @test d_zero_partial == d_f64
    @test hash(d_zero_partial) == hash(d_f64)

    nan1 = LeafFlux(ExchangeRates(NaN, 1.0, (; co2 = 0.0)))
    nan2 = LeafFlux(ExchangeRates(NaN, 1.0, (; co2 = 0.0)))
    @test nan1 != nan2
    @test isequal(nan1, nan2)
    @test hash(nan1) == hash(nan2)
    @test LeafFlux(ExchangeRates(1.0, 2.0, (; a = 0.0))) != LeafFlux(ExchangeRates(1.0, 2.0, (; b = 0.0)))

    # Base.promote_rule concrete AD container promotion vs incompatible non-promotion
    D64 = typeof(mkdual(1.0, 1.0))
    TD = NamedTuple{(:co2,),Tuple{D64}}
    r_dual = ExchangeRates(mkdual(1.0, 1.0), 2.0, (; co2 = -1.0))
    @test eltype([r64, r_dual]) === ExchangeRates{D64,TD}
    for W in (LeafFlux, StemFlux, FaceFlux, VolumeSource)
        v_prom = [W(r64), W(r_dual)]
        @test eltype(v_prom) === W{ExchangeRates{D64,TD}}
        @test isbitstype(eltype(v_prom))
        @test v_prom[1] == W(r64)
    end
    m64 = MoistAirTendencies(1.0e-5, 1.0e-5, 80.0, (; co2 = -1.0e-6))
    m_dual = MoistAirTendencies(mkdual(1.0e-5, 1.0), 1.0e-5, 80.0, (; co2 = -1.0e-6))
    @test eltype([m64, m_dual]) === MoistAirTendencies{D64,TD}
    @test eltype([LayerGeometry(1.5, 0.8), LayerGeometry(mkdual(1.5, 1.0), 0.8)]) ===
          LayerGeometry{D64}
    @test eltype([pg, PlantLayerGeometry(mkdual(2.0, 1.0), 1.2, 0.3)]) ===
          PlantLayerGeometry{D64}
    @test Base.promote_rule(typeof(r64), typeof(r32)) === Union{}
    @test Base.promote_rule(typeof(r64), typeof(r64_o2)) === Union{}
    @test Base.promote_rule(typeof(LeafFlux(r64)), typeof(LeafFlux(r32))) === Union{}
    @test Base.promote_rule(LayerGeometry{Float64}, LayerGeometry{Float32}) === Union{}

    # Mixed primal + Dual column accumulation via concrete container promotion
    mixed_faces = accumulate_steady_fluxes(
        FaceFlux(r64),
        [VolumeSource(r64), VolumeSource(ExchangeRates(mkdual(20.0, 1.0), 0.0, (; co2 = 0.0)))],
        [LayerGeometry(1.5, 0.8), LayerGeometry(2.0, mkdual(0.6, 0.5))],
    )
    @test eltype(mixed_faces) === FaceFlux{ExchangeRates{D64,TD}}
    @test mixed_faces[end].rates.energy == mkdual(100.0 + 1.5 * 100.0 + 2.0 * 20.0, 2.0)
end

@testset "ForwardDiff 1st/2nd-order AD, nonfinite partial guards, and BigFloat Duals" begin
    mkdual(v, p) = ForwardDiff.Dual{Nothing}(v, p)
    bad_dual = mkdual(1.0, NaN)
    inf_dual = mkdual(1.0, Inf)
    @test !CanopyExchange.all_finite(bad_dual)
    @test !CanopyExchange.all_finite(inf_dual)
    @test CanopyExchange.all_finite(mkdual(1.0, 2.0))
    @test CanopyExchange.exact_zero(mkdual(-0.0, 0.0))
    @test !CanopyExchange.exact_zero(mkdual(0.0, 1.0))
    @test CanopyExchange.primal(mkdual(mkdual(2.5, 1.0), mkdual(0.0, 1.0))) === 2.5
    @test CanopyExchange.primal_type(typeof(mkdual(mkdual(2.5f0, 1.0f0), mkdual(0.0f0, 1.0f0)))) ===
          Float32
    @test !valid_exchange(ExchangeRates(bad_dual, mkdual(0.0, 0.0), (; co2 = mkdual(0.0, 0.0))))
    @test !valid_exchange(
        ExchangeRates(mkdual(1.0, 0.0), mkdual(0.0, 0.0), (; co2 = mkdual(0.0, 0.0), o2 = bad_dual)),
    )
    @test !valid_exchange(
        MoistAirTendencies(bad_dual, mkdual(0.0, 0.0), mkdual(0.0, 0.0), (; co2 = mkdual(0.0, 0.0))),
    )
    @test_throws ArgumentError LayerGeometry(bad_dual, mkdual(1.0, 0.0))
    @test_throws ArgumentError PlantLayerGeometry(mkdual(1.0, 0.0), mkdual(1.0, 0.0), inf_dual)
    @test_throws ArgumentError moist_air_tendencies(
        VolumeSource(ExchangeRates(10.0, 1.0e-5, (; co2 = 0.0))),
        bad_dual,
    )

    # Empty-layer NaN guard vs right-tangent sensitivity at lad = 0
    @test area_source(mkdual(NaN, NaN), mkdual(0.0, 0.0)) === mkdual(0.0, 0.0)
    @test ForwardDiff.derivative(a -> area_source(25.0, a), 0.0) === 25.0
    @test area_source(mkdual(BigFloat(NaN), BigFloat(0)), mkdual(BigFloat(0), BigFloat(0))) ==
          mkdual(BigFloat(0), BigFloat(0))
    @test ForwardDiff.derivative(a -> area_source(BigFloat(25), a), BigFloat(0)) == BigFloat(25)

    # End-to-end 1st and 2nd derivatives of coupled plant + column + moist-air pipeline
    function pipeline_energy_tendency(x)
        P = CanopyExchange.primal_type(eltype(x))
        T2 = Tracers(:co2, :o2)
        sh_leaf, evap_leaf, hv, sh_stem, lad, sad, dz, eps_gas = x
        lf = sensible_water_flux(sh_leaf, evap_leaf, hv, T2(P(-1.2e-5), P(1.3e-5)))
        sf = sensible_stem_flux(sh_stem, P(0.2) * evap_leaf, hv, T2(P(1.0e-6), P(-1.0e-6)))
        geom = PlantLayerGeometry(dz, lad, sad)
        srcs = plant_air_sources(lf, sf, geom)
        bot = surface_static_energy_flux(P(30.0), P(70.0), P(2.8e-5), T2(P(2.0e-6), P(-2.0e-6)))
        bot_s = shift_water_energy_datum(bot, P(0.0))
        top = steady_top_flux(bot_s, P(0.5) * srcs.air, geom)
        tend = layer_tendency(bot_s, top, srcs.air, geom)
        moist = moist_air_tendencies(tend, eps_gas)
        return moist.rho_energy + P(1.0e6) * moist.rho_water + P(1.0e5) * moist.tracers.o2
    end

    x0 = [95.0, 4.2e-5, 2.54e6, 28.0, 1.4, 0.25, 1.8, 0.92]
    g_ad = ForwardDiff.gradient(pipeline_energy_tendency, x0)
    H_ad = ForwardDiff.hessian(pipeline_energy_tendency, x0)

    xb = BigFloat.(x0)
    g_ad_bf = ForwardDiff.gradient(pipeline_energy_tendency, xb)
    @test isapprox(g_ad, Float64.(g_ad_bf); rtol = 1e-12, atol = 1e-12)
    h_step = BigFloat(10)^(-12)
    for k in eachindex(x0)
        xp = copy(xb); xp[k] += h_step
        xm = copy(xb); xm[k] -= h_step
        df_bf = Float64((pipeline_energy_tendency(xp) - pipeline_energy_tendency(xm)) / (2 * h_step))
        @test isapprox(g_ad[k], df_bf; rtol = 1e-11, atol = 1e-12)
        d2f_bf = Float64(
            (
                pipeline_energy_tendency(xp) -
                2 * pipeline_energy_tendency(xb) +
                pipeline_energy_tendency(xm)
            ) / (h_step^2),
        )
        @test isapprox(H_ad[k, k], d2f_bf; rtol = 1e-9, atol = 1e-10)
    end
end

@testset "CliMA contracts, multi-layer geopotential telescoping, OffsetArrays, and zero allocs" begin
    for FT in (Float32, Float64)
        rtol = FT === Float32 ? 5.0f-5 : 1.0e-12
        geoms = [
            PlantLayerGeometry(FT(1.2), FT(0.0), FT(0.2)),
            PlantLayerGeometry(FT(1.8), FT(1.5), FT(0.3)),
            PlantLayerGeometry(FT(2.5), FT(1.1), FT(0.1)),
        ]
        z_faces = FT[0.0, 1.2, 3.0, 5.5]
        grav = FT(9.80665)
        hv0 = FT(2.501e6)

        leaf_fluxes = [
            sensible_water_flux(FT(NaN), FT(NaN), hv0, (; co2 = FT(NaN))),
            sensible_water_flux(FT(85), FT(4.0e-5), hv0, (; co2 = FT(-1.4e-5))),
            sensible_water_flux(FT(120), FT(6.0e-5), hv0, (; co2 = FT(-2.1e-5))),
        ]
        stem_fluxes = [
            sensible_stem_flux(FT(22), (; co2 = zero(FT))),
            sensible_stem_flux(FT(18), (; co2 = zero(FT))),
            sensible_stem_flux(FT(15), (; co2 = zero(FT))),
        ]
        paired = plant_air_sources.(leaf_fluxes, stem_fluxes, geoms)
        air_srcs = getproperty.(paired, :air)
        leaf_srcs = getproperty.(paired, :leaf)
        stem_srcs = getproperty.(paired, :stem)

        # Donor conservation air + (leaf + stem) == 0
        for k in eachindex(geoms)
            @test valid_exchange(air_srcs[k])
            @test air_srcs[k] + (leaf_srcs[k] + stem_srcs[k]) == zero(air_srcs[k])
            sum_left = (air_srcs[k] + leaf_srcs[k]) + stem_srcs[k]
            @test abs(sum_left.rates.energy) <= eps(FT) * FT(256)
            @test abs(sum_left.rates.water) <= eps(FT)
            @test abs(sum_left.rates.tracers.co2) <= eps(FT)
        end

        # Steady-air column accumulation & exact zero-tendency identity
        soil = surface_static_energy_flux(FT(45), FT(95), FT(3.8e-5), (; co2 = FT(3.0e-6)))
        faces = accumulate_steady_fluxes(soil, air_srcs, geoms)
        @test length(faces) == 4
        expected_top =
            soil + sum(layer_flux_increment(air_srcs[k], geoms[k]) for k in eachindex(geoms))
        @test isapprox(faces[end].rates.energy, expected_top.rates.energy; rtol)
        @test isapprox(faces[end].rates.water, expected_top.rates.water; rtol)
        @test isapprox(faces[end].rates.tracers.co2, expected_top.rates.tracers.co2; rtol)

        for k in eachindex(geoms)
            t_k = layer_tendency(faces[k], faces[k + 1], air_srcs[k], geoms[k])
            @test abs(t_k.rates.energy) <= rtol * abs(expected_top.rates.energy)
            @test abs(t_k.rates.water) <= rtol * abs(expected_top.rates.water)
        end

        # Spatially varying geopotential shift applied to SHARED streams BEFORE layer_tendency
        # must telescope interior faces identically across the column!
        z_mid = FT(0.5) .* (z_faces[1:(end - 1)] .+ z_faces[2:end])
        faces_geo = shift_water_energy_datum.(faces, grav .* z_faces)
        srcs_geo = shift_water_energy_datum.(air_srcs, grav .* z_mid)
        tends_geo = layer_tendency.(faces_geo[1:(end - 1)], faces_geo[2:end], srcs_geo, geoms)
        col_energy = sum(geoms[k].thickness * tends_geo[k].rates.energy for k in eachindex(geoms))
        expected_geo =
            faces_geo[1].rates.energy - faces_geo[end].rates.energy +
            sum(geoms[k].thickness * srcs_geo[k].rates.energy for k in eachindex(geoms))
        @test isapprox(col_energy, expected_geo; rtol, atol = eps(FT) * FT(1000))

        # OffsetArray, @view, and 0-layer column support
        off_srcs = OffsetArray(air_srcs, -1:1)
        off_geoms = OffsetArray(geoms, -1:1)
        off_faces = OffsetArray(similar(faces), 10:13)
        accumulate_steady_fluxes!(off_faces, soil, off_srcs, off_geoms)
        @test collect(off_faces) == faces

        view_faces = Vector{typeof(soil)}(undef, 3)
        accumulate_steady_fluxes!(view_faces, soil, @view(air_srcs[2:3]), @view(geoms[2:3]))
        @test view_faces[end] ==
              soil +
              layer_flux_increment(air_srcs[2], geoms[2]) +
              layer_flux_increment(air_srcs[3], geoms[3])

        empty_faces = accumulate_steady_fluxes(soil, eltype(air_srcs)[], eltype(geoms)[])
        @test empty_faces == [soil]

        @test_throws DimensionMismatch accumulate_steady_fluxes!(
            Vector{typeof(soil)}(undef, 2),
            soil,
            air_srcs,
            geoms,
        )
        @test_throws DimensionMismatch accumulate_steady_fluxes(
            soil,
            OffsetArray(air_srcs, 0:2),
            off_geoms,
        )
        other = FT === Float32 ? Float64 : Float32
        @test_throws ArgumentError accumulate_steady_fluxes(
            face_flux(one(other), zero(other), (; co2 = zero(other))),
            air_srcs,
            geoms,
        )

        # Gas-fraction moist-air tendency identity and validation
        m1 = moist_air_tendencies(air_srcs[2])
        @test m1.rho_energy === air_srcs[2].rates.energy
        @test m1.rho_total === air_srcs[2].rates.water
        @test m1.rho_water === air_srcs[2].rates.water
        @test m1.tracers.co2 === air_srcs[2].rates.tracers.co2
        for bad_eps in (zero(FT), -FT(0.1), FT(1.01), FT(NaN), FT(Inf))
            @test_throws ArgumentError moist_air_tendencies(air_srcs[2], bad_eps)
        end
    end

    # Zero-allocation check on scalar kernels
    function alloc_kernel(lf, sf, pg, bot)
        p = plant_air_sources(lf, sf, pg)
        top = steady_top_flux(bot, p.air, pg)
        t = layer_tendency(bot, top, p.air, pg)
        return moist_air_tendencies(shift_water_energy_datum(t, 100.0), 0.95)
    end
    measure_alloc_kernel(lf, sf, pg, bot) = @allocated alloc_kernel(lf, sf, pg, bot)
    T2 = Tracers(:co2, :o2)
    lf0 = sensible_water_flux(90.0, 3.0e-5, 2.5e6, T2(-1.0e-5, 1.1e-5))
    sf0 = sensible_stem_flux(25.0, T2(0.0, 0.0))
    pg0 = PlantLayerGeometry(2.0, 1.2, 0.3)
    bot0 = surface_static_energy_flux(40.0, 80.0, 3.2e-5, T2(1.0e-6, -1.0e-6))
    measure_alloc_kernel(lf0, sf0, pg0, bot0)
    @test measure_alloc_kernel(lf0, sf0, pg0, bot0) == 0

    measure_accum_alloc(dest, s, srcs, geoms) =
        @allocated accumulate_steady_fluxes!(dest, s, srcs, geoms)
    srcs_alloc = [plant_air_sources(lf0, sf0, pg0).air, plant_air_sources(lf0, sf0, pg0).air]
    geoms_alloc = [pg0, pg0]
    faces_alloc = Vector{typeof(bot0)}(undef, 3)
    measure_accum_alloc(faces_alloc, bot0, srcs_alloc, geoms_alloc)
    @test measure_accum_alloc(faces_alloc, bot0, srcs_alloc, geoms_alloc) == 0
end
