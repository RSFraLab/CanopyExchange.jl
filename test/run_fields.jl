# Optional check of broadcast construction and kernels on CuArray and on ClimaCore
# fields (CPU and CUDA), Float32 and Float64; not part of Pkg.test(). It builds its
# own temporary environment with CUDA, ClimaCore and this checkout of CanopyExchange:
#     julia --startup-file=no test/run_fields.jl
using Pkg
Pkg.activate(; temp = true)
Pkg.develop(PackageSpec(path = dirname(@__DIR__)); io = devnull)
Pkg.add(["CUDA", "ClimaCore", "ClimaComms"]; io = devnull)

using Test, CUDA
import ClimaCore as CC
import ClimaComms
using CanopyExchange

CUDA.functional() || error("CUDA is required for this verification")
CUDA.allowscalar(false)
println("CUDA ", pkgversion(CUDA), ", ClimaCore ", pkgversion(CC), ", Julia ", VERSION)
const T2 = Tracers(:co2, :o2)
air_of(s) = s.air

# Construction by broadcast goes through functions and the Tracers callable; GPU
# field broadcasts reject a type or a keyword argument in the function position.
function pipeline(h, w, hv, c, o, fb, ft, dz, lad, sad)
    leaf = sensible_water_flux.(h, w, hv, T2.(c, o))
    stem = sensible_stem_flux.(h, T2.(o, c)) .+ sensible_water_flux.(StemArea(), h, w, hv, T2.(c, o))
    g = layer_geometry.(dz, lad)
    pg = layer_geometry.(dz, lad, sad)
    src = leaf_air_sources.(leaf, g)
    air = air_of.(src)
    plant_air = air_of.(plant_air_sources.(leaf, stem, pg))
    bottom = face_flux.(fb, w, T2.(c, o))
    top = surface_static_energy_flux.(ft, ft, w, T2.(c, o))
    tend = layer_tendency.(bottom, top, air, g)
    steady = steady_top_flux.(shift_water_energy_datum.(bottom, hv), plant_air, pg)
    moist = moist_air_tendencies.(tend, eltype(dz)(0.96))
    return (; leaf, tend, steady, moist, ok = valid_exchange.(tend))
end

@testset "CuArray $FT" for FT in (Float32, Float64)
    n = 64
    cpu = (
        FT.(range(10, 90; length = n)), fill(FT(4e-5), n), fill(FT(2.5e6), n),
        fill(FT(-1e-5), n), fill(FT(1e-5), n), fill(FT(20), n), fill(FT(40), n),
        fill(FT(1.5), n), FT.(range(0, 1; length = n)), fill(FT(0.2), n),
    )
    ref = pipeline(cpu...)
    gpu = pipeline(map(CuArray, cpu)...)
    @test gpu.tend isa CuArray && isbitstype(eltype(gpu.tend)) && eltype(eltype(gpu.tend)) === FT
    # GPU arithmetic may contract a*b + c, so sums are compared to rounding tolerance
    # and pure products exactly.
    comps(r) = (r.energy, r.water, r.tracers.co2, r.tracers.o2)
    close(a, b) = all(isapprox.(a, b; rtol = 8 * eps(FT), atol = 8 * eps(FT) * FT(100)))
    @test all(close.(comps.(getproperty.(Array(gpu.tend), :rates)), comps.(getproperty.(ref.tend, :rates))))
    @test all(close.(comps.(getproperty.(Array(gpu.steady), :rates)), comps.(getproperty.(ref.steady, :rates))))
    @test all(isapprox.(getproperty.(Array(gpu.moist), :rho_energy), getproperty.(ref.moist, :rho_energy);
        rtol = 8 * eps(FT), atol = 800 * eps(FT)))
    @test getproperty.(getproperty.(Array(gpu.tend), :rates), :water) ==
          getproperty.(getproperty.(ref.tend, :rates), :water)
    @test all(Array(gpu.ok))
    @test sum(gpu.tend) isa VolumeSource
end

@testset "ClimaCore fields $(nameof(typeof(device))) $FT" for device in
        (ClimaComms.CPUSingleThreaded(), ClimaComms.CUDADevice()), FT in (Float32, Float64)
    height = FT[0, 1, 3, 6, 10]
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(first(height)),
        CC.Geometry.ZPoint(last(height));
        boundary_names = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, CC.Geometry.ZPoint.(height))
    space = CC.Spaces.CenterFiniteDifferenceSpace(device, mesh)
    fld(v) = (f = CC.Fields.Field(FT, space); parent(f) .= FT(v); f)
    out = pipeline(fld(80), fld(4e-5), fld(2.5e6), fld(-1e-5), fld(1e-5), fld(20), fld(40),
        fld(1.5), fld(0.8), fld(0.2))
    @test out.tend isa CC.Fields.Field
    e = Array(parent(out.tend.rates.energy))
    expected = (FT(20) - FT(2) * FT(40)) / FT(1.5) + (FT(80) + FT(4e-5) * FT(2.5e6)) * FT(0.8)
    @test all(x -> isapprox(x, expected; rtol = 8 * eps(FT)), e)
    @test all(Array(parent(out.ok)))
    @test eltype(parent(out.moist.tracers.o2)) === FT
    @test eltype(parent(out.steady.rates.tracers.co2)) === FT
end
