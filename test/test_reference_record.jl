using Test, TOML, CanopyExchange
@isdefined(exchange_fixture) || include("exchange_fixture.jl")
include("reference_record.jl")

# Bit-for-bit comparison with the outputs frozen from v0.1.1, commit c35ad8d. Provenance,
# keys and units of the record: test/reference/README.md and write_reference_record.jl.
@testset "Frozen v0.1.1 reference record (leaf, stem, datum, steady, moist air)" begin
    path = joinpath(@__DIR__, "reference", "exchange_v0_1_1.toml")
    frozen = TOML.parsefile(path)["profiles"]
    computed = reference_profiles()
    @test length(frozen) == length(computed)
    for profile in computed
        matches = filter(
            p -> p["precision"] == profile["precision"] && p["oxygen"] == profile["oxygen"],
            frozen,
        )
        @test length(matches) == 1
        stored = only(matches)
        @test sort(collect(keys(stored))) == sort(collect(keys(profile)))
        for key in keys(profile)
            @test stored[key] == profile[key]
        end
    end
end
