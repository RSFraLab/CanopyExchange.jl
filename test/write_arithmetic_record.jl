# Compare with the old actual package in a separate process/environment.
using TOML
length(ARGS) == 2 || error("Usage: write_arithmetic_record.jl legacy|canonical output.toml")
mode, output = ARGS
if mode == "legacy"
    @eval using CanopyColumn.ConservativeExchange
elseif mode == "canonical"
    @eval using CanopyExchange
else
    error("Unknown module owner")
end
include("exchange_fixture.jl")
records = []
for FT in (Float32, Float64), enabled in (false, true)
    case = exchange_fixture(FT, enabled)
    paired = leaf_air_sources.(case.leaves, case.geometry)
    tendency = layer_tendency.(case.boundaries[1:end-1], case.boundaries[2:end],
                              getproperty.(paired, :air), case.geometry)
    push!(records, Dict("precision" => string(FT), "oxygen" => enabled,
        "leaf_flux" => [bitstring.(collect(exchange_components(x.rates))) for x in case.leaves],
        "air_source" => [bitstring.(collect(exchange_components(x.air.rates))) for x in paired],
        "leaf_source" => [bitstring.(collect(exchange_components(x.leaf.rates))) for x in paired],
        "tendency" => [bitstring.(collect(exchange_components(x.rates))) for x in tendency]))
end
open(io -> TOML.print(io, Dict("profiles" => records); sorted=true), output, "w")
println(output)
