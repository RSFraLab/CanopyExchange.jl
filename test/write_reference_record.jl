# Writes the frozen reference record of CanopyExchange outputs as IEEE bit strings.
#
# Provenance: test/reference/exchange_v0_1_1.toml was written by this script from
# CanopyExchange v0.1.1 (commit c35ad8d) with Julia 1.12.6:
#     julia --startup-file=no --project=. test/write_reference_record.jl test/reference/exchange_v0_1_1.toml
# Do NOT regenerate that file from refactored code: it is the fixed reference that
# later versions are compared against (test/test_reference_record.jl). Write any
# comparison output to a different path.
using TOML
using CanopyExchange
length(ARGS) == 1 || error("Usage: write_reference_record.jl output.toml")
include(joinpath(@__DIR__, "exchange_fixture.jl"))
include(joinpath(@__DIR__, "reference_record.jl"))
open(io -> TOML.print(io, Dict("profiles" => reference_profiles()); sorted = true), ARGS[1], "w")
println(ARGS[1])
