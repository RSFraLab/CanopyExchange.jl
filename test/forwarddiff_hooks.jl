# Scalar hooks for ForwardDiff.Dual, defined here only while the package ships no
# ForwardDiff extension (the extension defines exactly these four methods).
using ForwardDiff
if Base.get_extension(CanopyExchange, :CanopyExchangeForwardDiffExt) === nothing
    CanopyExchange.primal_type(::Type{<:ForwardDiff.Dual{T,V}}) where {T,V} =
        CanopyExchange.primal_type(V)
    CanopyExchange.primal(d::ForwardDiff.Dual) = CanopyExchange.primal(ForwardDiff.value(d))
    CanopyExchange.all_finite(d::ForwardDiff.Dual) =
        CanopyExchange.all_finite(ForwardDiff.value(d)) &&
        all(CanopyExchange.all_finite, ForwardDiff.partials(d).values)
    CanopyExchange.exact_zero(d::ForwardDiff.Dual) =
        CanopyExchange.exact_zero(ForwardDiff.value(d)) &&
        all(CanopyExchange.exact_zero, ForwardDiff.partials(d).values)
end
