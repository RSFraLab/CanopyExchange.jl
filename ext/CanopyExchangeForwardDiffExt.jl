# Scalar hooks for ForwardDiff dual numbers. The checks of CanopyExchange compare the
# value of a dual; finiteness and exact zero include every derivative component. This
# behaves identically on ForwardDiff 0.10 and 1.x, whose `iszero` and `==` differ.
module CanopyExchangeForwardDiffExt

using ForwardDiff: Dual, value, partials
import CanopyExchange: primal_type, primal, all_finite, exact_zero

primal_type(::Type{<:Dual{T,V}}) where {T,V} = primal_type(V)
primal(d::Dual) = primal(value(d))
all_finite(d::Dual) = all_finite(value(d)) && all(all_finite, partials(d).values)
exact_zero(d::Dual) = exact_zero(value(d)) && all(exact_zero, partials(d).values)

end
