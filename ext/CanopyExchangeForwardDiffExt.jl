# Scalar hooks for ForwardDiff dual numbers, through ForwardDiff's public interface
# (`value`, and `partials` as an indexable collection). The checks of CanopyExchange
# compare the value of a dual; finiteness, exact zero and equality include every
# derivative component. This behaves identically on ForwardDiff 0.10 and 1.x, whose
# `iszero` and `==` differ.
module CanopyExchangeForwardDiffExt

using ForwardDiff: Dual, value, partials
import CanopyExchange: primal_type, primal, all_finite, exact_zero, scalar_equal

primal_type(::Type{<:Dual{T,V}}) where {T,V} = primal_type(V)
primal(d::Dual) = primal(value(d))
all_finite(d::Dual) = all_finite(value(d)) && all(all_finite, _partials(d))
exact_zero(d::Dual) = exact_zero(value(d)) && all(exact_zero, _partials(d))

# Derivative components as a tuple, by indexing `partials(d)`; the length is the type
# parameter `N` of `Dual{T,V,N}`, so the tuple is unrolled and allocation free.
@inline _partials(d::Dual{T,V,N}) where {T,V,N} = ntuple(i -> partials(d)[i], Val(N))

# Equal values and equal derivative components; a plain number has zero derivatives.
# Duals of different tags (nested differentiation) are compared as constants.
scalar_equal(compare, a::Dual{T}, b::Dual{T}) where {T} =
    scalar_equal(compare, value(a), value(b)) &&
    _pairs_equal(compare, _partials(a), _partials(b))
scalar_equal(compare, a::Dual, b::Dual) =
    _constant_equal(compare, a, b) && _constant_equal(compare, b, a)
scalar_equal(compare, a::Dual, b::Real) = _constant_equal(compare, a, b)
scalar_equal(compare, a::Real, b::Dual) = _constant_equal(compare, b, a)

_constant_equal(compare, a::Dual, b) =
    scalar_equal(compare, value(a), b) && _all_zero(compare, _partials(a))
@inline _pairs_equal(compare, ::Tuple{}, ::Tuple{}) = true
@inline _pairs_equal(compare, ::Tuple, ::Tuple) = false
@inline _pairs_equal(compare, p::Tuple{Any,Vararg{Any,N}}, q::Tuple{Any,Vararg{Any,N}}) where {N} =
    scalar_equal(compare, first(p), first(q)) && _pairs_equal(compare, Base.tail(p), Base.tail(q))
@inline _all_zero(compare, ::Tuple{}) = true
@inline _all_zero(compare, p::Tuple) =
    scalar_equal(compare, first(p), zero(first(p))) && _all_zero(compare, Base.tail(p))

end
