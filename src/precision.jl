# ============================================================================
# Scalar hooks and the precision policy
# ============================================================================
#
# The four hooks are exact for plain floating-point numbers. A number type that
# carries derivative information (ForwardDiff.Dual through the package extension
# CanopyExchangeForwardDiffExt) overloads them; nothing else in the package
# inspects the internals of a scalar type.

"""
    CanopyExchange.primal_type(::Type{T}) -> Type

Floating-point type that carries the value of scalar type `T`. The default returns
`T`; a derivative-carrying number type returns the floating-point type of its value
(for `ForwardDiff.Dual{Tag,Float64,N}` this is `Float64`).

The precision policy compares `primal_type` of every scalar that meets in a record
or kernel: they must agree and be an `AbstractFloat`.

# Examples
```jldoctest
julia> CanopyExchange.primal_type(Float32)
Float32
```
"""
primal_type(::Type{T}) where {T} = T

"""
    CanopyExchange.primal(x)

Value of scalar `x` without derivative information, in the units of `x`. Used only
for comparisons in validity checks (layer thickness, area densities, gas fraction),
so a check never depends on how a number type compares derivative components.

# Examples
```jldoctest
julia> CanopyExchange.primal(2.5f0)
2.5f0
```
"""
primal(x) = x

"""
    CanopyExchange.all_finite(x) -> Bool

Whether the value and every derivative component of scalar `x` are finite. Used by
[`valid_exchange`](@ref) and by the geometry and gas-fraction checks.

# Examples
```jldoctest
julia> (CanopyExchange.all_finite(1.5), CanopyExchange.all_finite(NaN))
(true, false)
```
"""
all_finite(x) = isfinite(x)

"""
    CanopyExchange.exact_zero(x) -> Bool

Whether the value and every derivative component of scalar `x` are zero (either sign).
[`area_source`](@ref) uses it for the empty-layer rule: an area density that is an
exact zero contributes an exact zero, while a zero value with a nonzero derivative
component keeps the right derivative `rate`.

# Examples
```jldoctest
julia> (CanopyExchange.exact_zero(-0.0), CanopyExchange.exact_zero(1.0e-12))
(true, false)
```
"""
exact_zero(x) = iszero(x)

# Precision policy, in one place: every scalar that meets in a record or a kernel
# must carry the same floating-point type. Float64 with a Float64-valued dual is
# accepted (and promotes); Float32 with Float64 is not; integers, rationals,
# irrationals and booleans are not exchange scalars. The checks depend only on
# types, so they fold away in compiled kernels; the messages are static.
@inline _float_type(x) = primal_type(eltype(typeof(x)))

@inline function same_precision(a, rest...)
    P = _float_type(a)
    P <: AbstractFloat ||
        throw(ArgumentError("canopy exchange scalars must be floating point"))
    _same_precision(P, rest...)
    return nothing
end
@inline _same_precision(::Type) = nothing
@inline function _same_precision(::Type{P}, b, rest...) where {P}
    _float_type(b) === P ||
        throw(ArgumentError("mixed floating-point precision in canopy exchange"))
    return _same_precision(P, rest...)
end

# Whether two scalar types differ at most by derivative information (container promotion).
@inline _promotable(::Type{A}, ::Type{B}) where {A,B} =
    primal_type(A) === primal_type(B) && primal_type(A) <: AbstractFloat

# Validity checks on primal values; derivative components only need to be finite.
@inline function _checked(x, bound_ok, message)
    all_finite(x) && bound_ok(primal(x)) || throw(ArgumentError(message))
    return x
end
@inline _positive(x) = x > zero(x)
@inline _nonnegative(x) = x >= zero(x)
@inline _unit_interval(x) = zero(x) < x <= one(x)
