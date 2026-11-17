# ============================================================================
# Shared behaviour of every record: value equality, hashing, broadcast scalar
# ============================================================================

# Internal supertype of every record this package defines.
abstract type ScalarRecord end

Base.Broadcast.broadcastable(x::ScalarRecord) = Ref(x)

@inline _fields(x::T) where {T} = ntuple(i -> getfield(x, i), Val(fieldcount(T)))
_measure(_) = nothing
@inline _same_kind(a, b) =
    nameof(typeof(a)) === nameof(typeof(b)) && _measure(a) === _measure(b)

# Componentwise comparison through the `scalar_equal` hook: value and derivative
# components, identical on every ForwardDiff version. Tuple recursion keeps the
# result a `Bool` and the call free of allocations.
@inline _all_equal(compare, ::Tuple{}, ::Tuple{}) = true
@inline _all_equal(compare, ::Tuple, ::Tuple) = false
@inline _all_equal(compare, a::Tuple{Any,Vararg{Any,N}}, b::Tuple{Any,Vararg{Any,N}}) where {N} =
    _value_equal(compare, first(a), first(b)) &&
    _all_equal(compare, Base.tail(a), Base.tail(b))
@inline _value_equal(compare, x::NamedTuple, y::NamedTuple) =
    keys(x) === keys(y) && _all_equal(compare, values(x), values(y))
@inline _value_equal(compare, x::ScalarRecord, y::ScalarRecord) = _record_equal(compare, x, y)
@inline _value_equal(compare, x, y) = scalar_equal(compare, x, y)
@inline _record_equal(compare, a, b) =
    _same_kind(a, b) && _all_equal(compare, _fields(a), _fields(b))

Base.:(==)(a::ScalarRecord, b::ScalarRecord) = _record_equal(==, a, b)
Base.isequal(a::ScalarRecord, b::ScalarRecord) = _record_equal(isequal, a, b)

# Hash the primal value with the zero sign normalised, so records that compare equal
# (`-0.0 == 0.0`, or a dual with zero derivative and a plain float) hash equal.
Base.hash(a::ScalarRecord, h::UInt) =
    _hash_values(_fields(a), hash(_measure(a), hash(nameof(typeof(a)), h)))
@inline _hash_values(::Tuple{}, h::UInt) = h
@inline _hash_values(t::Tuple, h::UInt) = _hash_values(Base.tail(t), _hash_value(first(t), h))
@inline _hash_value(x::ScalarRecord, h::UInt) = hash(x, h)
@inline _hash_value(x::NamedTuple, h::UInt) = _hash_values(values(x), hash(keys(x), h))
@inline _hash_value(x, h::UInt) = (p = primal(x); hash(ifelse(iszero(p), zero(p), p), h))
