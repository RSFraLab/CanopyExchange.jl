const UnsupportedDiscreteScalar = Union{Bool,Integer,Rational,AbstractIrrational}

"""
    primal_float_type(::Type{T}) -> Type

Return the underlying leaf floating-point type of continuous scalar type `T`.

For `T <: AbstractFloat`, return `T` directly. For concrete composite `Real`
types whose first field stores the primal value (such as `ForwardDiff.Dual`),
recurse until reaching the underlying `AbstractFloat`.
"""
@inline primal_float_type(::Type{T}) where {T<:AbstractFloat} = T
@inline function primal_float_type(::Type{T}) where {T<:Real}
    isconcretetype(T) || return T
    T <: UnsupportedDiscreteScalar && return T
    if fieldcount(T) > 0
        F1 = fieldtype(T, 1)
        if F1 <: Real && F1 !== T
            return primal_float_type(F1)
        end
    end
    return T
end
@inline primal_float_type(::Type{T}) where {T} = T

"""
    primal_scalar(x::Real) -> AbstractFloat

Extract the underlying primal floating-point value from continuous scalar `x`.
"""
@inline primal_scalar(x::AbstractFloat) = x
@inline function primal_scalar(x::T) where {T<:Real}
    if isconcretetype(T) && fieldcount(T) > 0
        v = getfield(x, 1)
        if v isa Real && typeof(v) !== T
            return primal_scalar(v)
        end
    end
    return x
end

"""
    is_continuous_scalar(::Type{T}) -> Bool

Return `true` if `T` is a concrete continuous `Real` scalar backed by an `AbstractFloat`.

Reject abstract types (`Real`, `AbstractFloat`), discrete scalar types (`Bool`,
`Integer`, `Rational`, `AbstractIrrational`), and composite duals built over
discrete scalar types.
"""
@inline function is_continuous_scalar(::Type{T}) where {T}
    isconcretetype(T) || return false
    T <: Real || return false
    T <: UnsupportedDiscreteScalar && return false
    P = primal_float_type(T)
    return isconcretetype(P) && (P <: AbstractFloat)
end

@inline function require_continuous_scalar(::Type{T}, ctor, args) where {T}
    is_continuous_scalar(T) || throw(MethodError(ctor, args))
    return nothing
end

"""
    compatible_scalar(::Type{A}, ::Type{B}) -> Bool

Return `true` if `A` and `B` are continuous scalars sharing the identical primal float type.

This permits automatic differentiation promotion (`Float64` with
`ForwardDiff.Dual{Tag,Float64,N}`) while rejecting mixed floating-point
precisions (`Float32` with `Float64`) and discrete scalar types.
"""
@inline function compatible_scalar(::Type{A}, ::Type{B}) where {A,B}
    return is_continuous_scalar(A) &&
           is_continuous_scalar(B) &&
           (primal_float_type(A) === primal_float_type(B))
end

@inline function promote_compatible_scalar(ctor, x::A) where {A}
    require_continuous_scalar(A, ctor, (x,))
    return A
end

@inline function promote_compatible_scalar(ctor, x::A, y::B) where {A,B}
    compatible_scalar(A, B) || throw(MethodError(ctor, (x, y)))
    FT = promote_type(A, B)
    is_continuous_scalar(FT) || throw(MethodError(ctor, (x, y)))
    return FT
end

@inline function promote_compatible_scalar(ctor, x::A, y::B, rest...) where {A,B}
    all_args = (x, y, rest...)
    compatible_scalar(A, B) || throw(MethodError(ctor, all_args))
    FT12 = promote_type(A, B)
    is_continuous_scalar(FT12) || throw(MethodError(ctor, all_args))
    return _promote_compatible_tail(ctor, all_args, FT12, rest...)
end

@inline _promote_compatible_tail(ctor, all_args, ::Type{FT}) where {FT} = FT
@inline function _promote_compatible_tail(ctor, all_args, ::Type{FT}, z::C, rest...) where {FT,C}
    compatible_scalar(FT, C) || throw(MethodError(ctor, all_args))
    FT_next = promote_type(FT, C)
    is_continuous_scalar(FT_next) || throw(MethodError(ctor, all_args))
    return _promote_compatible_tail(ctor, all_args, FT_next, rest...)
end

"""
    finite_scalar(x::Real) -> Bool

Return `true` if every leaf floating-point component of `x` (both primal and partials) is finite.

Unlike `Base.isfinite` on `ForwardDiff.Dual`, which inspects only the primal
value, `finite_scalar` unrolls all composite struct and tuple fields at compile
time in zero allocations and rejects `NaN` or `Inf` in derivative partials.
"""
@inline finite_scalar(x::AbstractFloat) = isfinite(x)
@generated function finite_scalar(x::T) where {T}
    T <: AbstractFloat && return :(isfinite(x))
    if T <: Tuple
        n = fieldcount(T)
        n == 0 && return :(true)
        checks = [:(finite_scalar(getfield(x, $i))) for i in 1:n]
        return Expr(:&&, checks...)
    elseif isconcretetype(T) && fieldcount(T) > 0
        checks = [:(finite_scalar(getfield(x, $i))) for i in 1:fieldcount(T)]
        return Expr(:&&, checks...)
    elseif T <: Real
        return :(isfinite(x))
    else
        return :(throw(MethodError(finite_scalar, (x,))))
    end
end

"""
    exact_zero_scalar(x::Real) -> Bool

Return `true` if every leaf floating-point component of `x` (both primal and partials) is zero.

Evaluate `iszero` structurally across all leaf `AbstractFloat` components in zero
allocations, handling IEEE `±0.0`, `BigFloat(0)`, and nested `ForwardDiff.Dual`
numbers without pointer-identity comparisons (`===`).
"""
@inline exact_zero_scalar(x::AbstractFloat) = iszero(x)
@generated function exact_zero_scalar(x::T) where {T}
    T <: AbstractFloat && return :(iszero(x))
    if T <: Tuple
        n = fieldcount(T)
        n == 0 && return :(true)
        checks = [:(exact_zero_scalar(getfield(x, $i))) for i in 1:n]
        return Expr(:&&, checks...)
    elseif isconcretetype(T) && fieldcount(T) > 0
        checks = [:(exact_zero_scalar(getfield(x, $i))) for i in 1:fieldcount(T)]
        return Expr(:&&, checks...)
    elseif T <: Real
        return :(iszero(x))
    else
        return :(throw(MethodError(exact_zero_scalar, (x,))))
    end
end

@inline function _hash_scalar(x::Real, h::UInt)
    p = primal_scalar(x)
    return hash(ifelse(iszero(p), zero(p), p), h)
end
