module SphericalFunctionsDoubleFloatsExt

# The conversions of a `HalfOddInteger` to the DoubleFloats types.  `DoubleFloats` defines
# `Double64(x::T) where {T<:Real}`, which would reach for a `promote_rule` that
# `HalfOddInteger` deliberately lacks, so each type is given a method of its own.  A method
# for the parametric `DoubleFloat{T}` would not do: it is exactly as specific as
# `DoubleFloats`' own, and the two would be ambiguous.  The value is converted as the equal
# `Rational`, which `DoubleFloats` converts exactly whenever the value is representable; its
# conversion of an `Int` passes through the float type of the high word, and so rounds a
# numerator wider than that type's significand.

import SphericalFunctions: HalfOddInteger
import DoubleFloats: Double16, Double32, Double64

(::Type{Double16})(a::HalfOddInteger) = Double16(Rational(a))
(::Type{Double32})(a::HalfOddInteger) = Double32(Rational(a))
(::Type{Double64})(a::HalfOddInteger) = Double64(Rational(a))

end # module SphericalFunctionsDoubleFloatsExt
