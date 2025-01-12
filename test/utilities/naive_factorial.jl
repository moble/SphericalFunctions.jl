"""Write factorials in the obvious way, for testing purposes.

This snippet lets us write things like `5❗` to get `120`, for example.  This should
probably only be used in tests, because factorials are not very efficient.  Therefore,
to ensure accuracy and avoid overflow, the argument is first converted to a `BigInt`,
making it even more inefficient.  This "big" conversion should be preserved as long as
possible to ensure accurate cancellations in the final result.

Use this snippet by including the following in your test file:

    include("../utilities/naive_factorial.jl")
    import .NaiveFactorials: ❗
"""
module NaiveFactorials
    struct Factorial end
    Base.:*(n::Integer, ::Factorial) = factorial(big(n))
    function Base.:*(n::Rational, ::Factorial) where {Rational}
        if denominator(n) == 1
            return factorial(big(numerator(n)))
        else
            throw(ArgumentError("Cannot compute factorial of a non-integer rational"))
        end
    end
    const ❗ = Factorial()
end
