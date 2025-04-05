@testmodule ConventionsUtilities begin
    import FastDifferentiation

    const 𝒾 = im

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

    function dʲsin²ᵏθdcosθʲ(;j, k, θ)
        if j < 0
            throw(ArgumentError("j=$j must be non-negative"))
        end
        if j == 0
            return sin(θ)^(2k)
        end
        x = FastDifferentiation.make_variables(:x)[1]
        ∂ₓʲfᵏ = FastDifferentiation.derivative((1 - x^2)^k, (x for _ ∈ 1:j)...)
        return FastDifferentiation.make_function([∂ₓʲfᵏ,], [x,])(cos(θ))[1]
    end

end


@testitem "dʲsin²ᵏθdcosθʲ" setup=[ConventionsUtilities, Utilities] begin
    # dʲsin²ᵏθdcosθʲ is intended to represent the jth derivative of sin(θ)^(2k) with respect
    # to cos(θ).  We can compare it to some actual derivatives of sin(θ)^(2k) to verify its
    # correctness.
    import .ConventionsUtilities: dʲsin²ᵏθdcosθʲ
    for θ ∈ βrange(Float64, 15)
        @test dʲsin²ᵏθdcosθʲ(j=0, k=0, θ=θ) ≈ 1
        @test dʲsin²ᵏθdcosθʲ(j=0, k=1, θ=θ) ≈ sin(θ)^2
        @test dʲsin²ᵏθdcosθʲ(j=1, k=1, θ=θ) ≈ -2cos(θ)
        @test dʲsin²ᵏθdcosθʲ(j=2, k=1, θ=θ) ≈ -2
        @test dʲsin²ᵏθdcosθʲ(j=3, k=1, θ=θ) ≈ 0
        @test dʲsin²ᵏθdcosθʲ(j=0, k=2, θ=θ) ≈ sin(θ)^4
        @test dʲsin²ᵏθdcosθʲ(j=1, k=2, θ=θ) ≈ -4 * cos(θ) * sin(θ)^2 atol=4eps()
        @test dʲsin²ᵏθdcosθʲ(j=2, k=2, θ=θ) ≈ -4 + 12cos(θ)^2
        @test dʲsin²ᵏθdcosθʲ(j=3, k=2, θ=θ) ≈ 24cos(θ)
        @test dʲsin²ᵏθdcosθʲ(j=4, k=2, θ=θ) ≈ 24
        @test dʲsin²ᵏθdcosθʲ(j=5, k=2, θ=θ) ≈ 0
        @test dʲsin²ᵏθdcosθʲ(j=0, k=3, θ=θ) ≈ sin(θ)^6
        @test dʲsin²ᵏθdcosθʲ(j=1, k=3, θ=θ) ≈ -6 * cos(θ) * sin(θ)^4 atol=4eps()
        @test dʲsin²ᵏθdcosθʲ(j=2, k=3, θ=θ) ≈ -6 * sin(θ)^4 + 24cos(θ)^2 * sin(θ)^2 atol=100eps()
    end
end
