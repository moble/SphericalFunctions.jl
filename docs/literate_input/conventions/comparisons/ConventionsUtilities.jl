@testmodule ConventionsUtilities begin
    import FastDifferentiation
    import SphericalFunctions
    import SphericalFunctions.Deprecated
    using Quaternionic: Rotor

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

    # `∂ⁿ(f, n)` returns the `n`th derivative of the single-variable function `f`, as a
    # function.  The derivative is computed symbolically with `FastDifferentiation`, so that
    # formulas quoted from the literature in terms of derivatives (e.g., Rodrigues' formula
    # for the Legendre polynomials) can be transcribed literally.  The returned function
    # accepts either a number or a `FastDifferentiation` variable, so that these derivatives
    # can be nested (e.g., ``d^m/dx^m`` of ``d^ℓ/dx^ℓ``).
    function ∂ⁿ(f, n)
        n < 0 && throw(ArgumentError("n=$n must be non-negative"))
        function (x)
            if x isa FastDifferentiation.Node
                n == 0 ? f(x) : FastDifferentiation.derivative(f(x), (x for _ ∈ 1:n)...)
            else
                v = FastDifferentiation.make_variables(:x)[1]
                expr = n == 0 ? f(v) : FastDifferentiation.derivative(f(v), (v for _ ∈ 1:n)...)
                FastDifferentiation.make_function([expr,], [v,])(x)[1]
            end
        end
    end

    # Reference implementations of *this package's settled conventions*, as documented on the
    # conventions "Summary" page.  Every comparison page compares the literature against these
    # functions, so that the pages state relations to the settled conventions rather than to
    # whatever the code happens to compute today.
    #
    # The scalar `Deprecated.Y` and `Deprecated.d` already agree with the settled conventions.
    # `Deprecated.D`, however, still carries the pre-3.0 convention, which is the complex
    # conjugate of the settled one, 𝔇ˡₘ′ₘ(α, β, γ) = exp(-𝒾 m′ α) dˡₘ′ₘ(β) exp(-𝒾 m γ).  The
    # `conj` below compensates for that.
    #
    # TODO: Remove the `conj` calls in `D` when `Deprecated.D!` is flipped to the settled
    # convention (item 4 of the implementation checklist in `src/redesign/README.md`).  The
    # `@test_broken` in the "ConventionsUtilities reference conventions" testitem below will
    # start passing at that point, as a reminder.
    Y(s, ℓ, m, θ, ϕ) = Deprecated.Y(s, ℓ, m, θ, ϕ)
    Y(ℓ, m, θ, ϕ) = Deprecated.Y(ℓ, m, θ, ϕ)
    d(ℓ, m′, m, β) = Deprecated.d(ℓ, m′, m, β)
    D(ℓ, m′, m, α, β, γ) = conj(Deprecated.D(ℓ, m′, m, α, β, γ))
    D(ℓ, m′, m, R::Rotor) = conj(Deprecated.D_matrices(R, ℓ)[Deprecated.WignerDindex(ℓ, m′, m)])

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


@testitem "∂ⁿ" setup=[ConventionsUtilities] begin
    # `∂ⁿ` should reproduce ordinary derivatives, including when nested.
    import .ConventionsUtilities: ∂ⁿ
    for x ∈ (-0.9, -0.3, 0.0, 0.4, 0.8)
        @test ∂ⁿ(x -> x^4, 0)(x) ≈ x^4
        @test ∂ⁿ(x -> x^4, 1)(x) ≈ 4x^3
        @test ∂ⁿ(x -> x^4, 3)(x) ≈ 24x
        @test ∂ⁿ(x -> (x^2 - 1)^2, 2)(x) ≈ 12x^2 - 4
        @test ∂ⁿ(∂ⁿ(x -> (x^2 - 1)^2, 1), 1)(x) ≈ 12x^2 - 4
        @test ∂ⁿ(x -> sin(x), 2)(x) ≈ -sin(x)
    end
end


@testitem "ConventionsUtilities reference conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin
    # The reference functions `D`, `d`, and `Y` in `ConventionsUtilities` are the oracle for
    # every comparison page.  Here we pin them to the explicit formulas on the conventions
    # "Summary" page, so that the pages are guaranteed to compare against the *documented*
    # conventions.
    import .ConventionsUtilities: D, d, Y, 𝒾
    import SphericalFunctions.Deprecated
    import Quaternionic: from_euler_angles

    ϵₐ = 100eps()
    ϵᵣ = 100eps()
    ℓₘₐₓ = 4
    n = 5

    # 𝔇ˡₘ′ₘ(α, β, γ) = exp(-𝒾 m′ α) dˡₘ′ₘ(β) exp(-𝒾 m γ)
    for (α, β, γ) ∈ αβγrange(Float64, n)
        for (ℓ, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
            @test D(ℓ, m′, m, α, β, γ) ≈ cis(-m′*α) * d(ℓ, m′, m, β) * cis(-m*γ) atol=ϵₐ rtol=ϵᵣ
        end
    end

    # The rotor form of `D` agrees with the Euler-angle form
    for (α, β, γ) ∈ αβγrange(Float64, n)
        R = from_euler_angles(α, β, γ)
        for (ℓ, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
            @test D(ℓ, m′, m, R) ≈ D(ℓ, m′, m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
        end
    end

    # Explicit d⁽¹⁾ matrix from the "Details" page (rows m′, columns m, decreasing from 1 to -1)
    for β ∈ βrange(Float64, n)
        @test d(1,  1,  1, β) ≈ (1 + cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
        @test d(1,  1,  0, β) ≈ -sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d(1,  1, -1, β) ≈ (1 - cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
        @test d(1,  0,  1, β) ≈ sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d(1,  0,  0, β) ≈ cos(β) atol=ϵₐ rtol=ϵᵣ
        @test d(1,  0, -1, β) ≈ -sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d(1, -1,  1, β) ≈ (1 - cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
        @test d(1, -1,  0, β) ≈ sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d(1, -1, -1, β) ≈ (1 + cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
    end

    # ₛYₗₘ(θ, ϕ) agrees with the explicit sum on the "Summary" page (implemented in `sYlm`
    # from the `Utilities` snippet), and with (-1)^s √((2ℓ+1)/4π) conj(𝔇ˡₘ,₋ₛ(ϕ, θ, 0))
    for (θ, ϕ) ∈ θϕrange(Float64, n)
        for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, 2)
            @test Y(s, ℓ, m, θ, ϕ) ≈ sYlm(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Y(s, ℓ, m, θ, ϕ) ≈ (-1)^s * √((2ℓ+1)/(4π)) * conj(D(ℓ, m, -s, ϕ, θ, 0)) atol=ϵₐ rtol=ϵᵣ
        end
        for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
            @test Y(ℓ, m, θ, ϕ) ≈ Y(0, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Y(ℓ, m, θ, ϕ) ≈ √((2ℓ+1)/(4π)) * conj(D(ℓ, m, 0, ϕ, θ, 0)) atol=ϵₐ rtol=ϵᵣ
        end
    end

    # Reminder: `Deprecated.D` still carries the pre-3.0 (conjugate) convention.  When it is
    # flipped, this test will start passing, and the `conj` in `ConventionsUtilities.D` should
    # be removed.
    @test_broken Deprecated.D(1, 1, 0, 0.3, 0.4, 0.5) ≈ D(1, 1, 0, 0.3, 0.4, 0.5)
end
