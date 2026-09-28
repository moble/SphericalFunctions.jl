@testmodule ConventionsUtilities begin
    import FastDifferentiation

    const 𝒾 = im

    struct Factorial end
    Base.:*(n::Integer, ::Factorial) = factorial(big(n))
    function Base.:*(n::Rational, ::Factorial)
        if denominator(n) == 1
            return factorial(big(numerator(n)))
        else
            throw(ArgumentError("Cannot compute factorial of a non-integer rational"))
        end
    end
    const ❗ = Factorial()

    # `dʲsin²ᵏθdcosθʲ(; j, k, θ)` is the `j`th derivative of sin²ᵏθ with respect to cos θ,
    # evaluated at `θ`.  The derivative is computed symbolically and compiled the first time
    # each `(j, k)` is needed.  The compiled function is kept in `dʲsin²ᵏθdcosθʲ_functions`
    # and reused by every later evaluation, because the differentiation and compilation cost
    # several thousand times as much as an evaluation.
    const dʲsin²ᵏθdcosθʲ_functions = Dict{Tuple{Int, Int}, Any}()
    function dʲsin²ᵏθdcosθʲ(;j, k, θ)
        if j < 0
            throw(ArgumentError("j=$j must be non-negative"))
        end
        if j == 0
            return sin(θ)^(2k)
        end
        ∂ₓʲfᵏ = get!(dʲsin²ᵏθdcosθʲ_functions, (j, k)) do
            x = FastDifferentiation.make_variables(:x)[1]
            expr = FastDifferentiation.derivative((1 - x^2)^k, (x for _ ∈ 1:j)...)
            FastDifferentiation.make_function([expr,], [x,])
        end
        return ∂ₓʲfᵏ(cos(θ))[1]
    end

    # `∂ⁿ(f, n)` returns the `n`th derivative of the single-variable function `f`, as a
    # function.  The derivative is computed symbolically with `FastDifferentiation`, so that
    # formulas quoted from the literature in terms of derivatives (e.g., Rodrigues' formula
    # for the Legendre polynomials) can be transcribed literally.  The returned function
    # accepts either a number or a `FastDifferentiation` variable, so that these derivatives
    # can be nested (e.g., ``d^m/dx^m`` of ``d^ℓ/dx^ℓ``).
    #
    # As for `dʲsin²ᵏθdcosθʲ`, the compiled numerical function is kept in `∂ⁿ_functions`,
    # keyed on `(f, n)`.  The callers build a new closure such as `x -> (x^2 - 1)^ℓ` on
    # every evaluation, but closures that capture only bits values (here, `ℓ`) compare and
    # hash by content, so each distinct derivative is still compiled only once.  A closure
    # that captured a mutable value would compare by identity instead, and would add an
    # entry on every call.
    const ∂ⁿ_functions = Dict{Any, Any}()
    function ∂ⁿ(f, n)
        n < 0 && throw(ArgumentError("n=$n must be non-negative"))
        function (x)
            if x isa FastDifferentiation.Node
                n == 0 ? f(x) : FastDifferentiation.derivative(f(x), (x for _ ∈ 1:n)...)
            else
                ∂ⁿf = get!(∂ⁿ_functions, (f, n)) do
                    v = FastDifferentiation.make_variables(:x)[1]
                    expr = n == 0 ? f(v) : FastDifferentiation.derivative(f(v), (v for _ ∈ 1:n)...)
                    FastDifferentiation.make_function([expr,], [v,])
                end
                ∂ⁿf(x)[1]
            end
        end
    end

end


@testitem "dʲsin²ᵏθdcosθʲ" setup=[ConventionsUtilities, Utilities] begin
    # dʲsin²ᵏθdcosθʲ is intended to represent the jth derivative of sin(θ)^(2k) with respect
    # to cos(θ).  We can compare it to some actual derivatives of sin(θ)^(2k) to verify its
    # correctness.
    import .ConventionsUtilities: dʲsin²ᵏθdcosθʲ
    import .Utilities: βrange
    using Random
    rng = Random.Xoshiro(1234)
    for θ ∈ βrange(rng, Float64, 15)
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
    # Every comparison page compares the literature with this package's own functions, called
    # as a user would call them: a calculator — or, where blocks at several points are needed
    # at once, the full series — at each sample point, with its blocks compared element by
    # element.  Here we pin those functions to the explicit formulas on the conventions
    # "Summary" page, so that the pages are guaranteed to compare against the *documented*
    # conventions.
    #
    # The blocks of a half-integer calculator are labelled with `HalfOddInteger`s, which are
    # for indexing and deliberately refuse arithmetic with floats and `Rational`s.  The
    # formulas here do that arithmetic, so the labels are converted with `Rational` first.
    import .ConventionsUtilities: 𝒾, ❗
    import SphericalFunctions: D, d, DCalculator, dCalculator, sYlmCalculator, YlmCalculator
    import .Utilities: αβγrange, βrange, θϕrange, sYlm_closed_form

    ϵₐ = 100eps()
    ϵᵣ = 100eps()
    ℓₘₐₓ = 4
    n = 5

    # Wigner's formula for d from the "Summary" page.  The factorials are computed exactly
    # and the sum is evaluated in `BigFloat`, because the alternating sum loses accuracy in
    # `Float64` as ℓ grows.  Every factorial and exponent has an integer argument, even for
    # half-integer indices.
    function d_Wigner(ℓ, m′, m, β)
        β = big(β)
        sum(
            (-1)^Int(k - m + m′)
            * √((ℓ+m)❗ * (ℓ-m)❗ * (ℓ+m′)❗ * (ℓ-m′)❗)
            / ((ℓ+m-k)❗ * (k)❗ * (ℓ-m′-k)❗ * (k-m+m′)❗)
            * cos(β/2)^Int(2ℓ+m-m′-2k)
            * sin(β/2)^Int(2k-m+m′)
            for k ∈ max(0, m-m′):min(ℓ+m, ℓ-m′)
        )
    end

    # Every element of d agrees with Wigner's formula, for integer ℓ ≤ 8 and half-integer
    # ℓ ≤ 15/2
    for β ∈ βrange(rng, Float64, n)
        for (ℓ, dˡ) ∈ dCalculator(β, 8), m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
            @test dˡ[m′, m] ≈ d_Wigner(ℓ, m′, m, β) atol=ϵₐ rtol=ϵᵣ
        end
        for (ℓ, dˡ) ∈ dCalculator(β, 15//2)
            ℓ = Rational(ℓ)
            for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                @test dˡ[m′, m] ≈ d_Wigner(ℓ, m′, m, β) atol=ϵₐ rtol=ϵᵣ
            end
        end
    end

    # 𝔇ˡₘ′ₘ(α, β, γ) = exp(-𝒾 m′ α) dˡₘ′ₘ(β) exp(-𝒾 m γ), for integer and half-integer ℓ
    for (α, β, γ) ∈ αβγrange(rng, Float64, n)
        𝔇, 𝒹 = DCalculator(α, β, γ, ℓₘₐₓ), dCalculator(β, ℓₘₐₓ)
        for ((ℓ, 𝔇ˡ), (_, dˡ)) ∈ zip(𝔇, 𝒹), m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
            @test 𝔇ˡ[m′, m] ≈ cis(-m′*α) * dˡ[m′, m] * cis(-m*γ) atol=ϵₐ rtol=ϵᵣ
        end
        𝔇, 𝒹 = DCalculator(α, β, γ, 7//2), dCalculator(β, 7//2)
        for ((ℓ, 𝔇ˡ), (_, dˡ)) ∈ zip(𝔇, 𝒹)
            ℓ = Rational(ℓ)
            for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                @test 𝔇ˡ[m′, m] ≈ cis(-m′*α) * dˡ[m′, m] * cis(-m*γ) atol=ϵₐ rtol=ϵᵣ
            end
        end
    end

    # Explicit d⁽¹⁾ matrix from the "Details" page (rows m′, columns m, decreasing from 1 to -1)
    for β ∈ βrange(rng, Float64, n)
        d¹ = d(β, 1)[1]
        @test d¹[ 1,  1] ≈ (1 + cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
        @test d¹[ 1,  0] ≈ -sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d¹[ 1, -1] ≈ (1 - cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
        @test d¹[ 0,  1] ≈ sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d¹[ 0,  0] ≈ cos(β) atol=ϵₐ rtol=ϵᵣ
        @test d¹[ 0, -1] ≈ -sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d¹[-1,  1] ≈ (1 - cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
        @test d¹[-1,  0] ≈ sin(β) / √2 atol=ϵₐ rtol=ϵᵣ
        @test d¹[-1, -1] ≈ (1 + cos(β)) / 2 atol=ϵₐ rtol=ϵᵣ
    end

    # ₛYₗₘ(θ, ϕ) agrees with the explicit sum on the "Summary" page (implemented in
    # `sYlm_closed_form` from the `Utilities` module), and with (-1)^s √((2ℓ+1)/4π)
    # conj(𝔇ˡₘ,₋ₛ(ϕ, θ, 0)).  A calculator for the range of spin weights -2:2 gives blocks
    # indexed [s, m], of which the rows with |s| ≤ ℓ are the harmonics.
    for (θ, ϕ) ∈ θϕrange(rng, Float64, n)
        𝔇 = D(ϕ, θ, 0, ℓₘₐₓ)
        for (ℓ, Yˡ) ∈ sYlmCalculator(θ, ϕ, ℓₘₐₓ, -2:2), s ∈ -min(ℓ, 2):min(ℓ, 2), m ∈ -ℓ:ℓ
            @test Yˡ[s, m] ≈ sYlm_closed_form(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Yˡ[s, m] ≈ (-1)^s * √((2ℓ+1)/(4π)) * conj(𝔇[ℓ][m, -s]) atol=ϵₐ rtol=ϵᵣ
        end
        for (ℓ, Yˡ) ∈ YlmCalculator(θ, ϕ, ℓₘₐₓ), m ∈ -ℓ:ℓ
            @test Yˡ[m] ≈ sYlm_closed_form(0, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Yˡ[m] ≈ √((2ℓ+1)/(4π)) * conj(𝔇[ℓ][m, 0]) atol=ϵₐ rtol=ϵᵣ
        end
        # For half-integer s, the (-1)^s of the definition is e^{iπs}
        𝔇 = D(ϕ, θ, 0, 7//2)
        for (ℓ, Yˡ) ∈ sYlmCalculator(θ, ϕ, 7//2, -7//2:7//2)
            ℓ = Rational(ℓ)
            for s ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                @test Yˡ[s, m] ≈ cis(π*s) * √((2ℓ+1)/(4π)) * conj(𝔇[ℓ][m, -s]) atol=ϵₐ rtol=ϵᵣ
            end
        end
    end
end
