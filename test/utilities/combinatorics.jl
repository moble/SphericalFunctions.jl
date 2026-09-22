# Tests of the combinatorial helpers in `src/utilities/utils.jl`.
#
# `sqrtbinomial` is published on `docs/src/20-interface/05-utilities.md`, and §9 of the v3 design
# memo directs callers to it ("use `logbinomial`/`sqrtbinomial` from `utils.jl`, never
# `binomial`").  Its only test used to be "Preliminaries: sqrtbinomial" in the v2
# `test/deprecated/ssht.jl`, which went with the `Deprecated` module in 3.0; this item
# replaces it, and extends it to `logbinomial` and to the overflow regime that is the whole
# reason the function exists.

@testitem "Combinatorics: sqrtbinomial and logbinomial" begin
    import SphericalFunctions: sqrtbinomial, logbinomial
    using DoubleFloats: Double64

    # Against exact `BigInt` binomials.  `binomial(::Int, ::Int)` overflows above n ≈ 66,
    # which is exactly what this function is for, so the reference is always `big`.
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        for ℓ ∈ (1, 2, 3, 4, 5, 13, 64, 65, 66, 67, 1025)
            for s ∈ -2:2
                a = sqrtbinomial(2ℓ, ℓ - s, T)
                b = T(√big(binomial(big(2ℓ), big(ℓ - s))))
                @test a isa T
                @test a ≈ b rtol=√eps(T)
            end
        end
    end

    # The edge cases of `logbinomial`, which `sqrtbinomial` reaches through its `k == 0`,
    # `k == n` and `k == 1` branches, and the k > n÷2 reflection.
    for n ∈ (0, 1, 2, 7, 100, 1025)
        @test logbinomial(n, 0) == 0
        @test logbinomial(n, n) == 0
        @test sqrtbinomial(n, 0) == 1
        @test sqrtbinomial(n, n) == 1
        if n ≥ 1
            @test logbinomial(n, 1) ≈ log(n)
            @test logbinomial(n, n - 1) ≈ log(n)
            @test sqrtbinomial(n, 1) ≈ √n
        end
        for k ∈ 0:min(n, 20)
            # binomial(n, k) == binomial(n, n-k), through the reflection branch
            @test logbinomial(n, k) ≈ logbinomial(n, n - k) atol=1e-12
            @test logbinomial(n, k) ≈ log(big(binomial(big(n), big(k)))) rtol=1e-12
        end
    end

    # It really does beat `binomial`: n = 2050 overflows an `Int`, and the value is far
    # outside `Float64`'s range, yet its square root is not.
    @test_throws OverflowError binomial(2050, 1025)
    @test sqrtbinomial(2050, 1025) ≈ Float64(√big(binomial(big(2050), big(1025)))) rtol=1e-10
    @test isfinite(sqrtbinomial(2050, 1025))

    # The computation type is honored, and defaults to Float64.
    @test sqrtbinomial(10, 5) isa Float64
    @test sqrtbinomial(10, 5, BigFloat) isa BigFloat
    @test Float64(sqrtbinomial(10, 5, BigFloat)) ≈ sqrtbinomial(10, 5)
end
