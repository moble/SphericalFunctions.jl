# Tests of the helpers in `src/utilities/utils.jl`: the combinatorial functions, and the
# functions that read the element type and the number of rotors from rotor data.
#
# `sqrtbinomial` is published on `docs/src/20-interface/05-utilities.md` as the way to form
# the square root of a binomial coefficient too large for `binomial`, and this item tests it
# and `logbinomial` in the overflow regime that is the whole reason the functions exist.

@testitem "Combinatorics: sqrtbinomial and logbinomial" begin
    import SphericalFunctions: sqrtbinomial, logbinomial
    using DoubleFloats: Double64

    # Against exact `BigInt` binomials.  `binomial(::Int, ::Int)` overflows above n ≈ 66,
    # which is exactly what this function is for, so the reference is always `big`.  The
    # value is the exponential of a logarithm, whose rounding error it magnifies by the size
    # of the logarithm, so the tolerance is proportional to that; the largest error measured
    # in any of these types is under 4 eps(T) times the logarithm.
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        for ℓ ∈ (1, 2, 3, 4, 5, 13, 64, 65, 66, 67, 1025)
            for s ∈ -2:2
                a = sqrtbinomial(2ℓ, ℓ - s, T)
                b = T(√big(binomial(big(2ℓ), big(ℓ - s))))
                @test a isa T
                @test a ≈ b rtol=8max(1, log(abs(b)))*eps(T)
            end
        end
    end

    # The edge cases of `logbinomial`, which `sqrtbinomial` reaches through its `k == 0`, `k
    # == n` and `k == 1` branches, and the k > n÷2 reflection.
    for n ∈ (0, 1, 2, 7, 100, 1025)
        @test logbinomial(n, 0) == 0
        @test logbinomial(n, n) == 0
        @test sqrtbinomial(n, 0) == 1
        @test sqrtbinomial(n, n) == 1
        if n ≥ 1
            @test logbinomial(n, 1) == log(n)
            @test logbinomial(n, n - 1) == log(n)
            @test sqrtbinomial(n, 1) ≈ √n rtol=2eps()
        end
        for k ∈ 0:min(n, 20)
            # binomial(n, k) == binomial(n, n-k), through the reflection branch
            @test logbinomial(n, k) == logbinomial(n, n - k)
            # The largest error measured here is under 4 eps
            @test logbinomial(n, k) ≈ log(big(binomial(big(n), big(k)))) rtol=8eps() atol=8eps()
        end
    end

    # As for `binomial`, a `k` outside 0:n gives zero
    for n ∈ (0, 2, 7), k ∈ (-2, -1, n + 1, n + 3)
        @test binomial(n, k) == 0
        @test sqrtbinomial(n, k) == 0
        @test sqrtbinomial(n, k, BigFloat) == 0
    end

    # It really does beat `binomial`: n = 2050 overflows an `Int`, and the value is far
    # outside `Float64`'s range, yet its square root is not.  The error measured is 67 eps.
    @test_throws OverflowError binomial(2050, 1025)
    @test sqrtbinomial(2050, 1025) ≈ Float64(√big(binomial(big(2050), big(1025)))) rtol=200eps()
    @test isfinite(sqrtbinomial(2050, 1025))
    @test !isfinite(sqrtbinomial(2100, 1050))

    # The computation type is honored, and defaults to Float64.
    @test sqrtbinomial(10, 5) isa Float64
    @test sqrtbinomial(10, 5, BigFloat) isa BigFloat
    @test sqrtbinomial(10, 5, Float32) isa Float32
    @test Float64(sqrtbinomial(10, 5, BigFloat)) ≈ sqrtbinomial(10, 5) rtol=4eps()
    @test logbinomial(10, 5) isa Float64
    @test logbinomial(Int32(10), Int32(5)) isa Float64
    @test logbinomial(10, 5, Double64) isa Double64
end


@testitem "Utilities: the element type and the number of rotors of rotor data" begin
    import SphericalFunctions: rotor_basetype, nrotors, check_rotor_type, DCalculator
    using Quaternionic: Rotor, Quaternion, QuatVec
    using DoubleFloats: Double64

    # The element type is the `float` of the component type of the data, for a single rotor,
    # angle or phase, and for a vector of any one of those
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        R = Rotor{T}(1, 0, 0, 0)
        for data ∈ (R, [R, R], T(0.3), T[0.3, 0.4], cis(T(0.3)), [cis(T(0.3))])
            @test rotor_basetype(data) === T
        end
    end
    @test rotor_basetype(3) === Float64
    @test rotor_basetype([1, 2]) === Float64
    @test rotor_basetype(1 + 0im) === Float64

    # Data whose component type is abstract does not say what type to work in, and is
    # refused rather than answered with a guess, which would be `float(Real) === Float64`
    for data ∈ (
        Complex{Real}(1, 2.0), Complex{Real}[1 + 0im, 2.0 + 0im], Complex{AbstractFloat}[1.0 + 0im],
        Rotor{Real}(1.0, 0, 0, 0), Rotor{Real}[Rotor(1.0)], Real[1.0, 2.0],
        Union{Float64, Float32}[1.0, 2.0f0],
    )
        @test_throws ArgumentError rotor_basetype(data)
        @test_throws "which does not say what floating-point type to work in" rotor_basetype(data)
    end
    @test_throws "has components of type Real" rotor_basetype(Rotor{Real}[Rotor(1.0)])
    @test_throws ArgumentError DCalculator(Rotor{Real}[Rotor(1.0)], 2)

    # Data that is not rotor data at all is refused with the forms that are accepted
    for data ∈ (Any[1.0, 2.0], Rotor[Rotor(1.0)], "β", [1.0 2.0; 3.0 4.0], nothing)
        @test_throws ArgumentError rotor_basetype(data)
        @test_throws "Cannot build a calculator from rotor data of type" rotor_basetype(data)
        @test_throws "the accepted forms are a Rotor" rotor_basetype(data)
    end
    # ... and a quaternion that is not a `Rotor` is told how to make one
    for data ∈ (Quaternion(1.0, 0, 0, 0), [Quaternion(1.0, 0, 0, 0)], QuatVec(0.0, 1, 0, 0))
        @test_throws ArgumentError rotor_basetype(data)
        @test_throws "Rotations are taken as `Rotor`s" rotor_basetype(data)
    end

    # The number of rotors is one for a single one, and the length of a non-empty vector
    @test nrotors(Rotor(1.0)) == nrotors(0.3) == nrotors(cis(0.3)) == 1
    @test nrotors([0.3]) == 1
    @test nrotors(fill(Rotor(1.0), 5)) == 5
    @test_throws ArgumentError nrotors(Float64[])
    @test_throws "needs at least one rotor, but got an empty Vector{Float64}" nrotors(Float64[])
    @test_throws ArgumentError nrotors("β")
    @test_throws "Cannot build a calculator from rotor data of type String" nrotors("β")

    # A calculator's element type is fixed by the data it was built from
    calc = DCalculator(Rotor(1.0), 2)
    @test check_rotor_type(calc, Rotor(0.0, 1.0, 0.0, 0.0)) === nothing
    @test check_rotor_type(calc, [0.3, 0.4]) === nothing
    @test_throws ArgumentError check_rotor_type(calc, Rotor{Float32}(1, 0, 0, 0))
    @test_throws "This calculator works in Float64, but the given data would give Float32" check_rotor_type(calc, Float32(0.3))

    # The element type is known from the type of the data alone
    for data ∈ (Rotor(1.0), [Rotor(1.0f0)], 0.3, Double64[0.3], cis(0.3), [cis(big(0.3))])
        @test Base.return_types(rotor_basetype, (typeof(data),)) == [Type{rotor_basetype(data)}]
    end
end
