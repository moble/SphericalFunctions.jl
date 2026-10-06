@testitem "Complex powers" setup=[Utilities] begin
    import .Utilities: array_equal

    # Each power zᵏ is compared with the exact `BigFloat` power, element by element.  The
    # recurrence accumulates error linearly in k; the worst error measured over these cases
    # is 0.37(k+1) eps(T).  In `Float16`, the bound 2(k+1) eps(T) exceeds 1 beyond k ≈ 500,
    # where it would accept any value on the unit circle, so m = 1000 is omitted there.
    for T in [Float64, Float32, Float16]
        nozpowers = Vector{Complex{T}}(undef, 0)
        fudge = one(T) + 2 * sqrt(eps(T))
        for k in 0:25
            z = cis(k*big(π)/10)
            for m in (T === Float16 ? [0, 1, 2, 3, 4, 100] : [0, 1, 2, 3, 4, 1_000])
                mine = complex_powers(Complex{T}(z), m)
                theirs = z.^collect(0:m)
                @test all(abs.(mine .- theirs) .≤ 2(1:m+1) .* eps(T))
                @test_throws DomainError complex_powers(Complex{T}(z*fudge), m)
                @test_throws DomainError complex_powers(Complex{T}(z/fudge), m)
                inplace = zeros(Complex{T}, size(mine))
                complex_powers!(inplace, Complex{T}(z))
                @test array_equal(mine, inplace)
            end
            @test length(complex_powers!(nozpowers, Complex{T}(z))) == 0
        end
    end

end


@testitem "Complex powers: accuracy of the fused modulus" begin
    using SphericalFunctions: complex_powers

    # `complex_powers!` computes `modulus` as `√(fma(z.re, z.re, z.im*z.im))` rather than
    # `√(abs2(z))`.  A one-ulp error there feeds the cancellation-sensitive `dc`, and the
    # recurrence amplifies it linearly in `m`; the spelling therefore matters far more than
    # it looks like it should.  Measured at ϕ = 0.3: with the fused form the error is
    # 5.5e-16 / 2.6e-15 / 7.7e-15 at m = 128 / 1024 / 4096, against 2.9e-15 / 2.2e-14 /
    # 9.1e-14 for `√(abs2(z))`.  The thresholds below sit between the two, so this test
    # discriminates, on every machine: `fma` rounds once whether or not the hardware has an
    # FMA instruction, where `muladd` would round twice on a machine that does not fuse it.
    setprecision(BigFloat, 512) do
        ϕ = 0.3
        for (m, tol) in ((128, 1.0e-15), (1024, 6.0e-15), (4096, 2.0e-14))
            computed = complex_powers(cis(ϕ), m)
            exact = [Complex{BigFloat}(cis(BigFloat(ϕ) * k)) for k in 0:m]
            @test maximum(abs.(Complex{BigFloat}.(computed) .- exact)) < tol
        end
    end
end


@testitem "Complex powers: the arguments of complex_powers" begin
    using SphericalFunctions: complex_powers

    z = cis(0.3)
    # The largest power may be of any integer type, and is refused if it is negative
    @test complex_powers(z, Int32(3)) == complex_powers(z, 3)
    @test complex_powers(z, big(3)) == complex_powers(z, 3)
    @test complex_powers(z, 0) == [one(z)]
    @test_throws ArgumentError complex_powers(z, -1)
    @test_throws "must be non-negative; got m=-2" complex_powers(z, -2)

    # A real phase, or one with integer or `Bool` components, is converted to a complex
    # floating-point number, so that the powers are computed and returned in that type
    @test complex_powers(1.0, 3) == ones(ComplexF64, 4)
    @test eltype(complex_powers(1.0, 3)) === ComplexF64
    @test eltype(complex_powers(1.0f0, 3)) === ComplexF32
    @test complex_powers(-1, 3) == ComplexF64[1, -1, 1, -1]
    @test eltype(complex_powers(1 + 0im, 3)) === ComplexF64
    @test complex_powers(im, 4) == ComplexF64[1, im, -1, -im, 1]
    @test eltype(complex_powers(im, 4)) === ComplexF64

    # A number that is not a phase is refused, because the recurrence would compute the
    # powers of a different number
    for w ∈ (2cis(0.3), 0.0im, 0.5)
        @test_throws DomainError complex_powers(w, 3)
    end
    @test_throws "complex amplitude approximately 1" complex_powers(2cis(0.3), 3)
end


@testitem "Complex powers: ComplexPowers runs the recurrence of complex_powers!" begin
    using SphericalFunctions: complex_powers!, ComplexPowers
    using DoubleFloats: Double64

    # The iterator and the in-place kernel share the start and the step of one recurrence,
    # so they agree to the last bit, for every phase, including those on the diagonals,
    # which are the edges of the sector -π/4 < arg z ≤ π/4 into which the recurrence rotates
    # `z`, and in every precision.
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        for θ ∈ (
            range(-T(π), T(π), length=41)..., T(π)/2, -T(π)/2, zero(T),
            T(π)/4, -T(π)/4, 3T(π)/4, -3T(π)/4,
        )
            z = cis(θ)
            @test isequal(first(ComplexPowers(z), 30), complex_powers!(zeros(Complex{T}, 30), z))
        end
        diagonals = (Complex{T}(σ₁, σ₂) / √T(2) for σ₁ ∈ (1, -1) for σ₂ ∈ (1, -1))
        for z ∈ (Complex{T}.((1, -1, im, -im))..., diagonals...)
            @test isequal(first(ComplexPowers(z), 30), complex_powers!(zeros(Complex{T}, 30), z))
        end
    end

    # It is an infinite iterator of the complex type of its argument, to which a real or
    # integer phase is converted
    let p = ComplexPowers(cis(0.1))
        @test Base.IteratorSize(typeof(p)) === Base.IsInfinite()
        @test eltype(p) === ComplexF64
        @test eltype(ComplexPowers(cis(0.1f0))) === ComplexF32
        @test first(ComplexPowers(1.0), 3) == ones(ComplexF64, 3)
        @test first(ComplexPowers(im), 5) == ComplexF64[1, im, -1, -im, 1]
    end

    # ... and, like `complex_powers`, it refuses a number that is not a phase
    for w ∈ (2cis(0.3), 0.0im, 0.5)
        @test_throws DomainError ComplexPowers(w)
    end
    @test_throws "complex amplitude approximately 1" ComplexPowers(0.0im)
end


@testitem "Complex powers: a quarter turn is the division by i" begin
    using SphericalFunctions: quarter_turn
    using DoubleFloats: Double64
    using Random: Xoshiro

    # The recurrence rotates its phase into the sector -π/4 < arg z ≤ π/4 by quarter turns,
    # which exchange the components where that is exactly the division by i, and divide
    # otherwise.  So the two agree to the last bit, signs of zeros included, for every pair
    # of these components: zeros of either sign, infinities, NaN, the smallest positive
    # number (a subnormal, except in `BigFloat`), the extremes of each type, the edges 1/2
    # and 2 of the range in which the components are exchanged and the numbers just outside
    # it, the components ±1/√2 of the phases on the diagonals, and numbers far from the unit
    # circle.  They agree for random phases too.
    rng = Xoshiro(20261001)
    for T ∈ (Float16, Float32, Float64, BigFloat, Double64)
        tiny = nextfloat(zero(T))
        special = T[
            0, -0.0, 1, -1, 0.5, -0.5, 2, -2, prevfloat(T(1)/2), nextfloat(T(2)),
            Inf, -Inf, NaN, tiny, -tiny, floatmin(T), -floatmin(T), floatmax(T),
            -floatmax(T), eps(T), -eps(T), 1/√T(2), -1/√T(2), T(1e-300), T(1e300),
        ]
        @test all(
            isequal(quarter_turn(Complex(x, y)), Complex(x, y) / 1im)
            for x ∈ special, y ∈ special
        )
        @test all(
            isequal(quarter_turn(z), z / 1im)
            for z ∈ (cis(2T(π) * rand(rng, T)) for _ ∈ 1:10_000)
        )
    end
end


@testitem "Complex powers: several phases at once are separate calls" begin
    using SphericalFunctions: complex_powers!, phase_powers!
    using DoubleFloats: Double64
    using Random: Xoshiro

    # The calculators fill the powers of the phases of two rotors in one call of
    # `phase_powers!`, and the powers of each phase are exactly those that `complex_powers!`
    # gives for that phase alone: in every precision, at every length, for one to four
    # phases at once, and for the phases on the axes and the diagonals, where the rotation
    # into the sector of the recurrence is decided.  Every phase of the list takes its turn
    # in each position of a call.
    rng = Xoshiro(20261001)
    for T ∈ (Float16, Float32, Float64, BigFloat, Double64)
        s = 1/√T(2)
        phases = Complex{T}[
            1, im, -1, -im, complex(1, -zero(T)), complex(-1, -zero(T)),
            complex(-zero(T), 1), complex(-zero(T), -1),
            complex(s, s), complex(-s, s), complex(s, -s), complex(-s, -s),
            (cis(2T(π) * rand(rng, T)) for _ ∈ 1:12)...,
        ]
        for n ∈ 1:4, M ∈ (0, 1, 2, 3, 4, 17, 65)
            @test all(1:length(phases)) do k
                z = ntuple(i -> phases[mod1(k + i - 1, length(phases))], n)
                together = phase_powers!(ntuple(_ -> fill(Complex{T}(NaN), M), n), z)
                alone = map(zᵢ -> complex_powers!(fill(Complex{T}(NaN), M), zᵢ), z)
                all(map(isequal, together, alone))
            end
        end
    end

    # A phase in a lower precision than its vector gives the powers that `complex_powers!`
    # gives for it, whatever the other phases of the call are
    for M ∈ (3, 17, 65)
        z = (cis(0.3f0), cis(2.9), cis(-1.3f0), ComplexF32(-0.0, 1))
        together = phase_powers!(ntuple(_ -> zeros(ComplexF64, M), 4), z)
        alone = map(zᵢ -> complex_powers!(zeros(ComplexF64, M), zᵢ), z)
        @test all(map(isequal, together, alone))
    end

    # The rows of a matrix, which is how the calculators pass their tables, are filled as
    # vectors of their own are
    let Z = fill(ComplexF64(NaN), 4, 17), z = (cis(0.3), cis(-2.0), 1.0im, cis(π/4))
        phase_powers!(ntuple(i -> view(Z, i, :), 4), z)
        alone = map(zᵢ -> complex_powers!(zeros(ComplexF64, 17), zᵢ), z)
        @test all(i -> isequal(Z[i, :], alone[i]), 1:4)
    end

    # Vectors of unequal lengths are refused before anything is written
    let a = fill(ComplexF64(NaN), 3), b = fill(ComplexF64(NaN), 4)
        @test_throws DimensionMismatch phase_powers!((a, b), (cis(0.1), cis(0.2)))
        @test_throws "lengths (3, 4)" phase_powers!((a, b), (cis(0.1), cis(0.2)))
        @test all(isnan, a) && all(isnan, b)
    end
end


@testitem "ComplexPowers" begin
    using SphericalFunctions: ComplexPowers

    # The error in zᵐ, measured here against the powers of the exact phase and so including
    # the rounding of z itself, grows linearly in m, and so does the bound.  Over this grid
    # the worst ratio of the error to (m+1) eps is 0.58 (at θ = 2.5, m = 10), so the bound
    # of twice that leaves a margin of more than 3.  The errors are accumulated and asserted
    # once per phase, rather than once per power.
    mₘₐₓ = 10_000
    for θ ∈ BigFloat(0):big(1//10):2big(π)
        z¹ = cis(θ)
        zᵐexact = one(z¹)
        worst = 0.0
        for (i, zᵐ) in enumerate(ComplexPowers(ComplexF64(z¹)))
            m = i-1
            worst = max(worst, Float64(abs(zᵐ - zᵐexact)) / (2(m+1) * eps(Float64)))
            zᵐexact *= z¹
            m == mₘₐₓ && break
        end
        @test worst < 1
    end
end


@testitem "Complex powers: a phase near an axis is as accurate on either side of it" begin
    using SphericalFunctions: complex_powers

    # The recurrence runs on `z` rotated by a power of i into the sector -π/4 < arg z ≤ π/4,
    # so a phase at a distance δ from an axis runs it on e^{±iδ}, near 1, from whichever
    # side it approaches the axis.  Measured against the exact powers of the rounded `z`,
    # over m ∈ 256:1024, the error in zᵐ at these points is at most 0.26 m eps in either
    # precision, and the bound is 0.4 m eps.  Were a phase just below an axis rotated to π/2
    # - δ instead, where Re z is near 0 and the recurrence is least accurate, the errors at
    # these points would reach 0.91 m eps, so the bound discriminates between the two.  The
    # smallest powers are left out, because their errors are a few roundings, which dominate
    # the ratio at small m whatever the rotation.
    setprecision(BigFloat, 128) do
        for T ∈ (Float64, Float32), k ∈ 0:3, σ ∈ (-1, 1)
            for δ ∈ (1e-4, 1e-3, 3e-3, 0.0175, 0.025, 0.035, 0.2)
                z = Complex{T}(cis(k * big(π) / 2 + σ * big(δ)))
                computed = complex_powers(z, 1024)
                exact = Complex{BigFloat}(z)
                zᵐ = one(exact)
                worst = 0.0
                for m ∈ 1:1024
                    zᵐ *= exact
                    if m ≥ 256
                        worst = max(worst, Float64(abs(computed[m+1] - zᵐ)) / (m * eps(T)))
                    end
                end
                @test worst < 0.4
            end
        end
    end
end


@testitem "Complex powers: a phase near 1 is as accurate on either side of it" begin
    using SphericalFunctions: complex_powers

    # A phase e^{iθ} with small |θ| is not rotated, whatever the sign of θ, and the recurrence
    # is then at its most accurate, since t ≈ -θ² is tiny.  At these phases, measured at m =
    # 4096 against the exact power of the rounded `z`, the error is at most 0.035 m eps in
    # Float64 and 6e-4 m eps in Float32, the same for either sign; the bounds are 0.1 m eps
    # and 0.01 m eps.  Were a small negative θ rotated to π/2 - |θ| instead, far from 1, the
    # errors would be 0.14 to 0.5 m eps for most of these phases.
    setprecision(BigFloat, 128) do
        M = 4096
        for (T, θs, bound) ∈ (
            (Float64, (1e-3, 1e-5, 1e-7, 1e-9, 1e-11, 1e-13, 1e-15), 0.1),
            (Float32, (1e-5, 1e-7), 0.01),
        )
            for θ ∈ θs, σ ∈ (-1, 1)
                z = Complex{T}(cis(σ * big(θ)))
                zᴹ = complex_powers(z, M)[M+1]
                @test Float64(abs(zᴹ - Complex{BigFloat}(z)^M)) / (M * eps(T)) < bound
            end
        end
    end
end


@testitem "Complex powers: the powers of -z and i z are (-1)ᵐ and iᵐ times those of z" begin
    using SphericalFunctions: complex_powers, ComplexPowers
    using DoubleFloats: Double64

    # The sector -π/4 < arg z ≤ π/4 into which the recurrence rotates `z` is half open, so
    # z, i z, -z and -i z are all rotated to the same phase in it, and differ only in the
    # power of i that is factored out and restored exactly.  That includes the phases on the
    # diagonals, which are the edges of the sector, and those on the axes.  (The comparison
    # is `==`, because the signs of zeros may differ.)  This is what makes 𝔇(-R) = 𝔇(R)
    # exact for integer ℓ, even when a phase of R lies on a diagonal.
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        diagonals = (Complex{T}(σ₁, σ₂) / √T(2) for σ₁ ∈ (1, -1) for σ₂ ∈ (1, -1))
        on_axes = Complex{T}.((1, -1, im, -im))
        for z ∈ (cis.(range(-T(π), T(π), length=101))..., diagonals..., on_axes...)
            zᵐ = complex_powers(z, 200)
            @test complex_powers(-z, 200) == [(-1)^m * zᵐ[m+1] for m ∈ 0:200]
            @test complex_powers(im * z, 200) == [(1im)^mod(m, 4) * zᵐ[m+1] for m ∈ 0:200]
            @test first(ComplexPowers(-z), 50) == [(-1)^m * zᵐ[m+1] for m ∈ 0:49]
        end
    end
end


@testitem "Complex powers: the powers of conj(z) are the conjugates of the powers of z" begin
    using SphericalFunctions: complex_powers, ComplexPowers

    # Away from the diagonals, the rotation that brings conj(z) into the sector -π/4 < arg z
    # ≤ π/4 is the conjugate of the one that brings `z`, so every step for conj(z) is
    # exactly the conjugate of the same step for `z`.  None of the phases of the grid lies
    # on a diagonal.  (The comparison is `==`, because the signs of zeros may differ.)  On a
    # diagonal, where the sector is half open, `z` and conj(z) are rotated to the same
    # phase, (1 + i)/√2, and their powers are conjugates only to within the error of the
    # recurrence: at most 0.71 m eps in these three precisions, and the bound is 2 m eps.
    for T ∈ (Float16, Float32, Float64)
        for θ ∈ range(-T(π), T(π), length=101)
            z = cis(θ)
            @test complex_powers(conj(z), 200) == conj.(complex_powers(z, 200))
            @test first(ComplexPowers(conj(z)), 50) == conj.(first(ComplexPowers(z), 50))
        end
        for σ₁ ∈ (1, -1), σ₂ ∈ (1, -1)
            z = Complex{T}(σ₁, σ₂) / √T(2)
            Δ = complex_powers(conj(z), 200) - conj.(complex_powers(z, 200))
            @test all(abs(Δ[m+1]) ≤ 2m * eps(T) for m ∈ 0:200)
        end
    end
end
