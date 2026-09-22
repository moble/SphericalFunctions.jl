# Tests for the v3 spin-weighted spherical-harmonic transforms: the `SSHT` constructor and the
# three concrete types `SSHTRS`, `SSHTMinimal` and `SSHTMatrix`, their `pixels`/`rotors`,
# synthesis (`𝒯 * f̃`, `mul!`), analysis (`𝒯 \ f`, `ldiv!`), in-place semantics, `ModeWeights`
# input and output, the `SSHTRS` quadrature options, and use from several tasks.
#
# The oracle is the closed-form ₛYₗₘ of the `Utilities` snippet — the explicit sum over
# factorials from the conventions pages, which owes nothing to this package.  The items here
# use `sYlm_pixels`, also from that snippet: the same formula with the pixel-independent
# factorials hoisted out of the loop (≈ 40 times faster, and checked against `sYlm` itself in
# "SSHT synthesis").  They compare synthesis against it pixel by pixel, analysis against the
# mode it came from, and round trips against their input.  The pixelizations are compared against the formulas in their own
# docstrings.  Cross-checks that stay inside the package — `sYlm_matrix`, the evaluation
# `w(R)` of a `ModeWeights`, and the agreement of the three algorithms on a shared grid — are
# labelled as such where they appear.
#
# Tolerances are set from measurement; the numbers quoted in the comments are the largest
# errors seen over every method, T, ℓₘₐₓ and s that the item covers, in units of eps(T).
#
# This file mirrors the intent of the v2 `test/deprecated/ssht.jl` (deleted in 3.0) against
# the v3 API.  The item names here all begin with "SSHT", which is what kept them from
# colliding with the v2 items while both existed.


@testitem "SSHT construction" begin
    import SphericalFunctions: SSHT, SSHTRS, SSHTMinimal, SSHTMatrix, pixels, rotors, spin, Ysize
    import SphericalFunctions: nmodes, npixels  # unexported
    import SphericalFunctions: golden_ratio_spiral_rotors
    using DoubleFloats: Double64
    using Logging: NullLogger, with_logger
    using Quaternionic: Rotor

    types = (("RS", SSHTRS), ("Minimal", SSHTMinimal), ("Matrix", SSHTMatrix))

    # Every method and type, including the edge cases s = ℓₘₐₓ and ℓₘₐₓ = 1
    for T in (Float64, Double64, Float32), (s, ℓₘₐₓ) in ((-2, 6), (0, 3), (1, 1), (3, 3), (-3, 4))
        for (method, Type) in types
            𝒯 = SSHT(s, ℓₘₐₓ; method, T)
            @test 𝒯 isa Type{T}
            @test 𝒯 isa SSHT{T}
            @test typeof(Type(s, ℓₘₐₓ; T)) === typeof(𝒯)  # the concrete constructor agrees
            @test spin(𝒯) == s
            @test SphericalFunctions.ℓₘₐₓ(𝒯) == ℓₘₐₓ
            @test SphericalFunctions.ℓₘᵢₙ(𝒯) == abs(s)
            @test nmodes(𝒯) == Ysize(abs(s), ℓₘₐₓ) == (ℓₘₐₓ + 1)^2 - s^2
            @test npixels(𝒯) == length(pixels(𝒯)) == length(rotors(𝒯))
            @test npixels(𝒯) ≥ nmodes(𝒯)
            @test eltype(rotors(𝒯)) === Rotor{T}

            # `show` (both forms), `repr` and `string` agree, and name the type, T, s and ℓₘₐₓ
            str = repr(𝒯)
            @test str == string(𝒯)
            @test str == sprint(show, 𝒯)
            @test str == sprint(show, MIME("text/plain"), 𝒯)
            @test occursin("SSHT$method", str)
            @test occursin(string(T), str)
            @test occursin("s=$s", str)
            @test occursin("ℓₘₐₓ=$ℓₘₐₓ", str)
        end
        # Default pixelizations: RS has (2ℓₘₐₓ+1) rings of (2ℓₘₐₓ+1) points; the others are
        # optimal, with exactly as many points as modes
        @test npixels(SSHT(s, ℓₘₐₓ; method="RS", T)) == (2ℓₘₐₓ + 1)^2
        @test npixels(SSHT(s, ℓₘₐₓ; method="Minimal", T)) == Ysize(abs(s), ℓₘₐₓ)
        @test npixels(SSHT(s, ℓₘₐₓ; method="Matrix", T)) == Ysize(abs(s), ℓₘₐₓ)
    end

    # Defaults: method "RS", T Float64
    @test SSHT(1, 4) isa SSHTRS{Float64}
    @test SSHT(1, 4; method="Minimal") isa SSHTMinimal{Float64}
    @test SSHT(1, 4; method="Matrix") isa SSHTMatrix{Float64}

    # The `inplace` option is part of the type of the in-place-capable methods (default true)
    @test SSHT(1, 4; method="Minimal") isa SSHTMinimal{Float64, true}
    @test SSHT(1, 4; method="Minimal", inplace=false) isa SSHTMinimal{Float64, false}
    @test SSHT(1, 4; method="Matrix") isa SSHTMatrix{Float64, true}
    @test SSHT(1, 4; method="Matrix", inplace=false) isa SSHTMatrix{Float64, false}
    @test SSHTMinimal(1, 4; T=Float32, inplace=false) isa SSHTMinimal{Float32, false}
    @test SSHTMatrix(1, 4; T=Float32, inplace=false) isa SSHTMatrix{Float32, false}

    # "Direct" is a deprecated alias of "Matrix"
    @test_deprecated r"renamed \"Matrix\"" SSHT(0, 3; method="Direct")
    if Base.JLOptions().depwarn != 2  # with --depwarn=error the call throws instead
        𝒯 = with_logger(NullLogger()) do
            SSHT(0, 3; method="Direct")
        end
        @test 𝒯 isa SSHTMatrix{Float64, true}
        @test rotors(𝒯) == rotors(SSHT(0, 3; method="Matrix"))
    end

    # A nonsense method is an error naming the offending string; the names are case-sensitive
    @test_throws ErrorException SSHT(-2, 8; method="NonsenseGarbage")
    @test_throws "NonsenseGarbage" SSHT(-2, 8; method="NonsenseGarbage")
    @test_throws "Unrecognized" SSHT(-2, 8; method="NonsenseGarbage")
    @test_throws ErrorException SSHT(-2, 8; method="rs")
    @test_throws ErrorException SSHT(-2, 8; method="matrix")

    # |s| > ℓₘₐₓ is an error for every method (there are no such modes) — including inside
    # the default `Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, T)` of the "Matrix" method, which
    # is evaluated before the constructor body runs its own check.
    for method in ("RS", "Minimal", "Matrix"), (s, ℓₘₐₓ) in ((3, 2), (-3, 2), (1, 0))
        @test_throws "exceeds ℓₘₐₓ" SSHT(s, ℓₘₐₓ; method)
    end
    # The Matrix check itself is reached when the rotors are supplied explicitly
    @test_throws "exceeds ℓₘₐₓ" SSHTMatrix(3, 2; Rθϕ=golden_ratio_spiral_rotors(0, 3))
    @test_throws "exceeds ℓₘₐₓ" SSHTMinimal(3, 2; θ=Float64[])
    @test_throws "exceeds ℓₘₐₓ" SSHTRS(3, 2)

    # "Matrix" warns when the dense matrix gets large (Ysize² > 65⁴, i.e. ℓₘₐₓ ≥ 65 for s = 0)
    @test_logs (:warn, r"\"Matrix\" method for s-SHT is only") match_mode=:any SSHTMatrix(0, 65)
    # ... and is quiet otherwise, as are the other methods
    @test_logs SSHTMatrix(0, 8)
    @test_logs SSHT(0, 8; method="RS")
    @test_logs SSHT(0, 8; method="Minimal")
end


@testitem "SSHT pixels and rotors" begin
    import SphericalFunctions: SSHT, pixels, rotors
    import SphericalFunctions: npixels  # unexported
    import SphericalFunctions: sorted_rings, sorted_ring_pixels, sorted_ring_rotors
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors, fejer1_rings
    using DoubleFloats: Double64
    using Quaternionic: Rotor, from_spherical_coordinates, to_spherical_coordinates
    using StaticArrays: SVector
    using Random

    rng = Random.Xoshiro(1729)

    for T in (Float64, Double64, Float32), ℓₘₐₓ in (3, 4, 5, 8, 13), s in -2:2
        # Minimal: ring j ∈ |s|:ℓₘₐₓ at colatitude `sorted_rings(s, ℓₘₐₓ, T)[j]` has 2j+1
        # equally spaced points starting at ϕ = 0; this is `sorted_ring_pixels` (up to the
        # rounding of ϕ, which the two compute differently)
        𝒯 = SSHT(s, ℓₘₐₓ; method="Minimal", T)
        p = pixels(𝒯)
        @test eltype(p) === SVector{2, T}
        @test eltype(rotors(𝒯)) === Rotor{T}
        @test p ≈ sorted_ring_pixels(s, ℓₘₐₓ, T)
        @test rotors(𝒯) ≈ sorted_ring_rotors(s, ℓₘₐₓ, T)
        @test rotors(𝒯) == from_spherical_coordinates.(p)
        θs = sorted_rings(s, ℓₘₐₓ, T)
        expected = [
            SVector(θⱼ, k * 2T(π) / (2j + 1))
            for (j, θⱼ) in zip(abs(s):ℓₘₐₓ, θs) for k in 0:2j
        ]
        @test p == expected
        @test length(p) == npixels(𝒯) == (ℓₘₐₓ + 1)^2 - s^2

        # ... and those rings are what `sorted_rings` documents: the interior points of an
        # equally spaced grid on [0, π] (so never a pole), ordered so that each successive
        # ring — which has one more pair of points than the last — lies at least as close
        # to the equator as its predecessor.  Measured: the set matches the grid exactly, and
        # the ordering is violated by at most 1 eps(T), among rings that are exactly equally
        # far from the equator and whose order is therefore arbitrary.
        @test sort(θs) == collect(LinRange{T}(0, T(π), ℓₘₐₓ - abs(s) + 3))[begin+1:end-1]
        @test all(0 < θⱼ < T(π) for θⱼ in θs)
        let d = abs.(θs .- T(π) / 2)
            @test all(d[j+1] ≤ d[j] + 4eps(T) for j in 1:length(d)-1)
        end

        # RS: ring-major, ϕₖ = 2πk/Nϕ fastest, on the Fejér-1 rings by default
        𝒯 = SSHT(s, ℓₘₐₓ; method="RS", T)
        p = pixels(𝒯)
        @test eltype(p) === SVector{2, T}
        @test eltype(rotors(𝒯)) === Rotor{T}
        @test length(p) == sum(𝒯.Nϕ) == npixels(𝒯)
        @test 𝒯.Nϕ == fill(2ℓₘₐₓ + 1, 2ℓₘₐₓ + 1)
        θs = fejer1_rings(2ℓₘₐₓ + 1, T)
        @test 𝒯.θ == θs
        expected = [SVector(θⱼ, k * 2T(π) / (2ℓₘₐₓ + 1)) for θⱼ in θs for k in 0:2ℓₘₐₓ]
        @test p == expected
        @test rotors(𝒯) == from_spherical_coordinates.(p)
        # Within a ring θ is constant and ϕ increases; across rings θ increases
        @test all(p[i][1] == p[i + 1][1] && p[i][2] < p[i + 1][2] for i in 1:2ℓₘₐₓ)
        @test issorted([p[1 + (2ℓₘₐₓ + 1) * r][1] for r in 0:2ℓₘₐₓ])

        # Matrix: the golden-ratio spiral by default, or exactly the rotors given
        R = golden_ratio_spiral_rotors(s, ℓₘₐₓ, T)
        𝒯 = SSHT(s, ℓₘₐₓ; method="Matrix", T)
        @test rotors(𝒯) == R
        @test eltype(rotors(𝒯)) === Rotor{T}
        @test pixels(𝒯) == to_spherical_coordinates.(R)
        Rshuffled = R[randperm(rng, length(R))]
        𝒯 = SSHT(s, ℓₘₐₓ; method="Matrix", T, Rθϕ=Rshuffled)
        @test rotors(𝒯) == Rshuffled
        @test pixels(𝒯) == to_spherical_coordinates.(Rshuffled)

        # ... and that default spiral is the one its docstring describes: N = (ℓₘₐₓ+1)² - s²
        # points, successive azimuths separated by exactly Δϕ = 2π(2-φ), and cos θ at the
        # midpoints of N equal subintervals of [-1, 1] — "uniformly distributed in cos θ",
        # with no point on either pole.  Measured: ϕ is reproduced exactly (0 eps(T)) and
        # cos θ to within 2 eps(T), for every (T, ℓₘₐₓ, s) here.
        gp = golden_ratio_spiral_pixels(s, ℓₘₐₓ, T)
        @test R == from_spherical_coordinates.(gp)
        let N = (ℓₘₐₓ + 1)^2 - s^2, Δϕ = 2T(π) * (2 - T(MathConstants.φ))
            @test length(gp) == N
            @test maximum(abs, [gp[i+1][2] - i * Δϕ for i in 0:N-1]) ≤ 4eps(T) * N
            @test maximum(abs, [cos(gp[i+1][1]) - (1 - T(2i + 1) / N) for i in 0:N-1]) ≤ 4eps(T)
            @test all(0 < θϕ[1] < T(π) for θϕ in gp)
        end
    end
end


# The five pixelizations that depend on the spin weight — the golden-ratio spiral and the
# sorted rings, as pixels and as rotors — accept half-odd-integer `s` and `ℓₘₐₓ`, spelled as
# `Rational`s with denominator 2 or as `HalfOddInteger`s.  The items below check the same
# geometric properties the integer item above checks, now with half-odd indices, and then
# that the integer path is exactly what it was before the indices were widened.

@testitem "SSHT pixelizations: half-integer pixel counts and spiral geometry" begin
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
    import SphericalFunctions: sorted_ring_pixels, sorted_ring_rotors
    import SphericalFunctions: Ysize
    using DoubleFloats: Double64
    using Quaternionic: Rotor
    using StaticArrays: SVector

    for T in (Float64, Double64, Float32), s in -7//2:7//2, ℓₘₐₓ in abs(s):abs(s)+6
        N = Ysize(abs(s), ℓₘₐₓ)
        @test N == (ℓₘₐₓ + 1)^2 - s^2  # the closed form holds for half-odd indices too
        for f in (golden_ratio_spiral_pixels, sorted_ring_pixels)
            p = f(s, ℓₘₐₓ, T)
            @test length(p) == N
            @test eltype(p) === SVector{2, T}
        end
        for f in (golden_ratio_spiral_rotors, sorted_ring_rotors)
            R = f(s, ℓₘₐₓ, T)
            @test length(R) == N
            @test eltype(R) === Rotor{T}
        end

        # The spiral is what its docstring describes, with N = Ysize(|s|, ℓₘₐₓ) points:
        # successive azimuths separated by exactly Δϕ = 2π(2-φ), and cos θ at the midpoints
        # of N equal subintervals of [-1, 1], with no point on either pole.  Measured: ϕ is
        # reproduced exactly and cos θ to within 1.5 eps(T), over every (T, s, ℓₘₐₓ) here.
        gp = golden_ratio_spiral_pixels(s, ℓₘₐₓ, T)
        let Δϕ = 2T(π) * (2 - T(MathConstants.φ))
            @test maximum(abs, [gp[i+1][2] - i * Δϕ for i in 0:N-1]) ≤ 4eps(T) * N
            @test maximum(abs, [cos(gp[i+1][1]) - (1 - T(2i + 1) / N) for i in 0:N-1]) ≤ 4eps(T)
            @test all(0 < θϕ[1] < T(π) for θϕ in gp)
        end
    end

    # A spin weight larger than ℓₘₐₓ is refused before any pixel is placed, by both
    # families and with one message.
    @test_throws "exceeds ℓₘₐₓ" golden_ratio_spiral_pixels(3//2, 1//2)
    @test_throws "exceeds ℓₘₐₓ" golden_ratio_spiral_rotors(-5//2, 3//2, Float32)
    @test_throws "exceeds ℓₘₐₓ" sorted_ring_pixels(3//2, 1//2)
    @test_throws "exceeds ℓₘₐₓ" sorted_ring_rotors(-5//2, 3//2, Float32)
    @test_throws "|s|=5//2 exceeds ℓₘₐₓ=3//2; there are no such modes." golden_ratio_spiral_pixels(-5//2, 3//2)
    @test_throws "|s|=5//2 exceeds ℓₘₐₓ=3//2; there are no such modes." sorted_ring_pixels(-5//2, 3//2)
end

@testitem "SSHT pixelizations: half-integer sorted rings" begin
    import SphericalFunctions: sorted_rings
    using DoubleFloats: Double64

    for T in (Float64, Double64, Float32), s in -7//2:7//2, ℓₘₐₓ in abs(s):abs(s)+9
        θs = sorted_rings(s, ℓₘₐₓ, T)
        n = Int(ℓₘₐₓ - abs(s) + 1)  # one ring per j ∈ |s|:ℓₘₐₓ
        @test length(θs) == n
        @test eltype(θs) === T
        @test all(0 < θⱼ < T(π) for θⱼ in θs)
        # The rings are the interior points of an equally spaced grid on [0, π], ordered so
        # that each successive ring lies at least as close to the equator as its predecessor.
        # Measured: the set matches the grid exactly, and the ordering is violated by at most
        # 1.5 eps(T), among rings that are equally far from the equator up to rounding and
        # whose order is settled by the ulp shift `ceil(s)`.  That shift is one ulp for
        # s = 1/2, which the rounding of the grid can match, so which side of the equator the
        # second ring falls on is not asserted here; the integer item does not assert it
        # either.
        @test sort(θs) == collect(LinRange{T}(0, T(π), n + 2))[begin+1:end-1]
        let d = abs.(θs .- T(π) / 2)
            @test all(d[j+1] ≤ d[j] + 4eps(T) for j in 1:length(d)-1)
        end
    end

    # A spin weight larger than ℓₘₐₓ describes no modes and is refused, with the message the
    # golden-ratio spiral uses; the smallest admissible ℓₘₐₓ is |s|, which gives one ring.
    @test_throws "exceeds ℓₘₐₓ" sorted_rings(3//2, 1//2)
    @test_throws "|s|=3//2 exceeds ℓₘₐₓ=1//2; there are no such modes." sorted_rings(3//2, 1//2)
    @test_throws "exceeds ℓₘₐₓ" sorted_rings(-7//2, 5//2, Float32)
    @test sorted_rings(3//2, 3//2) == [Float64(π) / 2]
end

@testitem "SSHT pixelizations: half-integer ring pixels and rotors" begin
    import SphericalFunctions: sorted_rings, sorted_ring_pixels, sorted_ring_rotors
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
    using DoubleFloats: Double64
    using Quaternionic: from_spherical_coordinates
    using StaticArrays: SVector

    for T in (Float64, Double64, Float32), s in -5//2:5//2, ℓₘₐₓ in abs(s):abs(s)+5
        # Ring j ∈ |s|:ℓₘₐₓ at colatitude `sorted_rings(s, ℓₘₐₓ, T)[j]` has 2j+1 equally
        # spaced points starting at ϕ = 0 — an even number of points for half-odd j.  The
        # colatitudes are reproduced exactly and every ring starts at exactly ϕ = 0; the
        # remaining azimuths come from a `LinRange`, whose rounding differs from that of
        # k·2π/(2j+1) by up to 4 eps(T) over every (T, s, ℓₘₐₓ) here (the exact equality the
        # integer item asserts holds only for an odd number of points).
        θs = sorted_rings(s, ℓₘₐₓ, T)
        p = sorted_ring_pixels(s, ℓₘₐₓ, T)
        expected = [
            SVector(θⱼ, k * 2T(π) / n)
            for (j, θⱼ) in zip(abs(s):ℓₘₐₓ, θs) for n in (Int(2j + 1),) for k in 0:n-1
        ]
        @test length(p) == length(expected)
        @test all(θϕ[1] == e[1] for (θϕ, e) in zip(p, expected))
        @test all(e[2] != 0 || θϕ[2] == 0 for (θϕ, e) in zip(p, expected))
        @test maximum(abs(θϕ[2] - e[2]) for (θϕ, e) in zip(p, expected)) ≤ 8eps(T)
        @test [count(θϕ -> θϕ[1] == θⱼ, p) for θⱼ in θs] == [Int(2j + 1) for j in abs(s):ℓₘₐₓ]

        # The rotors are the lifts of the pixels, for both families.
        @test sorted_ring_rotors(s, ℓₘₐₓ, T) == from_spherical_coordinates.(p)
        @test golden_ratio_spiral_rotors(s, ℓₘₐₓ, T) ==
            from_spherical_coordinates.(golden_ratio_spiral_pixels(s, ℓₘₐₓ, T))
    end
end

@testitem "SSHT pixelizations: half-integer spellings and mixed kinds" begin
    import SphericalFunctions: sorted_rings, sorted_ring_pixels, sorted_ring_rotors
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
    using SphericalFunctions: HalfOddInteger
    using DoubleFloats: Double64

    pixelizations = (
        golden_ratio_spiral_pixels, golden_ratio_spiral_rotors,
        sorted_rings, sorted_ring_pixels, sorted_ring_rotors
    )

    # The `Rational` spelling and the `HalfOddInteger` spelling are the same call.
    for T in (Float64, Double64, Float32), f in pixelizations, (s, ℓₘₐₓ) in ((1//2, 7//2), (-3//2, 9//2), (5//2, 5//2))
        r = f(s, ℓₘₐₓ, T)
        @test r == f(HalfOddInteger(s), HalfOddInteger(ℓₘₐₓ), T)
        @test r == f(HalfOddInteger(s), ℓₘₐₓ, T)
        @test r == f(s, HalfOddInteger(ℓₘₐₓ), T)
        @test typeof(r) === typeof(f(HalfOddInteger(s), HalfOddInteger(ℓₘₐₓ), T))
    end
    @test sorted_rings(1//2, 7//2) == sorted_rings(1//2, 7//2, Float64)
    @test golden_ratio_spiral_pixels(1//2, 7//2) == golden_ratio_spiral_pixels(1//2, 7//2, Float64)

    # Mixing the two kinds of index is refused with a message naming both spellings, for
    # every one of the five functions and in either order.
    for f in pixelizations
        @test_throws "all be integers, like 3, or all be half-odd-integers, like 7//2" f(1//2, 3)
        @test_throws "all be integers, like 3, or all be half-odd-integers, like 7//2" f(0, 7//2)
        @test_throws "all be integers, like 3, or all be half-odd-integers, like 7//2" f(HalfOddInteger(1//2), 3, Float32)
        @test_throws "all be integers, like 3, or all be half-odd-integers, like 7//2" f(2, HalfOddInteger(7//2), Float32)
        # A `Rational` with any other denominator is not an index at all.
        @test_throws "must have denominator 2" f(1//3, 7//2)
        @test_throws "must have denominator 2" f(2//1, 6//1)
    end
end

@testitem "SSHT pixelizations: half-integer support leaves the integer path unchanged" begin
    import SphericalFunctions: sorted_rings, sorted_ring_pixels, sorted_ring_rotors
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
    using StaticArrays: SVector

    # These literals were recorded from the functions as they stood before half-odd indices
    # were admitted, in a fresh session, and the comparisons are exact.
    @test sorted_rings(2, 6) == [
        2.6179938779914944, 0.5235987755982988, 2.0943951023931953, 1.0471975511965976,
        1.5707963267948966
    ]
    @test sorted_rings(-1, 3) == [0.7853981633974483, 2.356194490192345, 1.5707963267948966]
    @test sorted_rings(0, 3, Float32) == Float32[0.62831855, 2.5132742, 1.2566371, 1.8849556]
    @test eltype(sorted_rings(2, 6)) === Float64

    g = golden_ratio_spiral_pixels(2, 6)
    @test length(g) == 45
    @test g[1] == SVector(0.21121088035917845, 0.0)
    @test g[23] == SVector(1.5707963267948966, 52.799191054030366)
    @test g[end] == SVector(2.9303817732306148, 105.59838210806073)

    # ... and every pixel agrees exactly with the former inline formula, N = (ℓₘₐₓ+1)² - s²,
    # evaluated exactly as the former body evaluated it.
    for T in (Float64, Float32), (s, ℓₘₐₓ) in ((2, 6), (-1, 3), (0, 5), (3, 3))
        expected = let π = T(π), φ = T(MathConstants.φ)
            N = (ℓₘₐₓ+1)^2 - s^2
            Δϕ = 2π * (2 - φ)
            ϕ = (0:N-1) * Δϕ
            cosθ = LinRange{T}(1, -1, N+1)[begin:end-1] .- 1/T(N)
            [SVector(acos(cosθ), ϕ) for (cosθ, ϕ) in zip(cosθ, ϕ)]
        end
        @test golden_ratio_spiral_pixels(s, ℓₘₐₓ, T) == expected
        # ... as do the rings, with the spin weight itself as the count of ulps.
        expected_rings = let πo2 = prevfloat(T(π)/2, s)
            sort(
                collect(LinRange{T}(0, π, 2+ℓₘₐₓ-abs(s)+1))[begin+1:end-1],
                lt=(x,y)->(abs(x-πo2)<abs(y-πo2)),
                rev=true
            )
        end
        @test sorted_rings(s, ℓₘₐₓ, T) == expected_rings
    end

    # For |s| > ℓₘₐₓ, where the former bodies threw an obscure `LinRange` error or returned an
    # empty vector, both families now refuse the call with one clear message.  No input with
    # |s| ≤ ℓₘₐₓ is affected; ℓₘₐₓ = |s| gives one ring and 2|s|+1 pixels.
    for f in (golden_ratio_spiral_pixels, golden_ratio_spiral_rotors, sorted_rings, sorted_ring_pixels, sorted_ring_rotors)
        @test_throws "|s|=3 exceeds ℓₘₐₓ=2; there are no such modes." f(3, 2)
        @test_throws "|s|=4 exceeds ℓₘₐₓ=2; there are no such modes." f(-4, 2, Float32)
        @test length(f(2, 2)) == (f === sorted_rings ? 1 : 5)
    end

    # Integer indices of a narrower type are promoted, and give the same result.
    @test sorted_rings(Int8(2), 6) == sorted_rings(2, 6)
    @test golden_ratio_spiral_pixels(Int8(2), 6) == g
    @test sorted_ring_pixels(Int8(-1), Int16(3)) == sorted_ring_pixels(-1, 3)
    @test sorted_ring_rotors(2, 6) == sorted_ring_rotors(Int8(2), 6)
    @test golden_ratio_spiral_rotors(2, 6) == golden_ratio_spiral_rotors(Int8(2), 6)
end

@testitem "SSHT synthesis" setup=[Utilities] begin
    import SphericalFunctions: SSHT, pixels, rotors, Ysize, Yindex, sYlm_matrix, ModeWeights
    import SphericalFunctions: npixels  # unexported
    using DoubleFloats: Double64
    using LinearAlgebra: norm
    using StaticArrays: SVector
    using Random

    # `sYlm_pixels` — the `Utilities` closed form with the pixel-independent factorials
    # hoisted out of the loop, about 40 times faster — really is that closed form.  This is
    # the one place the shared helper is checked; the other items that use it say so.  The
    # two are compared over the whole ℓ range the file uses (not just the low ℓ where
    # cancellation in the alternating sum is mild), at four hand-picked points *and* on a
    # real pixelization, since they accumulate that sum in different precisions (BigFloat
    # in `sYlm`, `T` in `sYlm_pixels`).  Measured: ≤ 1.05 eps(T) at the four points and
    # ≤ 8.3 eps(T) on the pixelization, so 100 eps(T) leaves a factor of ≳ 12.
    for T in (Float64, Double64, Float32)
        q = [SVector{2, T}(θ, ϕ) for (θ, ϕ) in ((0.3, 1.1), (2.7, 5.2), (0.0, 0.0), (3.0, 2.0))]
        for s in -2:2
            qs = vcat(q, pixels(SSHT(s, 6; method="RS", T)))
            for ℓ in abs(s):6, m in -ℓ:ℓ
                @test maximum(
                    abs, sYlm_pixels(s, ℓ, m, qs) .- [sYlm(s, ℓ, m, θϕ[1], θϕ[2]) for θϕ in qs]
                ) ≤ 100eps(T)
            end
        end
    end

    rng = Random.Xoshiro(3141)

    for (method, T) in Iterators.product(("RS", "Minimal", "Matrix"), (Float64, Double64, Float32))
        kw = method == "RS" ? (;) : (; inplace=false)  # RS has no `inplace` option
        for ℓₘₐₓ in 3:6
            # Measured against the closed form over every method, T, s ∈ -2:2 and single
            # mode, the synthesis error is at most 26 eps(T) — smallest at ℓₘₐₓ = 3 (2 eps)
            # and growing roughly linearly with ℓₘₐₓ — so this leaves a factor of ≳ 11.
            ϵ = 50ℓₘₐₓ * eps(T)
            for s in -2:2
                𝒯 = SSHT(s, ℓₘₐₓ; method, T, kw...)
                n = Ysize(abs(s), ℓₘₐₓ)
                p = pixels(𝒯)
                # Column i of `Y` is the closed-form ₛYₗₘ of mode i sampled on the pixels
                Y = reduce(hcat, [sYlm_pixels(s, ℓ, m, p) for ℓ in abs(s):ℓₘₐₓ for m in -ℓ:ℓ])

                # Each single mode synthesizes to the closed-form harmonic on the pixels
                for ℓ in abs(s):ℓₘₐₓ, m in -ℓ:ℓ
                    f̃ = zeros(Complex{T}, n)
                    f̃[Yindex(ℓ, m, abs(s))] = one(T)
                    f = 𝒯 * f̃
                    @test f isa Vector{Complex{T}}
                    @test length(f) == npixels(𝒯)
                    @test f ≈ Y[:, Yindex(ℓ, m, abs(s))] atol=ϵ rtol=ϵ
                    # A ModeWeights with ℓₘᵢₙ = |s| is accepted and gives the same values
                    @test 𝒯 * ModeWeights(f̃, s) == f
                end

                # Random mode weights: linearity in the same closed-form harmonics.  Whole
                # blocks are compared at once, so a broken transform costs one failure rather
                # than thousands.  (Measured: ≤ 10.2 eps(T) ‖f̃‖ — worst at Matrix/Float64/
                # ℓₘₐₓ = 6 — so this leaves a factor of 9.8.)
                f̃ = randn(rng, Complex{T}, n)
                f = 𝒯 * f̃
                @test maximum(abs, f .- Y * f̃) ≤ 100eps(T) * norm(f̃)

                # Agreement with the package's own dense matrix of harmonics, and with the
                # pointwise evaluation of a ModeWeights.  These are internal cross-checks
                # rather than independent references; the closed form above is the reference.
                @test f ≈ sYlm_matrix(rotors(𝒯), ℓₘₐₓ, s) * f̃ rtol=100ℓₘₐₓ * eps(T)
                w = ModeWeights(f̃, s)
                R = rotors(𝒯)
                scale = maximum(abs, f)
                for k in 1:max(1, npixels(𝒯) ÷ 6):npixels(𝒯)
                    @test f[k] ≈ w(R[k]) atol=100ℓₘₐₓ * eps(T) * scale
                end
            end
        end
    end
end


@testitem "SSHT analysis" setup=[Utilities] begin
    import SphericalFunctions: SSHT, SSHTMatrix, pixels, rotors, Ysize, Yindex, ModeWeights, spin
    import SphericalFunctions: nmodes, npixels  # unexported
    using DoubleFloats: Double64
    using LinearAlgebra: mul!, ldiv!
    using StaticArrays: SVector
    using Random

    # `sYlm_pixels` — the closed-form ₛYₗₘ of the `Utilities` snippet on a list of pixels,
    # with the pixel-independent factorials hoisted out of the loop — comes from that
    # snippet; see "SSHT synthesis", where it is checked against `sYlm` itself.

    # Analysis tolerance.  Measured over every T, s ∈ -2:2 and single mode, in units of
    # eps(T): "RS" ≤ 8, "Matrix" ≤ 18, and round trips of random weights ≤ 37, all growing
    # only slowly with ℓₘₐₓ, so 100ℓₘₐₓ eps(T) leaves a factor of ≳ 13.  "Minimal" is a
    # different story: its linear system is increasingly ill-conditioned, and the error grows
    # by roughly a factor of 20 per unit ℓₘₐₓ — 47, 593, 7925 and 155955 eps(T) at ℓₘₐₓ = 3,
    # 4, 5 and 6.  The formula below tracks that growth with a factor of 30–150 to spare.
    #
    # The `min` is essential, not cosmetic.  Unclamped, "Minimal" at ℓₘₐₓ = 6 asks for
    # 2.4e7 eps(T), which is 2.86 in Float32 — larger than the unit coefficients being
    # tested, so an all-zero, sign-flipped or wrong-mode result would pass every assertion
    # below.  Capping at 0.05 keeps every assertion falsifiable (an all-zero result is off
    # by 1 and a sign flip by 2, i.e. 20× and 40× above the cap) while still admitting
    # everything measured for the combinations actually run: worst case Float32 "Minimal"
    # at ℓₘₐₓ = 5, single mode 9.4e-4 (margin 53×) and the three-algorithm cross-check
    # 2.6e-3 (margin 20×).
    tolerance(method, ℓₘₐₓ, ::Type{T}) where {T} = min(
        T(0.05),
        method == "Minimal" ? 500ℓₘₐₓ * 20.0^(ℓₘₐₓ - 3) * eps(T) : 100ℓₘₐₓ * eps(T)
    )

    rng = Random.Xoshiro(2718)

    for (method, T) in Iterators.product(("RS", "Minimal", "Matrix"), (Float64, Double64, Float32))
        kw = method == "RS" ? (;) : (; inplace=false)
        for ℓₘₐₓ in 3:6
            # The one combination no honest tolerance can accommodate: Float32 "Minimal" at
            # ℓₘₐₓ = 6 has a measured single-mode error of 1.9e-2 on a unit coefficient, so
            # even the capped tolerance would clear it by only 2.7×, and any tolerance that
            # admits it comfortably also admits a badly wrong answer.  (The conditioning is
            # intrinsic, not a regression: the v2 implementation — the `Deprecated` module,
            # removed in 3.0 — measured 1.97e5 eps there against v3's 1.56e5.)
            method == "Minimal" && T === Float32 && ℓₘₐₓ > 5 && continue
            ϵ = tolerance(method, ℓₘₐₓ, T)
            for s in -2:2
                𝒯 = SSHT(s, ℓₘₐₓ; method, T, kw...)
                n = nmodes(𝒯)
                p = pixels(𝒯)

                # Single modes evaluated from the closed form come back as unit vectors, as
                # ModeWeights
                for ℓ in abs(s):ℓₘₐₓ, m in -ℓ:ℓ
                    f = sYlm_pixels(s, ℓ, m, p)
                    f̃ = 𝒯 \ f
                    @test f̃ isa ModeWeights{Complex{T}}
                    @test spin(f̃) == s
                    @test SphericalFunctions.ℓₘᵢₙ(f̃) == abs(s)
                    @test SphericalFunctions.ℓₘₐₓ(f̃) == ℓₘₐₓ
                    @test length(f̃) == n
                    expected = zeros(Complex{T}, n)
                    expected[Yindex(ℓ, m, abs(s))] = one(T)
                    @test f̃ ≈ expected atol=ϵ rtol=ϵ
                    @test f̃[ℓ, m] ≈ one(T) atol=ϵ
                end

                # Round trip of random weights
                f̃ = randn(rng, Complex{T}, n)
                f = 𝒯 * f̃
                f̃′ = 𝒯 \ f
                @test f̃′ isa ModeWeights{Complex{T}}
                @test f̃′ ≈ f̃ atol=ϵ rtol=ϵ

                # ModeWeights in, ModeWeights out
                w = ModeWeights(f̃, s)
                # `ModeWeights` wraps its argument by reference (`parent(w) === f̃`), so
                # `w == f̃` could never fail; compare against a snapshot instead.
                f̃₀ = copy(f̃)
                fw = 𝒯 * w
                @test fw == f
                w′ = 𝒯 \ fw
                @test w′ isa ModeWeights{Complex{T}}
                @test spin(w′) == s
                @test SphericalFunctions.ℓₘᵢₙ(w′) == abs(s)
                @test SphericalFunctions.ℓₘₐₓ(w′) == ℓₘₐₓ
                @test w′ ≈ w atol=ϵ rtol=ϵ
                @test f̃ == f̃₀  # the input was not touched

                # Several columns at once: plain arrays in and out, each column independent
                F̃ = randn(rng, Complex{T}, n, 3)
                F = 𝒯 * F̃
                @test F isa Matrix{Complex{T}}
                @test size(F) == (npixels(𝒯), 3)
                @test all(F[:, j] ≈ 𝒯 * F̃[:, j] for j in 1:3)
                F̃′ = 𝒯 \ F
                @test F̃′ isa Matrix{Complex{T}}
                @test size(F̃′) == (n, 3)
                @test F̃′ ≈ F̃ atol=ϵ rtol=ϵ
                @test all(F̃′[:, j] ≈ 𝒯 \ F[:, j] for j in 1:3)
            end
        end
    end

    # The three methods are genuinely different algorithms — "RS" takes an FFT along each
    # ring and applies a quadrature rule across rings, "Minimal" solves one small system per
    # ring, "Matrix" factors one dense matrix — and they share no code path beyond the
    # harmonics themselves.  Handing the "RS" and "Minimal" sample points to "Matrix" puts
    # all three on a common grid, where they must agree.  (This is a cross-check inside the
    # package, not an independent reference; the closed-form comparisons above are that.)
    for T in (Float64, Double64, Float32), ℓₘₐₓ in (3, 5), s in (-2, 0, 1)
        n = Ysize(abs(s), ℓₘₐₓ)
        f̃ = randn(rng, Complex{T}, n)
        for method in ("RS", "Minimal")
            kw = method == "RS" ? (;) : (; inplace=false)
            𝒯 = SSHT(s, ℓₘₐₓ; method, T, kw...)
            𝒯ᴹ = SSHTMatrix(s, ℓₘₐₓ; T, Rθϕ=rotors(𝒯), inplace=false)
            @test npixels(𝒯ᴹ) == npixels(𝒯)
            # Measured over exactly these cases, in units of eps(T): synthesis ≤ 43,
            # analysis ≤ 20 against "RS", and — with the same ill-conditioning as above —
            # ≤ 8900 against "Minimal".  Margins are 12× and better.
            ϵ = tolerance(method, ℓₘₐₓ, T)
            f = 𝒯 * f̃
            @test f ≈ 𝒯ᴹ * f̃ atol=ϵ rtol=ϵ
            @test 𝒯 \ copy(f) ≈ 𝒯ᴹ \ copy(f) atol=ϵ rtol=ϵ
        end
    end

    # Inputs of the wrong length, or ModeWeights with the wrong ℓ range, are errors for every
    # method, in `*`, `\`, `mul!` and `ldiv!`
    for method in ("RS", "Minimal", "Matrix"), inplace in (true, false)
        s, ℓₘₐₓ = -1, 4
        kw = method == "RS" ? (;) : (; inplace)
        𝒯 = SSHT(s, ℓₘₐₓ; method, kw...)
        n, N = nmodes(𝒯), npixels(𝒯)
        @test_throws "first dimension of the mode weights" 𝒯 * zeros(ComplexF64, n + 1)
        @test_throws ErrorException 𝒯 * zeros(ComplexF64, n - 1)
        @test_throws ErrorException 𝒯 * zeros(ComplexF64, n + 1, 2)
        if N != n  # for "Minimal" and "Matrix" the two lengths agree, so this is legal
            @test_throws ErrorException 𝒯 * zeros(ComplexF64, N)
        end
        @test_throws "first dimension of the function values" 𝒯 \ zeros(ComplexF64, N + 1)
        @test_throws ErrorException 𝒯 \ zeros(ComplexF64, N - 1)
        @test_throws ErrorException 𝒯 \ zeros(ComplexF64, N - 1, 2)
        # ModeWeights with ℓₘᵢₙ = 0 (length Ysize(0, ℓₘₐₓ)) or the wrong ℓₘₐₓ
        @test_throws "ModeWeights have ℓ ∈ 0:4" 𝒯 * ModeWeights(zeros(ComplexF64, Ysize(0, ℓₘₐₓ)), s; ℓₘᵢₙ=0)
        @test_throws "ModeWeights have ℓ ∈ 1:5" 𝒯 * ModeWeights(zeros(ComplexF64, Ysize(1, ℓₘₐₓ + 1)), s)
        @test_throws ErrorException 𝒯 * ModeWeights(zeros(ComplexF64, Ysize(1, ℓₘₐₓ - 1)), s)
        @test_throws ErrorException mul!(zeros(ComplexF64, N + 1), 𝒯, zeros(ComplexF64, n))
        @test_throws ErrorException mul!(zeros(ComplexF64, N), 𝒯, zeros(ComplexF64, n + 1))
        @test_throws ErrorException ldiv!(zeros(ComplexF64, n + 1), 𝒯, zeros(ComplexF64, N))
        @test_throws ErrorException ldiv!(zeros(ComplexF64, n), 𝒯, zeros(ComplexF64, N + 1))
        @test_throws ErrorException ldiv!(ModeWeights(zeros(ComplexF64, Ysize(0, ℓₘₐₓ)), s; ℓₘᵢₙ=0), 𝒯, zeros(ComplexF64, N))
    end
end


@testitem "SSHT in-place semantics" begin
    import SphericalFunctions: SSHT, SSHTRS, SSHTMinimal, SSHTMatrix, ModeWeights, Ysize, spin
    import SphericalFunctions: nmodes, npixels  # unexported
    using DoubleFloats: Double64
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(1618)

    for T in (Float64, Double64), (s, ℓₘₐₓ) in ((-1, 4), (2, 5), (0, 3))
        n = Ysize(abs(s), ℓₘₐₓ)
        ϵ = 500ℓₘₐₓ^3 * eps(T) * 50
        f̃0 = randn(rng, Complex{T}, n)

        for method in ("Minimal", "Matrix")
            𝒯 = SSHT(s, ℓₘₐₓ; method, T)  # inplace=true is the default
            𝒯n = SSHT(s, ℓₘₐₓ; method, T, inplace=false)
            InplaceType = method == "Minimal" ? SSHTMinimal : SSHTMatrix
            @test 𝒯 isa InplaceType{T, true}
            @test 𝒯n isa InplaceType{T, false}

            # Reference results from the non-in-place object, which leaves its inputs alone
            f̃ = copy(f̃0)
            fref = 𝒯n * f̃
            @test f̃ == f̃0
            @test fref !== f̃
            @test fref isa Vector{Complex{T}}
            g = copy(fref)
            f̃ref = 𝒯n \ g
            @test g == fref
            @test f̃ref isa ModeWeights{Complex{T}}
            @test spin(f̃ref) == s
            @test f̃ref ≈ f̃0 atol=ϵ rtol=ϵ

            # In-place synthesis.  SSHTMinimal overwrites and returns its input.  SSHTMatrix's
            # `*` is a plain matrix-vector product, so it allocates even for the in-place type
            # (v2's "Direct" did the same); only its `\` acts in place.
            f̃ = copy(f̃0)
            f = 𝒯 * f̃
            if method == "Minimal"
                @test f === f̃
                @test f̃ == fref
            else
                @test f !== f̃
                @test f̃ == f̃0
                @test f == fref
            end

            # In-place analysis: the input array is overwritten with the mode weights and
            # returned as is — a plain Vector, not a ModeWeights
            g = copy(fref)
            g̃ = 𝒯 \ g
            @test g̃ === g
            @test g̃ isa Vector{Complex{T}}
            @test g == parent(f̃ref)

            # A ModeWeights input is overwritten and returned the same way
            w = ModeWeights(copy(f̃0), s)
            fw = 𝒯 * w
            if method == "Minimal"
                @test fw === w
                @test parent(w) == fref
            else
                @test fw == fref
                @test parent(w) == f̃0
            end
            wf = ModeWeights(copy(fref), s)  # function values can be stored in one, as N == n
            w̃ = 𝒯 \ wf
            @test w̃ === wf
            @test parent(wf) == parent(f̃ref)

            # Two-argument `ldiv!(𝒯, x)` (and `mul!(𝒯, x)` for Minimal) act in place on the
            # non-in-place type as well
            g = copy(fref)
            @test ldiv!(𝒯n, g) === g
            @test g == parent(f̃ref)
            if method == "Minimal"
                g̃ = copy(f̃0)
                @test mul!(𝒯n, g̃) === g̃
                @test g̃ == fref
            end

            # Several columns, in place
            F̃ = hcat(f̃0, -f̃0)
            F = 𝒯 * copy(F̃)
            @test F ≈ hcat(fref, -fref) atol=ϵ rtol=ϵ
            G = copy(F)
            @test (𝒯 \ G) === G
            @test G ≈ F̃ atol=ϵ rtol=ϵ
        end

        # Three-argument `mul!`/`ldiv!` with explicit outputs work for every type (in place
        # or not) and never touch their inputs
        for method in ("RS", "Minimal", "Matrix")
            objects = method == "RS" ? (SSHT(s, ℓₘₐₓ; method, T),) :
                (SSHT(s, ℓₘₐₓ; method, T), SSHT(s, ℓₘₐₓ; method, T, inplace=false))
            for 𝒯 in objects
                N = npixels(𝒯)
                fref = SSHT(s, ℓₘₐₓ; method, T, (method == "RS" ? (;) : (; inplace=false))...) * f̃0
                f̃ = copy(f̃0)
                f = zeros(Complex{T}, N)
                @test mul!(f, 𝒯, f̃) === f
                @test f̃ == f̃0
                @test f == fref
                g = copy(fref)
                g̃ = zeros(Complex{T}, n)
                @test ldiv!(g̃, 𝒯, g) === g̃
                @test g == fref
                @test g̃ ≈ f̃0 atol=ϵ rtol=ϵ
                # ... including into a ModeWeights output
                w̃ = ModeWeights{Complex{T}}(undef, s, ℓₘₐₓ)
                @test ldiv!(w̃, 𝒯, g) === w̃
                @test parent(w̃) == g̃
                @test g == fref
                # ... and from a ModeWeights input
                w = ModeWeights(copy(f̃0), s)
                fill!(f, zero(Complex{T}))
                @test mul!(f, 𝒯, w) === f
                @test parent(w) == f̃0
                @test f == fref
                # ... and with several columns
                F̃ = hcat(f̃0, 2f̃0, -f̃0)
                F = zeros(Complex{T}, N, 3)
                @test mul!(F, 𝒯, F̃) === F
                @test F̃ == hcat(f̃0, 2f̃0, -f̃0)
                # (Not bitwise: a matrix-matrix product need not round like a
                # matrix-vector product.)
                @test F[:, 1] ≈ fref atol=ϵ rtol=ϵ
                @test F ≈ hcat(fref, 2fref, -fref) atol=ϵ rtol=ϵ
                G̃ = zeros(Complex{T}, n, 3)
                @test ldiv!(G̃, 𝒯, F) === G̃
                @test G̃ ≈ F̃ atol=ϵ rtol=ϵ
            end
        end
    end
end


@testitem "SSHT trailing dimensions" begin
    import SphericalFunctions: SSHT, Ysize
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(4242)

    # The docstring promises that any dimensions after the first are broadcast over.  Two
    # trailing dimensions exercise the reshaping that a matrix input does not.
    for method in ("RS", "Minimal", "Matrix"), (s, ℓₘₐₓ) in ((1, 4), (-2, 3))
        ϵ = 500ℓₘₐₓ^3 * eps(Float64) * (method == "Minimal" ? 50 : 1)
        kw = method == "RS" ? (;) : (; inplace=false)
        𝒯 = SSHT(s, ℓₘₐₓ; method, kw...)
        n, N = nmodes(𝒯), npixels(𝒯)
        F̃ = randn(rng, ComplexF64, n, 3, 2)
        columns = [𝒯 * F̃[:, j, k] for j in 1:3, k in 1:2]

        F = 𝒯 * F̃
        @test size(F) == (N, 3, 2)
        @test all(F[:, j, k] ≈ columns[j, k] for j in 1:3, k in 1:2)
        F̃′ = 𝒯 \ F
        @test size(F̃′) == (n, 3, 2)
        @test F̃′ ≈ F̃ atol=ϵ rtol=ϵ

        G = similar(F)
        @test mul!(G, 𝒯, F̃) === G
        @test G ≈ F
        G̃ = similar(F̃)
        @test ldiv!(G̃, 𝒯, F) === G̃
        @test G̃ ≈ F̃ atol=ϵ rtol=ϵ

        if method != "RS"  # the in-place types too
            𝒯i = SSHT(s, ℓₘₐₓ; method)
            H = 𝒯i * copy(F̃)
            @test size(H) == (N, 3, 2)
            @test H ≈ F
            H̃ = 𝒯i \ copy(F)
            @test size(H̃) == (n, 3, 2)
            @test H̃ ≈ F̃ atol=ϵ rtol=ϵ
        end
    end
end


@testitem "SSHTRS options" setup=[Utilities] begin
    import SphericalFunctions: SSHT, SSHTRS, pixels, rotors, Ysize, Yindex, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    import SphericalFunctions: fejer1_rings, fejer2_rings, clenshaw_curtis_rings
    import SphericalFunctions: fejer1, fejer2, clenshaw_curtis
    using DoubleFloats: Double64
    using FFTW: FFTW
    using Logging: NullLogger, with_logger
    using Quaternionic: from_spherical_coordinates
    using StaticArrays: SVector
    using Random

    # `sYlm_pixels` — the closed-form ₛYₗₘ of the `Utilities` snippet on a list of pixels,
    # with the pixel-independent factorials hoisted out of the loop — comes from that
    # snippet; see "SSHT synthesis", where it is checked against `sYlm` itself.

    rng = Random.Xoshiro(577)

    for T in (Float64, Double64, Float32), ℓₘₐₓ in (3, 5), s in -2:2
        # Measured against the closed form over every quadrature rule, Nϕ, T, ℓₘₐₓ and s
        # below, in units of eps(T): synthesis ≤ 13, analysis of a single mode ≤ 7, round
        # trip of random weights ≤ 14.  This leaves a factor of ≳ 16.
        ϵ = 50ℓₘₐₓ * eps(T)
        n = Ysize(abs(s), ℓₘₐₓ)
        f̃ = randn(rng, Complex{T}, n)
        N = 2ℓₘₐₓ + 1
        quadratures = (
            ("fejer1", fejer1_rings(N, T), fejer1(N, T)),
            ("fejer2", fejer2_rings(N, T), fejer2(N, T)),
            ("clenshaw_curtis", clenshaw_curtis_rings(N, T), clenshaw_curtis(N, T)),
            ("fejer1 with extra rings", fejer1_rings(N + 3, T), fejer1(N + 3, T)),
        )
        for (label, θ, w) in quadratures
            Nθ = length(θ)
            # Nϕ as an Int (2ℓₘₐₓ+1, and 2ℓₘₐₓ+2 for a more FFT-friendly count) or a Vector
            # (uniform, or varying from ring to ring)
            for Nϕ in (N, N + 1, fill(N, Nθ), [N + k for k in 0:Nθ-1])
                𝒯 = SSHT(s, ℓₘₐₓ; method="RS", T, θ, quadrature_weights=w, Nϕ)
                @test 𝒯 isa SSHTRS{T}
                Nϕs = Nϕ isa Integer ? fill(Nϕ, Nθ) : Nϕ
                @test 𝒯.Nϕ == Nϕs
                @test npixels(𝒯) == sum(Nϕs)
                p = pixels(𝒯)
                @test p == [
                    SVector(θⱼ, k * 2T(π) / Nϕⱼ)
                    for (θⱼ, Nϕⱼ) in zip(θ, Nϕs) for k in 0:Nϕⱼ-1
                ]
                @test rotors(𝒯) == from_spherical_coordinates.(p)

                # Exact quadrature for band-limited functions: the round trip is exact
                f = 𝒯 * f̃
                f̃′ = 𝒯 \ f
                @test f̃′ isa ModeWeights{Complex{T}}
                @test f̃′ ≈ f̃ atol=ϵ rtol=ϵ

                # Single modes against the closed-form harmonics evaluated on this grid, and
                # analysis inverting that synthesis.  `fill(N, Nθ)` describes the very same
                # grid as the integer `N`, so it is skipped rather than measured twice.
                if !(Nϕ isa AbstractVector && allequal(Nϕ) && Nϕ[begin] == N)
                    for ℓ in abs(s):ℓₘₐₓ, m in -ℓ:ℓ
                        g̃ = zeros(Complex{T}, n)
                        g̃[Yindex(ℓ, m, abs(s))] = one(T)
                        g = 𝒯 * g̃
                        @test g ≈ sYlm_pixels(s, ℓ, m, p) atol=ϵ rtol=ϵ
                        @test 𝒯 \ g ≈ g̃ atol=ϵ rtol=ϵ
                    end
                end
            end
        end
    end

    # Rings and weights given in Float64 are converted to T
    𝒯 = SSHTRS(1, 5; T=Float32, θ=fejer1_rings(11), quadrature_weights=fejer1(11))
    @test 𝒯 isa SSHTRS{Float32}
    @test eltype(pixels(𝒯)) === SVector{2, Float32}
    @test 𝒯.θ == Float32.(fejer1_rings(11))

    # FFTW planner options are accepted and do not change the results
    𝒯 = SSHTRS(1, 5)
    𝒯m = SSHTRS(1, 5; plan_fft_flags=FFTW.MEASURE, plan_fft_timelimit=1.0)
    f̃ = randn(rng, ComplexF64, nmodes(𝒯))
    @test 𝒯m * f̃ ≈ 𝒯 * f̃ rtol=100eps()
    @test 𝒯m \ (𝒯m * f̃) ≈ f̃ rtol=500 * 5^3 * eps()

    # Too few points on a ring: a warning at construction, and aliasing of the largest |m|
    s, ℓₘₐₓ = 1, 5
    @test_logs (:warn, r"fewer than 2ℓₘₐₓ\+1") SSHTRS(s, ℓₘₐₓ; Nϕ=2ℓₘₐₓ)
    @test_logs (:warn, r"alias") SSHTRS(s, ℓₘₐₓ; Nϕ=[2ℓₘₐₓ + 1 - (k == 3) for k in 1:2ℓₘₐₓ+1])
    @test_logs SSHTRS(s, ℓₘₐₓ; Nϕ=2ℓₘₐₓ + 1)  # no warning at the band limit
    𝒯 = with_logger(NullLogger()) do
        SSHTRS(s, ℓₘₐₓ; Nϕ=2ℓₘₐₓ)
    end
    @test npixels(𝒯) == (2ℓₘₐₓ) * (2ℓₘₐₓ + 1)
    f̃ = randn(rng, ComplexF64, nmodes(𝒯))
    f̃[Yindex(ℓₘₐₓ, -ℓₘₐₓ, abs(s))] = 0
    f̃[Yindex(ℓₘₐₓ, ℓₘₐₓ, abs(s))] = 0
    @test 𝒯 \ (𝒯 * f̃) ≈ f̃ atol=500ℓₘₐₓ^3 * eps() rtol=500ℓₘₐₓ^3 * eps()  # |m| < ℓₘₐₓ is fine
    f̃[Yindex(ℓₘₐₓ, ℓₘₐₓ, abs(s))] = 1
    @test !isapprox(𝒯 \ (𝒯 * f̃), f̃; atol=1e-3)  # m = ±ℓₘₐₓ alias onto each other

    # Inconsistent options are errors
    @test_throws "same length" SSHTRS(s, ℓₘₐₓ; θ=fejer1_rings(11), quadrature_weights=fejer1(12))
    @test_throws "same length" SSHTRS(s, ℓₘₐₓ; θ=fejer1_rings(12), quadrature_weights=fejer1(11))
    @test_throws "Nϕ must be a single number or have the same length" SSHTRS(s, ℓₘₐₓ; Nϕ=fill(11, 10))
    @test_throws ErrorException SSHTRS(s, ℓₘₐₓ; Nϕ=fill(11, 12))
    @test_throws "at least one point" SSHTRS(s, ℓₘₐₓ; Nϕ=0)
    @test_throws ErrorException SSHTRS(s, ℓₘₐₓ; Nϕ=[11, 11, 11, 11, 11, 0, 11, 11, 11, 11, 11])
end


@testitem "SSHTMatrix options" begin
    import SphericalFunctions: SSHT, SSHTMatrix, pixels, rotors, Ysize, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    import SphericalFunctions: golden_ratio_spiral_rotors, sorted_ring_rotors
    using DoubleFloats: Double64
    using LinearAlgebra: LinearAlgebra, lu, qr
    using Random

    rng = Random.Xoshiro(8128)

    for T in (Float64, Double64), (s, ℓₘₐₓ) in ((-2, 5), (1, 4))
        ϵ = 500ℓₘₐₓ^3 * eps(T)
        n = Ysize(abs(s), ℓₘₐₓ)
        f̃ = randn(rng, Complex{T}, n)

        # Any set of rotors may be used; with exactly n of them the default is LU and in place
        R = sorted_ring_rotors(s, ℓₘₐₓ, T)
        𝒯 = SSHTMatrix(s, ℓₘₐₓ; T, Rθϕ=R)
        @test 𝒯 isa SSHTMatrix{T, true, <:LinearAlgebra.LU}
        @test rotors(𝒯) == R
        @test npixels(𝒯) == n
        f = 𝒯 * f̃
        @test 𝒯 \ copy(f) ≈ f̃ atol=ϵ rtol=ϵ

        # An explicit decomposition, e.g. QR, must also solve the problem
        𝒯qr = SSHTMatrix(s, ℓₘₐₓ; T, decomposition=qr, inplace=false)
        @test 𝒯qr isa SSHTMatrix{T, false}
        @test !(𝒯qr isa SSHTMatrix{T, false, <:LinearAlgebra.LU})
        @test 𝒯qr * f̃ == SSHTMatrix(s, ℓₘₐₓ; T, decomposition=lu, inplace=false) * f̃
        @test 𝒯qr \ (𝒯qr * f̃) ≈ f̃ atol=ϵ rtol=ϵ

        # More points than modes: a least-squares analysis, QR by default, never in place
        Rmore = [golden_ratio_spiral_rotors(s, ℓₘₐₓ, T); sorted_ring_rotors(0, ℓₘₐₓ + 1, T)]
        𝒯ls = SSHTMatrix(s, ℓₘₐₓ; T, Rθϕ=Rmore)
        @test 𝒯ls isa SSHTMatrix{T, false}
        @test !(𝒯ls isa SSHTMatrix{T, false, <:LinearAlgebra.LU})
        @test npixels(𝒯ls) == length(Rmore) > nmodes(𝒯ls) == n
        fls = 𝒯ls * f̃
        @test length(fls) == length(Rmore)
        f̃′ = 𝒯ls \ fls
        @test f̃′ isa ModeWeights{Complex{T}}
        @test f̃′ ≈ f̃ atol=ϵ rtol=ϵ
        # ... which is the use case of the docstring: sample on the grid appropriate to the
        # *lowest* |s| (the most modes, here s = 0) and reuse those points for this spin
        if s != 0
            𝒯s = SSHTMatrix(s, ℓₘₐₓ; T, Rθϕ=golden_ratio_spiral_rotors(0, ℓₘₐₓ, T))
            @test 𝒯s isa SSHTMatrix{T, false}
            @test npixels(𝒯s) == Ysize(0, ℓₘₐₓ) > nmodes(𝒯s) == n
            g̃ = randn(rng, Complex{T}, nmodes(𝒯s))
            @test 𝒯s \ (𝒯s * g̃) ≈ g̃ atol=ϵ rtol=ϵ
        end

        # Fewer points than modes is underdetermined, and in place needs exactly n points
        @test_throws "underdetermined" SSHTMatrix(s, ℓₘₐₓ; T, Rθϕ=R[1:end-1])
        @test_throws "exactly as many sample points as modes" SSHTMatrix(s, ℓₘₐₓ; T, Rθϕ=Rmore, inplace=true)
    end
end


@testitem "SSHT thread safety via separate objects" begin
    import SphericalFunctions: SSHT, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(8675309)

    # Each SSHT holds its own workspace, so one object per task gives results identical to a
    # serial computation (the algorithms are deterministic: no cross-thread reductions).
    for method in ("RS", "Minimal", "Matrix"), T in (Float64, Float32)
        s, ℓₘₐₓ = -2, 8
        kw = method == "RS" ? (;) : (; inplace=false)

        # Two objects of the same parameters on two tasks, as in the assignment
        𝒯s = [SSHT(s, ℓₘₐₓ; method, T, kw...) for _ in 1:2]
        F̃s = [randn(rng, Complex{T}, nmodes(𝒯s[1]), 3) for _ in 1:2]
        serial = map(F̃s) do F̃
            F = 𝒯s[1] * F̃
            (F, 𝒯s[1] \ F)
        end
        tasks = map(zip(𝒯s, F̃s)) do (𝒯, F̃)
            Threads.@spawn begin
                F = 𝒯 * F̃
                (F, 𝒯 \ F)
            end
        end
        parallel = fetch.(tasks)
        @test parallel == serial

        # Eight objects on eight tasks, each transforming many inputs (including in-place
        # `mul!`/`ldiv!` and ModeWeights), to raise the odds of overlapping execution
        𝒯s = [SSHT(s, ℓₘₐₓ; method, T, kw...) for _ in 1:8]
        inputs = [[randn(rng, Complex{T}, nmodes(𝒯s[1])) for _ in 1:12] for _ in 1:8]
        function work(𝒯, f̃s)
            out = Vector{Vector{Complex{T}}}(undef, 3length(f̃s))
            f = zeros(Complex{T}, npixels(𝒯))
            g̃ = ModeWeights{Complex{T}}(undef, s, ℓₘₐₓ)
            for (i, f̃) in enumerate(f̃s)
                mul!(f, 𝒯, f̃)
                out[3i - 2] = copy(f)
                ldiv!(g̃, 𝒯, f)
                out[3i - 1] = copy(parent(g̃))
                out[3i] = collect(𝒯 \ (𝒯 * ModeWeights(f̃, s)))
            end
            out
        end
        serial = [work(𝒯s[1], f̃s) for f̃s in inputs]
        parallel = fetch.([Threads.@spawn work(𝒯, f̃s) for (𝒯, f̃s) in zip(𝒯s, inputs)])
        @test parallel == serial
    end
end

@testitem "SSHT half-integer construction" begin
    using SphericalFunctions
    using SphericalFunctions: HalfOddInteger, nmodes, npixels

    𝒯 = SSHT(1//2, 7//2)
    @test 𝒯 isa SSHTRS{Float64}
    @test spin(𝒯) === HalfOddInteger(1//2)
    @test SphericalFunctions.ℓₘₐₓ(𝒯) === HalfOddInteger(7//2)
    @test SphericalFunctions.ℓₘᵢₙ(𝒯) === HalfOddInteger(1//2)
    @test nmodes(𝒯) == Ysize(1//2, 7//2) == 20
    @test npixels(𝒯) == 8 * 8  # 2ℓₘₐₓ+1 rings of 2ℓₘₐₓ+1 points, both even numbers
    @test occursin("s=1//2", sprint(show, 𝒯)) && occursin("ℓₘₐₓ=7//2", sprint(show, 𝒯))
    𝒯m = SSHT(-3//2, 7//2; method="Matrix")
    @test 𝒯m isa SSHTMatrix{Float64, true}
    @test npixels(𝒯m) == nmodes(𝒯m) == Ysize(3//2, 7//2)
    @test SSHTMatrix(1//2, 5//2; T=Float32) isa SSHTMatrix{Float32}
    @test SSHTRS(1//2, 5//2; T=Float32) isa SSHTRS{Float32}

    # The spellings agree, and the two kinds of index may not be mixed
    @test typeof(SSHT(HalfOddInteger(1//2), HalfOddInteger(7//2))) === typeof(𝒯)
    @test_throws "must all be integers, like 3, or all be half-odd-integers" SSHT(1//2, 4)
    @test_throws "must all be integers, like 3, or all be half-odd-integers" SSHTMatrix(1, 7//2)
    @test_throws "must all be integers, like 3, or all be half-odd-integers" SSHTRS(1//2, 4)
    @test_throws "must all be integers, like 3, or all be half-odd-integers" map2salm(zeros(ComplexF64, 8, 8), 1//2, 3)
    @test_throws "exceeds ℓₘₐₓ" SSHT(5//2, 3//2)

    # The Minimal method is integer-only (decision E3), and names the methods that are not
    @test_throws "\"RS\"" SSHT(1//2, 7//2; method="Minimal")
    @test_throws "\"Matrix\"" SSHTMinimal(1//2, 7//2)

    # The integer path stores `Int`, as it always did
    @test SSHT(Int8(1), Int8(4)).s === 1
    @test SSHT(1, 4; method="Minimal").s === 1 && SSHT(1, 4; method="Matrix").ℓₘₐₓ === 4
end

@testitem "SSHT half-integer round trips and cross-method agreement" begin
    using SphericalFunctions
    using SphericalFunctions: HalfOddInteger
    using Quaternionic
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = MersenneTwister(2026_09_17)
    for s ∈ (1//2, -1//2, 3//2, -5//2), ℓₘₐₓ ∈ (abs(s), 7//2, 11//2)
        ℓₘₐₓ ≥ abs(s) || continue
        f̃ = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
        𝒯r = SSHT(s, ℓₘₐₓ)
        𝒯m = SSHT(s, ℓₘₐₓ; method="Matrix")

        # Synthesis by the ring algorithm is the direct evaluation of the harmonics at its own
        # rotors, which is the strongest check available: the two computations share nothing
        # but the calculator
        fr = 𝒯r * f̃
        @test fr ≈ sYlm_matrix(rotors(𝒯r), ℓₘₐₓ, s) * parent(f̃) rtol=1e-12
        fm = 𝒯m * copy(f̃)
        @test fm ≈ sYlm_matrix(rotors(𝒯m), ℓₘₐₓ, s) * parent(f̃) rtol=1e-12

        # Analysis inverts synthesis, for both methods, and labels its result
        g̃r = 𝒯r \ fr
        @test g̃r isa ModeWeights
        @test spin(g̃r) === HalfOddInteger(s) && SphericalFunctions.ℓₘₐₓ(g̃r) === HalfOddInteger(ℓₘₐₓ)
        @test parent(g̃r) ≈ parent(f̃) rtol=1e-11
        @test parent(𝒯m \ copy(fm)) ≈ parent(f̃) rtol=1e-11

        # The Matrix method on the ring grid agrees with the ring algorithm
        𝒯mr = SSHTMatrix(s, ℓₘₐₓ; Rθϕ=rotors(𝒯r))
        @test 𝒯mr * copy(f̃) ≈ fr rtol=1e-12
        @test parent(𝒯mr \ copy(fr)) ≈ parent(f̃) rtol=1e-10

        # The in-place forms compute the same values
        f2 = similar(fr)
        mul!(f2, 𝒯r, f̃)
        @test f2 == fr
        g2 = similar(f̃)
        ldiv!(g2, 𝒯r, fr)
        @test g2 == g̃r

        # Trailing dimensions are transformed column by column
        F̃ = hcat(parent(f̃), 2 .* parent(f̃))
        F = 𝒯r * F̃
        @test F[:, 1] ≈ fr && F[:, 2] ≈ 2fr
        @test 𝒯r \ F ≈ F̃ rtol=1e-11
    end
end

@testitem "SSHT half-integer values live on the double cover" begin
    using SphericalFunctions
    using Quaternionic

    s, ℓₘₐₓ = 1//2, 7//2
    𝒯 = SSHT(s, ℓₘₐₓ)
    f̃ = ModeWeights(ComplexF64.(1:Ysize(1//2, 7//2)), s)
    f = 𝒯 * f̃
    # The function values are those at `rotors(𝒯)`.  A full circuit of the azimuth reaches the
    # antipodal rotor, at which every harmonic — and so the function — changes sign.
    R₊ = [from_spherical_coordinates(θ, ϕ + 2π) for (θ, ϕ) ∈ pixels(𝒯)]
    @test sYlm_matrix(R₊, ℓₘₐₓ, s) * parent(f̃) ≈ -f rtol=1e-12
    @test sYlm_matrix(rotors(𝒯), ℓₘₐₓ, s) * parent(f̃) ≈ f rtol=1e-12
end
