# Tests for the spin-weighted spherical-harmonic transforms: the `SSHT` constructor and the
# three concrete types `SSHTRS`, `SSHTMinimal` and `SSHTMatrix`, their `pixels`/`rotors`,
# synthesis (`𝒯 * f̃`, `mul!`), analysis (`𝒯 \ f`, `ldiv!`), in-place semantics,
# `ModeWeights` input and output, the `SSHTRS` quadrature options, and use from several
# tasks, with `copy`, `deepcopy` and serialization.
#
# The oracle is the closed-form ₛYₗₘ of the `Utilities` module — the explicit sum over
# factorials from the conventions pages, which owes nothing to this package.  The items here
# use `sYlm_closed_form_pixels`, also from that module: the same formula with the
# pixel-independent factorials hoisted out of the loop (≈ 40 times faster, and checked
# against `sYlm_closed_form` itself in "SSHT synthesis").  They compare synthesis against it
# pixel by pixel, analysis against the mode it came from, and round trips against their
# input.  The pixelizations are compared against the formulas in their own docstrings.
# Cross-checks that stay inside the package — `sYlm_matrix`, the evaluation `w(R)` of a
# `ModeWeights`, and the agreement of the three algorithms on a shared grid — are labelled
# as such where they appear.
#
# Tolerances are set from measurement; the numbers quoted in the comments are the largest
# errors seen over every method, T, ℓₘₐₓ and s that the item covers, in units of eps(T).
#
# The item names all begin with "SSHT", so that they can be picked out together — by a human
# reading a results list, and by the name filters of `juliati` and the MCP runner.


@testitem "SSHT construction" begin
    import SphericalFunctions: SSHT, SSHTRS, SSHTMinimal, SSHTMatrix, pixels, rotors, spin, Ysize
    import SphericalFunctions: nmodes, npixels, floattype  # unexported
    import SphericalFunctions: golden_ratio_spiral_rotors
    using DoubleFloats: Double64
    using LinearAlgebra: lu
    using Logging: NullLogger, with_logger
    using Quaternionic: Rotor

    types = (("RS", SSHTRS), ("Minimal", SSHTMinimal), ("Matrix", SSHTMatrix))

    # Every method and type, including the edge cases s = ℓₘₐₓ, ℓₘₐₓ = 1 and ℓₘₐₓ = 0
    for T in (Float64, Double64, Float32), (s, ℓₘₐₓ) in ((-2, 6), (0, 3), (1, 1), (3, 3), (-3, 4), (0, 0))
        for (method, Type) in types
            𝒯 = SSHT(s, ℓₘₐₓ, T; method)
            @test 𝒯 isa Type{T}
            @test 𝒯 isa SSHT{T}
            @test typeof(Type(s, ℓₘₐₓ, T)) === typeof(𝒯)  # the concrete constructor agrees
            @test spin(𝒯) == s
            @test floattype(𝒯) === T
            @test floattype(typeof(𝒯)) === T
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
        @test npixels(SSHT(s, ℓₘₐₓ, T; method="RS")) == (2ℓₘₐₓ + 1)^2
        @test npixels(SSHT(s, ℓₘₐₓ, T; method="Minimal")) == Ysize(abs(s), ℓₘₐₓ)
        @test npixels(SSHT(s, ℓₘₐₓ, T; method="Matrix")) == Ysize(abs(s), ℓₘₐₓ)
    end

    # Defaults: method "RS", T Float64
    @test SSHT(1, 4) isa SSHTRS{Float64}
    @test SSHT(1, 4; method="Minimal") isa SSHTMinimal{Float64}
    @test SSHT(1, 4; method="Matrix") isa SSHTMatrix{Float64}

    # The element type is the third positional argument, of every constructor; there is no
    # keyword for it
    for (method, Type) in types
        @test SSHT(1, 4, Float32; method) isa Type{Float32}
        @test_throws MethodError SSHT(1, 4; method, T=Float32)
        @test_throws MethodError Type(1, 4; T=Float32)
        # ... and it must be a floating-point type, which it computes in
        @test_throws ArgumentError SSHT(1, 4, Int; method)
        @test_throws "T=$Int is not one" Type(1, 4, Int)
    end

    # The `inplace` option is part of the type of the in-place-capable methods (default
    # true), and so must be a `Bool`
    @test SSHT(1, 4; method="Minimal") isa SSHTMinimal{Float64, true}
    @test SSHT(1, 4; method="Minimal", inplace=false) isa SSHTMinimal{Float64, false}
    @test SSHT(1, 4; method="Matrix") isa SSHTMatrix{Float64, true}
    @test SSHT(1, 4; method="Matrix", inplace=false) isa SSHTMatrix{Float64, false}
    @test SSHTMinimal(1, 4, Float32; inplace=false) isa SSHTMinimal{Float32, false}
    @test SSHTMatrix(1, 4, Float32; inplace=false) isa SSHTMatrix{Float32, false}
    for Type in (SSHTMinimal, SSHTMatrix)
        @test_throws ArgumentError Type(0, 2; inplace=1)
        @test_throws "inplace=1 must be `true` or `false`" Type(0, 2; inplace=1)
    end

    # "Direct" is a deprecated alias of "Matrix"
    @test_deprecated r"renamed \"Matrix\"" SSHT(0, 3; method="Direct")
    if Base.JLOptions().depwarn != 2  # with --depwarn=error the call throws instead
        𝒯 = with_logger(NullLogger()) do
            SSHT(0, 3; method="Direct")
        end
        @test 𝒯 isa SSHTMatrix{Float64, true}
        @test rotors(𝒯) == rotors(SSHT(0, 3; method="Matrix"))
    end

    # The method may be a `Symbol` as well as a string
    for (method, Type) in types
        @test SSHT(1, 4; method=Symbol(method)) isa Type{Float64}
    end
    # A nonsense method is an error naming the offending value as it was written; the names
    # are case-sensitive
    @test_throws ArgumentError SSHT(-2, 8; method="NonsenseGarbage")
    @test_throws "Unrecognized s-SHT method \"NonsenseGarbage\"" SSHT(-2, 8; method="NonsenseGarbage")
    @test_throws "Unrecognized s-SHT method :rs;" SSHT(-2, 8; method=:rs)
    @test_throws ArgumentError SSHT(-2, 8; method="rs")
    @test_throws ArgumentError SSHT(-2, 8; method="matrix")

    # |s| > ℓₘₐₓ is an error for every method (there are no such modes), including every s
    # when ℓₘₐₓ is negative, and it is reached before the default `Rθϕ=leja_rotors(s, ℓₘₐₓ,
    # T)` of the "Matrix" method is computed.
    for method in ("RS", "Minimal", "Matrix"), (s, ℓₘₐₓ) in ((3, 2), (-3, 2), (1, 0), (0, -1))
        @test_throws ArgumentError SSHT(s, ℓₘₐₓ; method)
        @test_throws "|s|=$(abs(s)) exceeds ℓₘₐₓ=$ℓₘₐₓ" SSHT(s, ℓₘₐₓ; method)
    end
    @test_throws "exceeds ℓₘₐₓ" SSHTMatrix(3, 2; Rθϕ=golden_ratio_spiral_rotors(0, 3))
    @test_throws "exceeds ℓₘₐₓ" SSHTMinimal(3, 2; θ=Float64[])
    @test_throws "exceeds ℓₘₐₓ" SSHTRS(3, 2)

    # "Minimal" refuses two rings at one colatitude, which would make its system singular,
    # including a zero of either sign and colatitudes given as integers
    for (s, θ) ∈ ((0, [1.0, 1.0]), (1, [0.5, 2.0, 0.5]), (-2, [0.0, 1.0, -0.0]), (0, [2, 1, 2]))
        ℓₘₐₓ = abs(s) + length(θ) - 1
        @test_throws ArgumentError SSHTMinimal(s, ℓₘₐₓ; θ)
        @test_throws "θ must be distinct" SSHT(s, ℓₘₐₓ; method="Minimal", θ)
    end
    # ... and explains any other layout that makes it singular, such as a ring of more than
    # one point at a pole
    for (s, θ) ∈ ((0, [2.4, 0.0, 1.6]), (1, [0.0, 1.6, 2.4]), (-2, [0.8, 1.6, 0.0]))
        ℓₘₐₓ = abs(s) + length(θ) - 1
        @test_throws ArgumentError SSHTMinimal(s, ℓₘₐₓ; θ)
        @test_throws "their system is singular" SSHT(s, ℓₘₐₓ; method="Minimal", θ)
    end

    # "Matrix" warns when the dense matrix gets large (Ysize² > 65⁴, i.e. ℓₘₐₓ ≥ 65 for s =
    # 0).  The warning comes before the count of the points is checked, so a single point
    # shows it without the cost of building the transform.
    @test_logs (:warn, r"\"Matrix\" method for s-SHT is only") match_mode=:any (
        @test_throws "underdetermined" SSHTMatrix(0, 65; Rθϕ=[Rotor(1.0, 0.0, 0.0, 0.0)])
    )
    # ... and is quiet otherwise, as are the other methods
    @test_logs SSHTMatrix(0, 8)
    @test_logs SSHT(0, 8; method="RS")
    @test_logs SSHT(0, 8; method="Minimal")
    # Too few points are refused, and so is an LU decomposition of more points than modes,
    # which could not solve the least-squares problem
    @test_throws DimensionMismatch SSHTMatrix(0, 2; Rθϕ=golden_ratio_spiral_rotors(0, 1))
    @test_throws ArgumentError SSHTMatrix(0, 2; Rθϕ=golden_ratio_spiral_rotors(0, 3), decomposition=lu, inplace=false)
    @test_throws "LU decomposition cannot solve the least-squares problem" SSHTMatrix(
        0, 2; Rθϕ=golden_ratio_spiral_rotors(0, 3), decomposition=lu, inplace=false
    )

    # "Minimal" and "Matrix" warn when their sample points are too badly conditioned for
    # half the digits of T to survive a round trip, which for "Minimal" happens at large
    # ℓₘₐₓ — sooner for s ≠ 0 — and is quiet below.  (With the rings of `sorted_rings`
    # rather than `minimal_rings`, s = 2 warned already at ℓₘₐₓ = 10, and was garbage by
    # ℓₘₐₓ = 14.  At s = 0, measured: 1.3e-8 at ℓₘₐₓ = 64, 1.7e-8 at 66, and 1.6e-6 at 80.)
    @test_logs (:warn, r"\"Minimal\" s-SHT with s=2, ℓₘₐₓ=32 and T=Float64 is inaccurate") SSHT(2, 32; method="Minimal")
    @test_logs (:warn, r"is inaccurate") SSHT(0, 80; method="Minimal", inplace=false)
    @test_logs SSHT(0, 64; method="Minimal", inplace=false)
    @test_logs (:warn, r"T=Float32 is inaccurate") SSHT(2, 24, Float32; method="Minimal")
    @test_logs SSHT(2, 16; method="Minimal")
    @test_logs SSHT(-2, 24; method="Minimal")
    @test_logs SSHT(0, 32; method="Minimal")
    @test_logs SSHT(2, 16, Float32; method="Minimal")
    # For "Matrix" the points are the caller's to choose, and the golden-ratio spiral, which
    # spreads its points evenly, is badly conditioned by ℓₘₐₓ = 48 (a round trip measured
    # 3.2e-5), while the default Leja points are not (5.6e-13; and 1.4e-4 in Float32 at ℓₘₐₓ
    # = 32, against √eps(Float32) = 3.5e-4)
    @test_logs (:warn, r"\"Matrix\" s-SHT with s=0, ℓₘₐₓ=48 and T=Float64 is inaccurate") SSHTMatrix(
        0, 48; Rθϕ=golden_ratio_spiral_rotors(0, 48), inplace=false
    )
    @test_logs SSHTMatrix(0, 48; inplace=false)
    @test_logs SSHTMatrix(0, 32, Float32; inplace=false)
end


@testitem "SSHT pixels and rotors" begin
    import SphericalFunctions: SSHT, pixels, rotors
    import SphericalFunctions: npixels  # unexported
    import SphericalFunctions: sorted_rings, sorted_ring_pixels, sorted_ring_rotors, minimal_rings
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors, fejer1_rings
    import SphericalFunctions: leja_rotors
    using DoubleFloats: Double64
    using Quaternionic: Rotor, from_spherical_coordinates, to_spherical_coordinates
    using StaticArrays: SVector
    using Random

    rng = Random.Xoshiro(1729)

    for T in (Float64, Double64, Float32), ℓₘₐₓ in (3, 4, 5, 8, 13), s in -2:2
        # Minimal: the rings of `minimal_rings`, in its order, each of Nϕ equally spaced points
        # starting at ϕ = 0
        𝒯 = SSHT(s, ℓₘₐₓ, T; method="Minimal")
        p = pixels(𝒯)
        @test eltype(p) === SVector{2, T}
        @test eltype(rotors(𝒯)) === Rotor{T}
        @test rotors(𝒯) == from_spherical_coordinates.(p)
        rings = minimal_rings(s, ℓₘₐₓ, T)
        @test 𝒯.θ == rings.θ && 𝒯.Nϕ == rings.Nϕ
        expected = [SVector(θᵣ, k * 2T(π) / N) for (θᵣ, N) in zip(rings.θ, rings.Nϕ) for k in 0:N-1]
        @test p == expected
        @test length(p) == npixels(𝒯) == (ℓₘₐₓ + 1)^2 - s^2
        @test length(rings.θ) == ℓₘₐₓ - abs(s) + 1 && all(isodd, rings.Nϕ)
        @test issorted(rings.Nϕ)

        # ... and those rings are what `minimal_rings` documents: the interior points of an
        # equally spaced grid on [0, π] (so never a pole); every m has as many rings whose
        # windows include it as there are modes with that m; a ring centered on the side of
        # -s lies in the northern hemisphere and one on the side of +s in the southern.  For
        # s = 0 they are the rings of `sorted_rings` and `sorted_ring_pixels`.
        @test sort(rings.θ) == collect(LinRange{T}(0, T(π), ℓₘₐₓ - abs(s) + 3))[begin+1:end-1]
        for m in -ℓₘₐₓ:ℓₘₐₓ
            covering = count(abs(m - c) ≤ N ÷ 2 for (N, c) in zip(rings.Nϕ, rings.centers))
            @test covering == ℓₘₐₓ - max(abs(m), abs(s)) + 1
        end
        @test all(c * s ≤ 0 || θᵣ > T(π) / 2 for (θᵣ, c) in zip(rings.θ, rings.centers))
        @test all(c * s ≥ 0 || θᵣ < T(π) / 2 for (θᵣ, c) in zip(rings.θ, rings.centers))
        if s == 0
            @test rings.θ == sorted_rings(s, ℓₘₐₓ, T)
            @test rings.Nϕ == [2j + 1 for j in 0:ℓₘₐₓ]
            # These differ only through the rounding of the azimuths (≤ 8 eps, as measured
            # below)
            @test length(p) == length(sorted_ring_pixels(s, ℓₘₐₓ, T))
            @test all(isapprox(a, b; atol=8eps(T)) for (a, b) in zip(p, sorted_ring_pixels(s, ℓₘₐₓ, T)))
            @test all(isapprox(a, b; atol=8eps(T)) for (a, b) in zip(rotors(𝒯), sorted_ring_rotors(s, ℓₘₐₓ, T)))
        end

        # `sorted_rings` itself: the same interior grid, ordered so that each successive
        # ring — which has one more pair of points than the last — lies at least as close to
        # the equator as its predecessor.  Measured: the ordering is violated by at most 1
        # eps(T), among rings that are exactly equally far from the equator and whose order
        # is therefore arbitrary.
        θs = sorted_rings(s, ℓₘₐₓ, T)
        @test sort(θs) == collect(LinRange{T}(0, T(π), ℓₘₐₓ - abs(s) + 3))[begin+1:end-1]
        @test all(0 < θⱼ < T(π) for θⱼ in θs)
        let d = abs.(θs .- T(π) / 2)
            @test all(d[j+1] ≤ d[j] + 4eps(T) for j in 1:length(d)-1)
        end

        # RS: ring-major, ϕₖ = 2πk/Nϕ fastest, on the Fejér-1 rings by default
        𝒯 = SSHT(s, ℓₘₐₓ, T; method="RS")
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

        # Matrix: the Leja points by default, or exactly the rotors given
        R = leja_rotors(s, ℓₘₐₓ, T)
        𝒯 = SSHT(s, ℓₘₐₓ, T; method="Matrix")
        @test rotors(𝒯) == R
        @test eltype(rotors(𝒯)) === Rotor{T}
        @test pixels(𝒯) == to_spherical_coordinates.(R)
        Rshuffled = R[randperm(rng, length(R))]
        𝒯 = SSHT(s, ℓₘₐₓ, T; method="Matrix", Rθϕ=Rshuffled)
        @test rotors(𝒯) == Rshuffled
        @test pixels(𝒯) == to_spherical_coordinates.(Rshuffled)
        # ... and the vector returned is the caller's to modify, as it is for the other
        # methods, which build theirs anew: the transform's own points are not touched
        let Rᵥ = rotors(𝒯)
            Rᵥ[1] = Rotor{T}(1, 0, 0, 0)
            @test rotors(𝒯) == Rshuffled
            @test pixels(𝒯) == to_spherical_coordinates.(Rshuffled)
        end

        # ... and the golden-ratio spiral, the previous default, is the one its docstring
        # describes: N = (ℓₘₐₓ+1)² - s² points, successive azimuths separated by Δϕ =
        # 2π(2-φ) and reduced to [0, 2π), and cos θ at the midpoints of N equal subintervals
        # of [-1, 1] — "uniformly distributed in cos θ", with no point on either pole.
        # Measured against a BigFloat reference, for every (T, ℓₘₐₓ, s) here: ϕ to within 5
        # eps(T) (Double64; half an ulp in the others), and cos θ to within 2 eps(T).
        gp = golden_ratio_spiral_pixels(s, ℓₘₐₓ, T)
        @test golden_ratio_spiral_rotors(s, ℓₘₐₓ, T) == from_spherical_coordinates.(gp)
        let N = (ℓₘₐₓ + 1)^2 - s^2, ϕ(i) = 2big(π) * mod(i * (2 - big(MathConstants.φ)), 1)
            @test length(gp) == N
            @test all(0 ≤ gp[i][2] < 2T(π) for i in 1:N)
            @test maximum(i -> abs(gp[i+1][2] - ϕ(i)), 0:N-1) ≤ 8eps(T)
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
        # successive azimuths separated by Δϕ = 2π(2-φ) and reduced to [0, 2π), compared with
        # a BigFloat reference as in the integer item above, and cos θ at the midpoints of N
        # equal subintervals of [-1, 1], with no point on either pole.  Measured: cos θ to
        # within 1.5 eps(T), over every (T, s, ℓₘₐₓ) here.
        gp = golden_ratio_spiral_pixels(s, ℓₘₐₓ, T)
        let ϕ(i) = 2big(π) * mod(i * (2 - big(MathConstants.φ)), 1)
            @test all(0 ≤ gp[i][2] < 2T(π) for i in 1:N)
            @test maximum(i -> abs(gp[i+1][2] - ϕ(i)), 0:N-1) ≤ 8eps(T)
            @test maximum(abs, [cos(gp[i+1][1]) - (1 - T(2i + 1) / N) for i in 0:N-1]) ≤ 4eps(T)
            @test all(0 < θϕ[1] < T(π) for θϕ in gp)
        end
    end

    # A spin weight larger than ℓₘₐₓ is refused before any pixel is placed, by both families
    # and with one message.
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
        # that each successive ring lies at least as close to the equator as its
        # predecessor.  Measured: the set matches the grid exactly, and the ordering is
        # violated by at most 1.5 eps(T), among rings that are equally far from the equator
        # up to rounding and whose order is settled by the ulp shift, s rounded away from
        # zero.
        @test sort(θs) == collect(LinRange{T}(0, T(π), n + 2))[begin+1:end-1]
        let d = abs.(θs .- T(π) / 2)
            @test all(d[j+1] ≤ d[j] + 4eps(T) for j in 1:length(d)-1)
        end
        # The order for -s mirrors that for s, ring by ring, as it does for integer s, in
        # every T: the side of each ring is decided by the position of its slot in the grid,
        # and not by rounding.  (The equator itself is a slot when the number of rings is
        # odd, and holds the largest ring for either sign.)
        if s > 0
            mirrored = sorted_rings(-s, ℓₘₐₓ, T)
            @test [θ == T(π) / 2 ? 0 : sign(θ - T(π) / 2) for θ in mirrored] ==
                [θ == T(π) / 2 ? 0 : -sign(θ - T(π) / 2) for θ in θs]
        end
    end
    @test sorted_rings(-1//2, 5//2) == Float64(π) .- sorted_rings(1//2, 5//2)

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
    mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
    for f in pixelizations
        @test_throws ArgumentError f(1//2, 3)
        @test_throws mixed f(1//2, 3)
        @test_throws "mixes integers (s) with half-odd-integers (ℓₘₐₓ)" f(0, 7//2)
        @test_throws mixed f(HalfOddInteger(1//2), 3, Float32)
        @test_throws mixed f(2, HalfOddInteger(7//2), Float32)
        # A `Rational` with any other denominator is not an index at all, and one with
        # denominator 1 is to be written as an integer.
        @test_throws ArgumentError f(1//3, 7//2)
        @test_throws "1//3 is neither an integer nor a half-odd-integer" f(1//3, 7//2)
        @test_throws "2//1 is a whole number; write it as the integer 2" f(2//1, 6//1)
    end
end

@testitem "SSHT pixelizations: pinned integer values and refused integer types" begin
    import SphericalFunctions: sorted_rings, sorted_ring_pixels, sorted_ring_rotors
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
    using StaticArrays: SVector

    # These literals pin the points, so that a change of them has to be deliberate; for s =
    # 0, of each pair of mirror-image slots the larger ring is in the north.
    @test sorted_rings(2, 6) == [
        2.6179938779914944, 0.5235987755982988, 2.0943951023931953, 1.0471975511965976,
        1.5707963267948966
    ]
    @test sorted_rings(-1, 3) == [0.7853981633974483, 2.356194490192345, 1.5707963267948966]
    @test sorted_rings(0, 3, Float32) == Float32[2.5132742, 0.62831855, 1.8849556, 1.2566371]
    @test eltype(sorted_rings(2, 6)) === Float64

    g = golden_ratio_spiral_pixels(2, 6)
    @test length(g) == 45
    @test g[1] == SVector(0.2112108803591783, 0.0)
    @test g[23] == SVector(1.5707963267948966, 2.5337085965936814)
    @test g[end] == SVector(2.9303817732306148, 5.067417193187363)

    # ... and every pixel agrees with the formulas of the docstrings, N = (ℓₘₐₓ+1)² - s²
    # points with cos θ at the midpoints of N equal subintervals of [-1, 1] and the azimuth
    # 2π frac(k(2-φ)).  The azimuth is correctly rounded, and so is the colatitude in
    # `Float32`, which is computed in `Float64`; in `Float64` the colatitude is within one
    # ulp.
    for T in (Float64, Float32), (s, ℓₘₐₓ) in ((2, 6), (-1, 3), (0, 5), (3, 3))
        N = (ℓₘₐₓ+1)^2 - s^2
        θ = [T(acos(1 - big(2k-1)/N)) for k in 1:N]
        ϕ = [T(2big(π) * mod(k * (2 - big(MathConstants.φ)), 1)) for k in 0:N-1]
        pixels = golden_ratio_spiral_pixels(s, ℓₘₐₓ, T)
        @test last.(pixels) == ϕ
        @test all(abs.(first.(pixels) .- θ) .≤ (T === Float64 ? eps.(θ) : 0))
        # ... as do the rings: slot i of n, counted from the north pole, is |2i - (n+1)| half
        # spacings from the equator, and the rings fill the slots from the farthest in, the
        # larger of each mirror-image pair to the north for s ≥ 0 and to the south for s <
        # 0.  (This expression is the implementation's own, so it pins the current ordering
        # bit for bit; the independent checks of the documented properties are in
        # `test/utilities/pixelizations.jl`.)
        expected_rings = let n = ℓₘₐₓ - abs(s) + 1
            slots = collect(LinRange{T}(0, π, n + 2))[begin+1:end-1]
            slots[sort(1:n, by=i -> (-abs(2i - (n + 1)), s ≥ 0 ? -i : i))]
        end
        @test sorted_rings(s, ℓₘₐₓ, T) == expected_rings
    end

    # For |s| > ℓₘₐₓ, both families refuse the call with one clear message.  No input with
    # |s| ≤ ℓₘₐₓ is affected; ℓₘₐₓ = |s| gives one ring and 2|s|+1 pixels.
    for f in (golden_ratio_spiral_pixels, golden_ratio_spiral_rotors, sorted_rings, sorted_ring_pixels, sorted_ring_rotors)
        @test_throws "|s|=3 exceeds ℓₘₐₓ=2; there are no such modes." f(3, 2)
        @test_throws "|s|=4 exceeds ℓₘₐₓ=2; there are no such modes." f(-4, 2, Float32)
        @test length(f(2, 2)) == (f === sorted_rings ? 1 : 5)

        # Integer indices of another type than `Int` are refused, with the reason
        @test_throws ArgumentError f(Int8(2), 6)
        @test_throws "`Int8` is narrower than `Int`" f(Int8(2), 6)
        @test_throws "`Int16` is narrower than `Int`" f(Int8(-1), Int16(3))
        @test_throws "is unsigned" f(0, UInt(3))
    end
end

@testitem "SSHT synthesis" setup=[Utilities] begin
    import SphericalFunctions: SSHT, pixels, rotors, Ysize, Yindex, sYlm_matrix, ModeWeights
    import SphericalFunctions: npixels  # unexported
    import .Utilities: sYlm_closed_form, sYlm_closed_form_pixels
    using DoubleFloats: Double64
    using LinearAlgebra: norm
    using StaticArrays: SVector
    using Random

    # `sYlm_closed_form_pixels` — the `Utilities` closed form with the pixel-independent
    # factorials hoisted out of the loop, about 40 times faster — really is that closed
    # form.  This is the one place the shared helper is checked; the other items that use it
    # say so.  The two are compared over the whole ℓ range the file uses (not just the low ℓ
    # where cancellation in the alternating sum is mild), at four hand-picked points *and*
    # on a real pixelization, since they accumulate that sum in different precisions
    # (BigFloat in `sYlm_closed_form`, `T` in `sYlm_closed_form_pixels`).  Measured: ≤ 1.05
    # eps(T) at the four points and ≤ 8.3 eps(T) on the pixelization, so 100 eps(T) leaves a
    # factor of ≳ 12.
    for T in (Float64, Double64, Float32)
        q = [SVector{2, T}(θ, ϕ) for (θ, ϕ) in ((0.3, 1.1), (2.7, 5.2), (0.0, 0.0), (3.0, 2.0))]
        for s in -2:2
            qs = vcat(q, pixels(SSHT(s, 6, T; method="RS")))
            for ℓ in abs(s):6, m in -ℓ:ℓ
                @test maximum(
                    abs,
                    sYlm_closed_form_pixels(s, ℓ, m, qs)
                    .- [sYlm_closed_form(s, ℓ, m, θϕ[1], θϕ[2]) for θϕ in qs]
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
                𝒯 = SSHT(s, ℓₘₐₓ, T; method, kw...)
                n = Ysize(abs(s), ℓₘₐₓ)
                p = pixels(𝒯)
                # Column i of `Y` is the closed-form ₛYₗₘ of mode i sampled on the pixels
                Y = reduce(
                    hcat, [sYlm_closed_form_pixels(s, ℓ, m, p) for ℓ in abs(s):ℓₘₐₓ for m in -ℓ:ℓ]
                )

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
    import SphericalFunctions: salm2map, map2salm_plan
    import SphericalFunctions: nmodes, npixels  # unexported
    import .Utilities: sYlm_closed_form_pixels
    using DoubleFloats: Double64
    using LinearAlgebra: mul!, ldiv!
    using StaticArrays: SVector
    using Random

    # `sYlm_closed_form_pixels` is the closed-form ₛYₗₘ of the `Utilities` module on a list
    # of pixels, with the pixel-independent factorials hoisted out of the loop; see "SSHT
    # synthesis", where it is checked against `sYlm_closed_form` itself.

    # Analysis tolerance.  Measured over every T, s ∈ -2:2 and single mode, in units of
    # eps(T): "RS" ≤ 8, "Matrix" ≤ 18, "Minimal" ≤ 50, and round trips of random weights ≤
    # 37, all growing only slowly with ℓₘₐₓ at these sizes, so 100ℓₘₐₓ eps(T) leaves a
    # factor of ≳ 6.  (On the rings of `sorted_rings` rather than `minimal_rings`, "Minimal"
    # would be a different story: the error would grow by a factor of about 20 per unit
    # ℓₘₐₓ, to 1.6e5 eps(T) at ℓₘₐₓ = 6, and Float32 could not be tested beyond ℓₘₐₓ = 5.)
    tolerance(method, ℓₘₐₓ, ::Type{T}) where {T} = 100ℓₘₐₓ * eps(T)

    rng = Random.Xoshiro(2718)

    for (method, T) in Iterators.product(("RS", "Minimal", "Matrix"), (Float64, Double64, Float32))
        kw = method == "RS" ? (;) : (; inplace=false)
        for ℓₘₐₓ in 3:6
            ϵ = tolerance(method, ℓₘₐₓ, T)
            for s in -2:2
                𝒯 = SSHT(s, ℓₘₐₓ, T; method, kw...)
                n = nmodes(𝒯)
                p = pixels(𝒯)

                # Single modes evaluated from the closed form come back as unit vectors, as
                # ModeWeights
                for ℓ in abs(s):ℓₘₐₓ, m in -ℓ:ℓ
                    f = sYlm_closed_form_pixels(s, ℓ, m, p)
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
                @test all(isapprox(F[:, j], 𝒯 * F̃[:, j]; atol=ϵ, rtol=ϵ) for j in 1:3)
                F̃′ = 𝒯 \ F
                @test F̃′ isa Matrix{Complex{T}}
                @test size(F̃′) == (n, 3)
                @test F̃′ ≈ F̃ atol=ϵ rtol=ϵ
                @test all(isapprox(F̃′[:, j], 𝒯 \ F[:, j]; atol=ϵ, rtol=ϵ) for j in 1:3)
            end
        end
    end

    # The three methods are very different algorithms — "RS" takes an FFT along each ring
    # and applies a quadrature rule across rings, "Minimal" solves a sequence of small
    # systems for groups of m values, "Matrix" factors one dense matrix — and they share no
    # code path beyond the harmonics themselves.  Handing the "RS" and "Minimal" sample
    # points to "Matrix" puts all three on a common grid, where they must agree.  (This is a
    # cross-check inside the package, not an independent reference; the closed-form
    # comparisons above are that.)
    for T in (Float64, Double64, Float32), ℓₘₐₓ in (3, 5), s in (-2, 0, 1)
        n = Ysize(abs(s), ℓₘₐₓ)
        f̃ = randn(rng, Complex{T}, n)
        for method in ("RS", "Minimal")
            kw = method == "RS" ? (;) : (; inplace=false)
            𝒯 = SSHT(s, ℓₘₐₓ, T; method, kw...)
            𝒯ᴹ = SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=rotors(𝒯), inplace=false)
            @test npixels(𝒯ᴹ) == npixels(𝒯)
            # Measured over exactly these cases, in units of eps(T): synthesis ≤ 43, and
            # analysis ≤ 20 against "RS" and ≤ 30 against "Minimal".
            ϵ = tolerance(method, ℓₘₐₓ, T)
            f = 𝒯 * f̃
            @test f ≈ 𝒯ᴹ * f̃ atol=ϵ rtol=ϵ
            @test 𝒯 \ copy(f) ≈ 𝒯ᴹ \ copy(f) atol=ϵ rtol=ϵ
        end
    end

    # Inputs of the wrong length are errors for every method, in `*`, `\`, `mul!` and
    # `ldiv!`, and so are ModeWeights with the wrong labels.  ("RS" has no `inplace`
    # option.)
    for method in ("RS", "Minimal", "Matrix"), inplace in (method == "RS" ? (false,) : (true, false))
        s, ℓₘₐₓ = -1, 4
        kw = method == "RS" ? (;) : (; inplace)
        𝒯 = SSHT(s, ℓₘₐₓ; method, kw...)
        n, N = nmodes(𝒯), npixels(𝒯)
        @test_throws DimensionMismatch 𝒯 * zeros(ComplexF64, n + 1)
        @test_throws "first dimension of the mode weights has length $(n + 1)" 𝒯 * zeros(ComplexF64, n + 1)
        @test_throws DimensionMismatch 𝒯 * zeros(ComplexF64, n - 1)
        @test_throws DimensionMismatch 𝒯 * zeros(ComplexF64, n + 1, 2)
        if N != n  # for "Minimal" and "Matrix" the two lengths agree, so this is legal
            @test_throws DimensionMismatch 𝒯 * zeros(ComplexF64, N)
        end
        @test_throws DimensionMismatch 𝒯 \ zeros(ComplexF64, N + 1)
        @test_throws "first dimension of the function values" 𝒯 \ zeros(ComplexF64, N + 1)
        @test_throws DimensionMismatch 𝒯 \ zeros(ComplexF64, N - 1)
        @test_throws DimensionMismatch 𝒯 \ zeros(ComplexF64, N - 1, 2)
        @test_throws DimensionMismatch mul!(zeros(ComplexF64, N + 1), 𝒯, zeros(ComplexF64, n))
        @test_throws DimensionMismatch mul!(zeros(ComplexF64, N), 𝒯, zeros(ComplexF64, n + 1))
        # ModeWeights with ℓₘᵢₙ = 0 (length Ysize(0, ℓₘₐₓ)) or another ℓₘₐₓ.  Synthesis
        # takes any range up to the transform's ℓₘₐₓ (see "SSHT synthesis of mode weights of
        # another range of ℓ"), but not beyond it; `mul!` and the output of analysis need
        # the transform's range exactly.  Each refusal names the constructor that re-ranges
        # them.
        rerange = "`ModeWeights(w; ℓₘᵢₙ=1, ℓₘₐₓ=4)` copies mode weights"
        w₀ = ModeWeights(zeros(ComplexF64, Ysize(0, ℓₘₐₓ)), s; ℓₘᵢₙ=0)
        @test 𝒯 * w₀ == zeros(ComplexF64, N)
        @test_throws ArgumentError 𝒯 * ModeWeights(zeros(ComplexF64, Ysize(1, ℓₘₐₓ + 1)), s)
        @test_throws "ModeWeights have ℓ ∈ 1:5, but the transform synthesizes ℓ ≤ 4" 𝒯 * ModeWeights(zeros(ComplexF64, Ysize(1, ℓₘₐₓ + 1)), s)
        @test_throws rerange 𝒯 * ModeWeights(zeros(ComplexF64, Ysize(1, ℓₘₐₓ + 1)), s)
        @test_throws ArgumentError mul!(zeros(ComplexF64, N), 𝒯, w₀)
        @test_throws "ModeWeights have ℓ ∈ 0:4, but the transform requires ℓ ∈ 1:4" mul!(zeros(ComplexF64, N), 𝒯, w₀)
        @test_throws rerange mul!(zeros(ComplexF64, N), 𝒯, ModeWeights(zeros(ComplexF64, Ysize(1, ℓₘₐₓ - 1)), s))
        # (An output vector *longer* than the modes is accepted, and labelled over its first n
        # entries; a shorter one is not.)
        let out = zeros(ComplexF64, n + 1)
            w = ldiv!(out, 𝒯, zeros(ComplexF64, N))
            @test w isa ModeWeights && parent(array_view(w)) === out && length(array_view(w)) == n
        end
        @test_throws "at least" ldiv!(zeros(ComplexF64, n - 1), 𝒯, zeros(ComplexF64, N))
        @test_throws DimensionMismatch ldiv!(zeros(ComplexF64, n), 𝒯, zeros(ComplexF64, N + 1))
        @test_throws ArgumentError ldiv!(w₀, 𝒯, zeros(ComplexF64, N))
        @test_throws rerange ldiv!(w₀, 𝒯, zeros(ComplexF64, N))
        # ModeWeights of the wrong spin weight but the right length: -s, or 0 with ℓₘᵢₙ = |s|.
        # Both synthesis and the output of `ldiv!` refuse them, rather than treating them as
        # (or filling them with) weights of spin s.
        for w in (
            ModeWeights(ones(ComplexF64, n), -s),
            ModeWeights(ones(ComplexF64, n), 0; ℓₘᵢₙ=abs(s)),
        )
            @test_throws ArgumentError 𝒯 * w
            @test_throws "ModeWeights have spin weight s=$(spin(w)), but the transform is for s=$s" 𝒯 * w
            @test_throws "ModeWeights have spin weight" mul!(zeros(ComplexF64, N), 𝒯, w)
            @test_throws "ModeWeights have spin weight" ldiv!(w, 𝒯, zeros(ComplexF64, N))
            if method == "RS"  # salm2map needs its own Clenshaw–Curtis plan
                plan = map2salm_plan(zeros(ComplexF64, 2ℓₘₐₓ+1, 2ℓₘₐₓ+1), s, ℓₘₐₓ)
                @test_throws "ModeWeights have spin weight" salm2map(w, plan)
            end
        end
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
        # Measured over these cases: round trips ≤ 12 eps(T), and the products of several
        # columns at once within 4 eps(T) of the single products.  This leaves a factor of ≳ 25.
        ϵ = 100ℓₘₐₓ * eps(T)
        f̃0 = randn(rng, Complex{T}, n)

        for method in ("Minimal", "Matrix")
            𝒯 = SSHT(s, ℓₘₐₓ, T; method)  # inplace=true is the default
            𝒯n = SSHT(s, ℓₘₐₓ, T; method, inplace=false)
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
            # `*` is a plain matrix-vector product, so it allocates even for the in-place type;
            # only its `\` acts in place.
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

            # In-place analysis: the input array is overwritten with the mode weights, which
            # come back as a ModeWeights wrapping that same storage — so that they are indexed
            # by (ℓ, m), as from every other analysis.  (Returned as the bare Vector, they
            # would let `(𝒯 \ f)[ℓ, m]` silently read the linear index instead.)
            g = copy(fref)
            g̃ = 𝒯 \ g
            @test g̃ isa ModeWeights{Complex{T}}
            @test parent(g̃) === g
            @test spin(g̃) == s
            @test g == parent(f̃ref)
            @test all(g̃[ℓ, m] == f̃ref[ℓ, m] for ℓ in abs(s):ℓₘₐₓ for m in -ℓ:ℓ)

            # A ModeWeights input is overwritten the same way.  What comes back from synthesis
            # is its storage, holding function values, as a plain Vector — not the ModeWeights,
            # whose labels no longer describe it — and from analysis a ModeWeights wrapping
            # the storage.
            w = ModeWeights(copy(f̃0), s)
            fw = 𝒯 * w
            if method == "Minimal"
                @test fw === parent(w)
                @test fw isa Vector{Complex{T}}
                @test parent(w) == fref
            else
                @test fw == fref
                @test parent(w) == f̃0
            end
            wf = ModeWeights(copy(fref), s)  # function values can be stored in one, as N == n
            w̃ = 𝒯 \ wf
            @test w̃ isa ModeWeights{Complex{T}}
            @test parent(w̃) === parent(wf)
            @test parent(wf) == parent(f̃ref)

            # Two-argument `ldiv!(𝒯, x)` (and `mul!(𝒯, x)` for Minimal) act in place on the
            # non-in-place type as well — the first returning the mode weights as a
            # ModeWeights over the argument's storage, and the second the function values as
            # that storage, a plain array, even when it came in a ModeWeights
            g = copy(fref)
            g̃ = ldiv!(𝒯n, g)
            @test g̃ isa ModeWeights{Complex{T}} && parent(g̃) === g && spin(g̃) == s
            @test g == parent(f̃ref)
            if method == "Minimal"
                g̃ = copy(f̃0)
                @test mul!(𝒯n, g̃) === g̃
                @test g̃ == fref
                wg = ModeWeights(copy(f̃0), s)
                @test mul!(𝒯n, wg) === parent(wg)
                @test parent(wg) == fref
            else
                @test_throws MethodError mul!(𝒯n, copy(f̃0))
            end

            # Acting in place needs storage that can hold the complex results: real storage
            # is refused, naming the alternative, and the transform that is not in place
            # accepts real data
            r = real.(fref)
            @test_throws ArgumentError 𝒯 \ copy(r)
            @test_throws "use `𝒯 \\ f` with a transform constructed with `inplace=false`" 𝒯 \ copy(r)
            @test_throws "whose element type is $T" ldiv!(𝒯n, copy(r))
            @test array_view(𝒯n \ r) ≈ array_view(𝒯n \ complex.(r)) atol=ϵ rtol=ϵ
            if method == "Minimal"
                @test_throws "use `𝒯 * f̃` with a transform constructed with `inplace=false`" 𝒯 * real.(f̃0)
                @test_throws ArgumentError mul!(𝒯n, real.(f̃0))
            end
            @test 𝒯n * real.(f̃0) ≈ 𝒯n * complex.(real.(f̃0)) atol=ϵ rtol=ϵ

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
            objects = method == "RS" ? (SSHT(s, ℓₘₐₓ, T; method),) :
                (SSHT(s, ℓₘₐₓ, T; method), SSHT(s, ℓₘₐₓ, T; method, inplace=false))
            for 𝒯 in objects
                N = npixels(𝒯)
                fref = SSHT(s, ℓₘₐₓ, T; method, (method == "RS" ? (;) : (; inplace=false))...) * f̃0
                f̃ = copy(f̃0)
                f = zeros(Complex{T}, N)
                @test mul!(f, 𝒯, f̃) === f
                @test f̃ == f̃0
                @test f == fref
                g = copy(fref)
                g̃ = zeros(Complex{T}, n)
                w̃g = ldiv!(g̃, 𝒯, g)  # a bare output comes back labelled, over its storage
                @test w̃g isa ModeWeights{Complex{T}} && parent(array_view(w̃g)) === g̃
                @test spin(w̃g) == s
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
        # Measured: round trips ≤ 10 eps, and the columns of a product of several at once
        # within 4 eps of the single products (a matrix-matrix product need not round like a
        # matrix-vector product).  The three-argument forms and the in-place types run the
        # same kernels as `*` and `\` on the same data, so they agree exactly.
        ϵ = 100ℓₘₐₓ * eps(Float64)
        kw = method == "RS" ? (;) : (; inplace=false)
        𝒯 = SSHT(s, ℓₘₐₓ; method, kw...)
        n, N = nmodes(𝒯), npixels(𝒯)
        F̃ = randn(rng, ComplexF64, n, 3, 2)
        columns = [𝒯 * F̃[:, j, k] for j in 1:3, k in 1:2]

        F = 𝒯 * F̃
        @test size(F) == (N, 3, 2)
        @test all(isapprox(F[:, j, k], columns[j, k]; atol=ϵ, rtol=ϵ) for j in 1:3, k in 1:2)
        F̃′ = 𝒯 \ F
        @test size(F̃′) == (n, 3, 2)
        @test F̃′ ≈ F̃ atol=ϵ rtol=ϵ

        G = similar(F)
        @test mul!(G, 𝒯, F̃) === G
        @test G == F
        G̃ = similar(F̃)
        @test ldiv!(G̃, 𝒯, F) === G̃
        @test G̃ == F̃′

        if method != "RS"  # the in-place types too
            𝒯i = SSHT(s, ℓₘₐₓ; method)
            H = 𝒯i * copy(F̃)
            @test size(H) == (N, 3, 2)
            @test H == F
            H̃ = 𝒯i \ copy(F)
            @test size(H̃) == (n, 3, 2)
            @test H̃ == F̃′
        end
    end
end


@testitem "SSHTRS options" setup=[Utilities] begin
    import SphericalFunctions: SSHT, SSHTRS, pixels, rotors, Ysize, Yindex, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    import SphericalFunctions: fejer1_rings, fejer2_rings, clenshaw_curtis_rings
    import SphericalFunctions: fejer1, fejer2, clenshaw_curtis
    import .Utilities: sYlm_closed_form_pixels
    using DoubleFloats: Double64
    using FFTW: FFTW
    using Logging: NullLogger, with_logger
    using Quaternionic: from_spherical_coordinates
    using StaticArrays: SVector
    using Random

    # `sYlm_closed_form_pixels` is the closed-form ₛYₗₘ of the `Utilities` module on a list of
    # pixels, with the pixel-independent factorials hoisted out of the loop; see "SSHT
    # synthesis", where it is checked against `sYlm_closed_form` itself.

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
                𝒯 = SSHT(s, ℓₘₐₓ, T; method="RS", θ, quadrature_weights=w, Nϕ)
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
                        @test g ≈ sYlm_closed_form_pixels(s, ℓ, m, p) atol=ϵ rtol=ϵ
                        @test 𝒯 \ g ≈ g̃ atol=ϵ rtol=ϵ
                    end
                end
            end
        end
    end

    # Rings and weights must be given in T; Float64 ones are refused for T=Float32 rather
    # than rounded, and integer ones, which convert exactly, are accepted for every T, while
    # rational ones, which do not, are refused
    @test_throws ArgumentError SSHTRS(
        1, 5, Float32; θ=fejer1_rings(11), quadrature_weights=fejer1(11, Float32)
    )
    @test_throws "must be a vector of Float32" SSHTRS(
        1, 5, Float32; θ=fejer1_rings(11), quadrature_weights=fejer1(11, Float32)
    )
    @test_throws "`quadrature_weights` must be a vector of Float32" SSHTRS(
        1, 5, Float32; θ=fejer1_rings(11, Float32), quadrature_weights=fejer1(11)
    )
    𝒯 = SSHTRS(1, 5, Float32; θ=fejer1_rings(11, Float32), quadrature_weights=fejer1(11, Float32))
    @test 𝒯 isa SSHTRS{Float32}
    @test eltype(pixels(𝒯)) === SVector{2, Float32}
    @test 𝒯.θ == fejer1_rings(11, Float32)
    @test SSHTRS(0, 0; θ=[1], quadrature_weights=[2]).θ == [1.0]
    for T in (Float32, Double64, BigFloat)
        local 𝒯 = SSHTRS(0, 0, T; θ=[1], quadrature_weights=[2])
        @test 𝒯.θ == [1] && 𝒯.θ isa Vector{T}
        @test_throws "`θ` must be a vector of $T" SSHTRS(0, 0, T; θ=[1//2], quadrature_weights=[2])
    end
    @test_throws "`θ` must be a vector of Float64" SSHTRS(0, 0; θ=[1//2], quadrature_weights=[2])

    # A single ring is still a batch of rings, as any vector of rotor data is
    𝒯 = SSHTRS(0, 0)
    @test length(𝒯.θ) == 1
    @test 𝒯 \ (𝒯 * [0.3 + 0.4im]) ≈ [0.3 + 0.4im] rtol=4eps()

    # FFTW planner options are accepted and do not change the results
    𝒯 = SSHTRS(1, 5)
    𝒯m = SSHTRS(1, 5; plan_fft_flags=FFTW.MEASURE, plan_fft_timelimit=1.0)
    f̃ = randn(rng, ComplexF64, nmodes(𝒯))
    @test 𝒯m * f̃ ≈ 𝒯 * f̃ rtol=100eps()
    @test 𝒯m \ (𝒯m * f̃) ≈ f̃ rtol=100eps()  # measured: ≤ 3.5 eps over 20 draws

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
    # |m| < ℓₘₐₓ is fine (measured: ≤ 9.2 eps over 20 draws)
    @test 𝒯 \ (𝒯 * f̃) ≈ f̃ atol=500eps() rtol=500eps()
    f̃[Yindex(ℓₘₐₓ, ℓₘₐₓ, abs(s))] = 1
    @test !isapprox(𝒯 \ (𝒯 * f̃), f̃; atol=1e-3)  # m = ±ℓₘₐₓ alias onto each other

    # Inconsistent options are errors: lengths that disagree, and rings of no points
    @test_throws DimensionMismatch SSHTRS(s, ℓₘₐₓ; θ=fejer1_rings(11), quadrature_weights=fejer1(12))
    @test_throws "same length" SSHTRS(s, ℓₘₐₓ; θ=fejer1_rings(11), quadrature_weights=fejer1(12))
    @test_throws "same length" SSHTRS(s, ℓₘₐₓ; θ=fejer1_rings(12), quadrature_weights=fejer1(11))
    @test_throws "Nϕ must be a single number or have the same length" SSHTRS(s, ℓₘₐₓ; Nϕ=fill(11, 10))
    @test_throws DimensionMismatch SSHTRS(s, ℓₘₐₓ; Nϕ=fill(11, 12))
    @test_throws ArgumentError SSHTRS(s, ℓₘₐₓ; Nϕ=0)
    @test_throws "at least one point" SSHTRS(s, ℓₘₐₓ; Nϕ=0)
    @test_throws "at least one point" SSHTRS(s, ℓₘₐₓ; Nϕ=[11, 11, 11, 11, 11, 0, 11, 11, 11, 11, 11])
end


@testitem "SSHTMatrix options" begin
    import SphericalFunctions: SSHT, SSHTMatrix, pixels, rotors, Ysize, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    import SphericalFunctions: golden_ratio_spiral_rotors, sorted_ring_rotors, array_view
    import SphericalFunctions: sYlm_matrix
    using DoubleFloats: Double64
    using LinearAlgebra: LinearAlgebra, lu, qr, ldiv!, norm
    using Random

    rng = Random.Xoshiro(8128)

    for T in (Float64, Double64), (s, ℓₘₐₓ) in ((-2, 5), (1, 4))
        # Measured over these cases: round trips ≤ 18 eps(T) on the default points and with
        # QR, and least-squares solutions ≤ 6 eps(T).  This leaves a factor of ≳ 20.  On the
        # sorted rings, which are badly conditioned for s ≠ 0, a round trip measured up to
        # 7200 eps(T) (s = -2, ℓₘₐₓ = 5).
        ϵ = 100ℓₘₐₓ * eps(T)
        ϵ_rings = 3e4 * eps(T)
        n = Ysize(abs(s), ℓₘₐₓ)
        f̃ = randn(rng, Complex{T}, n)

        # Any set of rotors may be used; with exactly n of them the default is LU and in place
        R = sorted_ring_rotors(s, ℓₘₐₓ, T)
        𝒯 = SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=R)
        @test 𝒯 isa SSHTMatrix{T, true, <:LinearAlgebra.LU}
        @test rotors(𝒯) == R
        @test npixels(𝒯) == n
        f = 𝒯 * f̃
        @test 𝒯 \ copy(f) ≈ f̃ atol=ϵ_rings rtol=ϵ_rings

        # An explicit decomposition, e.g. QR, must also solve the problem
        𝒯qr = SSHTMatrix(s, ℓₘₐₓ, T; decomposition=qr, inplace=false)
        @test 𝒯qr isa SSHTMatrix{T, false}
        @test !(𝒯qr isa SSHTMatrix{T, false, <:LinearAlgebra.LU})
        @test 𝒯qr * f̃ == SSHTMatrix(s, ℓₘₐₓ, T; decomposition=lu, inplace=false) * f̃
        @test 𝒯qr \ (𝒯qr * f̃) ≈ f̃ atol=ϵ rtol=ϵ
        # ... but need only solve it with `\`: a decomposition with no in-place `ldiv!`, such
        # as the matrix itself, serves a transform that does not act in place, and the
        # constructor's measurement of the round trip does not need one either
        𝒯id = SSHTMatrix(s, ℓₘₐₓ, T; decomposition=identity, inplace=false)
        @test 𝒯id isa SSHTMatrix{T, false, Matrix{Complex{T}}}
        @test 𝒯id \ (𝒯id * f̃) ≈ f̃ atol=ϵ rtol=ϵ

        # More points than modes: a least-squares analysis, QR by default, never in place
        Rmore = [golden_ratio_spiral_rotors(s, ℓₘₐₓ, T); sorted_ring_rotors(0, ℓₘₐₓ + 1, T)]
        𝒯ls = SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=Rmore)
        @test 𝒯ls isa SSHTMatrix{T, false}
        @test !(𝒯ls isa SSHTMatrix{T, false, <:LinearAlgebra.LU})
        @test npixels(𝒯ls) == length(Rmore) > nmodes(𝒯ls) == n
        fls = 𝒯ls * f̃
        @test length(fls) == length(Rmore)
        f̃′ = 𝒯ls \ fls
        @test f̃′ isa ModeWeights{Complex{T}}
        @test f̃′ ≈ f̃ atol=ϵ rtol=ϵ
        # The two-argument `ldiv!` solves in place, leaving the solution in the first n entries,
        # and labels just those — rather than returning the whole array, residual and all
        g = copy(fls)
        wg = ldiv!(𝒯ls, g)
        @test wg isa ModeWeights{Complex{T}}
        @test parent(array_view(wg)) === g && length(array_view(wg)) == n
        @test wg ≈ f̃ atol=ϵ rtol=ϵ
        # Data that no set of modes fits exactly exercises the least-squares property itself:
        # the residual is orthogonal to every column of the harmonic matrix (the normal
        # equations), which an interpolation through any n of the points would not achieve
        Y = sYlm_matrix(rotors(𝒯ls), ℓₘₐₓ, s)
        fnoisy = fls + randn(rng, Complex{T}, length(fls))
        f̃ls = array_view(𝒯ls \ fnoisy)
        r = fnoisy - Y * f̃ls
        @test norm(r) > norm(fnoisy) / 100  # the data really do not fit
        @test norm(Y' * r) ≤ ϵ * norm(Y) * norm(fnoisy)
        f̃interp = Y[1:n, :] \ fnoisy[1:n]
        @test norm(Y' * (fnoisy - Y * f̃interp)) > 100ϵ * norm(Y) * norm(fnoisy)
        @test array_view(ldiv!(𝒯ls, copy(fnoisy))) ≈ f̃ls atol=ϵ rtol=ϵ
        # ... which is the use case of the docstring: sample on the grid appropriate to the
        # *lowest* |s| (the most modes, here s = 0) and reuse those points for this spin
        if s != 0
            𝒯s = SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=golden_ratio_spiral_rotors(0, ℓₘₐₓ, T))
            @test 𝒯s isa SSHTMatrix{T, false}
            @test npixels(𝒯s) == Ysize(0, ℓₘₐₓ) > nmodes(𝒯s) == n
            g̃ = randn(rng, Complex{T}, nmodes(𝒯s))
            @test 𝒯s \ (𝒯s * g̃) ≈ g̃ atol=ϵ rtol=ϵ
        end

        # Fewer points than modes is underdetermined, and in place needs exactly n points
        @test_throws DimensionMismatch SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=R[1:end-1])
        @test_throws "underdetermined" SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=R[1:end-1])
        @test_throws ArgumentError SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=Rmore, inplace=true)
        @test_throws "exactly as many sample points as modes" SSHTMatrix(s, ℓₘₐₓ, T; Rθϕ=Rmore, inplace=true)
    end
end


@testitem "SSHT sample data must be in the transform's type" begin
    import SphericalFunctions: SSHT, SSHTMatrix, SSHTRS, SSHTMinimal, rotors, leja_rotors,
        fejer1_rings, fejer1, minimal_rings
    using Quaternionic: Rotor, Quaternion, QuatVec

    # A `QuatVec` is not the rotor of its direction: converted to a rotor, it would become
    # 𝐢, whose pixel is the south pole, so it is refused.  A `Quaternion` is the rotation
    # of its normalization, which is the rotor that the transform keeps.
    @test_throws ArgumentError SSHTMatrix(0, 0; Rθϕ=[QuatVec(0.0, 1.0, 0.0, 0.0)])
    @test_throws "Rotations are taken as" SSHTMatrix(0, 0; Rθϕ=[QuatVec(0.0, 1.0, 0.0, 0.0)])
    @test rotors(SSHTMatrix(0, 0; Rθϕ=[2.0 * Quaternion(1.0, 0, 0, 0)])) == [Rotor(1.0)]
    @test_throws "Rotations are taken as" SSHT(0, 0; method="Matrix", Rθϕ=[QuatVec(0.0, 1.0, 0.0, 0.0)])

    # Data of another precision is refused rather than rounded or widened, for each method
    # and through the `SSHT` front end; asking for its type works
    R = leja_rotors(0, 2, BigFloat)
    @test_throws ArgumentError SSHTMatrix(0, 2; Rθϕ=R)
    @test_throws "must be a vector of `Rotor{Float64}`s" SSHTMatrix(0, 2; Rθϕ=R)
    @test_throws "must be a vector of `Rotor{Float64}`s" SSHT(0, 2; method="Matrix", Rθϕ=R)
    @test_throws "must be a vector of `Rotor{BigFloat}`s" SSHTMatrix(0, 2, BigFloat; Rθϕ=leja_rotors(0, 2))
    @test rotors(SSHTMatrix(0, 2, BigFloat; Rθϕ=R)) == R
    @test_throws "`θ` must be a vector of Float64" SSHTRS(
        0, 2; θ=fejer1_rings(5, BigFloat), quadrature_weights=fejer1(5)
    )
    @test_throws ArgumentError SSHTMinimal(1, 3; θ=minimal_rings(1, 3, BigFloat).θ)
    @test_throws "`θ` must be a vector of Float64" SSHTMinimal(1, 3; θ=minimal_rings(1, 3, BigFloat).θ)
    @test_throws "`θ` must be a vector of Float64" SSHT(1, 3; method="Minimal", θ=minimal_rings(1, 3, BigFloat).θ)
    @test_throws "`θ` must be a vector of BigFloat" SSHTMinimal(1, 3, BigFloat; θ=minimal_rings(1, 3).θ)
    @test SSHTMinimal(1, 3, BigFloat; θ=minimal_rings(1, 3, BigFloat).θ).θ == minimal_rings(1, 3, BigFloat).θ

    # A vector whose element type does not fix a precision is refused too
    @test_throws "must be a vector of `Rotor{Float64}`s" SSHTMatrix(0, 0; Rθϕ=Rotor[Rotor(1.0, 0, 0, 0)])
    @test_throws "`θ` must be a vector of Float64" SSHTMinimal(0, 0; θ=Real[0.5])
    # ... while integer colatitudes, which convert exactly, are accepted in every type
    @test SSHTMinimal(0, 0, BigFloat; θ=[1]).θ == [1]
    @test length(SSHTMinimal(0, 1, Float32; θ=[1, 2]).θ) == 2
    @test_throws DimensionMismatch SSHTMinimal(0, 1; θ=[1, 2, 3])
end

@testitem "SSHT thread safety via separate objects" begin
    import SphericalFunctions: SSHT, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(8675309)

    # Each SSHT holds its own workspace, so one object per task gives results identical to a
    # serial computation (each transform runs on its caller's thread, in a fixed order).
    # With a single thread the tasks run one after another, since a transform never yields,
    # so these checks can fail only when Julia has several threads, as it does in CI; the
    # structural checks of "SSHT copies share their tables and have workspace of their own"
    # guard against shared workspace regardless.
    for method in ("RS", "Minimal", "Matrix"), T in (Float64, Float32)
        s, ℓₘₐₓ = -2, 8
        kw = method == "RS" ? (;) : (; inplace=false)

        # Two objects of the same parameters on two tasks, as in the assignment
        𝒯s = [SSHT(s, ℓₘₐₓ, T; method, kw...) for _ in 1:2]
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
        𝒯s = [SSHT(s, ℓₘₐₓ, T; method, kw...) for _ in 1:8]
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

@testitem "SSHT copies share their tables and have workspace of their own" begin
    import SphericalFunctions: SSHT
    import SphericalFunctions: nmodes  # unexported
    using Random

    rng = Random.Xoshiro(20260926)

    # `copy` of an "RS" or "Minimal" transform shares what no transform modifies — the rings,
    # weights, FFT plans, phase factors, tables and decompositions — and allocates new
    # workspace, so that the copy and the original may be used by two tasks at once.  Its
    # results are those of the original to the last bit.
    shared = Dict(
        "RS" => (
            :θ, :quadrature_weights, :Nϕ, :ring_ranges, :synthesis_phases, :analysis_phases
        ),
        "Minimal" => (:θ, :Nϕ, :centers, :ring_ranges, :mode_m, :Λ, :blocks),
    )
    workspace = Dict("RS" => (:F, :G), "Minimal" => (:F, :f̃, :rhs))
    for (method, s, ℓₘₐₓ, kw) in (
        ("RS", -2, 8, (;)), ("RS", 3//2, 13//2, (;)),
        ("Minimal", 2, 8, (;)), ("Minimal", 0, 6, (; inplace=false)),
    )
        𝒯 = SSHT(s, ℓₘₐₓ; method, kw...)
        𝒯′ = copy(𝒯)
        @test typeof(𝒯′) === typeof(𝒯)
        @test all(getfield(𝒯′, name) === getfield(𝒯, name) for name in shared[method])
        @test 𝒯′.plans.forward === 𝒯.plans.forward && 𝒯′.plans.backward === 𝒯.plans.backward
        for name in workspace[method]
            a, b = getfield(𝒯′, name), getfield(𝒯, name)
            @test a !== b
            if eltype(a) <: AbstractArray  # a vector of buffers, each of them new
                @test all(x !== y for (x, y) in zip(a, b))
            end
        end
        if method == "RS"
            @test 𝒯′.λ !== 𝒯.λ && 𝒯′.λ.Yˡ !== 𝒯.λ.Yˡ
        end
        F̃ = randn(rng, ComplexF64, nmodes(𝒯), 3)
        F = 𝒯 * copy(F̃)
        @test 𝒯′ * copy(F̃) == F
        @test 𝒯′ \ copy(F) == 𝒯 \ copy(F)
    end

    # "Matrix" holds no workspace, so its copy is the object itself
    𝒯 = SSHT(1, 4; method="Matrix")
    @test copy(𝒯) === 𝒯
end


@testitem "SSHT use from several tasks: a copy for each task, and one shared SSHTMatrix" begin
    import SphericalFunctions: SSHT, SSHTMatrix, golden_ratio_spiral_rotors
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: mul!, ldiv!, qr
    using Random

    rng = Random.Xoshiro(20260927)

    # The pattern the `SSHT` docstring recommends: the data divided into chunks, and one task
    # spawned for each chunk, with its own copy of the transform.  Every result must be that of
    # the serial computation, to the last bit.  (Without separate workspace, tasks that share
    # an "RS" or "Minimal" object produce mode weights wrong by O(1).)  With a single thread
    # the tasks run one after another, since a transform never yields, so this check would pass
    # even if the copies shared their workspace; there, the structural checks of "SSHT copies
    # share their tables and have workspace of their own" are what guard against that.
    roundtrip(𝒯, f̃) = collect(𝒯 \ (𝒯 * copy(f̃)))
    function in_tasks(𝒯, inputs, ntasks)
        chunks = Iterators.partition(inputs, cld(length(inputs), ntasks))
        tasks = [
            Threads.@spawn(let 𝒯ₖ = copy(𝒯); [roundtrip(𝒯ₖ, f̃) for f̃ ∈ chunk]; end)
            for chunk ∈ chunks
        ]
        reduce(vcat, fetch.(tasks))
    end
    for (method, s, ℓₘₐₓ, kw) in (
        ("RS", -2, 12, (;)), ("RS", 1//2, 15//2, (;)),
        ("Minimal", 0, 10, (;)), ("Minimal", 2, 8, (; inplace=false)),
    )
        𝒯 = SSHT(s, ℓₘₐₓ; method, kw...)
        inputs = [randn(rng, ComplexF64, nmodes(𝒯)) for _ in 1:48]
        serial = [roundtrip(𝒯, f̃) for f̃ in inputs]
        @test in_tasks(𝒯, inputs, 8) == serial
    end

    # An SSHTMatrix holds no workspace, so one object may serve any number of tasks at once,
    # both in place and not, and with the QR decomposition of a least-squares problem
    for 𝒯 in (
        SSHTMatrix(2, 10), SSHTMatrix(2, 10; inplace=false),
        SSHTMatrix(1, 6; Rθϕ=golden_ratio_spiral_rotors(0, 7)),  # more points than modes
        SSHTMatrix(1//2, 11//2; decomposition=qr, inplace=false),
    )
        inputs = [randn(rng, ComplexF64, nmodes(𝒯)) for _ in 1:12]
        function work(f̃s)
            f = zeros(ComplexF64, npixels(𝒯))
            g̃ = zeros(ComplexF64, nmodes(𝒯))
            map(f̃s) do f̃
                mul!(f, 𝒯, f̃)
                ldiv!(g̃, 𝒯, f)
                (copy(f), copy(g̃), roundtrip(𝒯, f̃))
            end
        end
        serial = work(inputs)
        @test all(==(serial), fetch.([Threads.@spawn(work(inputs)) for _ in 1:16]))
    end
end


@testitem "SSHT deepcopy and serialization" begin
    import SphericalFunctions
    import SphericalFunctions: SSHT
    import SphericalFunctions: nmodes  # unexported
    using DoubleFloats: Double64
    using FFTW: FFTW
    using Serialization: Serialization, serialize, deserialize
    using Random

    rng = Random.Xoshiro(20260928)

    # An FFTW plan object wraps a pointer to a plan that belongs to the process that made it.
    # A deep copy of a transform shares the plans of the original: a copied plan object would
    # wrap the pointer without owning it, and the copy would execute freed memory once the
    # original had been collected.  A deserialized transform makes its plans again, since
    # Serialization writes every pointer as a null pointer.  Either failure would crash the
    # process rather than fail a test, so each round trip runs only after the property that
    # makes it safe has been checked.
    function reference(𝒯, T)
        f̃ = randn(rng, Complex{T}, nmodes(𝒯), 2)
        f = 𝒯 * copy(f̃)
        (f̃, f, 𝒯 \ copy(f))
    end
    cases = (
        ("RS", 2, 12, Float64, (;)), ("RS", 1//2, 7//2, Float32, (;)),
        ("RS", 1, 4, BigFloat, (;)), ("RS", -1, 6, Float64, (; plan_fft_flags=FFTW.MEASURE)),
        ("Minimal", 2, 8, Float64, (; inplace=false)),
        ("Minimal", 0, 5, Double64, (; inplace=false)),
        ("Matrix", 2, 8, Float64, (;)), ("Matrix", 1//2, 7//2, Float32, (;)),
    )

    # A deep copy, of which the original is collected before it is used
    make_deep_copy(s, ℓₘₐₓ, T, method, kw) = deepcopy(SSHT(s, ℓₘₐₓ, T; method, kw...))
    for (method, s, ℓₘₐₓ, T, kw) in cases
        𝒯 = SSHT(s, ℓₘₐₓ, T; method, kw...)
        if method != "Matrix"
            d = deepcopy(𝒯)
            safe = d.plans.forward === 𝒯.plans.forward && d.F !== 𝒯.F
            @test safe
            safe || continue
        end
        f̃, f, g̃ = reference(𝒯, T)
        𝒯₂ = make_deep_copy(s, ℓₘₐₓ, T, method, kw)
        for _ in 1:3
            GC.gc(true)
        end
        @test 𝒯₂ * copy(f̃) == f
        @test 𝒯₂ \ copy(f) == g̃
    end

    # A round trip through Serialization, as a transform sent to another process makes.  (The
    # plans are made again with the same options, so for the planner's default of
    # `FFTW.ESTIMATE` they are the same plans, and the results are identical.)
    roundtrip(x) = (io = IOBuffer(); serialize(io, x); deserialize(seekstart(io)))
    for (method, s, ℓₘₐₓ, T, kw) in cases
        𝒯 = SSHT(s, ℓₘₐₓ, T; method, kw...)
        if method != "Matrix"
            m = which(
                Serialization.deserialize,
                Tuple{Serialization.Serializer{IOBuffer}, Type{typeof(𝒯.plans)}}
            )
            safe = m.module === SphericalFunctions
            @test safe
            safe || continue
        end
        f̃, f, g̃ = reference(𝒯, T)
        𝒯₂ = roundtrip(𝒯)
        @test typeof(𝒯₂) === typeof(𝒯)
        if method != "Matrix"
            @test 𝒯₂.plans.forward[1] !== 𝒯.plans.forward[1]  # made again
            @test (𝒯₂.plans.flags, 𝒯₂.plans.timelimit) == (𝒯.plans.flags, 𝒯.plans.timelimit)
        end
        if get(kw, :plan_fft_flags, FFTW.ESTIMATE) == FFTW.ESTIMATE
            @test 𝒯₂ * copy(f̃) == f
            @test 𝒯₂ \ copy(f) == g̃
        else  # a measured plan may be another algorithm, which rounds differently
            ϵ = 100ℓₘₐₓ * eps(T)
            @test 𝒯₂ * copy(f̃) ≈ f atol=ϵ rtol=ϵ
            @test 𝒯₂ \ copy(f) ≈ g̃ atol=ϵ rtol=ϵ
        end
    end
end


@testitem "SSHT transforms run on the caller's thread, with one plan per ring size" begin
    import SphericalFunctions: SSHT, SSHTRS, SSHTMinimal, Ysize
    import SphericalFunctions: nmodes, npixels  # unexported
    using FFTW: FFTW
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(20260929)

    # A transform runs on its caller's thread.  A `@threads` loop inside it, or an FFTW plan
    # made for several threads, allocates for every task it spawns, so the allocation of a
    # transform that writes into given storage bounds both.  (Measured: 176–224 bytes, the
    # headers of the reshaped arrays and views.  A single `@threads` loop allocates about
    # 1.5 kB, and each execution of a threaded plan for a ring of 33 points about 8 kB.)
    function allocation(f, args...)
        f(args...)
        @allocated f(args...)
    end
    for (method, s, ℓₘₐₓ) in (
        ("RS", -2, 16), ("RS", 1//2, 31//2), ("Minimal", 0, 12), ("Minimal", 2, 10)
    )
        𝒯 = SSHT(s, ℓₘₐₓ; method)
        f̃ = randn(rng, ComplexF64, nmodes(𝒯))
        f = zeros(ComplexF64, npixels(𝒯))
        g̃ = zeros(ComplexF64, nmodes(𝒯))
        @test allocation(mul!, f, 𝒯, f̃) < 1000
        @test allocation(ldiv!, g̃, 𝒯, f) < 1000
        if method == "Minimal"
            @test allocation(mul!, 𝒯, copy(f̃)) < 1000
            @test allocation(ldiv!, 𝒯, copy(f)) < 1000
        end
    end
    # ... and so does one of several columns, once the workspace of a chunk of columns exists
    let 𝒯 = SSHT(-2, 16)
        F̃ = randn(rng, ComplexF64, nmodes(𝒯), 11)
        F = zeros(ComplexF64, npixels(𝒯), 11)
        @test allocation(mul!, F, 𝒯, F̃) < 1000
        @test allocation(ldiv!, similar(F̃), 𝒯, F) < 1000
    end

    # One pair of FFT plans is made for each distinct number of points on a ring — a single pair
    # for the default grid, whose rings are all alike — and each is made for a single thread.
    # (FFTW describes a plan for several threads as a "thr" solver.)
    for T in (Float64, Float32)
        𝒯 = SSHTRS(2, 32, T)
        @test length(𝒯.plans.forward) == length(𝒯.plans.backward) == 1
        @test 𝒯.plans.sizes == [65] && all(==(1), 𝒯.plans.index)
        @test !any(p -> occursin("thr", string(p)), [𝒯.plans.forward; 𝒯.plans.backward])
    end
    let Nϕ = [15 + k % 3 for k in 1:15], 𝒯 = SSHTRS(1, 7; Nϕ)
        @test 𝒯.plans.sizes == unique(Nϕ)
        @test 𝒯.plans.sizes[𝒯.plans.index] == Nϕ
    end
    # The rings of "Minimal" for s ≠ 0 come in pairs of equal size
    let 𝒯 = SSHTMinimal(2, 10)
        @test 𝒯.plans.sizes == unique(𝒯.Nϕ) && length(𝒯.plans.sizes) < length(𝒯.Nϕ)
        @test 𝒯.plans.sizes[𝒯.plans.index] == 𝒯.Nϕ
    end
end


@testitem "SSHTRS transforms several columns a chunk at a time, each column exactly" begin
    import SphericalFunctions: SSHT, SSHTRS, ModeWeights, sYlm_matrix, rotors
    import SphericalFunctions: nmodes, npixels  # unexported
    using DoubleFloats: Double64
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(20260930)

    # The columns of batched data are transformed a chunk of up to 8 at a time, with the
    # harmonics of each ℓ computed once for the chunk.  The arithmetic of each column is that of
    # a transform of the column alone, in the same order, so the results agree to the last bit —
    # for 11 columns (a chunk of 8 and one of 3), in three dimensions, and for a half-integer
    # spin weight, whose ring phase factors are tabulated for each size of ring.
    for (s, ℓₘₐₓ, T) in (
        (-2, 16, Float64), (0, 6, Float32), (3//2, 23//2, Float64), (1, 5, Double64)
    )
        𝒯 = SSHT(s, ℓₘₐₓ, T)
        n, N = nmodes(𝒯), npixels(𝒯)
        F̃ = randn(rng, Complex{T}, n, 11)
        F = 𝒯 * F̃
        @test all(F[:, k] == 𝒯 * F̃[:, k] for k in 1:11)
        G̃ = 𝒯 \ F
        @test all(G̃[:, k] == parent(𝒯 \ F[:, k]) for k in 1:11)
        F₃ = reshape(F, N, 1, 11)
        @test 𝒯 \ F₃ == reshape(G̃, n, 1, 11)
        # No columns at all is an empty result
        @test size(𝒯 * zeros(Complex{T}, n, 0)) == (N, 0)
        @test size(𝒯 \ zeros(Complex{T}, N, 0)) == (n, 0)
    end

    # Rings of several different sizes, for a half-integer spin weight: each ring uses the
    # plans and phase factors of its own size.  Synthesis is the evaluation of the harmonics at
    # the transform's rotors, and the analysis inverts it (the rings number 2ℓₘₐₓ+1, of at least
    # 2ℓₘₐₓ+1 points each).  Measured over 20 draws: synthesis within 14 eps of the largest
    # value, and round trips within 6 eps.
    for (s, ℓₘₐₓ) in ((1//2, 11//2), (-5//2, 9//2))
        Nθ = Int(2ℓₘₐₓ + 1)
        𝒯 = SSHTRS(s, ℓₘₐₓ; Nϕ=[Nθ + k % 3 for k in 1:Nθ])
        @test length(𝒯.plans.sizes) == 3
        ϵ = 100ℓₘₐₓ * eps()
        F̃ = randn(rng, ComplexF64, nmodes(𝒯), 2)
        F = 𝒯 * F̃
        @test F ≈ sYlm_matrix(rotors(𝒯), ℓₘₐₓ, s) * F̃ atol=ϵ rtol=ϵ
        @test 𝒯 \ F ≈ F̃ atol=ϵ rtol=ϵ
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
    @test SSHTMatrix(1//2, 5//2, Float32) isa SSHTMatrix{Float32}
    @test SSHTRS(1//2, 5//2, Float32) isa SSHTRS{Float32}

    # The spellings agree, and the two kinds of index may not be mixed
    @test typeof(SSHT(HalfOddInteger(1//2), HalfOddInteger(7//2))) === typeof(𝒯)
    mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
    @test_throws ArgumentError SSHT(1//2, 4)
    @test_throws mixed SSHT(1//2, 4)
    @test_throws mixed SSHTMatrix(1, 7//2)
    @test_throws mixed SSHTRS(1//2, 4)
    @test_throws mixed map2salm(zeros(ComplexF64, 8, 8), 1//2, 3)
    @test_throws mixed salm2map(zeros(ComplexF64, 8), 1//2, 3, 8, 8)
    @test_throws mixed SphericalFunctions.map2salm_plan(zeros(ComplexF64, 8, 8), 1, 7//2)
    @test_throws "exceeds ℓₘₐₓ" SSHT(5//2, 3//2)

    # The Minimal method is integer-only (decision E3), and names the methods that are not
    for f in (() -> SSHT(1//2, 7//2; method="Minimal"), () -> SSHTMinimal(1//2, 7//2))
        @test_throws ArgumentError f()
        @test_throws "does not accept half-odd-integers" f()
        @test_throws "The \"RS\" method (the default) and the \"Matrix\" method both accept" f()
    end
    # ... and so are its rings
    @test_throws "`sorted_rings` and `sorted_ring_pixels` accept half-integer ones" SphericalFunctions.minimal_rings(1//2, 5//2)

    # The integer path stores `Int`, and integers of other types are refused
    @test SSHT(1, 4).s === 1
    @test SSHT(1, 4; method="Minimal").s === 1 && SSHT(1, 4; method="Matrix").ℓₘₐₓ === 4
    for method in ("RS", "Minimal", "Matrix")
        @test_throws ArgumentError SSHT(Int8(1), Int8(4); method)
        @test_throws "`Int8` is narrower than `Int`" SSHT(Int8(1), 4; method)
    end
end

@testitem "SSHT half-integer round trips and cross-method agreement" begin
    using SphericalFunctions
    using SphericalFunctions: HalfOddInteger
    using Quaternionic
    using LinearAlgebra: mul!, ldiv!
    using Random

    # Measured over these cases, in the norm of the error relative to that of the result: at
    # most 3.5ℓₘₐₓ eps, for every comparison below, so 100ℓₘₐₓ eps leaves a factor of ≳ 28.
    rng = MersenneTwister(2026_09_17)
    for s ∈ (1//2, -1//2, 3//2, -5//2), ℓₘₐₓ ∈ (abs(s), 7//2, 11//2)
        ϵ = 100ℓₘₐₓ * eps()
        f̃ = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
        𝒯r = SSHT(s, ℓₘₐₓ)
        𝒯m = SSHT(s, ℓₘₐₓ; method="Matrix")

        # Synthesis by the ring algorithm is the direct evaluation of the harmonics at its own
        # rotors, which is the strongest check available: the two computations share nothing
        # but the calculator
        fr = 𝒯r * f̃
        @test fr ≈ sYlm_matrix(rotors(𝒯r), ℓₘₐₓ, s) * parent(f̃) rtol=ϵ
        fm = 𝒯m * copy(f̃)
        @test fm ≈ sYlm_matrix(rotors(𝒯m), ℓₘₐₓ, s) * parent(f̃) rtol=ϵ

        # Analysis inverts synthesis, for both methods, and labels its result
        g̃r = 𝒯r \ fr
        @test g̃r isa ModeWeights
        @test spin(g̃r) === HalfOddInteger(s) && SphericalFunctions.ℓₘₐₓ(g̃r) === HalfOddInteger(ℓₘₐₓ)
        @test parent(g̃r) ≈ parent(f̃) rtol=ϵ
        @test parent(𝒯m \ copy(fm)) ≈ parent(f̃) rtol=ϵ

        # The Matrix method on the ring grid agrees with the ring algorithm
        𝒯mr = SSHTMatrix(s, ℓₘₐₓ; Rθϕ=rotors(𝒯r))
        @test 𝒯mr * copy(f̃) ≈ fr rtol=ϵ
        @test parent(𝒯mr \ copy(fr)) ≈ parent(f̃) rtol=ϵ

        # The in-place forms compute the same values
        f2 = similar(fr)
        mul!(f2, 𝒯r, f̃)
        @test f2 == fr
        g2 = similar(f̃)
        ldiv!(g2, 𝒯r, fr)
        @test g2 == g̃r

        # Trailing dimensions are transformed column by column, each exactly as it would be
        # alone (and doubling the weights doubles every rounded result exactly)
        F̃ = hcat(parent(f̃), 2 .* parent(f̃))
        F = 𝒯r * F̃
        @test F[:, 1] == fr && F[:, 2] == 2fr
        @test 𝒯r \ F ≈ F̃ rtol=ϵ
        @test (𝒯r \ F)[:, 1] == parent(g̃r)
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
    # (Measured: within 2.1ℓₘₐₓ eps, in the norm of the error relative to that of `f`.)
    ϵ = 100ℓₘₐₓ * eps()
    R₊ = [from_spherical_coordinates(θ, ϕ + 2π) for (θ, ϕ) ∈ pixels(𝒯)]
    @test sYlm_matrix(R₊, ℓₘₐₓ, s) * parent(f̃) ≈ -f rtol=ϵ
    @test sYlm_matrix(rotors(𝒯), ℓₘₐₓ, s) * parent(f̃) ≈ f rtol=ϵ
end

@testitem "SSHTMatrix mul! and ldiv! with aliased arguments" begin
    import SphericalFunctions: SSHTMatrix, ModeWeights, Ysize
    using LinearAlgebra: mul!, ldiv!, qr
    using Random

    # With as many points as modes the input and output of the three-argument forms can be
    # the same array, which is the obvious way to act in place with them.  A BLAS product
    # cannot write over its own input — `mul!(x, 𝒯, x)` would return zeros — so the input is
    # copied whenever it shares memory with the output.
    rng = Random.Xoshiro(20260923)
    for (s, ℓₘₐₓ) ∈ ((0, 4), (2, 5)), inplace ∈ (true, false), decomposition ∈ (nothing, qr)
        kw = decomposition === nothing ? (; inplace) : (; inplace, decomposition)
        𝒯 = SSHTMatrix(s, ℓₘₐₓ; kw...)
        n = Ysize(abs(s), ℓₘₐₓ)
        f̃ = randn(rng, ComplexF64, n)
        f = 𝒯.Y * f̃
        # (Measured: the products exact, or within 2.4ℓₘₐₓ eps, and the solves 6.5ℓₘₐₓ eps.)
        ϵ = 100ℓₘₐₓ * eps()

        x = copy(f̃)
        @test mul!(x, 𝒯, x) === x
        @test x ≈ f atol=ϵ
        y = copy(f)
        @test parent(array_view(ldiv!(y, 𝒯, y))) === y
        @test y ≈ f̃ atol=ϵ

        # ... and when the aliasing is through a ModeWeights wrapping the output's storage
        x = copy(f̃)
        @test mul!(x, 𝒯, ModeWeights(x, s)) === x
        @test x ≈ f atol=ϵ
        y = copy(f)
        w = ModeWeights(y, s)
        @test ldiv!(w, 𝒯, y) === w
        @test parent(w) ≈ f̃ atol=ϵ

        # Several columns
        X = hcat(f̃, 2f̃)
        @test mul!(X, 𝒯, X) === X
        @test X ≈ hcat(f, 2f) atol=2ϵ
    end
end


@testitem "SSHTRS: rings without their weights are refused, and an inexact rule warns" begin
    import SphericalFunctions: SSHT, SSHTRS
    import SphericalFunctions: fejer1_rings, fejer2_rings, clenshaw_curtis_rings
    import SphericalFunctions: fejer1, fejer2, clenshaw_curtis
    using DoubleFloats: Double64

    # The weights belong to the rule that placed the rings, and only the caller knows which
    # rule that was
    @test_throws ArgumentError SSHTRS(1, 4; θ=clenshaw_curtis_rings(9))
    @test_throws ArgumentError SSHTRS(1, 4; θ=fejer1_rings(11))
    @test_throws ArgumentError SSHT(1, 4; θ=clenshaw_curtis_rings(9))
    @test_throws ArgumentError SSHTRS(1//2, 7//2; θ=fejer2_rings(8))
    @test_throws "quadrature_weights" SSHTRS(1, 4; θ=clenshaw_curtis_rings(9))

    # The analysis is exact only when the rule integrates polynomials of degree 2ℓₘₐₓ in cos θ,
    # which too few rings cannot do, and which one rule's weights on another's rings do not
    @test_logs (:warn, r"exact") SSHTRS(1, 4; θ=fejer1_rings(5), quadrature_weights=fejer1(5))
    @test_logs (:warn, r"exact") SSHTRS(
        1, 4; θ=clenshaw_curtis_rings(9), quadrature_weights=fejer1(9)
    )
    @test_logs (:warn, r"exact") SSHTRS(
        1//2, 7//2; θ=fejer1_rings(6), quadrature_weights=fejer1(6)
    )
    # For a half-odd ℓₘₐₓ the degree 2ℓₘₐₓ is odd, and a rule symmetric about the equator
    # integrates odd degrees exactly, so 2ℓₘₐₓ of its rings suffice (measured: a round trip at
    # ℓₘₐₓ = 7/2 on 7 Fejér rings is exact to 1e-15, and wrong by 4e-2 on 6)
    @test_logs SSHTRS(1//2, 7//2; θ=fejer1_rings(7), quadrature_weights=fejer1(7))
    @test_logs SSHTRS(1//2, 7//2; θ=fejer2_rings(7), quadrature_weights=fejer2(7))

    # Exact rules raise no warning, at any size and in any of the types the transforms are
    # tested in
    for T ∈ (Float64, Float32, Double64), (s, ℓₘₐₓ) ∈ ((0, 0), (1, 4), (-2, 12), (1//2, 7//2))
        N = Int(2ℓₘₐₓ + 1)
        @test_logs SSHTRS(s, ℓₘₐₓ, T)
        @test_logs SSHTRS(s, ℓₘₐₓ, T; θ=fejer1_rings(N, T), quadrature_weights=fejer1(N, T))
        @test_logs SSHTRS(s, ℓₘₐₓ, T; θ=fejer2_rings(N, T), quadrature_weights=fejer2(N, T))
        @test_logs SSHTRS(
            s, ℓₘₐₓ, T; θ=fejer1_rings(N + 3, T), quadrature_weights=fejer1(N + 3, T)
        )
        if N ≥ 2
            @test_logs SSHTRS(
                s, ℓₘₐₓ, T; θ=clenshaw_curtis_rings(N, T), quadrature_weights=clenshaw_curtis(N, T)
            )
        end
    end
end


@testitem "SSHT: three-argument mul! and ldiv! refuse mismatched trailing dimensions" begin
    import SphericalFunctions: SSHT
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(2718)

    # Every method transforms each column of the trailing dimensions separately, so the input
    # and the output must agree in them; broadcasting would otherwise copy one input column
    # into every column of the output, or fail with a message that says nothing of the
    # transform.
    for method ∈ ("RS", "Minimal", "Matrix"), (s, ℓₘₐₓ) ∈ ((1, 4), (-2, 3))
        kw = method == "RS" ? (;) : (; inplace=false)
        𝒯 = SSHT(s, ℓₘₐₓ; method, kw...)
        n, N = nmodes(𝒯), npixels(𝒯)
        calls = (
            () -> mul!(zeros(ComplexF64, N, 3), 𝒯, randn(rng, ComplexF64, n)),
            () -> mul!(zeros(ComplexF64, N), 𝒯, randn(rng, ComplexF64, n, 3)),
            () -> mul!(zeros(ComplexF64, N, 2), 𝒯, randn(rng, ComplexF64, n, 3)),
            () -> ldiv!(zeros(ComplexF64, n, 3), 𝒯, randn(rng, ComplexF64, N)),
            () -> ldiv!(zeros(ComplexF64, n), 𝒯, randn(rng, ComplexF64, N, 3)),
            () -> ldiv!(zeros(ComplexF64, n, 2), 𝒯, randn(rng, ComplexF64, N, 3)),
        )
        for call ∈ calls
            @test_throws DimensionMismatch call()
            @test_throws "Trailing dimensions" call()
        end
    end
end


@testitem "SSHT synthesis of mode weights of another range of ℓ" begin
    import SphericalFunctions: SSHT, ModeWeights, Ysize, Yindex, spin, ð, ð̄, salm2map
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: mul!
    using Random

    rng = Random.Xoshiro(20260924)

    # Synthesis accepts mode weights of the transform's spin weight with any range of ℓ up to
    # its ℓₘₐₓ: the modes they lack are zero, and their entries with ℓ < |s|, which belong to
    # no function of spin weight s, are ignored.  The reference is the same weights copied by
    # hand into the transform's range.
    function padded(w, ℓₘₐₓ)
        s = spin(w)
        v = zeros(eltype(w), Ysize(abs(s), ℓₘₐₓ))
        for ℓ in max(abs(s), SphericalFunctions.ℓₘᵢₙ(w)):SphericalFunctions.ℓₘₐₓ(w), m in -ℓ:ℓ
            v[Yindex(ℓ, m, abs(s))] = w[ℓ, m]
        end
        ModeWeights(v, s)
    end
    for method in ("RS", "Minimal", "Matrix"), (s, ℓₘₐₓ) in ((1, 6), (-2, 6), (1//2, 11//2))
        (method == "Minimal" && !(s isa Integer)) && continue
        𝒯 = SSHT(s, ℓₘₐₓ; method)
        # Ranges inside the transform's, below |s| as well as inside it, and — for an integer
        # s — only below |s|, when the function is zero
        ℓ₀ = s isa Integer ? 0 : 1//2
        ranges = [(abs(s) + 1, ℓₘₐₓ - 1), (ℓ₀, ℓₘₐₓ), (ℓ₀, abs(s)), (abs(s), abs(s))]
        s isa Integer && push!(ranges, (0, abs(s) - 1))
        for (ℓₘᵢₙ, ℓₘₐₓ′) in ranges
            w = ModeWeights(randn(rng, ComplexF64, Ysize(ℓₘᵢₙ, ℓₘₐₓ′)), s, ℓₘᵢₙ, ℓₘₐₓ′)
            w₀ = copy(parent(w))
            reference = padded(w, ℓₘₐₓ)
            @test 𝒯 * w == 𝒯 * ModeWeights(copy(parent(reference)), s)
            # A copy in the transform's range is what synthesis uses, even for the in-place
            # "Minimal", which would otherwise overwrite the storage of `w`
            @test parent(w) == w₀
        end
        # Weights beyond the transform's ℓₘₐₓ are refused rather than dropped
        w = ModeWeights(zeros(ComplexF64, Ysize(abs(s), ℓₘₐₓ + 1)), s)
        @test_throws ArgumentError 𝒯 * w
        @test_throws "synthesizes ℓ ≤ $ℓₘₐₓ only" 𝒯 * w
        # `mul!` needs the transform's range exactly, so that it allocates nothing
        @test_throws ArgumentError mul!(
            zeros(ComplexF64, npixels(𝒯)), 𝒯, ModeWeights(zeros(ComplexF64, Ysize(abs(s), ℓₘₐₓ - 1)), s)
        )
    end

    # This is what makes a spin-changing operator's output synthesizable directly: ð and ð̄
    # keep the range of ℓ of their input, which then either includes ℓ below the new |s| or
    # lacks some of the new range.  The analysis → operator → synthesis chain gives exactly
    # what the same weights copied into the new range by hand give.
    for (s, op) in ((1, ð), (1, ð̄), (-2, ð), (-2, ð̄), (0, ð), (1//2, ð), (3//2, ð̄)), method in ("RS", "Matrix")
        ℓₘₐₓ = s isa Integer ? 6 : 11//2
        𝒯 = SSHT(s, ℓₘₐₓ; method)
        w = ModeWeights(randn(rng, ComplexF64, nmodes(𝒯)), s)
        dw = op * (𝒯 \ (𝒯 * w))
        𝒯′ = SSHT(spin(dw), ℓₘₐₓ; method)
        @test 𝒯′ * dw == 𝒯′ * padded(dw, ℓₘₐₓ)
        N = Int(2ℓₘₐₓ + 2)
        @test salm2map(dw, N, N) == salm2map(padded(dw, ℓₘₐₓ), N, N)
    end
end


@testitem "SSHT real data" begin
    import SphericalFunctions: SSHT, ModeWeights, Ysize, map2salm, salm2map
    import SphericalFunctions: nmodes, npixels  # unexported
    using DoubleFloats: Double64
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(314159)

    # Real mode weights and real function values — the ordinary case for spin weight 0 — are
    # accepted by every method that does not act in place, and give the complex results of the
    # same data made complex, in the transform's type.  (Measured: identical for "RS" and
    # "Minimal", and within 2 eps for "Matrix", whose BLAS product differs for real input.)
    for T in (Float64, Double64), method in ("RS", "Minimal", "Matrix")
        kw = method == "RS" ? (;) : (; inplace=false)
        𝒯 = SSHT(0, 4, T; method, kw...)
        ϵ = 20eps(T)
        f̃ = randn(rng, T, nmodes(𝒯))
        f = 𝒯 * f̃
        @test f isa Vector{Complex{T}}
        @test f ≈ 𝒯 * complex.(f̃) atol=ϵ rtol=ϵ
        @test 𝒯 * ModeWeights(f̃, 0) ≈ f atol=ϵ rtol=ϵ
        r = randn(rng, T, npixels(𝒯))
        w = 𝒯 \ r
        @test w isa ModeWeights{Complex{T}}
        @test parent(w) ≈ parent(𝒯 \ complex.(r)) atol=ϵ rtol=ϵ
        # The outputs of `mul!` and `ldiv!` must hold the complex results
        @test_throws ArgumentError mul!(zeros(T, npixels(𝒯)), 𝒯, f̃)
        @test_throws "must hold complex floating-point numbers" ldiv!(zeros(T, nmodes(𝒯)), 𝒯, r)
    end

    # Data in another precision is transformed in the transform's type, by "Minimal" as by the
    # others, rather than in the type of the data
    for method in ("RS", "Minimal", "Matrix")
        𝒯 = SSHT(0, 4; method, (method == "RS" ? (;) : (; inplace=false))...)
        @test 𝒯 * randn(rng, ComplexF32, nmodes(𝒯)) isa Vector{ComplexF64}
        @test 𝒯 \ randn(rng, ComplexF32, npixels(𝒯)) isa ModeWeights{ComplexF64}
    end

    # `map2salm` and `salm2map` take real data, working in its floating-point type
    let ℓₘₐₓ = 4, N = 2ℓₘₐₓ + 1
        m = randn(rng, N, N)
        @test map2salm(m, 0, ℓₘₐₓ) isa ModeWeights{ComplexF64}
        @test parent(map2salm(m, 0, ℓₘₐₓ)) == parent(map2salm(complex.(m), 0, ℓₘₐₓ))
        @test map2salm(Float32.(m), 0, ℓₘₐₓ) isa ModeWeights{ComplexF32}
        a = randn(rng, Ysize(0, ℓₘₐₓ))
        @test salm2map(a, 0, ℓₘₐₓ, N, N) == salm2map(complex.(a), 0, ℓₘₐₓ, N, N)
        @test salm2map(ModeWeights(a, 0), N, N) == salm2map(a, 0, ℓₘₐₓ, N, N)
    end
end


@testitem "SSHT in-place analysis of storage that is not contiguous" begin
    import SphericalFunctions: SSHT, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: ldiv!, qr
    using Random

    rng = Random.Xoshiro(2718281)

    # A strided view, or a row of a matrix, is analyzed in place by every in-place form, as by
    # those that are not, although LAPACK solves only in contiguous storage, and refuses other
    # storage with "matrix does not have contiguous columns".
    for (method, kw) in (("Matrix", (;)), ("Matrix", (; inplace=false)), ("Matrix", (; decomposition=qr)),
                         ("Minimal", (;)), ("Minimal", (; inplace=false)))
        𝒯 = SSHT(-1, 4; method, kw...)
        n = nmodes(𝒯)
        f̃ = randn(rng, ComplexF64, n)
        f = 𝒯 * copy(f̃)
        reference = collect(𝒯 \ copy(f))
        A = zeros(ComplexF64, 2n)
        A[1:2:end] .= f
        C = zeros(ComplexF64, 3, n)
        C[2, :] .= f
        # (The same numbers are solved for, in a copy, so the results agree exactly.)
        @test collect(𝒯 \ view(copy(A), 1:2:2n)) == reference
        @test collect(𝒯 \ view(copy(C), 2, :)) == reference
        x = view(copy(A), 1:2:2n)
        @test collect(ldiv!(𝒯, x)) == reference
        @test x == reference  # in place
    end
end


@testitem "SSHT function values held in a ModeWeights" begin
    import SphericalFunctions: SSHT, ModeWeights
    import SphericalFunctions: nmodes, npixels  # unexported
    using LinearAlgebra: mul!, ldiv!
    using Random

    rng = Random.Xoshiro(141421)

    # Function values may be held in a ModeWeights where there are as many of them as modes,
    # or by `array_view` in any case; every analysis path unwraps them, as every synthesis
    # path unwraps its mode weights.
    for (method, kw) in (("RS", (;)), ("Matrix", (;)), ("Matrix", (; inplace=false)),
                         ("Minimal", (;)), ("Minimal", (; inplace=false)))
        𝒯 = SSHT(0, 4; method, kw...)
        f̃ = randn(rng, ComplexF64, nmodes(𝒯))
        f = 𝒯 * copy(f̃)
        reference = parent(𝒯 \ copy(f))
        # (The same numbers go through the same kernels, so the results agree exactly.)
        wf = ModeWeights(copy(f), 0)
        @test parent(𝒯 \ wf) == reference
        @test parent(ldiv!(zeros(ComplexF64, nmodes(𝒯)), 𝒯, ModeWeights(copy(f), 0))) == reference
        # ... and a ModeWeights may receive function values from `mul!`, which returns it
        out = ModeWeights(zeros(ComplexF64, npixels(𝒯)), 0)
        @test mul!(out, 𝒯, f̃) === out
        @test parent(out) == f
    end
end


@testitem "SSHT inference" begin
    import SphericalFunctions: SSHT, SSHTRS, SSHTMinimal, SSHTMatrix, ModeWeights, Ysize
    import SphericalFunctions: map2salm, salm2map, map2salm_plan
    import SphericalFunctions: golden_ratio_spiral_pixels, leja_pixels, sorted_rings
    using Test: @inferred

    # The entry points infer concrete types from the types of their arguments alone, for
    # indices whose values are not known to the compiler.  The element type is positional,
    # and so part of the type of the call.  (Passed as a keyword, it defeated inference, and
    # `salm2map` inferred `Any`.)  The type of an "RS" transform follows from its arguments;
    # those of "Minimal" and "Matrix" also depend on the value of `inplace`, and that of
    # `SSHT` on the value of `method`, so these are not asserted.
    for (s, ℓₘₐₓ) in ((1, 4), (1//2, 7//2))
        @test @inferred(SSHTRS(s, ℓₘₐₓ)) isa SSHTRS{Float64}
        @test @inferred(SSHTRS(s, ℓₘₐₓ, Float32)) isa SSHTRS{Float32}
        N = Int(2ℓₘₐₓ + 1)
        m = randn(ComplexF64, N, N)
        @test @inferred(map2salm_plan(m, s, ℓₘₐₓ)) isa SSHTRS{Float64}
        w = @inferred map2salm(m, s, ℓₘₐₓ)
        @test w isa ModeWeights{ComplexF64}
        @test @inferred(map2salm(real.(m), s, ℓₘₐₓ)) isa ModeWeights{ComplexF64}
        @test @inferred(salm2map(w, s, ℓₘₐₓ, N, N)) isa Matrix{ComplexF64}
        @test @inferred(salm2map(parent(w), s, ℓₘₐₓ, N, N)) isa Matrix{ComplexF64}
        @test @inferred(salm2map(w, N, N)) isa Matrix{ComplexF64}
        𝒯 = SSHTRS(s, ℓₘₐₓ)
        @test @inferred(𝒯 * w) isa Vector{ComplexF64}
        @test @inferred(𝒯 \ (𝒯 * w)) isa ModeWeights{ComplexF64}
        @test @inferred(golden_ratio_spiral_pixels(s, ℓₘₐₓ, Float32)) isa Vector
        @test @inferred(sorted_rings(s, ℓₘₐₓ)) isa Vector{Float64}
    end
    # The in-place analysis of "Minimal" and of an in-place "Matrix" always labels the whole
    # of its argument, so its type follows from the transform's
    for 𝒯 in (SSHT(0, 3; method="Minimal"), SSHT(0, 3; method="Matrix"))
        @test @inferred(𝒯 \ randn(ComplexF64, 16)) isa ModeWeights{ComplexF64}
    end
    @test @inferred(leja_pixels(0, 3, Float32)) isa Vector
end


@testitem "SSHT round trips in BigFloat, at ℓₘₐₓ = 0, and of half-integer spin in other types" begin
    import SphericalFunctions: SSHT, ModeWeights, sYlm_matrix, rotors
    import SphericalFunctions: nmodes  # unexported
    using DoubleFloats: Double64
    using Random

    rng = Random.Xoshiro(20260925)

    # Measured, in units of eps(T): BigFloat at (s, ℓₘₐₓ) = (1, 4) ≤ 25 ("Minimal"; ≤ 8
    # otherwise), a single mode at ℓₘₐₓ = 0 ≤ 1.1, and half-integer spin in Float32 and Double64
    # ≤ 14.  100ℓₘₐₓ eps(T) (and 100 eps(T) at ℓₘₐₓ = 0) leaves a factor of ≳ 7.
    roundtrip(𝒯, f̃) = collect(𝒯 \ (𝒯 * copy(f̃)))
    for method in ("RS", "Minimal", "Matrix")
        kw = method == "RS" ? (;) : (; inplace=false)
        for (T, s, ℓₘₐₓ) in ((BigFloat, 1, 4), (Float64, 0, 0), (Float32, 0, 0), (Double64, 0, 0))
            𝒯 = SSHT(s, ℓₘₐₓ, T; method, kw...)
            ϵ = 100max(ℓₘₐₓ, 1) * eps(T)
            f̃ = randn(rng, Complex{T}, nmodes(𝒯))
            @test roundtrip(𝒯, f̃) ≈ f̃ atol=ϵ rtol=ϵ
            @test 𝒯 * f̃ isa Vector{Complex{T}}
        end
    end
    for method in ("RS", "Matrix"), T in (Float32, Double64), (s, ℓₘₐₓ) in ((1//2, 7//2), (-3//2, 11//2), (5//2, 5//2))
        𝒯 = SSHT(s, ℓₘₐₓ, T; method)
        ϵ = 100ℓₘₐₓ * eps(T)
        f̃ = randn(rng, Complex{T}, nmodes(𝒯))
        @test roundtrip(𝒯, f̃) ≈ f̃ atol=ϵ rtol=ϵ
        # Synthesis is the evaluation of the harmonics at the transform's own rotors
        @test 𝒯 * f̃ ≈ sYlm_matrix(rotors(𝒯), ℓₘₐₓ, s) * f̃ atol=ϵ rtol=ϵ
    end
end
