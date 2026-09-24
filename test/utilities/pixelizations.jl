# Tests of the pixelizations in `src/utilities/pixelizations.jl`.
#
# The golden-ratio spiral and the sorted rings are exercised thoroughly by the transform
# tests in `test/ssht/`, which is where they matter; what is left over — and what this file
# covers — is the part of the module no transform reaches.  That is the two equiangular
# grids, Driscoll–Healy and McEwen–Wiaux, and the Leja points, which are public but which
# nothing in the package defaults to, and the two-argument entry points whose only job is to
# supply `T=Float64`.
#
# Each grid is checked against the formula its docstring quotes from the paper it cites,
# rather than against a stored table, so a change of convention has to be deliberate.

@testitem "Pixelizations: Driscoll–Healy grid" begin
    using StaticArrays: SVector
    using Quaternionic: Rotor, from_spherical_coordinates
    import SphericalFunctions: driscoll_healy_pixels, driscoll_healy_rotors

    # Eq. quoted in the docstring: θᵢ = πi/2b and ϕⱼ = πj/b, for i and j both running over
    # 0:2b-1, where b = ℓₘₐₓ+1 is the band limit.
    for T ∈ (Float64, Float32), ℓₘₐₓ ∈ (0, 1, 2, 5)
        b = ℓₘₐₓ + 1
        p = driscoll_healy_pixels(0, ℓₘₐₓ, T)

        @test p isa Vector{<:SVector{2, T}}
        @test length(p) == (2b)^2
        @test p == [
            SVector{2, T}(T(π) * i / 2b, T(π) * j / b)
            for i ∈ 0:(2b - 1) for j ∈ 0:(2b - 1)
        ]

        # The grid starts at the north pole and never reaches the south pole, since the
        # largest θ is π(2b-1)/2b.
        @test first(p)[1] == 0
        @test maximum(q[1] for q ∈ p) < T(π)
        # ϕ covers [0, 2π) once per θ: the j index runs twice around the b-fold division.
        @test maximum(q[2] for q ∈ p) ≈ T(π) * (2b - 1) / b

        # The `s` argument is documented as ignored, and the one-argument form is the same
        # grid with `s` defaulted to 0.
        @test driscoll_healy_pixels(2, ℓₘₐₓ, T) == p
        @test driscoll_healy_pixels(-1//2, ℓₘₐₓ, T) == p
        @test driscoll_healy_pixels(ℓₘₐₓ, T) == p

        # The rotors are just the same points handed to the coordinate map
        R = driscoll_healy_rotors(0, ℓₘₐₓ, T)
        @test R isa Vector{<:Rotor{T}}
        @test length(R) == length(p)
        @test R == from_spherical_coordinates.(p)
        @test driscoll_healy_rotors(ℓₘₐₓ, T) == R
        @test driscoll_healy_rotors(2, ℓₘₐₓ, T) == R
    end

    # `Float64` is the default element type
    @test driscoll_healy_pixels(0, 3) == driscoll_healy_pixels(0, 3, Float64)
    @test driscoll_healy_pixels(3) == driscoll_healy_pixels(0, 3, Float64)
    @test driscoll_healy_rotors(3) == driscoll_healy_rotors(0, 3, Float64)
end

@testitem "Pixelizations: McEwen–Wiaux grid" begin
    using StaticArrays: SVector
    using Quaternionic: Rotor
    using LinearAlgebra: rank
    import SphericalFunctions: mcewen_wiaux_pixels, mcewen_wiaux_rotors, sYlm_matrix, Ysize

    # Eqs. (17) and (18) of McEwen & Wiaux (2011), with their band limit L defined by
    # ₛf_ℓm = 0 for ℓ ≥ L, so L = ℓₘₐₓ + 1: θₜ = π(2t+1)/(2L-1) for t ∈ 0:L-1, and
    # ϕₚ = 2πp/(2L-1) for p ∈ 0:2L-2.  Transcribed from the paper, not from the code.
    for T ∈ (Float64, Float32), ℓₘₐₓ ∈ (0, 1, 2, 5)
        L = ℓₘₐₓ + 1
        p = mcewen_wiaux_pixels(0, ℓₘₐₓ, T)

        @test p isa Vector{<:SVector{2, T}}
        @test length(p) == L * (2L - 1)
        @test p ≈ [
            SVector{2, T}(T(π) * (2t + 1) / (2L - 1), 2T(π) * q / (2L - 1))
            for t ∈ 0:(L - 1) for q ∈ 0:(2L - 2)
        ] rtol=2eps(T)

        # Every point is on the sphere: θ ∈ (0, π], with no sample at the north pole, and the
        # last ring — 2L-1 points — exactly at the south pole
        @test all(0 < q[1] ≤ T(π) for q ∈ p)
        @test count(q -> q[1] == T(π), p) == 2L - 1

        @test mcewen_wiaux_pixels(2, ℓₘₐₓ, T) == p
        @test mcewen_wiaux_pixels(ℓₘₐₓ, T) == p

        R = mcewen_wiaux_rotors(0, ℓₘₐₓ, T)
        @test R isa Vector{<:Rotor{T}}
        @test length(R) == length(p)
        @test mcewen_wiaux_rotors(ℓₘₐₓ, T) == R
    end

    # The grid determines every mode of the band limit it was built for.  (With L = ℓₘₐₓ,
    # as the grid was once built, the harmonic matrix at ℓₘₐₓ = 2 had rank 6 for 9 modes.)
    for ℓₘₐₓ ∈ (1, 2, 5, 8), s ∈ (0, 1, -2)
        abs(s) ≤ ℓₘₐₓ || continue
        @test rank(sYlm_matrix(mcewen_wiaux_rotors(ℓₘₐₓ), ℓₘₐₓ, s)) == Ysize(abs(s), ℓₘₐₓ)
    end

    @test mcewen_wiaux_pixels(0, 4) == mcewen_wiaux_pixels(0, 4, Float64)
    @test mcewen_wiaux_pixels(4) == mcewen_wiaux_pixels(0, 4, Float64)
    @test mcewen_wiaux_rotors(4) == mcewen_wiaux_rotors(0, 4, Float64)
end

@testitem "Pixelizations: quadrature ring sets" begin
    import SphericalFunctions: fejer1_rings, fejer2_rings, clenshaw_curtis_rings

    # Eqs. (12)-(14) of Reinecke and Seljebotn, as cited in the source.  Fejér's rules place
    # no sample at either pole; Clenshaw–Curtis places one at each.
    for T ∈ (Float64, Float32), N ∈ (2, 3, 8, 15)
        θ1 = fejer1_rings(N, T)
        θ2 = fejer2_rings(N, T)
        θcc = clenshaw_curtis_rings(N, T)

        @test θ1 isa Vector{T} && length(θ1) == N
        @test θ2 isa Vector{T} && length(θ2) == N
        @test θcc isa Vector{T} && length(θcc) == N

        @test θ1 == [(2n + 1) * T(π) / 2N for n ∈ 0:N-1]
        @test θ2 == [n * T(π) / (N + 1) for n ∈ 1:N]
        @test θcc == [n * T(π) / (N - 1) for n ∈ 0:N-1]

        # All three are strictly increasing and stay within [0, π]
        for θ ∈ (θ1, θ2, θcc)
            @test issorted(θ)
            @test allunique(θ)
            @test all(0 .≤ θ .≤ T(π))
        end

        # Neither Fejér rule samples a pole; Clenshaw–Curtis samples both
        @test 0 ∉ θ1 && T(π) ∉ θ1
        @test 0 ∉ θ2 && T(π) ∉ θ2
        @test first(θcc) == 0 && last(θcc) == T(π)

        # Fejér's first rule is symmetric about the equator by construction
        @test θ1 ≈ reverse(T(π) .- θ1)
        @test θ2 ≈ reverse(T(π) .- θ2)
    end

    # `Float64` is the default element type
    @test fejer1_rings(7) == fejer1_rings(7, Float64)
    @test fejer2_rings(7) == fejer2_rings(7, Float64)
    @test clenshaw_curtis_rings(7) == clenshaw_curtis_rings(7, Float64)
end

@testitem "Pixelizations: default element type of the ring pixelizations" begin
    import SphericalFunctions: sorted_ring_pixels, sorted_ring_rotors
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors

    # The two-argument forms exist only to supply `T=Float64` and re-dispatch, and the
    # transform tests always pass `T` explicitly, so they are checked here.
    for (s, ℓₘₐₓ) ∈ ((0, 4), (2, 5), (-2, 3), (1//2, 7//2), (-3//2, 5//2))
        @test sorted_ring_pixels(s, ℓₘₐₓ) == sorted_ring_pixels(s, ℓₘₐₓ, Float64)
        @test sorted_ring_rotors(s, ℓₘₐₓ) == sorted_ring_rotors(s, ℓₘₐₓ, Float64)
        @test golden_ratio_spiral_pixels(s, ℓₘₐₓ) == golden_ratio_spiral_pixels(s, ℓₘₐₓ, Float64)
        @test golden_ratio_spiral_rotors(s, ℓₘₐₓ) == golden_ratio_spiral_rotors(s, ℓₘₐₓ, Float64)
    end

    # A spin weight describing no modes is refused by both families, at the boundary, before
    # any point is placed.
    @test_throws "exceeds ℓₘₐₓ" sorted_ring_pixels(4, 3)
    @test_throws "exceeds ℓₘₐₓ" sorted_ring_rotors(4, 3)
    @test_throws "exceeds ℓₘₐₓ" golden_ratio_spiral_pixels(-4, 3)
    @test_throws "exceeds ℓₘₐₓ" golden_ratio_spiral_rotors(-4, 3)

    # Mixing the two kinds of index is refused with an explanation rather than a MethodError
    @test_throws ArgumentError sorted_ring_pixels(1//2, 3)
    @test_throws ArgumentError golden_ratio_spiral_pixels(1//2, 3)
end

@testitem "Pixelizations: Leja points" begin
    import SphericalFunctions: leja_pixels, leja_rotors, golden_ratio_spiral_pixels,
        golden_ratio_spiral_rotors, sYlm_matrix, Ysize, SSHT
    import SphericalFunctions
    using DoubleFloats: Double64
    using LinearAlgebra: cond
    using Quaternionic: Rotor, from_spherical_coordinates
    using StaticArrays: SVector
    using Random

    # The shape of the result: exactly as many points as modes, a subset of the golden-ratio
    # spiral of twice as many candidates, in the spiral's order (so strictly increasing θ),
    # for either kind of index and each element type
    for T ∈ (Float64, Double64, Float32), (s, ℓₘₐₓ) ∈ ((0, 5), (2, 6), (-1, 4), (1//2, 9//2))
        p = leja_pixels(s, ℓₘₐₓ, T)
        N = Ysize(abs(s), ℓₘₐₓ)
        @test p isa Vector{SVector{2, T}}
        @test length(p) == N
        @test issubset(p, SphericalFunctions.golden_ratio_spiral(2N, T))
        @test issorted(first.(p)) && allunique(p)
        @test leja_rotors(s, ℓₘₐₓ, T) == from_spherical_coordinates.(p)
        @test leja_rotors(s, ℓₘₐₓ, T) isa Vector{Rotor{T}}
    end
    @test leja_pixels(2, 6) == leja_pixels(2, 6, Float64)
    @test leja_rotors(2, 6) == leja_rotors(2, 6, Float64)

    # The point of them: the harmonics on them are well conditioned where those on the spiral
    # of the same number of points are not.  Measured condition numbers: 71 and 100 against
    # 9.0e4 and 5.8e4 at ℓₘₐₓ = 32 for s = 0 and 2, and 19 against 4540 at ℓₘₐₓ = 31/2 for
    # s = 1/2.
    for (s, ℓₘₐₓ) ∈ ((0, 32), (2, 32), (1//2, 31//2))
        @test cond(sYlm_matrix(leja_rotors(s, ℓₘₐₓ), ℓₘₐₓ, s)) < 250
        @test cond(sYlm_matrix(golden_ratio_spiral_rotors(s, ℓₘₐₓ), ℓₘₐₓ, s)) > 1000
    end
    # ... which is what the "Matrix" transform's accuracy depends on: a round trip on them at
    # ℓₘₐₓ = 32 measured 2.2e-13, against 1.9e-10 on its default spiral
    let s = 2, ℓₘₐₓ = 32
        𝒯 = SSHT(s, ℓₘₐₓ; method="Matrix", Rθϕ=leja_rotors(s, ℓₘₐₓ), inplace=false)
        f̃ = randn(Random.Xoshiro(3), ComplexF64, Ysize(abs(s), ℓₘₐₓ))
        @test collect(𝒯 \ (𝒯 * f̃)) ≈ f̃ atol=5e-12 rtol=0
    end

    # `oversampling`: with no extra candidates every one is chosen, and the spiral comes back;
    # more candidates still give the right number of points; fewer than the number of points
    # is refused
    @test leja_pixels(0, 8; oversampling=1) == golden_ratio_spiral_pixels(0, 8)
    @test length(leja_pixels(2, 6; oversampling=3.5)) == Ysize(2, 6)
    @test_throws "must be at least 1" leja_pixels(0, 4; oversampling=0.5)
    @test_throws "must be at least 1" leja_rotors(0, 4; oversampling=0.5)

    # The same refusals as the other pixelizations, at the boundary
    @test_throws "exceeds ℓₘₐₓ" leja_pixels(4, 3)
    @test_throws "exceeds ℓₘₐₓ" leja_rotors(-4, 3)
    @test_throws ArgumentError leja_pixels(1//2, 3)
end
