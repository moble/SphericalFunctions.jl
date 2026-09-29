# Tests of the pixelizations in `src/utilities/pixelizations.jl`.
#
# The geometry of the golden-ratio spiral and of the sorted rings is exercised by the transform
# tests in `test/ssht/`, which is where they matter.  This file covers the rest: the two
# equiangular grids, Driscoll–Healy and McEwen–Wiaux, which no transform uses; the Leja
# points, which are the default points of the "Matrix" transform; the two-argument entry
# points whose only job is to supply `T=Float64`; and the promise that every pixelization is
# the same set of points in every element type `T`, up to the rounding of each coordinate —
# the spiral's azimuths reduced modulo 2π exactly, the Leja points chosen in `Float64`, and
# the rings ordered by exact keys, so that the order for -s is the mirror image of that for s.
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

    # The grid determines every mode of the band limit it was built for.  (Built with
    # L = ℓₘₐₓ, it would give a harmonic matrix of rank 6 for the 9 modes at ℓₘₐₓ = 2.)
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

        # Both Fejér rules are symmetric about the equator by construction, up to the rounding
        # of π - θ
        @test all(abs.(θ1 .- reverse(T(π) .- θ1)) .≤ 2eps(T(π)))
        @test all(abs.(θ2 .- reverse(T(π) .- θ2)) .≤ 2eps(T(π)))
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
    for f ∈ (sorted_ring_pixels, sorted_ring_rotors, golden_ratio_spiral_pixels, golden_ratio_spiral_rotors)
        @test_throws ArgumentError f(4, 3)
        @test_throws "|s|=4 exceeds ℓₘₐₓ=3" f(-4, 3)
        # ... including every spin weight, when ℓₘₐₓ is negative
        @test_throws "|s|=0 exceeds ℓₘₐₓ=-1" f(0, -1)
    end

    # Mixing the two kinds of index is refused with an explanation rather than a MethodError
    mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
    for f ∈ (sorted_ring_pixels, golden_ratio_spiral_pixels)
        @test_throws ArgumentError f(1//2, 3)
        @test_throws mixed f(1//2, 3)
    end
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

    # The candidates are chosen in `Float64` for every T, so every type gets the same points:
    # the same positions in the spiral, each point computed in T.  (Were they chosen in each
    # type, with Julia's generic LU for the types that LAPACK does not handle, only 58 of the
    # 121 points would be shared by Double64 and Float64 at ℓₘₐₓ = 10, and Float32 would depart
    # from Float64 at ℓₘₐₓ = 32.)
    for (s, ℓₘₐₓ, types) ∈ (
        (0, 10, (Float32, Double64, BigFloat)), (2, 8, (Float32, Double64)),
        (1//2, 9//2, (Float32, Double64)), (0, 32, (Float32,)), (2, 32, (Float32,)),
    )
        N = Ysize(abs(s), ℓₘₐₓ)
        positions(T) = indexin(leja_pixels(s, ℓₘₐₓ, T), SphericalFunctions.golden_ratio_spiral(2N, T))
        reference = positions(Float64)
        @test !any(isnothing, reference)
        for T ∈ types
            @test positions(T) == reference
        end
    end

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
    # more candidates still give the right number of points; fewer than the number of points,
    # or infinitely many, is refused
    @test leja_pixels(0, 8; oversampling=1) == golden_ratio_spiral_pixels(0, 8)
    @test length(leja_pixels(2, 6; oversampling=3.5)) == Ysize(2, 6)
    for oversampling ∈ (0.5, Inf, NaN)
        @test_throws ArgumentError leja_pixels(0, 4; oversampling)
        @test_throws "must be at least 1, and finite" leja_rotors(0, 4; oversampling)
    end

    # The same refusals as the other pixelizations, at the boundary
    @test_throws ArgumentError leja_pixels(4, 3)
    @test_throws "exceeds ℓₘₐₓ" leja_rotors(-4, 3)
    @test_throws ArgumentError leja_pixels(1//2, 3)
    @test_throws "mixes integers (ℓₘₐₓ) with half-odd-integers (s)" leja_pixels(1//2, 3)
end

@testitem "Pixelizations: the golden-ratio spiral's azimuths are reduced exactly" begin
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors, Ysize
    using DoubleFloats: Double64
    using Quaternionic: from_spherical_coordinates

    # The azimuth of point k is 2π times the fractional part of k(2-φ), which lies in
    # [0, 2π), computed here from a 400-bit reference.  Measured: correctly rounded in Float64
    # at every point of ℓₘₐₓ = 64 (4225 points), within half an ulp in Float32 and Float16,
    # and within 1.2 ulps in BigFloat; in Double64 the absolute error is below 2.4e-31.  (As
    # k Δϕ in T, unreduced, the azimuths reached 10137 rad at ℓₘₐₓ = 64, and the Float32 and
    # Float64 spirals differed by 4.3e-4 rad at ℓₘₐₓ = 32.)
    exact(k, T) = setprecision(BigFloat, 400) do
        2big(π) * mod(k * (2 - big(MathConstants.φ)), 1)
    end
    for (T, ℓₘₐₓ, ϵ) ∈ (
        (Float64, 64, 0.5eps(2π)), (Float32, 40, 0.5eps(2Float32(π))),
        (Float16, 10, 0.5eps(2Float16(π))), (Double64, 20, 3e-31), (BigFloat, 10, 2eps(2BigFloat(π))),
    )
        p = golden_ratio_spiral_pixels(0, ℓₘₐₓ, T)
        @test length(p) == Ysize(0, ℓₘₐₓ)
        @test all(0 ≤ θϕ[2] < 2T(π) for θϕ ∈ p)
        @test p[1][2] == 0
        @test maximum(k -> abs(big(p[k+1][2]) - exact(k, T)), 0:length(p)-1) ≤ ϵ
        @test golden_ratio_spiral_rotors(0, ℓₘₐₓ, T) == from_spherical_coordinates.(p)
    end
    # The spiral of each type is therefore the Float64 spiral, rounded
    let p64 = golden_ratio_spiral_pixels(1, 32), p32 = golden_ratio_spiral_pixels(1, 32, Float32)
        @test maximum(i -> maximum(abs, p32[i] - Float32.(p64[i])), eachindex(p64)) ≤ 2eps(2Float32(π))
    end
end

@testitem "Pixelizations: the sorted rings are ordered by exact keys" begin
    import SphericalFunctions: sorted_rings, minimal_rings
    using DoubleFloats: Double64

    # The position of each ring among the n equally spaced interior slots of [0, π]
    function slots(θs, n, T)
        grid = collect(LinRange{T}(0, T(π), n + 2))[begin+1:end-1]
        indexin(θs, grid)
    end

    # The rings fill the slots from those farthest from the equator to the nearest, the
    # distance of slot i being |2i - (n+1)| half slot spacings, and of the two rings assigned
    # to a pair of mirror-image slots the larger goes north for s ≥ 0 and south for s < 0.
    # The order is decided by the slots, not by their rounded colatitudes, so it is the same
    # in every T.  (An order decided by the rounded colatitudes would differ between Float32
    # and Float64 at 196 of the ℓₘₐₓ in 1:200 for s = 0, and for s = ±1 would fail to be the
    # mirror image at 58 of them.)
    for s ∈ (0, 1, -1, 2, -3, 1//2, -1//2, 5//2), ℓₘₐₓ ∈ abs(s) .+ (0:40)
        n = Int(ℓₘₐₓ - abs(s) + 1)
        order = slots(sorted_rings(s, ℓₘₐₓ), n, Float64)
        distance = [abs(2i - (n + 1)) for i ∈ order]
        @test issorted(distance; rev=true)
        for q ∈ 2:n
            if distance[q] == distance[q-1]  # a mirror-image pair: the larger ring is second
                @test (order[q] < order[q-1]) == (s ≥ 0)
            end
        end
        for T ∈ (Float32, Double64)
            @test slots(sorted_rings(s, ℓₘₐₓ, T), n, T) == order
        end
    end

    # For s ≠ 0, the order for -s is exactly the mirror image of that for s, out to large
    # ℓₘₐₓ, where an order decided by the rounded colatitudes would differ (first at ℓₘₐₓ = 42
    # for s = 1, and 83/2 for s = 1/2)
    for s ∈ (1//2, 1, 3//2, 2), ℓₘₐₓ ∈ abs(s) .+ (0:100)
        n = Int(ℓₘₐₓ - abs(s) + 1)
        @test slots(sorted_rings(-s, ℓₘₐₓ), n, Float64) == (n + 1) .- slots(sorted_rings(s, ℓₘₐₓ), n, Float64)
    end

    # The rings of the "Minimal" method use the same order for the rings centered on 0, and
    # so are the same in every T; for s = 0 they are exactly the sorted rings.
    for s ∈ (0, 1, -2), ℓₘₐₓ ∈ abs(s) .+ (0:20)
        n = ℓₘₐₓ - abs(s) + 1
        rings = minimal_rings(s, ℓₘₐₓ)
        for T ∈ (Float32, Double64)
            ringsT = minimal_rings(s, ℓₘₐₓ, T)
            @test ringsT.Nϕ == rings.Nϕ && ringsT.centers == rings.centers
            @test slots(ringsT.θ, n, T) == slots(rings.θ, n, Float64)
        end
        s == 0 && @test rings.θ == sorted_rings(0, ℓₘₐₓ)
    end
end

@testitem "Pixelizations: the quadrature ring sets validate the number of rings" begin
    import SphericalFunctions: fejer1_rings, fejer2_rings, clenshaw_curtis_rings

    # As for the weights of the same rules: at least one ring, and at least two for the
    # Clenshaw–Curtis rule, whose rings include both poles
    for T ∈ (Float64, Float32, BigFloat)
        for N ∈ (0, -1)
            @test_throws ArgumentError fejer1_rings(N, T)
            @test_throws ArgumentError fejer2_rings(N, T)
        end
        for N ∈ (1, 0, -1)
            @test_throws ArgumentError clenshaw_curtis_rings(N, T)
        end
        @test fejer1_rings(1, T) ≈ [T(π) / 2] atol=eps(T(π)) rtol=0
        @test fejer2_rings(1, T) ≈ [T(π) / 2] atol=eps(T(π)) rtol=0
        @test clenshaw_curtis_rings(2, T) ≈ [0, T(π)] atol=eps(T(π)) rtol=0
    end
    @test_throws ArgumentError fejer1_rings(0)
    @test_throws ArgumentError fejer2_rings(0)
    @test_throws ArgumentError clenshaw_curtis_rings(1)
end

@testitem "Pixelizations: the equiangular grids refuse a band limit that is not a non-negative integer" begin
    import SphericalFunctions: driscoll_healy_pixels, driscoll_healy_rotors
    import SphericalFunctions: mcewen_wiaux_pixels, mcewen_wiaux_rotors, HalfOddInteger

    # Both grids are defined for integer indices only, the spin weight included although
    # neither uses it; spin-weighted functions of half-integer spin are sampled with the
    # golden-ratio, Leja or sorted-ring pixelizations, which the refusal names.
    integer_only = "does not accept half-odd-integers"
    for f ∈ (driscoll_healy_pixels, driscoll_healy_rotors, mcewen_wiaux_pixels, mcewen_wiaux_rotors)
        for ℓₘₐₓ ∈ (7//2, HalfOddInteger(7//2))
            @test_throws ArgumentError f(ℓₘₐₓ)
            @test_throws integer_only f(ℓₘₐₓ, Float32)
            @test_throws integer_only f(1//2, ℓₘₐₓ)
            @test_throws integer_only f(0, ℓₘₐₓ, Float32)
            # ... with an explanation that names the ones to use
            @test_throws r"may be sampled with `golden_ratio_spiral_(pixels|rotors)`" f(1//2, ℓₘₐₓ)
        end
        @test_throws integer_only f(-1//2, 3)
        @test_throws r"may be sampled with `golden_ratio_spiral_(pixels|rotors)`" f(-1//2, 3)
        # Integers of another type are refused as for every function of indices, and floats
        # are not indices at all
        @test_throws ArgumentError f(Int32(3))
        @test_throws "narrower than `Int`" f(0, Int32(3))
        @test_throws MethodError f(2.5)
        @test_throws MethodError f(0, 2.5)
        # A negative band limit is refused
        @test_throws ArgumentError f(-1)
        @test_throws ArgumentError f(0, -1)
        @test_throws "non-negative band limit" f(-1, Float32)

        # The spin weight is not used
        @test f(-1, 3) == f(3)
        @test f(2, 0) == f(0)
        @test !isempty(f(0))
    end
end
