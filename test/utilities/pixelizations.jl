# Tests of the pixelizations in `src/utilities/pixelizations.jl`.
#
# The golden-ratio spiral and the sorted rings are exercised thoroughly by the transform
# tests in `test/ssht/`, which is where they matter; what is left over — and what this file
# covers — is the part of the module no transform reaches.  That is the two equiangular
# grids, Driscoll–Healy and McEwen–Wiaux, which are public but which nothing in the package
# defaults to, and the two-argument entry points whose only job is to supply `T=Float64`.
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
    import SphericalFunctions: mcewen_wiaux_pixels, mcewen_wiaux_rotors

    # Eq. quoted in the docstring: θₜ = π(2t+1)/(2ℓₘₐₓ-1) for t ∈ 0:ℓₘₐₓ-1, extended to
    # t ∈ 0:2ℓₘₐₓ-1 so that θ covers [0, 2π), and ϕₚ = 2πp/(2ℓₘₐₓ-1) for p ∈ 0:2ℓₘₐₓ-2.
    for T ∈ (Float64, Float32), ℓₘₐₓ ∈ (1, 2, 5)
        p = mcewen_wiaux_pixels(0, ℓₘₐₓ, T)

        @test p isa Vector{<:SVector{2, T}}
        @test length(p) == 2ℓₘₐₓ * (2ℓₘₐₓ - 1)
        @test p == [
            SVector{2, T}(T(π) * (2t + 1) / (2ℓₘₐₓ - 1), 2T(π) * q / (2ℓₘₐₓ - 1))
            for t ∈ 0:(2ℓₘₐₓ - 1) for q ∈ 0:(2ℓₘₐₓ - 2)
        ]

        # The θ extension is what distinguishes this grid: no sample sits at θ=0, and θ runs
        # past π into the second half of its period.
        @test minimum(q[1] for q ∈ p) > 0
        if ℓₘₐₓ > 1
            @test maximum(q[1] for q ∈ p) > T(π)
        end

        @test mcewen_wiaux_pixels(2, ℓₘₐₓ, T) == p
        @test mcewen_wiaux_pixels(ℓₘₐₓ, T) == p

        R = mcewen_wiaux_rotors(0, ℓₘₐₓ, T)
        @test R isa Vector{<:Rotor{T}}
        @test length(R) == length(p)
        @test mcewen_wiaux_rotors(ℓₘₐₓ, T) == R
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
