# Tests of the equiangular-grid convenience functions `map2salm` and `salm2map` from
# `src/ssht/rs.jl`.
#
# The grid is the one documented for `map2salm`: an `Nϕ × Nθ` array with `ϕₖ = 2πk/Nϕ` along
# the first dimension and `θ = clenshaw_curtis_rings(Nθ)` (both poles included) along the
# second.  The oracle is `sYlm_closed_form(s, ℓ, m, θ, ϕ)` from the `Utilities` module
# sampled on that grid — the explicit sum over factorials from the conventions pages, which
# owes nothing to this package — applied to an isolated mode, to a random combination of
# modes (linearity), and to the round trip `map2salm ∘ salm2map`.  The colatitudes and
# azimuths of the grid itself are checked against the formulas the docstring names.
#
# `map2salm` and `salm2map` are documented as front ends that build an `SSHTRS` on the
# Clenshaw–Curtis rings, so the comparisons against an explicitly constructed `SSHTRS` below
# are metamorphic: they pin the wiring — grid, ordering, ℓ range, trailing dimensions — but
# they run the same numerics, and it is the closed form that is the independent reference.
#
# The closed-form comparisons are made in Float64, Float32 and BigFloat.  Float16, for which
# FFTW has no transforms and GenericFFT serves instead, is checked by a round trip in the
# last item.
#
# The item names use a "Transforms: " prefix so that they can be picked out together — by a
# human reading a results list, and by the name filters of `juliati` and the MCP runner.

@testitem "Transforms: map2salm" setup=[Utilities] begin
    import SphericalFunctions: map2salm, SSHTRS, ModeWeights, Yindex, Yrange, Ysize,
        clenshaw_curtis_rings, clenshaw_curtis, spin, ℓₘᵢₙ, ℓₘₐₓ
    import .Utilities: sYlm_closed_form
    using Random
    using LinearAlgebra: norm

    # The closed-form ₛYₗₘ sampled on the `map2salm` grid, as an Nϕ×Nθ array
    function sampled_sYlm(::Type{T}, s, ℓ, m, Nϕ, Nθ) where {T}
        θs = clenshaw_curtis_rings(Nθ, T)
        ϕs = [2T(π) * k / Nϕ for k ∈ 0:Nϕ-1]
        [sYlm_closed_form(s, ℓ, m, θ, ϕ) for ϕ ∈ ϕs, θ ∈ θs]
    end

    ℓmax = 7
    Nθ = 2ℓmax + 1
    Nϕ = 2ℓmax + 2
    rng = Random.Xoshiro(1234)

    # Measured over every T and s ∈ -2:2, in units of eps(T): single-mode analysis ≤ 10.9
    # (worst in BigFloat; 5.5 in Float64, 8.0 in Float32) and linearity ≤ 10.0 eps(T)‖a‖.
    # 100 eps(T) therefore leaves a factor of ≳ 9.  (A tolerance of 30 eps would leave only
    # 2.8 — the tightest margin anywhere in these files, and not enough to be stable across
    # platforms and BLAS versions.)
    for T ∈ (Float64, Float32, BigFloat)
        ϵ = 100eps(T)

        # The grid assumed by `sampled_sYlm` is the one the docstring specifies: the
        # Clenshaw–Curtis colatitudes θₙ = nπ/(Nθ-1) for n = 0, …, Nθ-1, which include both
        # poles, and the equally spaced azimuths ϕₖ = 2πk/Nϕ
        @test clenshaw_curtis_rings(Nθ, T) == [n * T(π) / (Nθ - 1) for n ∈ 0:Nθ-1]
        @test clenshaw_curtis_rings(Nθ, T)[begin] == 0
        @test clenshaw_curtis_rings(Nθ, T)[end] == T(π)

        for s ∈ -2:2
            nmodes = Ysize(abs(s), ℓmax)
            modelist = Yrange(abs(s), ℓmax)
            # The transform `map2salm` is documented to build for this grid, constructed
            # explicitly.  Comparing to it is metamorphic — the same code path — so it
            # checks the plumbing, not the numbers; the closed-form checks below are the
            # reference.
            𝒯 = SSHTRS(
                s, ℓmax, T; θ=clenshaw_curtis_rings(Nθ, T),
                quadrature_weights=clenshaw_curtis(Nθ, T), Nϕ
            )
            # Coefficients of a random band-limited function, accumulated below from the
            # single-mode samples
            a = randn(rng, Complex{T}, nmodes)
            combined = zeros(Complex{T}, Nϕ, Nθ)

            for (i, (ℓ, m)) ∈ enumerate(modelist)
                f = sampled_sYlm(T, s, ℓ, m, Nϕ, Nθ)
                combined .+= a[i] .* f
                computed = map2salm(f, s, ℓmax)

                # A single ₛYₗₘ decomposes to the unit vector at (ℓ, m), in the canonical
                # ordering that starts at ℓ = |s|
                expected = zeros(Complex{T}, nmodes)
                expected[Yindex(ℓ, m, abs(s))] = one(T)
                @test computed ≈ expected atol=ϵ rtol=ϵ
                @test computed[ℓ, m] ≈ one(T) atol=ϵ

                # The result is a ModeWeights containing the spin and the ℓ range
                @test computed isa ModeWeights{Complex{T}}
                @test spin(computed) == s
                @test ℓₘᵢₙ(computed) == abs(s)
                @test ℓₘₐₓ(computed) == ℓmax
                @test length(computed) == nmodes
            end

            # Linearity: the random combination decomposes to its coefficients, each to a
            # precision set by the size of the function
            c = map2salm(combined, s, ℓmax)
            @test c isa ModeWeights{Complex{T}}
            @test spin(c) == s
            @test maximum(abs, c .- a) ≤ ϵ * norm(a)
            # ... and it is the "RS" transform on this grid applied to the flattened map
            @test c == 𝒯 \ vec(copy(combined))

            # A map with a trailing dimension gives a plain array whose columns are the
            # transforms of the individual maps
            (ℓ₁, m₁) = rand(rng, modelist)
            f₁ = sampled_sYlm(T, s, ℓ₁, m₁, Nϕ, Nθ)
            expected₁ = zeros(Complex{T}, nmodes)
            expected₁[Yindex(ℓ₁, m₁, abs(s))] = one(T)
            F = cat(f₁, combined; dims=3)
            @test size(F) == (Nϕ, Nθ, 2)
            C = map2salm(F, s, ℓmax)
            @test C isa Matrix{Complex{T}}
            @test size(C) == (nmodes, 2)
            @test C[:, 1] ≈ expected₁ atol=ϵ rtol=ϵ
            @test C[:, 1] ≈ map2salm(f₁, s, ℓmax) atol=ϵ rtol=ϵ
            @test maximum(abs, C[:, 2] .- a) ≤ ϵ * norm(a)
            @test maximum(abs, C[:, 2] .- c) ≤ ϵ * norm(a)
            @test C == 𝒯 \ reshape(copy(F), Nϕ * Nθ, 2)
        end
    end
end


@testitem "Transforms: map2salm plan reuse" setup=[Utilities] begin
    import SphericalFunctions: map2salm, map2salm_plan, salm2map, SSHTRS, ModeWeights, Ysize,
        clenshaw_curtis_rings, fejer1, pixels, spin, ℓₘₐₓ
    import .Utilities: array_equal
    using Random

    ℓmax = 7
    Nθ = 2ℓmax + 1
    Nϕ = 2ℓmax + 2
    rng = Random.Xoshiro(2345)

    for T ∈ (Float64, Float32, BigFloat), s ∈ (-2, 0, 1)
        f = randn(rng, Complex{T}, Nϕ, Nθ)
        plan = map2salm_plan(f, s, ℓmax)

        # The plan is an RS transform on the Clenshaw–Curtis rings of the map, with the
        # pixels ordered ring by ring and ϕ varying fastest — the storage order of the map
        @test plan isa SSHTRS{T}
        @test spin(plan) == s
        @test ℓₘₐₓ(plan) == ℓmax
        @test SphericalFunctions.nmodes(plan) == Ysize(abs(s), ℓmax)
        @test SphericalFunctions.npixels(plan) == Nϕ * Nθ
        @test pixels(plan) == [
            [θ, 2T(π) * k / Nϕ] for θ ∈ clenshaw_curtis_rings(Nθ, T) for k ∈ 0:Nϕ-1
        ]

        # Using the plan gives exactly the result of the direct call ...
        direct = map2salm(f, s, ℓmax)
        @test direct isa ModeWeights{Complex{T}}
        @test array_equal(map2salm(f, plan), direct)
        # ... and the plan can be reused for other maps of the same shape
        g = randn(rng, Complex{T}, Nϕ, Nθ)
        @test array_equal(map2salm(g, plan), map2salm(g, s, ℓmax))
        @test array_equal(map2salm(f, plan), direct)
        # ... including maps with trailing dimensions
        F = cat(f, g; dims=3)
        @test array_equal(map2salm(F, plan), map2salm(F, s, ℓmax))

        # The spin and ℓmax are those of the plan
        plan′ = map2salm_plan(f, -s, ℓmax - 1)
        @test array_equal(map2salm(f, plan′), map2salm(f, -s, ℓmax - 1))

        # A plan for a different grid shape is rejected, in either direction
        for (nϕ, nθ) ∈ ((Nϕ + 1, Nθ), (Nϕ, Nθ + 1), (Nθ, Nϕ), (2Nϕ, 2Nθ))
            h = randn(rng, Complex{T}, nϕ, nθ)
            @test_throws DimensionMismatch map2salm(h, plan)
            @test_throws "planned for a different grid" map2salm(h, plan)
            @test_throws DimensionMismatch map2salm(f, map2salm_plan(h, s, ℓmax))
        end

        # So is a plan of the right shape on other rings, or with other quadrature weights.
        # The default `SSHT(s, ℓₘₐₓ)` is the trap: its Fejér grid of 2ℓₘₐₓ+1 rings with
        # 2ℓₘₐₓ+1 points has exactly the shape of the smallest Clenshaw–Curtis grid.
        cc = "works on the Clenshaw–Curtis grid"
        fejér = SSHTRS(s, ℓmax, T; Nϕ)
        @test size(f) == (only(unique(fejér.Nϕ)), length(fejér.θ))
        @test_throws ArgumentError map2salm(f, fejér)
        @test_throws cc map2salm(f, fejér)
        @test_throws cc salm2map(map2salm(f, plan), fejér)
        # (Its constructor warns that one rule's weights on another's rings are inexact.)
        θcc = clenshaw_curtis_rings(Nθ, T)
        wrong_weights = @test_logs (:warn, r"exact") SSHTRS(
            s, ℓmax, T; Nϕ, θ=θcc, quadrature_weights=fejer1(Nθ, T)
        )
        @test_throws cc map2salm(f, wrong_weights)
        @test_throws cc salm2map(map2salm(f, plan), wrong_weights)
    end
end


@testitem "Transforms: salm2map" setup=[Utilities] begin
    import SphericalFunctions: salm2map, map2salm, map2salm_plan, ModeWeights, Yindex, Yrange,
        Ysize, clenshaw_curtis_rings, spin, ℓₘᵢₙ, ℓₘₐₓ
    import .Utilities: array_equal, sYlm_closed_form
    using Random
    using LinearAlgebra: norm

    # The closed-form ₛYₗₘ sampled on the `map2salm` grid, as an Nϕ×Nθ array
    function sampled_sYlm(::Type{T}, s, ℓ, m, Nϕ, Nθ) where {T}
        θs = clenshaw_curtis_rings(Nθ, T)
        ϕs = [2T(π) * k / Nϕ for k ∈ 0:Nϕ-1]
        [sYlm_closed_form(s, ℓ, m, θ, ϕ) for ϕ ∈ ϕs, θ ∈ θs]
    end

    ℓmax = 7
    Nθ = 2ℓmax + 1
    Nϕ = 2ℓmax + 2
    rng = Random.Xoshiro(3456)

    for T ∈ (Float64, Float32, BigFloat)
        ϵ = 100ℓmax * eps(T)

        for s ∈ -2:2
            nmodes = Ysize(abs(s), ℓmax)
            modelist = Yrange(abs(s), ℓmax)
            # A random band-limited function as a ModeWeights, and its values on the grid,
            # accumulated below from the single-mode samples
            w = ModeWeights(randn(rng, Complex{T}, nmodes), s)
            @test spin(w) == s
            @test ℓₘᵢₙ(w) == abs(s)
            @test ℓₘₐₓ(w) == ℓmax
            w_sampled = zeros(Complex{T}, Nϕ, Nθ)

            for (i, (ℓ, m)) ∈ enumerate(modelist)
                Y = sampled_sYlm(T, s, ℓ, m, Nϕ, Nθ)
                w_sampled .+= w[i] .* Y
                salm = zeros(Complex{T}, nmodes)
                salm[Yindex(ℓ, m, abs(s))] = one(T)
                f = salm2map(salm, s, ℓmax, Nϕ, Nθ)
                @test f isa Matrix{Complex{T}}
                @test size(f) == (Nϕ, Nθ)
                # Synthesis of a single mode reproduces the closed-form ₛYₗₘ on the grid ...
                @test f ≈ Y atol=ϵ rtol=ϵ
                # ... and analysis recovers the mode
                @test map2salm(f, s, ℓmax) ≈ salm atol=ϵ rtol=ϵ
            end

            # A ModeWeights input is accepted, and gives the same values as its storage
            fw = salm2map(w, s, ℓmax, Nϕ, Nθ)
            @test fw isa Matrix{Complex{T}}
            @test size(fw) == (Nϕ, Nθ)
            @test array_equal(fw, salm2map(parent(w), s, ℓmax, Nϕ, Nθ))
            @test maximum(abs, fw .- w_sampled) ≤ ϵ * norm(w)
            # The round trip returns a ModeWeights with the same spin
            rt = map2salm(fw, s, ℓmax)
            @test rt isa ModeWeights{Complex{T}}
            @test spin(rt) == s
            @test ℓₘᵢₙ(rt) == abs(s)
            @test ℓₘₐₓ(rt) == ℓmax
            @test maximum(abs, rt .- w) ≤ ϵ * norm(w)
            # Mode weights of the wrong length are rejected.  A ModeWeights knows its range
            # of ℓ, and may cover less than the grid's ℓₘₐₓ, the modes it lacks being zero;
            # but not more.  (With enough points along each ring for the larger ℓmax, so
            # that only the mismatch is at issue.)
            @test_throws DimensionMismatch salm2map(parent(w), s, ℓmax + 1, 2ℓmax + 3, Nθ)
            @test_throws DimensionMismatch salm2map(zeros(Complex{T}, nmodes + 1), s, ℓmax, Nϕ, Nθ)
            let padded = ModeWeights([parent(w); zeros(Complex{T}, 2ℓmax + 3)], s)
                @test salm2map(w, s, ℓmax + 1, 2ℓmax + 3, Nθ) == salm2map(padded, s, ℓmax + 1, 2ℓmax + 3, Nθ)
            end
            @test_throws ArgumentError salm2map(w, s, ℓmax - 1, Nϕ, Nθ)
            @test_throws "synthesizes ℓ ≤ $(ℓmax - 1) only" salm2map(w, s, ℓmax - 1, Nϕ, Nθ)
            # ... and its spin weight and ℓₘₐₓ need not be restated
            @test salm2map(w, Nϕ, Nθ) == fw

            # Trailing dimensions are broadcast over
            salm2 = randn(rng, Complex{T}, nmodes, 2)
            f2 = salm2map(salm2, s, ℓmax, Nϕ, Nθ)
            @test f2 isa Array{Complex{T}, 3}
            @test size(f2) == (Nϕ, Nθ, 2)
            for j ∈ 1:2
                @test f2[:, :, j] ≈ salm2map(salm2[:, j], s, ℓmax, Nϕ, Nθ) atol=ϵ*norm(salm2[:, j]) rtol=ϵ
            end
            rt2 = map2salm(f2, s, ℓmax)
            @test rt2 isa Matrix{Complex{T}}
            @test size(rt2) == (nmodes, 2)
            @test maximum(abs, rt2 .- salm2) ≤ ϵ * norm(salm2)

            # The transform planned by `map2salm_plan` serves synthesis too
            plan = map2salm_plan(fw, s, ℓmax)
            @test array_equal(salm2map(w, plan), fw)
            @test array_equal(salm2map(salm2, plan), f2)
        end
    end

    # Fewer than 2ℓmax+1 points along a ring cannot resolve m = ±ℓmax, so analysis on such a
    # grid warns about the aliasing.  Synthesis is exact on any grid, folding the aliased
    # frequencies correctly, and stays quiet; its values are those of the closed form.
    # (Measured: ≤ 14 eps ‖salm‖ over ten draws.)
    let salm = randn(rng, ComplexF64, Ysize(ℓmax)), Nϕ = 2ℓmax
        f = @test_logs salm2map(salm, 0, ℓmax, Nϕ, Nθ)
        θs, ϕs = clenshaw_curtis_rings(Nθ), [2π * k / Nϕ for k ∈ 0:Nϕ-1]
        g = sum(
            salm[i] .* [sYlm_closed_form(0, ℓ, m, θ, ϕ) for ϕ ∈ ϕs, θ ∈ θs]
            for (i, (ℓ, m)) ∈ enumerate(Yrange(0, ℓmax))
        )
        @test maximum(abs, f - g) ≤ 100eps() * norm(salm)
    end
    @test_logs (:warn, r"fewer than 2ℓₘₐₓ\+1") map2salm(
        zeros(Complex{Float64}, 2ℓmax, Nθ), 0, ℓmax
    )

    # Too few rings make the analysis inexact, and are warned about too, exactly below the
    # number at which it becomes exact: 2⌊ℓₘₐₓ⌋+1 rings.  Measured at ℓₘₐₓ = 7 (for s = 0,
    # ±1 and ±2 alike): exact at Nθ = 15, wrong by 4e-3 at 14 and by 1.7 at 8.  Synthesis is
    # exact on any number of rings, and stays quiet.
    let L = 7
        @test_logs (:warn, r"Nθ=14 rings, but the Clenshaw–Curtis analysis needs at least 15") map2salm(
            zeros(ComplexF64, 2L + 1, 2L), 0, L
        )
        @test_logs map2salm(zeros(ComplexF64, 2L + 1, 2L + 1), 2, L)
        @test_logs salm2map(zeros(ComplexF64, Ysize(L)), 0, L, 2L + 1, 8)
    end
    # For a half-odd ℓₘₐₓ one ring fewer suffices (measured: exact at 2ℓₘₐₓ rings)
    @test_logs map2salm(zeros(ComplexF64, 8, 7), 1//2, 7//2)
    @test_logs (:warn, r"needs at least 7 for ℓₘₐₓ=7//2") map2salm(zeros(ComplexF64, 8, 6), 1//2, 7//2)

    # The grid has two dimensions, and at least two rings, since the Clenshaw–Curtis rings
    # include both poles; anything else is refused with an explanation, rather than by FFTW
    @test_throws DimensionMismatch map2salm(zeros(ComplexF64, 81), 0, 4)
    @test_throws "takes a map of size Nϕ × Nθ" map2salm(zeros(ComplexF64, 81), 0, 4)
    @test_throws DimensionMismatch map2salm(
        zeros(ComplexF64, 81), map2salm_plan(zeros(ComplexF64, 9, 9), 0, 4)
    )
    @test_throws "at least two rings" map2salm(zeros(ComplexF64, 9, 1), 0, 4)
    @test_throws "at least two rings" salm2map(zeros(ComplexF64, Ysize(4)), 0, 4, 9, 1)
    @test_throws ArgumentError salm2map(zeros(ComplexF64, Ysize(4)), 0, 4, 9, 0)
    @test_throws "at least one point" salm2map(zeros(ComplexF64, Ysize(4)), 0, 4, 0, 9)
end

@testitem "Transforms: map2salm and salm2map with half-integer spin weight" begin
    using SphericalFunctions
    using SphericalFunctions: HalfOddInteger
    using Random

    # Measured over 30 draws, in the norm of the error relative to that of the result: at
    # most 0.9ℓₘₐₓ eps, so 100ℓₘₐₓ eps leaves a factor of ≳ 100.
    rng = MersenneTwister(17)
    for (s, ℓₘₐₓ) ∈ ((1//2, 7//2), (-3//2, 9//2))
        ϵ = 100ℓₘₐₓ * eps()
        Nϕ, Nθ = Int(2ℓₘₐₓ + 3), Int(2ℓₘₐₓ + 2)  # at least 2ℓₘₐₓ+1 each
        f̃ = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
        map = salm2map(f̃, s, ℓₘₐₓ, Nϕ, Nθ)
        @test size(map) == (Nϕ, Nθ)
        g̃ = map2salm(map, s, ℓₘₐₓ)
        @test g̃ isa ModeWeights && spin(g̃) === HalfOddInteger(s)
        @test parent(g̃) ≈ parent(f̃) rtol=ϵ
        # The map holds the function's values at the plan's rotors
        plan = SphericalFunctions.map2salm_plan(map, s, ℓₘₐₓ)
        @test vec(map) ≈ sYlm_matrix(rotors(plan), ℓₘₐₓ, s) * parent(f̃) rtol=ϵ
        # A stack of maps is transformed map by map, each exactly as it would be alone (and
        # doubling a map doubles every rounded result exactly)
        G̃ = map2salm(cat(map, 2map; dims=3), s, ℓₘₐₓ)
        @test G̃[:, 1] == parent(g̃) && G̃[:, 2] == 2 .* parent(g̃)
    end
end


@testitem "Transforms: map2salm and salm2map in Float16" begin
    import SphericalFunctions: map2salm, salm2map, ModeWeights, Ysize, array_view, spin
    import Random

    # FFTW has no Float16 transforms, so the ring transforms go through GenericFFT, and the
    # quadrature weights, computed through FFTW in Float32, are returned in Float16.  A
    # round trip on a grid of 2ℓₘₐₓ+1 rings recovers the mode weights; the error measured
    # over these cases is at most 0.71ℓₘₐₓ eps.
    rng = Random.Xoshiro(20260924)
    for s ∈ (0, 1, -2), ℓₘₐₓ ∈ (2, 4, 8)
        f̃ = ModeWeights(ComplexF16.(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)) ./ 4), s)
        f = salm2map(f̃, 2ℓₘₐₓ + 1, 2ℓₘₐₓ + 1)
        @test f isa Matrix{ComplexF16}
        g̃ = map2salm(f, s, ℓₘₐₓ)
        @test eltype(g̃) === ComplexF16 && spin(g̃) == s
        @test maximum(abs, ComplexF64.(array_view(g̃)) .- ComplexF64.(array_view(f̃))) ≤
            10ℓₘₐₓ * eps(Float16)
    end
end
