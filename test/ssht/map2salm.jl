# Tests of the equiangular-grid convenience functions `map2salm` and `salm2map` from
# `src/ssht/rs.jl`, ported from the v2 `test/deprecated/map2salm.jl` (deleted in 3.0) to the
# v3 API.
#
# The grid is the one documented for `map2salm`: an `Nϕ × Nθ` array with `ϕₖ = 2πk/Nϕ` along
# the first dimension and `θ = clenshaw_curtis_rings(Nθ)` (both poles included) along the
# second.  The oracle is the closed-form `sYlm(s, ℓ, m, θ, ϕ)` from the `Utilities` snippet
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
# Float16 is not tested, because FFTW has no Float16 transforms (the v2 tests skipped it for
# the same reason).
#
# The item names use a "Transforms: " prefix so that they are distinguishable — to a human
# reading a results list, and to `runtests.jl`'s `occursin` filters — from the v2 items,
# which were also named `map2salm`.  Those items are gone with the `Deprecated` module, but
# the prefix is kept: the names are what people filter on.

@testitem "Transforms: map2salm" setup=[Utilities] begin
    import SphericalFunctions: map2salm, SSHTRS, ModeWeights, Yindex, Yrange, Ysize,
        clenshaw_curtis_rings, clenshaw_curtis, spin, ℓₘᵢₙ, ℓₘₐₓ
    using Random
    using LinearAlgebra: norm

    # The closed-form ₛYₗₘ sampled on the `map2salm` grid, as an Nϕ×Nθ array
    function sampled_sYlm(::Type{T}, s, ℓ, m, Nϕ, Nθ) where {T}
        θs = clenshaw_curtis_rings(Nθ, T)
        ϕs = [2T(π) * k / Nϕ for k ∈ 0:Nϕ-1]
        [sYlm(s, ℓ, m, θ, ϕ) for ϕ ∈ ϕs, θ ∈ θs]
    end

    ℓmax = 7
    Nθ = 2ℓmax + 1
    Nϕ = 2ℓmax + 2
    rng = Random.Xoshiro(1234)

    # Measured over every T and s ∈ -2:2, in units of eps(T): single-mode analysis ≤ 10.9
    # (worst in BigFloat; 5.5 in Float64, 8.0 in Float32) and linearity ≤ 10.0 eps(T)‖a‖.
    # 100 eps(T) therefore leaves a factor of ≳ 9.  (The 30 eps left over unexamined
    # from the v2 test left only 2.8 — the tightest margin anywhere in these files, and not
    # enough to be stable across platforms and BLAS versions.)
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
            # explicitly.  Comparing to it is metamorphic — the same code path — so it checks
            # the plumbing, not the numbers; the closed-form checks below are the reference.
            𝒯 = SSHTRS(
                s, ℓmax; T, θ=clenshaw_curtis_rings(Nθ, T),
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
    import SphericalFunctions: map2salm, map2salm_plan, SSHTRS, ModeWeights, Ysize,
        clenshaw_curtis_rings, pixels, spin, ℓₘₐₓ
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
            @test_throws ErrorException map2salm(h, plan)
            @test_throws "planned for a different grid" map2salm(h, plan)
            @test_throws ErrorException map2salm(f, map2salm_plan(h, s, ℓmax))
        end
    end
end


@testitem "Transforms: salm2map" setup=[Utilities] begin
    import SphericalFunctions: salm2map, map2salm, map2salm_plan, ModeWeights, Yindex, Yrange,
        Ysize, clenshaw_curtis_rings, spin, ℓₘᵢₙ, ℓₘₐₓ
    using Random
    using LinearAlgebra: norm

    # The closed-form ₛYₗₘ sampled on the `map2salm` grid, as an Nϕ×Nθ array
    function sampled_sYlm(::Type{T}, s, ℓ, m, Nϕ, Nθ) where {T}
        θs = clenshaw_curtis_rings(Nθ, T)
        ϕs = [2T(π) * k / Nϕ for k ∈ 0:Nϕ-1]
        [sYlm(s, ℓ, m, θ, ϕ) for ϕ ∈ ϕs, θ ∈ θs]
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
            # Mode weights of the wrong length or ℓ range are rejected (with enough points along
            # each ring for the larger ℓmax, so that only the mismatch is reported)
            @test_throws ErrorException salm2map(w, s, ℓmax + 1, 2ℓmax + 3, Nθ)
            @test_throws ErrorException salm2map(parent(w), s, ℓmax + 1, 2ℓmax + 3, Nθ)
            @test_throws ErrorException salm2map(zeros(Complex{T}, nmodes + 1), s, ℓmax, Nϕ, Nθ)

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

    # Fewer than 2ℓmax+1 points along a ring cannot resolve m = ±ℓmax, and both directions
    # warn about the aliasing
    @test_logs (:warn, r"fewer than 2ℓₘₐₓ\+1") salm2map(
        zeros(Complex{Float64}, Ysize(ℓmax)), 0, ℓmax, 2ℓmax, Nθ
    )
    @test_logs (:warn, r"fewer than 2ℓₘₐₓ\+1") map2salm(
        zeros(Complex{Float64}, 2ℓmax, Nθ), 0, ℓmax
    )
end

@testitem "Transforms: map2salm and salm2map with half-integer spin weight" begin
    using SphericalFunctions
    using SphericalFunctions: HalfOddInteger
    using Random

    rng = MersenneTwister(17)
    for (s, ℓₘₐₓ) ∈ ((1//2, 7//2), (-3//2, 9//2))
        Nϕ, Nθ = Int(2ℓₘₐₓ + 3), Int(2ℓₘₐₓ + 2)  # at least 2ℓₘₐₓ+1 each
        f̃ = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
        map = salm2map(f̃, s, ℓₘₐₓ, Nϕ, Nθ)
        @test size(map) == (Nϕ, Nθ)
        g̃ = map2salm(map, s, ℓₘₐₓ)
        @test g̃ isa ModeWeights && spin(g̃) === HalfOddInteger(s)
        @test parent(g̃) ≈ parent(f̃) rtol=1e-11
        # The map holds the function's values at the plan's rotors
        plan = SphericalFunctions.map2salm_plan(map, s, ℓₘₐₓ)
        @test vec(map) ≈ sYlm_matrix(rotors(plan), ℓₘₐₓ, s) * parent(f̃) rtol=1e-12
        # A stack of maps is transformed map by map
        G̃ = map2salm(cat(map, 2map; dims=3), s, ℓₘₐₓ)
        @test G̃[:, 1] ≈ parent(f̃) && G̃[:, 2] ≈ 2 .* parent(f̃)
    end
end
