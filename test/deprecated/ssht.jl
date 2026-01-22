@testsnippet SSHT begin
    using DoubleFloats
    FloatTypes = [Double64, Float64, Float32]
    methods = ["Direct", "Minimal", "RS"]
    inplacemethods = ["Direct", "Minimal"]
    cases = Iterators.product(methods, FloatTypes)
    inplacecases = Iterators.product(inplacemethods, FloatTypes)

    function explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
        if ≉(computed, expected, atol=ϵ, rtol=ϵ)
            @show method T ℓmax s ℓ m ϵ
            comp = copy(computed)
            @. comp[abs(comp)<ϵ]=0
            @show comp expected
            println("max_diff = ", maximum(abs, computed .- expected), ";")
            println()
            #error("")
        end
    end
end

# Preliminary check that `sqrtbinomial` works as expected
@testitem "Preliminaries: sqrtbinomial" begin
    import SphericalFunctions: Deprecated
    using DoubleFloats
    for T ∈ [Float16, Float32, Float64, Double64, BigFloat]
        for ℓ ∈ [1, 2, 3, 4, 5, 13, 64, 1025]
            for s ∈ -2:2
                # Note that `ℓ-abs(s)` is more relevant, but we test without `abs` here
                # to exercise more code paths
                a = Deprecated.sqrtbinomial(2ℓ, ℓ-s, T)
                b = T(√binomial(big(2ℓ), big(ℓ-s)))
                @test a ≈ b
            end
        end
    end
end

# Check that an error results from a nonsense method request
@testitem "Preliminaries: Nonsense method" begin
    let s=-2, ℓmax=8
        import SphericalFunctions: Deprecated
        @test_throws ErrorException Deprecated.SSHT(s, ℓmax; method="NonsenseGarbage")
    end
end

# Check what `show` looks like
@testitem "Preliminaries: SSHT show" begin
    let io=IOBuffer(), s=-2, ℓmax=8, T=Float64, method="Direct"
        import SphericalFunctions: Deprecated
        TD = "LinearAlgebra.LU{ComplexF64, Matrix{ComplexF64}, Vector{Int64}}"
        for inplace ∈ [true, false]
            expected = "SphericalFunctions.Deprecated.SSHT$method{$T, $inplace, $TD}($s, $ℓmax)"
            𝒯 = Deprecated.SSHT(s, ℓmax; T, method, inplace)
            Base.show(io, MIME("text/plain"), 𝒯)
            @test String(take!(io)) == expected
        end
    end
end

# Check that SSHTDirect warns if ℓₘₐₓ is too large
@testitem "Preliminaries: Direct ℓₘₐₓ" begin
    let s=0, ℓₘₐₓ=65
        import SphericalFunctions: Deprecated
        @test_warn """ "Direct" method for s-SHT is only """ Deprecated.SSHT(s, ℓₘₐₓ; method="Direct")
    end
end

# Check pixels and rotors of Minimal
@testitem "Preliminaries: Minimal pixels" setup=[SSHT] begin
    import SphericalFunctions: Deprecated
    for T ∈ FloatTypes
        for ℓmax ∈ [3, 4, 5, 13, 64]
            for s ∈ -min(2,abs(ℓmax)-1):min(2,abs(ℓmax)-1)
                𝒯 = Deprecated.SSHT(s, ℓmax; T=T, method="Minimal")
                @test Deprecated.pixels(𝒯) ≈ sorted_ring_pixels(s, ℓmax, T)
                @test Deprecated.rotors(𝒯) ≈ sorted_ring_rotors(s, ℓmax, T)
            end
        end
    end
end


# These test the ability of ssht to precisely reconstruct a pure `sYlm`.
@testitem "Synthesis" setup=[SSHT] begin
    import SphericalFunctions: Deprecated
    for (method, T) in cases

        for ℓmax ∈ 3:7

            # This was huge because we used to use NINJA expressions, which were
            # low-accuracy; we can probably reduce this now.
            ϵ = 500ℓmax^3 * eps(T)

            for s in -2:2
                𝒯 = Deprecated.SSHT(s, ℓmax; T=T, method=method)

                #for ℓmin in 0:abs(s)
                let ℓmin = abs(s)
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f = zeros(Complex{T}, Deprecated.Ysize(ℓmin, ℓmax))
                            f[Deprecated.Yindex(ℓ, m, ℓmin)] = one(T)
                            computed = 𝒯 * f
                            expected = Deprecated.Y.(s, ℓ, m, Deprecated.pixels(𝒯))
                            explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end
end


# These test the ability of ssht to precisely decompose the results of `sYlm`.
@testitem "Analysis" setup=[SSHT] begin
    import SphericalFunctions: Deprecated
    for (method, T) in cases

        for ℓmax ∈ 3:7

            # This was huge because we used to use NINJA expressions, which were
            # low-accuracy; we can probably reduce this now.
            ϵ = 500ℓmax^3 * eps(T)
            if method == "Minimal"
                ϵ *= 50
            end

            for s in -2:2
                𝒯 = Deprecated.SSHT(s, ℓmax; T=T, method=method)
                let ℓmin = abs(s)
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f = Deprecated.Y.(s, ℓ, m, Deprecated.pixels(𝒯))
                            computed = 𝒯 \ f
                            expected = zeros(Complex{T}, size(computed))
                            expected[Deprecated.Yindex(ℓ, m, ℓmin)] = one(T)
                            explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end
end

# These test the ability of ssht to precisely reconstruct a pure `sYlm`,
# and then reverse that process to find the pure mode again.
@testitem "A ∘ S" setup=[SSHT] begin
    import SphericalFunctions: Deprecated
    for (method, T) in cases
        # Note that the number of tests here scales as ℓmax^2, and
        # the time needed for each scales as (ℓmax log(ℓmax))^2,
        # so we don't bother going to very high ℓmax.
        @testset "$ℓmax" for ℓmax ∈ 3:7
            #ϵ = 20ℓmax^2 * eps(T)
            ϵ = 500ℓmax^3 * eps(T)
            if method == "Minimal"
                ϵ *= 50
            end
            for s in -2:2
                𝒯 = Deprecated.SSHT(s, ℓmax; T=T, method=method)
                let ℓmin = abs(s)
                    f = zeros(Complex{T}, Deprecated.Ysize(ℓmin, ℓmax))
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f[:] .= false
                            f[Deprecated.Yindex(ℓ, m, ℓmin)] = one(T)
                            expected = copy(f)
                            computed = 𝒯 \ (𝒯 * f)
                            explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end
end

# These test A ∘ S in the RS method when using different quadratures
@testitem "RS quadratures" setup=[SSHT] begin
    import SphericalFunctions: Deprecated
    using StaticArrays
    using Quaternionic
    for T in FloatTypes
        method = "RS"
        @testset "$ℓmax" for ℓmax ∈ 3:7
            #ϵ = 20ℓmax^2 * eps(T)
            ϵ = 500ℓmax^3 * eps(T)
            for s in -2:2
                for (θ, w) ∈ [
                    (fejer1_rings(2ℓmax+1, T), fejer1(2ℓmax+1, T)),
                    (fejer2_rings(2ℓmax+1, T), fejer2(2ℓmax+1, T)),
                    (clenshaw_curtis_rings(2ℓmax+1, T), clenshaw_curtis(2ℓmax+1, T))
                ]
                    𝒯 = Deprecated.SSHT(s, ℓmax; T=T, θ=θ, quadrature_weights=w, method="RS")
                    p1 = [
                        @SVector [θi, ϕi]
                        for θi ∈ θ
                        for ϕi ∈ LinRange(T(0), 2T(π), 2ℓmax+2)[begin:end-1]
                    ]
                    p2 = Deprecated.pixels(𝒯)
                    @test p1 ≈ p2
                    r1 = [from_spherical_coordinates(θϕ...) for θϕ ∈ p1]
                    r2 = Deprecated.rotors(𝒯)
                    @test r1 ≈ r2
                    let ℓmin = abs(s)
                        f = zeros(Complex{T}, Deprecated.Ysize(ℓmin, ℓmax))
                        for ℓ in abs(s):ℓmax
                            for m in -ℓ:ℓ
                                f[:] .= false
                                f[Deprecated.Yindex(ℓ, m, ℓmin)] = one(T)
                                expected = copy(f)
                                computed = 𝒯 \ (𝒯 * f)
                                explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                                @test computed ≈ expected atol=ϵ rtol=ϵ
                            end
                        end
                    end
                end
            end
        end
    end
end

# These test that the non-inplace versions of transformers that *can* work in place
# still work.
@testitem "Non-inplace" setup=[SSHT] begin
    import SphericalFunctions: Deprecated
    using LinearAlgebra
    @testset verbose=false "Non-inplace: $T $method" for (method, T) in inplacecases
        @testset "$ℓmax" for ℓmax ∈ [4,5]
            #ϵ = 20ℓmax^2 * eps(T)
            ϵ = 100ℓmax^3 * eps(T)
            if method == "Minimal"
                ϵ *= 50
            end
            for s in [-1, 1]
                𝒯 = Deprecated.SSHT(s, ℓmax; T=T, method=method, inplace=false)
                let ℓmin = abs(s)
                    f = zeros(Complex{T}, Deprecated.Ysize(ℓmin, ℓmax))
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f[:] .= false
                            f[Deprecated.Yindex(ℓ, m, ℓmin)] = one(T)
                            expected = f
                            f′ = similar(f)
                            LinearAlgebra.mul!(f′, 𝒯, f)
                            f′′ = 𝒯 * copy(f)
                            @test f′′ == f′
                            LinearAlgebra.ldiv!(f′′, 𝒯, f′)
                            computed = 𝒯 \ copy(f′)
                            @test f′′ == computed
                            explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end
end
