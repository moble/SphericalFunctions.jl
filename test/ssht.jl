@testset verbose=true "SSHT" begin

    FloatTypes = [Double64, Float64, Float32]
    methods = ["Direct", "Minimal", "RS"]
    inplacemethods = ["Direct", "Minimal"]
    cases = Iterators.product(methods, FloatTypes)
    inplacecases = Iterators.product(inplacemethods, FloatTypes)

    function sYlm(s, ℓ, m, θϕ)
        NINJA.sYlm(s, ℓ, m, θϕ[1], θϕ[2])
    end

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

    @testset "Preliminaries" begin

        # Preliminary check that `sqrtbinomial` works as expected
        @testset "sqrtbinomial" for T ∈ [Float16, Float32, Float64, Double64, BigFloat]
            for ℓ ∈ [1, 2, 3, 4, 5, 13, 64, 1025]
                for s ∈ -2:2
                    # Note that `ℓ-abs(s)` is more relevant, but we test without `abs` here
                    # to exercise more code paths
                    a = SphericalFunctions.sqrtbinomial(2ℓ, ℓ-s, T)
                    b = T(√binomial(big(2ℓ), big(ℓ-s)))
                    @test a ≈ b
                end
            end
        end

        # Check that an error results from a nonsense method request
        @testset "Nonsense method" begin
            let s=-2, ℓmax=8
                @test_throws ErrorException SSHT(s, ℓmax; method="NonsenseGarbage")
            end
        end

        # Check what `show` looks like
        @testset "SSHT show" begin
            let io=IOBuffer(), s=-2, ℓmax=8, T=Float64, method="Direct"
                for inplace ∈ [true, false]
                    expected = "SphericalFunctions.SSHT$method{$T, $inplace}($s, $ℓmax)"
                    𝒯 = SSHT(s, ℓmax; T=T, method=method, inplace=inplace)
                    Base.show(io, MIME("text/plain"), 𝒯)
                    @test String(take!(io)) == expected
                end
            end
        end

        # Check that SSHTDirect warns if ℓₘₐₓ is too large
        @testset "Direct ℓₘₐₓ" begin
            let s=0, ℓₘₐₓ=65
                @test_warn """ "Direct" method for s-SHT is only """ SSHT(s, ℓₘₐₓ; method="Direct")
            end
        end

        # Check that SSHTDirect warns if `check_blas_threads` is too low
        @testset "Direct threads" begin
            let cores=num_physical_cores(), blas_threads=LinearAlgebra.BLAS.get_num_threads()
                if cores > 1
                    LinearAlgebra.BLAS.set_num_threads(1)
                    try
                        @test_warn """ all available threads """ SSHT(0, 5; method="Direct")
                    finally
                        LinearAlgebra.BLAS.set_num_threads(blas_threads)
                    end
                end
            end
        end

        # Check pixels and rotors of Minimal
        @testset "Minimal pixels" for T ∈ FloatTypes
            for ℓmax ∈ [3, 4, 5, 13, 64]
                for s ∈ -min(2,abs(ℓmax)-1):min(2,abs(ℓmax)-1)
                    𝒯 = SSHT(s, ℓmax; T=T, method="Minimal")
                    @test pixels(𝒯) ≈ sorted_ring_pixels(s, ℓmax, T)
                    @test rotors(𝒯) ≈ sorted_ring_rotors(s, ℓmax, T)
                end
            end
        end

    end  # Preliminaries

    # These test the ability of ssht to precisely reconstruct a pure `sYlm`.
    @testset "Synthesis: $T $method" for (method, T) in cases

        # We can't go to very high ℓ, because NINJA.sYlm fails for low-precision numbers
        for ℓmax ∈ 3:7

            # We need ϵ to be huge, seemingly mostly due to the low-precision method
            # used for NINJA.sYlm; it is used because it is a simple reference method.
            ϵ = 500ℓmax^3 * eps(T)

            for s in -2:2
                𝒯 = SSHT(s, ℓmax; T=T, method=method)

                #for ℓmin in 0:abs(s)
                let ℓmin = abs(s)
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f = zeros(Complex{T}, SphericalFunctions.Ysize(ℓmin, ℓmax))
                            f[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
                            computed = 𝒯 * f
                            expected = sYlm.(s, ℓ, m, pixels(𝒯))
                            explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end  # Synthesis


    # These test the ability of ssht to precisely decompose the results of `sYlm`.
    @testset "Analysis: $T $method" for (method, T) in cases

        # We can't go to very high ℓ, because NINJA.sYlm fails for low-precision numbers
        for ℓmax ∈ 3:7

            # We need ϵ to be huge, seemingly mostly due to the low-precision method
            # used for NINJA.sYlm; it is used because it is a simple reference method.
            ϵ = 500ℓmax^3 * eps(T)
            if method == "Minimal"
                ϵ *= 50
            end

            for s in -2:2
                𝒯 = SSHT(s, ℓmax; T=T, method=method)
                let ℓmin = abs(s)
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f = sYlm.(s, ℓ, m, pixels(𝒯))
                            computed = 𝒯 \ f
                            expected = zeros(Complex{T}, size(computed))
                            expected[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
                            explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end  # Analysis

    # These test the ability of ssht to precisely reconstruct a pure `sYlm`,
    # and then reverse that process to find the pure mode again.
    @testset verbose=false "A ∘ S: $T $method" for (method, T) in cases
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
                𝒯 = SSHT(s, ℓmax; T=T, method=method)
                let ℓmin = abs(s)
                    f = zeros(Complex{T}, SphericalFunctions.Ysize(ℓmin, ℓmax))
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f[:] .= false
                            f[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
                            expected = copy(f)
                            computed = 𝒯 \ (𝒯 * f)
                            explain(computed, expected, method, T, ℓmax, s, ℓ, m, ϵ)
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end  # A ∘ S

    # These test A ∘ S in the RS method when using different quadratures
    @testset verbose=false "RS quadratures: $T" for T in FloatTypes
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
                    𝒯 = SSHT(s, ℓmax; T=T, θ=θ, quadrature_weights=w, method="RS")
                    p1 = [
                        @SVector [θi, ϕi]
                        for θi ∈ θ
                        for ϕi ∈ LinRange(T(0), 2T(π), 2ℓmax+2)[begin:end-1]
                    ]
                    p2 = pixels(𝒯)
                    @test p1 ≈ p2
                    r1 = [from_spherical_coordinates(θϕ...) for θϕ ∈ p1]
                    r2 = rotors(𝒯)
                    @test r1 ≈ r2
                    let ℓmin = abs(s)
                        f = zeros(Complex{T}, SphericalFunctions.Ysize(ℓmin, ℓmax))
                        for ℓ in abs(s):ℓmax
                            for m in -ℓ:ℓ
                                f[:] .= false
                                f[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
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
    end  # RS quadratures

    # These test that the non-inplace versions of transformers that *can* work in place
    # still work.
    @testset verbose=false "Non-inplace: $T $method" for (method, T) in inplacecases
        @testset "$ℓmax" for ℓmax ∈ [4,5]
            #ϵ = 20ℓmax^2 * eps(T)
            ϵ = 100ℓmax^3 * eps(T)
            if method == "Minimal"
                ϵ *= 50
            end
            for s in [-1, 1]
                𝒯 = SSHT(s, ℓmax; T=T, method=method, inplace=false)
                let ℓmin = abs(s)
                    f = zeros(Complex{T}, SphericalFunctions.Ysize(ℓmin, ℓmax))
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f[:] .= false
                            f[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
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
    end  # Non-inplace

end
