@testset verbose=true "map2salm" begin
    # NOTE: Float16 irfft returns Float32, which leads to type conflicts, so we just don't support
    # map2salm on Float16

    function sYlm(s::Int, ell::Int, m::Int, theta::T, phi::T) where {T<:Real}
        # Eqs. (II.7) and (II.8) of https://arxiv.org/abs/0709.0093v3
        # Note their weird definition w.r.t. `-s`
        k_min = max(0, m + s)
        k_max = min(ell + m, ell + s)
        sin_half_theta, cos_half_theta = sincos(theta / 2)
        return (-1)^(-s) * sqrt((2 * ell + 1) / (4 * T(π))) *
            T(sum(
                (-1) ^ (k)
                * sqrt(factorial(big(ell + m)) * factorial(big(ell - m)) * factorial(big(ell - s)) * factorial(big(ell + s)))
                * (cos_half_theta ^ (2 * ell + m + s - 2 * k))
                * (sin_half_theta ^ (2 * k - s - m))
                / (factorial(big(ell + m - k)) * factorial(big(ell + s - k)) * factorial(big(k)) * factorial(big(k - s - m)))
                for k in k_min:k_max
            )) *
            cis(m * phi)
    end

    m_m2Y2m = [
        (2, NINJA.m2Y22),
        (1, NINJA.m2Y21),
        (0, NINJA.m2Y20),
        (-1, NINJA.m2Y2m1),
        (-2, NINJA.m2Y2m2)
    ]

    @testset "Input expressions $T" for T in [BigFloat, Float64, Float32]
        # These are just internal consistency tests of the sYlm function
        # above, against the explicit expressions `mY2.`
        s = -2
        ℓ = 2
        Nϑ = 17
        Nφ = 18
        for (m, m2Y2m) in m_m2Y2m
            f1 = mapslices(ϕθ -> sYlm(s, ℓ, m, ϕθ[2], ϕθ[1]), phi_theta(Nφ, Nϑ, T), dims=[3])
            f2 = mapslices(ϕθ -> m2Y2m(ϕθ[2], ϕθ[1]), phi_theta(Nφ, Nϑ, T), dims=[3])
            @test f1 ≈ f2 atol=10eps(T) rtol=10eps(T)
        end
    end

    @testset "map2salm $T" for T in [BigFloat, Float64, Float32]
        # These test the ability of map2salm to precisely decompose the results of `sYlm`.
        ℓmax = 7
        Nϑ = 2ℓmax + 1
        Nφ = 2ℓmax + 2
        for s in -2:2
            #for ℓmin in 0:abs(s)
            let ℓmin = 0
                for ℓ in abs(s):ℓmax
                    for m in -ℓ:ℓ
                        f = mapslices(
                            ϕθ -> sYlm(s, ℓ, m, ϕθ[2], ϕθ[1]),
                            phi_theta(Nφ, Nϑ, T),
                            dims=[3]
                        )
                        computed = map2salm(f, s, ℓmax)
                        expected = zeros(Complex{T}, size(computed))
                        expected[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
                        if ≉(computed, expected, atol=30eps(T), rtol=30eps(T))
                            @show T
                            @show ℓmax
                            @show Nϑ
                            @show Nφ
                            @show s
                            @show ℓmin
                            @show ℓ
                            @show m
                            println("computed = $computed")
                            println("expected = $expected")
                            println("max_diff = ", maximum(abs, computed .- expected), ";")
                            println()
                        end
                        @test computed ≈ expected atol=30eps(T) rtol=30eps(T)

                        plan = plan_map2salm(f, s, ℓmax)
                        computed2 = map2salm(f, plan)
                        @test array_equal(computed, computed2)
                    end
                end
            end
        end
    end

end
