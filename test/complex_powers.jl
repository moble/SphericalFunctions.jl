@testitem "Complex powers" setup=[Utilities] begin

    complex_powers_comparison(z, m, T=Float64) = (
        complex_powers(Complex{T}(z), m),
        z.^collect(0:m),
        eps(T)*m
    )

    for T in [Float64, Float32, Float16]
        nozpowers = Vector{Complex{T}}(undef, 0)
        fudge = one(T) + 2 * sqrt(eps(T))
        for k in 0:25
            z = cis(k*big(π)/10)
            for m in [0, 1, 2, 3, 4, 1_000]
                mine, theirs, ϵ = complex_powers_comparison(z, m, T)
                @test mine ≈ theirs rtol=2ϵ
                @test_throws DomainError complex_powers(Complex{T}(z*fudge), m)
                @test_throws DomainError complex_powers(Complex{T}(z/fudge), m)
                inplace = zeros(Complex{T}, size(mine))
                complex_powers!(inplace, Complex{T}(z))
                @test array_equal(mine, inplace)
            end
            @test length(complex_powers!(nozpowers, Complex{T}(z))) == 0
        end
    end

end


@testitem "Complex powers: accuracy of the fused modulus" begin
    using SphericalFunctions: complex_powers

    # `complex_powers!` computes `modulus` as `√(muladd(z.re, z.re, z.im*z.im))` rather than
    # `√(abs2(z))`.  A one-ulp error there feeds the cancellation-sensitive `dc`, and the
    # recurrence amplifies it linearly in `m`; the spelling therefore matters far more than
    # it looks like it should.  Measured at ϕ = 0.3: with the fused form the error is
    # 5.5e-16 / 2.6e-15 / 7.7e-15 at m = 128 / 1024 / 4096, against 2.9e-15 / 2.2e-14 /
    # 9.1e-14 for `√(abs2(z))`.  The thresholds below sit between the two, so this test
    # genuinely discriminates.  (This is also why the function no longer uses `@fastmath`,
    # which used to supply the same accuracy non-deterministically.)
    setprecision(BigFloat, 512) do
        ϕ = 0.3
        for (m, tol) in ((128, 1.0e-15), (1024, 6.0e-15), (4096, 2.0e-14))
            computed = complex_powers(cis(ϕ), m)
            exact = [Complex{BigFloat}(cis(BigFloat(ϕ) * k)) for k in 0:m]
            @test maximum(abs.(Complex{BigFloat}.(computed) .- exact)) < tol
        end
    end
end
