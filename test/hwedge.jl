@testitem "HWedge" begin
    using SphericalFunctions.Redesign: HWedge, HWedge_size, ℓₘᵢₙ

    ## First, I just need a way of storing the various index combinations
    ## in a single integer — so that it can go into an HWedge, for example — for testing purposes.
    encoder(i::Int) = i
    encoder(i::Rational) = numerator(i)
    function encoder(iᵣ, m′, m)
        (
            encoder(abs(m))
            + (m < 0) * 10
            + encoder(abs(m′)) * 10^2
            + (m′ < 0) * 10^3
            + iᵣ * 10^4
        )
    end
    function decoder(i)
        d = digits(i)
        m = d[1] * (d[2] == 1 ? -1 : 1)
        m′ = d[3] * (d[4] == 1 ? -1 : 1)
        iᵣ = d[5]
        iᵣ, m′, m
    end
    # Test encoder/decoder
    for ℓₘₐₓ ∈ (5, 7//2)
        for iᵣ ∈ 1:7
            for m′ ∈ -ℓₘₐₓ:ℓₘₐₓ
                for m ∈ -ℓₘₐₓ:ℓₘₐₓ
                    code = encoder(iᵣ, m′, m)
                    (iᵣ₂, m′₂, m₂) = decoder(code)
                    @test iᵣ == iᵣ₂
                    if ℓₘₐₓ isa Int
                        @test m′ == m′₂
                        @test m == m₂
                    else
                        @test m′ == m′₂ // 2
                        @test m == m₂ // 2
                    end
                end
            end
        end
    end

    for ℓₘₐₓ ∈ (5, 9//2)
        for Nᵣ ∈ (1, 2, 3, 7)
            for m′ₘₐₓ ∈ ℓₘₐₓ:-1:ℓₘᵢₙ(ℓₘₐₓ)
                for m′ₘᵢₙ in -ℓₘₐₓ:-ℓₘᵢₙ(ℓₘₐₓ)
                    RT = Float64
                    H = HWedge(RT, Nᵣ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)

                    @test H.Nᵣ == Nᵣ
                    @test H.maxℓ == ℓₘₐₓ
                    @test H.maxm′ₘₐₓ == m′ₘₐₓ
                    @test H.minm′ₘᵢₙ == m′ₘᵢₙ

                    # When first created, the indices should have their smallest values
                    @test H.ℓ == ℓₘᵢₙ(ℓₘₐₓ)
                    @test H.m′ₘₐₓ == H.ℓ
                    @test H.m′ₘᵢₙ == -H.ℓ

                    # But the storage should be full size
                    expected_size = HWedge_size(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
                    @test size(H.parent) == (Nᵣ * expected_size, )

                    # Test changing ℓ
                    for new_ell in (ℓₘᵢₙ(ℓₘₐₓ):ℓₘₐₓ)
                        H.ℓ = new_ell
                        @test H.ℓ == new_ell
                        @test H.m′ₘₐₓ == min(new_ell, m′ₘₐₓ)
                        @test H.m′ₘᵢₙ == max(-new_ell, m′ₘᵢₙ)
                    end
                end
            end
        end
    end
end
