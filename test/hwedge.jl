@testitem "HWedge" setup=[EncodeDecode] begin
    using SphericalFunctions.Redesign: HWedge, HWedge_size, ℓₘᵢₙ
    using .EncodeDecode: encode, decode

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
