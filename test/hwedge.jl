@testitem "HWedge" setup=[EncodeDecode] begin
    using SphericalFunctions: HWedge, HWedge_size, Nᵣ, ℓ, ℓₘᵢₙ, m′ₘᵢₙ, m′ₘₐₓ
    using .EncodeDecode: encode, decode

    # We will fill the HWedge with integers that encode their indices.  By iterating over
    # valid indices in order, we can test that the storage layout is correct.  Specifically,
    # we want an inner loop over iᵣ, then m, then m′ because this is the order in which the
    # recurrence relations will fill the data, vectorizing over iᵣ, then iterating most
    # quickly over m.  We will check that we can fill that data both using 1D indexing and
    # 3D indexing, then verify that both methods give the same result.  We will also test
    # both methods for reading the data back out, which will verify that the indexing logic
    # is correct in both setindex! and getindex.
    function fill_1index!(w::HWedge{IT}) where {IT}
        let Nᵣ = Nᵣ(w), ℓ = ℓ(w), m′ₘₐₓ = m′ₘₐₓ(w), m′ₘᵢₙ = m′ₘᵢₙ(w)
            i = 1
            for m′ ∈ m′ₘᵢₙ:m′ₘₐₓ
                for m ∈ abs(m′):ℓ
                    for iᵣ ∈ 1:Nᵣ
                        w[i] = encode(iᵣ, m′, m)
                        i += 1
                    end
                end
            end
        end
        return w
    end
    function fill_3index!(w::HWedge{IT}) where {IT}
        let Nᵣ = Nᵣ(w), ℓ = ℓ(w), m′ₘₐₓ = m′ₘₐₓ(w), m′ₘᵢₙ = m′ₘᵢₙ(w)
            for m′ ∈ m′ₘᵢₙ:m′ₘₐₓ
                for m ∈ abs(m′):ℓ
                    for iᵣ ∈ 1:Nᵣ
                        w[iᵣ, m′, m] = encode(iᵣ, m′, m)
                    end
                end
            end
        end
        return w
    end
    function test_1index(w::HWedge{IT}) where {IT}
        let Nᵣ = Nᵣ(w), ℓ = ℓ(w), m′ₘₐₓ = m′ₘₐₓ(w), m′ₘᵢₙ = m′ₘᵢₙ(w)
            i = 1
            for m′ ∈ m′ₘᵢₙ:m′ₘₐₓ
                for m ∈ abs(m′):ℓ
                    for iᵣ ∈ 1:Nᵣ
                        @test decode(w[i]) == (iᵣ, numerator(m′), numerator(m))
                        i += 1
                    end
                end
            end
        end
    end
    function test_3index(w::HWedge{IT}) where {IT}
        let Nᵣ = Nᵣ(w), ℓ = ℓ(w), m′ₘₐₓ = m′ₘₐₓ(w), m′ₘᵢₙ = m′ₘᵢₙ(w)
            for m′ ∈ m′ₘᵢₙ:m′ₘₐₓ
                for m ∈ abs(m′):ℓ
                    for iᵣ ∈ 1:Nᵣ
                        @test decode(w[iᵣ, m′, m]) == (iᵣ, numerator(m′), numerator(m))
                    end
                end
            end
        end
    end

    for ℓₘₐₓ ∈ (5, 9//2)
        for Nᵣ ∈ (1, 2, 3, 7)
            for m′ₘₐₓ ∈ ℓₘᵢₙ(ℓₘₐₓ):ℓₘₐₓ
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
                        fill_1index!(H)
                        test_3index(H)
                        fill_3index!(H)
                        test_1index(H)
                    end
                end
            end
        end
    end
end
