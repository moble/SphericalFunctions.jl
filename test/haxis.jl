@testitem "HAxis" setup=[EncodeDecode] begin
    using SphericalFunctions: HAxis, Nᵣ, ℓ, ℓₘᵢₙ, maxℓ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, half_integer
    using .EncodeDecode: encode, decode

    # HAxis stores only the m′=ℓₘᵢₙ axis (0 or 1/2), with m ranging from ℓₘᵢₙ to ℓₘₐₓ.
    # The data layout is: [value for iᵣ ∈ 1:Nᵣ, m ∈ ℓₘᵢₙ:ℓₘₐₓ]
    # We want inner loop over iᵣ, outer loop over m for vectorization.
    
    function fill_1index!(h::HAxis{IT}) where {IT}
        let Nᵣ = Nᵣ(h), ℓ = ℓ(h), ℓₘᵢₙ = ℓₘᵢₙ(IT)
            i = 1
            for m ∈ ℓₘᵢₙ:ℓ
                for iᵣ ∈ 1:Nᵣ
                    h[i] = encode(iᵣ, ℓₘᵢₙ, m)
                    i += 1
                end
            end
        end
        return h
    end
    
    function fill_2index!(h::HAxis{IT}) where {IT}
        let Nᵣ = Nᵣ(h), ℓ = ℓ(h), ℓₘᵢₙ = ℓₘᵢₙ(IT)
            for m ∈ ℓₘᵢₙ:ℓ
                for iᵣ ∈ 1:Nᵣ
                    h[iᵣ, m] = encode(iᵣ, ℓₘᵢₙ, m)
                end
            end
        end
        return h
    end
    
    function fill_3index!(h::HAxis{IT}) where {IT}
        let Nᵣ = Nᵣ(h), ℓ = ℓ(h), ℓₘᵢₙ = ℓₘᵢₙ(IT)
            for m ∈ ℓₘᵢₙ:ℓ
                for iᵣ ∈ 1:Nᵣ
                    h[iᵣ, ℓₘᵢₙ, m] = encode(iᵣ, ℓₘᵢₙ, m)
                end
            end
        end
        return h
    end
    
    function test_1index(h::HAxis{IT}) where {IT}
        let Nᵣ = Nᵣ(h), ℓ = ℓ(h), ℓₘᵢₙ = ℓₘᵢₙ(IT)
            i = 1
            for m ∈ ℓₘᵢₙ:ℓ
                for iᵣ ∈ 1:Nᵣ
                    @test decode(h[i]) == (iᵣ, numerator(ℓₘᵢₙ), numerator(m))
                    i += 1
                end
            end
        end
    end
    
    function test_2index(h::HAxis{IT}) where {IT}
        let Nᵣ = Nᵣ(h), ℓ = ℓ(h), ℓₘᵢₙ = ℓₘᵢₙ(IT)
            for m ∈ ℓₘᵢₙ:ℓ
                for iᵣ ∈ 1:Nᵣ
                    @test decode(h[iᵣ, m]) == (iᵣ, numerator(ℓₘᵢₙ), numerator(m))
                end
            end
        end
    end
    
    function test_3index(h::HAxis{IT}) where {IT}
        let Nᵣ = Nᵣ(h), ℓ = ℓ(h), ℓₘᵢₙ = ℓₘᵢₙ(IT)
            for m ∈ ℓₘᵢₙ:ℓ
                for iᵣ ∈ 1:Nᵣ
                    @test decode(h[iᵣ, ℓₘᵢₙ, m]) == (iᵣ, numerator(ℓₘᵢₙ), numerator(m))
                end
            end
        end
    end

    # Test both integer and half-integer ℓ
    for ℓₘₐₓ ∈ (5, half_integer(9//2))  # an `Int` and a `HalfOddInteger`
        for Nᵣ ∈ (1, 2, 3, 7)
            IT = typeof(ℓₘₐₓ)
            RT = Float64
            h = HAxis(RT, Nᵣ, ℓₘₐₓ)

            # Check fields
            @test h.Nᵣ == Nᵣ
            @test h.maxℓ == ℓₘₐₓ

            # When first created, ℓ should be at its minimum value
            @test ℓ(h) == ℓₘᵢₙ(IT)

            # Check index ranges
            @test m′ₘᵢₙ(h) == ℓₘᵢₙ(IT)
            @test m′ₘₐₓ(h) == ℓₘᵢₙ(IT)
            @test mₘᵢₙ(h) == ℓₘᵢₙ(IT)
            @test mₘₐₓ(h) == ℓₘᵢₙ(IT)  # mₘₐₓ should equal current ℓ, not ℓₘₐₓ

            # Check storage size (allocated for maximum ℓₘₐₓ)
            expected_length = Nᵣ * (Int(ℓₘₐₓ - ℓₘᵢₙ(IT)) + 1)
            @test length(h.parent) == expected_length

            # Test changing ℓ
            for new_ell in (ℓₘᵢₙ(IT):ℓₘₐₓ)
                h.ℓ = new_ell
                @test h.ℓ == new_ell
                @test ℓ(h) == new_ell
                @test mₘₐₓ(h) == new_ell  # mₘₐₓ should track current ℓ

                # Test all three indexing methods (1D, 2D, 3D)
                fill_1index!(h)
                test_2index(h)
                test_3index(h)
                
                fill_2index!(h)
                test_1index(h)
                test_3index(h)
                
                fill_3index!(h)
                test_1index(h)
                test_2index(h)

                # Test bounds checking for current ℓ; check just for the string, because
                # some tests will throw a `FixedSizeArrays.BoundsErrorLight` instead of a
                # standard `BoundsError`.
                @test_throws "BoundsError" h[0]
                @test_throws "BoundsError" h[length(h.parent) + 1]
                @test_throws "BoundsError" h[0, ℓₘᵢₙ(IT)]
                @test_throws "BoundsError" h[Nᵣ + 1, ℓₘᵢₙ(IT)]
                @test_throws "BoundsError" h[1, ℓₘᵢₙ(IT) - 1]
                @test_throws "BoundsError" h[1, new_ell + 1]  # Beyond current ℓ

                # 3D indexing: m′ must equal ℓₘᵢₙ
                @test_throws "BoundsError" h[1, ℓₘᵢₙ(IT) - 1, ℓₘᵢₙ(IT)]
                @test_throws "BoundsError" h[1, ℓₘᵢₙ(IT) + 1, ℓₘᵢₙ(IT)]
                @test_throws "BoundsError" h[1, ℓₘᵢₙ(IT), ℓₘᵢₙ(IT) - 1]
                @test_throws "BoundsError" h[1, ℓₘᵢₙ(IT), new_ell + 1]  # Beyond current ℓ
            end

            # Test error conditions for changing ℓ
            @test_throws "greater than maxℓ" h.ℓ = ℓₘₐₓ + 1
            @test_throws "less than ℓₘᵢₙ" h.ℓ = ℓₘᵢₙ(IT) - 1
            
            # Test that we can't change other properties
            @test_throws "only `ℓ` is allowed to be changed" h.Nᵣ = 10
            @test_throws "only `ℓ` is allowed to be changed" h.maxℓ = ℓₘₐₓ + 1
        end
    end
end
