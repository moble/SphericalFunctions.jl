@testitem "HAxis" setup=[EncodeDecode, RefusalChecks] begin
    using SphericalFunctions: HAxis, Nᵣ, maxℓ, HalfOddInteger
    using .EncodeDecode: encode, decode

    # HAxis stores only the m′=0 axis of an integer order, with m ranging from 0 to the
    # current order.  The data layout is `[value for m ∈ 0:ℓ for iᵣ ∈ 1:Nᵣ]`, with the inner
    # loop over iᵣ for vectorization, and the recurrence reads and writes it by the linear
    # index.

    function fill_linear!(h::HAxis)
        let Nᵣ = Nᵣ(h), ℓ = h.ℓ
            i = 1
            for m ∈ 0:ℓ
                for iᵣ ∈ 1:Nᵣ
                    h[i] = encode(iᵣ, 0, m)
                    i += 1
                end
            end
        end
        return h
    end

    function test_linear(h::HAxis)
        let Nᵣ = Nᵣ(h), ℓ = h.ℓ
            for m ∈ 0:ℓ
                for iᵣ ∈ 1:Nᵣ
                    @test decode(h[iᵣ + Nᵣ * m]) == (iᵣ, 0, m)
                    @test decode(parent(h)[iᵣ + Nᵣ * m]) == (iᵣ, 0, m)
                end
            end
        end
    end

    for ℓₘₐₓ ∈ (0, 1, 5)
        for n ∈ (1, 2, 3, 7)
            RT = Float64
            h = HAxis(RT, n, ℓₘₐₓ)
            @test h isa HAxis{RT}

            # Check fields
            @test h.Nᵣ == n == Nᵣ(h)
            @test h.maxℓ == ℓₘₐₓ == maxℓ(h)

            # When first created, ℓ should be at its minimum value
            @test h.ℓ == 0

            # Check storage size (allocated for maximum ℓₘₐₓ)
            expected_length = n * (ℓₘₐₓ + 1)
            @test length(parent(h)) == expected_length
            @test parent(h) === h.parent

            # Test changing ℓ
            for new_ell in 0:ℓₘₐₓ
                h.ℓ = new_ell
                @test h.ℓ == new_ell

                # Linear indexing reads and writes the storage in place
                fill_linear!(h)
                test_linear(h)

                # Linear indexing runs over the whole allocation, so it is the storage that
                # sets the bounds, whatever the current ℓ; check just for the string,
                # because some tests will throw a `FixedSizeArrays.BoundsErrorLight` instead
                # of a standard `BoundsError`.
                @test checkbounds(Bool, h, 1) && checkbounds(Bool, h, expected_length)
                @test !checkbounds(Bool, h, 0) && !checkbounds(Bool, h, expected_length + 1)
                @test_throws "BoundsError" h[0]
                @test_throws "BoundsError" h[length(h.parent) + 1]
                @test_throws "BoundsError" h[0] = 1.0
                @test_throws "BoundsError" h[length(h.parent) + 1] = 1.0
            end

            # Test error conditions for changing ℓ
            @test refuses(() -> h.ℓ = ℓₘₐₓ + 1, ArgumentError, "greater than maxℓ")
            @test refuses(() -> h.ℓ = -1, ArgumentError, "less than ℓₘᵢₙ")
            @test refuses(
                () -> h.ℓ = HalfOddInteger(1//2), ArgumentError, "they must be the same"
            )

            # Test that we can't change other properties
            @test refuses(() -> h.Nᵣ = 10, ArgumentError, "only `ℓ` is allowed to be changed")
            @test refuses(
                () -> h.maxℓ = ℓₘₐₓ + 1, ArgumentError, "only `ℓ` is allowed to be changed"
            )
        end
    end

    # The constructor refuses an axis that could not hold its first order
    @test refuses(() -> HAxis(Float64, 0, 3), ArgumentError, "must be at least 1")
    @test refuses(() -> HAxis(Float64, 2, -1), ArgumentError, "must be at least ℓₘᵢₙ=0")
end
