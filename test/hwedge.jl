@testitem "HWedge" setup=[EncodeDecode] begin
    using SphericalFunctions: HWedge, HWedge_size, Nᵣ, ℓ, ℓₘᵢₙ, m′ₘᵢₙ, m′ₘₐₓ, half_integer
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

    for ℓₘₐₓ ∈ (5, half_integer(9//2))  # an `Int` and a `HalfOddInteger`
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

# The item above exercises the storage layout, which is what `recurrence!` depends on.  The
# items below cover the rest of the `HWedge` interface: the ways it refuses a bad request,
# the indexing paths no calculator takes, and the symmetry helper through which every read
# of an element outside the stored wedge has to pass.

@testitem "HWedge: construction and `ℓ` reassignment refuse bad arguments" begin
    using SphericalFunctions: HWedge, ℓₘᵢₙ, maxm′ₘₐₓ, minm′ₘᵢₙ, half_integer

    # At least one rotor, always
    @test_throws "must be at least 1" HWedge(Float64, 0, 5, 5, -5)
    @test_throws "must be at least 1" HWedge(Float64, -3, 5, 5, -5)

    for ℓₘₐₓ ∈ (5, half_integer(9//2))
        IT = typeof(ℓₘₐₓ)
        H = HWedge(Float64, 2, ℓₘₐₓ)

        # `ℓ` is the one property that may be reassigned, and only within its own type and
        # within the range the storage was allocated for.
        @test_throws "greater than maxℓ" H.ℓ = ℓₘₐₓ + 1
        @test_throws "less than ℓₘᵢₙ" H.ℓ = ℓₘᵢₙ(IT) - 1
        @test_throws "only `ℓ` is allowed to be changed" H.Nᵣ = 10
        @test_throws "only `ℓ` is allowed to be changed" H.maxℓ = ℓₘₐₓ + 1
        @test_throws "only `ℓ` is allowed to be changed" H.m′ₘₐₓ = ℓₘₐₓ

        # Assigning an index of the other kind is refused rather than silently converted,
        # which is what keeps an integer wedge from acquiring half-odd-integer indices.
        other = IT <: Integer ? half_integer(3//2) : 2
        @test_throws "they must be the same" H.ℓ = other

        # The assignment returns the new `ℓ`, and the `m′` bounds follow it
        for new_ℓ ∈ ℓₘᵢₙ(IT):ℓₘₐₓ
            @test (H.ℓ = new_ℓ) == new_ℓ
            @test H.m′ₘₐₓ == min(new_ℓ, maxm′ₘₐₓ(H))
            @test H.m′ₘᵢₙ == max(-new_ℓ, minm′ₘᵢₙ(H))
        end
    end

    # A wedge narrower than the full range clamps against the narrower bound, not against ℓ
    H = HWedge(Float64, 1, 6, 2, -3)
    H.ℓ = 6
    @test H.m′ₘₐₓ == 2
    @test H.m′ₘᵢₙ == -3
end

@testitem "HWedge: axes, bounds checking and display" begin
    using SphericalFunctions: HWedge, Nᵣ, ℓ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, half_integer

    for ℓₘₐₓ ∈ (4, half_integer(7//2))
        H = HWedge(Float64, 3, ℓₘₐₓ)
        H.ℓ = ℓₘₐₓ

        @test axes(H) == (1:Nᵣ(H), m′ₘᵢₙ(H):m′ₘₐₓ(H), mₘᵢₙ(H):mₘₐₓ(H))
        @test length(axes(H)) == 3
        @test first(axes(H)[1]) == 1 && last(axes(H)[1]) == Nᵣ(H)

        # Linear indexing runs over the whole allocation, so it is the parent that sets the
        # bounds.  These throw `FixedSizeArrays.BoundsErrorLight` in some cases, so match on
        # the name rather than the type, as the `HAxis` item does.
        @test_throws "BoundsError" H[0]
        @test_throws "BoundsError" H[length(parent(H)) + 1]
        H[1] = 17.0
        @test H[1] == 17.0
        @test_throws "BoundsError" H[0] = 1.0
        @test_throws "BoundsError" H[length(parent(H)) + 1] = 1.0

        # Three-index bounds: the rotor index, the stored wedge condition |m′| ≤ m ≤ ℓ, and
        # the m′ range are each checked.
        @test_throws "BoundsError" H[0, m′ₘᵢₙ(H), mₘₐₓ(H)]
        @test_throws "BoundsError" H[Nᵣ(H) + 1, m′ₘᵢₙ(H), mₘₐₓ(H)]
        @test_throws "BoundsError" H[1, m′ₘᵢₙ(H) - 1, mₘₐₓ(H)]
        @test_throws "BoundsError" H[1, m′ₘₐₓ(H) + 1, mₘₐₓ(H)]
        @test_throws "BoundsError" H[1, m′ₘₐₓ(H), ℓ(H) + 1]
        @test_throws "BoundsError" H[1, m′ₘₐₓ(H), m′ₘₐₓ(H) - 1]  # m < |m′| is not stored
        @test_throws "BoundsError" H[0, m′ₘᵢₙ(H), mₘₐₓ(H)] = 1.0
        @test_throws "BoundsError" H[1, m′ₘₐₓ(H), ℓ(H) + 1] = 1.0

        # `summary` names the type and the ranges currently in use; `show` uses it
        s = sprint(summary, H)
        @test occursin("HWedge", s)
        @test occursin("ℓ=$(ℓ(H))", s)
        @test occursin("iᵣ=1:$(Nᵣ(H))", s)
        @test sprint(show, H) == s

        # The three-argument `show` adds the storage and the portion currently in use
        s3 = sprint(show, MIME("text/plain"), H)
        @test occursin("HWedge", s3)
        @test occursin("stored in", s3)
        @test occursin("currently using", s3)
    end
end

@testitem "HWedge: `Rational` indices reach the half-odd-integer wedge" begin
    using SphericalFunctions: HWedge, Nᵣ, ℓ, m′ₘᵢₙ, m′ₘₐₓ, half_integer

    # Indexing a half-odd-integer wedge with `Rational`s is a documented convenience, so
    # `H[1, 1//2, 3//2]` has to land on the same element as the `HalfOddInteger` form.
    H = HWedge(Float64, 2, half_integer(7//2))
    H.ℓ = half_integer(7//2)

    for m′ ∈ m′ₘᵢₙ(H):m′ₘₐₓ(H), m ∈ abs(m′):ℓ(H), iᵣ ∈ 1:Nᵣ(H)
        q′, q = Rational(m′), Rational(m)
        H[iᵣ, m′, m] = 100iᵣ + 10Float64(m′) + Float64(m)
        @test H[iᵣ, q′, q] == H[iᵣ, m′, m]
        H[iᵣ, q′, q] = -1.0
        @test H[iᵣ, m′, m] == -1.0
    end
end

@testitem "HWedge: `wedge_source` is the only encoding of the H symmetries" begin
    using SphericalFunctions: wedge_source, wedge_source_error, transpose_sign, sgn
    using SphericalFunctions: half_integer

    # `sgn` is the Gumerov–Duraiswami convention, which differs from `Base.sign` at zero
    @test sgn(0) == 1
    @test sgn(3) == 1 && sgn(-3) == -1

    # For integer indices the transpose costs no sign at all; for half-odd-integer indices
    # it is sgn(m)sgn(m′).
    @test transpose_sign(2, 3) == 1
    @test transpose_sign(-2, 3) == 1
    @test transpose_sign(0, 0) == 1
    @test transpose_sign(half_integer(1//2), half_integer(3//2)) == 1
    @test transpose_sign(half_integer(-1//2), half_integer(3//2)) == -1
    @test transpose_sign(half_integer(-1//2), half_integer(-3//2)) == 1

    # Whatever (m′, m) is asked for, the source must lie in the stored wedge: b ≥ |a| and
    # |a| ≤ m′ₘₐₓ.  Check that for every element of a full matrix, at several m′ₘₐₓ.
    for ℓ ∈ 0:4, m′ₘₐₓ ∈ 0:ℓ
        for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
            if abs(m′) ≤ m′ₘₐₓ || abs(m) ≤ m′ₘₐₓ
                a, b, σ = wedge_source(m′, m, m′ₘₐₓ)
                @test b ≥ abs(a)
                @test abs(a) ≤ m′ₘₐₓ
                @test σ == 1              # integer indices never pick up a sign
                @test max(abs(a), abs(b)) ≤ ℓ
                # The map is an involution on the four symmetry images
                @test (a, b) ∈ ((m′, m), (-m, -m′), (m, m′), (-m′, -m))
            else
                # Neither index is small enough, so no stored element can supply it
                @test_throws ArgumentError wedge_source(m′, m, m′ₘₐₓ)
            end
        end
    end

    # Half-odd-integer indices do pick up a sign, and the identity H[m′,m] = σ H[a,b] has to
    # hold with that σ
    for twoℓ ∈ 1:2:7
        ℓ = half_integer(twoℓ//2)
        for twom′ ∈ -twoℓ:2:twoℓ, twom ∈ -twoℓ:2:twoℓ
            m′, m = half_integer(twom′//2), half_integer(twom//2)
            a, b, σ = wedge_source(m′, m, ℓ)
            @test b ≥ abs(a)
            @test σ ∈ (-1, 1)
            @test (a, b) ∈ ((m′, m), (-m, -m′), (m, m′), (-m′, -m))
        end
    end

    # The error is raised by a separate `@noinline` function so the hot path stays clean
    @test_throws ArgumentError wedge_source_error(3, 4, 2)
    @test_throws "both |m′| and |m| exceed m′ₘₐₓ" wedge_source_error(3, 4, 2)
    @test_throws "H[3, 4]" wedge_source_error(3, 4, 2)
end
