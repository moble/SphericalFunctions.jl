@testitem "HWedge" setup=[EncodeDecode] begin
    using SphericalFunctions: HWedge, HWedge_size, Nᵣ, ℓ, ℓₘᵢₙ, m′ₘᵢₙ, m′ₘₐₓ, minm′ₘᵢₙ,
        maxm′ₘₐₓ, HalfOddInteger
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

    for ℓₘₐₓ ∈ (5, HalfOddInteger(9//2))  # an `Int` and a `HalfOddInteger`
        for Nᵣ ∈ (1, 2, 3, 7)
            for m′ₘₐₓ ∈ ℓₘᵢₙ(ℓₘₐₓ):ℓₘₐₓ
                RT = Float64
                H = HWedge(RT, Nᵣ, ℓₘₐₓ, m′ₘₐₓ)

                @test H.Nᵣ == Nᵣ
                @test H.maxℓ == ℓₘₐₓ
                @test maxm′ₘₐₓ(H) == m′ₘₐₓ
                # The range of m′ is symmetric
                @test minm′ₘᵢₙ(H) == -m′ₘₐₓ

                # When first created, the indices should have their smallest values
                @test H.ℓ == ℓₘᵢₙ(ℓₘₐₓ)
                @test H.m′ₘₐₓ == H.ℓ
                @test m′ₘᵢₙ(H) == -H.ℓ

                # But the storage should be full size: every rotor's copy of every stored
                # row, m′ running over the symmetric range and m over |m′|:ℓₘₐₓ, counted
                # here directly rather than by the constructor's own formula
                expected_size = sum(Int(ℓₘₐₓ - abs(m′)) + 1 for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ)
                @test length(H.parent) == Nᵣ * expected_size

                # Test changing ℓ
                for new_ell in (ℓₘᵢₙ(ℓₘₐₓ):ℓₘₐₓ)
                    H.ℓ = new_ell
                    @test H.ℓ == new_ell
                    @test H.m′ₘₐₓ == min(new_ell, m′ₘₐₓ)
                    @test m′ₘᵢₙ(H) == -min(new_ell, m′ₘₐₓ)
                    fill_1index!(H)
                    test_3index(H)
                    fill_3index!(H)
                    test_1index(H)
                end
            end
        end
    end

    # `HWedge_size` counts the elements of a wedge for any range of m′ that brackets the
    # axis, not only a symmetric one; it is checked against the direct count here, for both
    # kinds of index.
    @test all(
        HWedge_size(L, mp, mn) == sum(L - abs(k) + 1 for k ∈ mn:mp)
        for L ∈ 0:12 for mp ∈ 0:L for mn ∈ -L:0
    )
    @test all(
        HWedge_size(L, mp, mn) == sum(Int(L - abs(k)) + 1 for k ∈ mn:mp)
        for L ∈ HalfOddInteger.(1//2:1:25//2)
        for mp ∈ HalfOddInteger(1//2):L for mn ∈ -L:HalfOddInteger(-1//2)
    )
end

# The item above exercises the storage layout, which is what `recurrence!` depends on.  The
# items below cover the rest of the `HWedge` interface: the ways it refuses a bad request,
# the indexing paths no calculator takes, and the symmetry helper through which every read
# of an element outside the stored wedge has to pass.

@testitem "HWedge: construction and `ℓ` reassignment refuse bad arguments" setup=[RefusalChecks] begin
    using SphericalFunctions: HWedge, ℓₘᵢₙ, m′ₘₐₓ, m′ₘᵢₙ, maxm′ₘₐₓ, minm′ₘᵢₙ,
        HalfOddInteger

    # At least one rotor, always
    @test refuses(() -> HWedge(Float64, 0, 5, 5), ArgumentError, "must be at least 1")
    @test refuses(() -> HWedge(Float64, -3, 5), ArgumentError, "must be at least 1")

    # The range of m′ is symmetric, so a lower limit is not an argument at all
    @test_throws MethodError HWedge(Float64, 1, 6, 2, -3)
    # ... and the limits are indices like any other: of one kind, and not of a narrow type
    @test refuses(() -> HWedge(Float64, 1, 7//2, 1), ArgumentError, "mixes integers")
    @test refuses(() -> HWedge(Float64, 1, Int8(4)), ArgumentError, "narrower than `Int`")
    @test refuses(() -> HWedge(2, 7//3), ArgumentError, "neither an integer nor")
    # (`validate_index_ranges` owns this message and its exception type)
    @test_throws "is too large for ℓₘₐₓ" HWedge(Float64, 1, 4, 5)

    # A half-odd-integer may be spelled as a `Rational`, which gives the very same wedge
    @test HWedge(Float64, 2, 7//2, 3//2) isa HWedge{HalfOddInteger, Float64}
    @test maxm′ₘₐₓ(HWedge(2, 7//2, 3//2)) === HalfOddInteger(3//2)
    @test length(parent(HWedge(2, 7//2, 3//2))) ==
        length(parent(HWedge(2, HalfOddInteger(7//2), HalfOddInteger(3//2))))

    for ℓₘₐₓ ∈ (5, HalfOddInteger(9//2))
        IT = typeof(ℓₘₐₓ)
        H = HWedge(Float64, 2, ℓₘₐₓ)

        # `ℓ` is the one property that may be reassigned, and only within its own type and
        # within the range the storage was allocated for.
        @test refuses(() -> H.ℓ = ℓₘₐₓ + 1, ArgumentError, "greater than maxℓ")
        @test refuses(() -> H.ℓ = ℓₘᵢₙ(IT) - 1, ArgumentError, "less than ℓₘᵢₙ")
        for property ∈ (:Nᵣ, :maxℓ, :m′ₘₐₓ)
            @test refuses(
                () -> setproperty!(H, property, ℓₘₐₓ), ArgumentError,
                "only `ℓ` is allowed to be changed"
            )
        end

        # Assigning an index of the other kind is refused rather than silently converted,
        # which is what keeps an integer wedge from acquiring half-odd-integer indices.
        other = IT <: Integer ? HalfOddInteger(3//2) : 2
        @test refuses(() -> H.ℓ = other, ArgumentError, "they must be the same")

        # The assignment returns the new `ℓ`, and the `m′` bounds follow it
        for new_ℓ ∈ ℓₘᵢₙ(IT):ℓₘₐₓ
            @test (H.ℓ = new_ℓ) == new_ℓ
            @test m′ₘₐₓ(H) == min(new_ℓ, maxm′ₘₐₓ(H))
            @test m′ₘᵢₙ(H) == max(-new_ℓ, minm′ₘᵢₙ(H))
        end
    end

    # A half-integer order may be written as a `Rational`, as it may everywhere else, but it
    # must be a half-odd-integer
    H = HWedge(Float64, 1, HalfOddInteger(9//2))
    H.ℓ = 5//2
    @test H.ℓ === HalfOddInteger(5//2)
    @test refuses(() -> H.ℓ = 2//1, ArgumentError, "so ℓ must be one too; got ℓ = 2//1")
    @test refuses(() -> H.m′ₘₐₓ = 1//2, ArgumentError, "only `ℓ` is allowed to be changed")

    # A wedge narrower than the full range clamps against the narrower bound, not against ℓ
    H = HWedge(Float64, 1, 6, 2)
    H.ℓ = 6
    @test m′ₘₐₓ(H) == 2
    @test m′ₘᵢₙ(H) == -2
    H.ℓ = 1
    @test m′ₘₐₓ(H) == 1
    @test m′ₘᵢₙ(H) == -1
end

@testitem "HWedge: axes, bounds checking and display" begin
    using SphericalFunctions: HWedge, Nᵣ, ℓ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, HalfOddInteger

    for ℓₘₐₓ ∈ (4, HalfOddInteger(7//2))
        H = HWedge(Float64, 3, ℓₘₐₓ)
        H.ℓ = ℓₘₐₓ

        @test axes(H) == (1:Nᵣ(H), m′ₘᵢₙ(H):m′ₘₐₓ(H), mₘᵢₙ(H):mₘₐₓ(H))
        @test length(axes(H)) == 3
        @test first(axes(H)[1]) == 1 && last(axes(H)[1]) == Nᵣ(H)
        # `ndims`, `axes(H, d)` and `size(H, d)` agree with those three indices, while
        # `length` and `size(H)` describe the flat storage that linear indexing runs over
        @test ndims(H) == ndims(typeof(H)) == 3
        @test axes(H, 3) == mₘᵢₙ(H):mₘₐₓ(H)
        @test axes(H, 4) == Base.OneTo(1)
        @test (size(H, 1), size(H, 2), size(H, 3), size(H, 4)) ==
            (Nᵣ(H), length(m′ₘᵢₙ(H):m′ₘₐₓ(H)), length(mₘᵢₙ(H):mₘₐₓ(H)), 1)
        @test size(H) == (length(parent(H)),) && length(H) == length(parent(H))

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

        # A copy is a wedge of its own, laid out for the same ℓ and holding the same numbers
        for i ∈ eachindex(parent(H))
            H[i] = i
        end
        Hc = copy(H)
        @test Hc isa typeof(H)
        @test Hc == H && ℓ(Hc) == ℓ(H) && m′ₘₐₓ(Hc) == m′ₘₐₓ(H)
        @test parent(Hc) !== parent(H) && parent(Hc) == parent(H)
        @test Hc.row_index !== H.row_index
        H[1, m′ₘₐₓ(H), mₘₐₓ(H)] = -1.0
        H.ℓ = ℓ(H) - 1
        @test Hc != H && ℓ(Hc) == ℓₘₐₓ && Hc[1, m′ₘₐₓ(Hc), mₘₐₓ(Hc)] != -1.0
    end
    # ... even where the storage holds entries that have never been written
    Hb = HWedge(BigFloat, 2, 3)
    @test copy(Hb) isa HWedge{Int, BigFloat}
end

@testitem "HWedge: `Rational` indices reach the half-odd-integer wedge" begin
    using SphericalFunctions: HWedge, Nᵣ, ℓ, m′ₘᵢₙ, m′ₘₐₓ, HalfOddInteger

    # Indexing a half-odd-integer wedge with `Rational`s is a documented convenience, so
    # `H[1, 1//2, 3//2]` has to land on the same element as the `HalfOddInteger` form.
    H = HWedge(Float64, 2, HalfOddInteger(7//2))
    H.ℓ = HalfOddInteger(7//2)

    for m′ ∈ m′ₘᵢₙ(H):m′ₘₐₓ(H), m ∈ abs(m′):ℓ(H), iᵣ ∈ 1:Nᵣ(H)
        q′, q = Rational(m′), Rational(m)
        H[iᵣ, m′, m] = 100iᵣ + 10Float64(m′) + Float64(m)
        @test H[iᵣ, q′, q] == H[iᵣ, m′, m]
        H[iᵣ, q′, q] = -1.0
        @test H[iᵣ, m′, m] == -1.0
    end
end

@testitem "HWedge: `wedge_source` supplies every element from the stored wedge" setup=[RefusalChecks] begin
    using SphericalFunctions: wedge_source, wedge_source_error, transpose_sign, sgn
    using SphericalFunctions: HalfOddInteger

    # `sgn` is the Gumerov–Duraiswami convention, which differs from `Base.sign` at zero
    @test sgn(0) == 1
    @test sgn(3) == 1 && sgn(-3) == -1

    # For integer indices the transpose costs no sign at all; for half-odd-integer indices
    # it is sgn(m)sgn(m′).
    @test transpose_sign(2, 3) == 1
    @test transpose_sign(-2, 3) == 1
    @test transpose_sign(0, 0) == 1
    @test transpose_sign(HalfOddInteger(1//2), HalfOddInteger(3//2)) == 1
    @test transpose_sign(HalfOddInteger(-1//2), HalfOddInteger(3//2)) == -1
    @test transpose_sign(HalfOddInteger(-1//2), HalfOddInteger(-3//2)) == 1

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
                # The source is one of the four symmetry images of (m′, m)
                @test (a, b) ∈ ((m′, m), (-m, -m′), (m, m′), (-m′, -m))
            else
                # Neither index is small enough, so no stored element can supply it
                @test_throws ArgumentError wedge_source(m′, m, m′ₘₐₓ)
            end
        end
    end

    # Half-odd-integer indices do pick up a sign, σ = sgn(m′) sgn(m) whenever the source is
    # a transposition; the values that the identity H[m′,m] = σ H[a,b] relates are compared
    # in "HCalculator: half-integer wedge_value vs closed-form H", since this item has none
    for twoℓ ∈ 1:2:7, twom′ₘₐₓ ∈ 1:2:twoℓ
        ℓ = HalfOddInteger(twoℓ//2)
        m′ₘₐₓ = HalfOddInteger(twom′ₘₐₓ//2)
        for twom′ ∈ -twoℓ:2:twoℓ, twom ∈ -twoℓ:2:twoℓ
            m′, m = HalfOddInteger(twom′//2), HalfOddInteger(twom//2)
            if abs(m′) ≤ m′ₘₐₓ || abs(m) ≤ m′ₘₐₓ
                a, b, σ = wedge_source(m′, m, m′ₘₐₓ)
                @test b ≥ abs(a)
                @test abs(a) ≤ m′ₘₐₓ
                @test (a, b) ∈ ((m′, m), (-m, -m′), (m, m′), (-m′, -m))
                # The images (m, m′) and (-m′, -m) are transpositions, and (-m, -m′) is not
                @test σ == ((a, b) ∈ ((m, m′), (-m′, -m)) && (a, b) != (m′, m) ?
                    transpose_sign(m′, m) : 1)
            else
                @test_throws ArgumentError wedge_source(m′, m, m′ₘₐₓ)
            end
        end
    end

    # The error is raised by a separate `@noinline` function so the hot path stays clean
    @test refuses(() -> wedge_source_error(3, 4, 2), ArgumentError, "both |m′| and |m| exceed m′ₘₐₓ")
    @test refuses(() -> wedge_source_error(3, 4, 2), ArgumentError, "H[3, 4]")
end

@testitem "HWedge, HAxis and WignerMatrixBatch refuse two indices, and compare element by element" begin
    import SphericalFunctions: HWedge, HAxis, HCalculator, DCalculator, recurrence!, Nᵣ
    using Quaternionic: from_euler_angles

    # Two indices mean `(m′, m)` only for a `WignerMatrix`.  These containers are indexed by
    # three (`[iᵣ, m′, m]`), and a two-index call that fell through to the `WignerMatrix`
    # method would silently read the wrong element: `H[0, 0]` would give 0.0066 where `H[1,
    # 0, 0]` is 0.869, and a batch's `wb[0, 0]` rotor 3's element `[0, -2]`.  So it is an
    # error.
    H = recurrence!(HCalculator(0.3, 3), 3)
    @test H isa HWedge
    @test_throws MethodError H[0, 0]
    @test_throws MethodError H[0, 0] = 1.0
    @test_throws MethodError Matrix(H)
    R = [from_euler_angles(0.1i, 0.2i, 0.3i) for i ∈ 1:3]
    wb = recurrence!(DCalculator(R, 2), 2)
    @test Nᵣ(wb) == 3
    @test_throws MethodError wb[0, 0]
    @test_throws MethodError wb[0, 0] = 1.0
    @test wb[3, 0, -2] == wb[3][0, -2]  # the three-index form is untouched

    # `==` compares every element a wedge holds, rather than a `Matrix` view that ignores m
    # (and so would stay true after stored elements had changed)
    H₁ = recurrence!(HCalculator(0.3, 3), 3)
    H₂ = recurrence!(HCalculator(0.3, 3), 3)
    @test H₁ == H₂
    H₂[1, 1, 2] += 1
    @test H₁ != H₂
    @test recurrence!(HCalculator(0.3, 3), 2) != H₁  # a different ℓ
    @test recurrence!(HCalculator([0.3, 0.4], 3), 3) != H₁  # a different number of rotors

    # ... and so does `==` for the m′ = 0 axis
    a₁, a₂ = HAxis(Float64, 2, 4), HAxis(Float64, 2, 4)
    a₁.ℓ = a₂.ℓ = 3
    for m ∈ 0:3, iᵣ ∈ 1:2
        a₁[iᵣ, m] = a₂[iᵣ, m] = 10iᵣ + m
    end
    @test a₁ == a₂
    a₂[2, 3] = 0.0
    @test a₁ != a₂
end
