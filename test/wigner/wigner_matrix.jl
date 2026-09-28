# Tests of the containers in `src/wigner/wigner_matrix.jl`.
#
# The calculators in `test/wigner/calculators.jl` and the harmonics in `test/sYlm/` build
# these containers constantly, so the parts that hold numbers are well covered already.
# What those tests never touch is the rest of the container interface: constructing one by
# hand from storage, the errors that refuses, `copy`/`similar`/`==`/`iterate`, conversion to
# a plain array, and display.  Each item below takes one container through that interface.
#
# `WignerRange` is the index range these containers report from `axes`, and it comes first
# because everything else depends on it behaving like a range over a possibly half-odd index.

@testitem "WignerRange: a unit range over either kind of index" begin
    import SphericalFunctions: WignerRange, HalfOddInteger

    # An integer range behaves as `UnitRange` does
    r = WignerRange(-2:3)
    @test first(r) == -2 && last(r) == 3
    @test length(r) == 6
    @test step(r) == 1
    @test collect(r) == -2:3
    @test r == -2:3
    @test firstindex(r) == 1 && lastindex(r) == length(r)
    @test r[1] == -2 && r[length(r)] == 3
    @test_throws BoundsError r[0]
    @test_throws BoundsError r[length(r) + 1]

    # A half-odd-integer range needs `step` and `length` given directly, because `Base`
    # forms them from `zero`/`oneunit`, which `HalfOddInteger` deliberately lacks.
    h = WignerRange(HalfOddInteger(-3//2):HalfOddInteger(5//2))
    @test first(h) == HalfOddInteger(-3//2) && last(h) == HalfOddInteger(5//2)
    @test step(h) == 1
    @test length(h) == 5
    @test firstindex(h) == 1 && lastindex(h) == 5
    @test h[1] == HalfOddInteger(-3//2)
    @test h[5] == HalfOddInteger(5//2)
    @test_throws BoundsError h[0]
    @test_throws BoundsError h[6]

    # An empty half-odd range has length 0 rather than a negative length
    @test length(WignerRange(HalfOddInteger(3//2):HalfOddInteger(1//2))) == 0

    # A `Bool` is not an index, even though it is an `Integer`
    @test_throws ArgumentError r[true]
    @test_throws ArgumentError h[true]

    # Membership is just the bracket, and a value of the other index kind is never a member.
    # Without the extra `in` methods, `1 ∈ r` is an ambiguity error rather than an answer.
    # These two are the regression test for the `in` ambiguity, which needs *two* guard
    # methods and not one.  The loose `in(::Integer, ::WignerRange{<:Integer})` settles
    # `in(::Real, ::WignerRange)` against `Base`'s `in(::Integer, ::AbstractUnitRange)`, and
    # Aqua fails if it goes missing.  Being loose, it is then ambiguous in turn with the
    # `IntegerHalf` method whenever the value and the range share one integer type, which
    # `in(x::T, ::WignerRange{T}) where {T<:Integer}` settles by being more specific than
    # both.  Aqua does not see that second pair — `Test.detect_ambiguities` reports nothing
    # for it on either Julia 1.12 or 1.13 — so these direct calls are its only guard.
    @test 0 ∈ r && -2 ∈ r && 3 ∈ r
    @test 4 ∉ r && -3 ∉ r
    @test HalfOddInteger(1//2) ∉ r
    @test 0.5 ∉ r
    @test HalfOddInteger(1//2) ∈ h && HalfOddInteger(-3//2) ∈ h
    @test HalfOddInteger(7//2) ∉ h
    @test 1 ∉ h
    # Any other value is a member exactly when it is `==` to one, as for `Base`'s ranges.  (A
    # test of the index type alone would answer `false` for every `Rational` and float, so that
    # `1//2 ∈ h` would be false while `1//2 ∈ collect(h)` is true.)
    @test 1//2 ∈ h && -3//2 ∈ h && 5//2 ∈ h && 1.5 ∈ h && Int8(3)//Int8(2) ∈ h
    @test 7//2 ∉ h && -5//2 ∉ h && 1//3 ∉ h && 1.0 ∉ h && 0.25 ∉ h && π ∉ h && 2//1 ∉ h
    @test all((x ∈ h) == (x ∈ collect(h)) for x ∈ (-5//2, -3//2, 1//2, 1, 5//2, 7//2))
    @test 1.0 ∈ r && 2//1 ∈ r && -2.0 ∈ r && Int8(3) ∈ r && big(0) ∈ r
    @test 4.0 ∉ r && 1//2 ∉ r && π ∉ r && HalfOddInteger(1//2) ∉ r
    # ... and the plain `UnitRange{HalfOddInteger}` that `keys` of a half-integer container is
    k = HalfOddInteger(1//2):HalfOddInteger(7//2)
    @test 3//2 ∈ k && 7//2 ∈ k && 3.5 ∈ k && HalfOddInteger(5//2) ∈ k
    @test 2 ∉ k && 9//2 ∉ k && 1//4 ∉ k
end

@testitem "Wigner containers: `validate_index_ranges` refuses every bad range" begin
    import SphericalFunctions: validate_index_ranges, HalfOddInteger

    # The five-argument form, used by the two-dimensional containers
    @test validate_index_ranges(3, 3, -3, 3, -3) === nothing
    @test_throws "must be non-negative" validate_index_ranges(-1, 0, 0, 0, 0)
    @test_throws "is less than" validate_index_ranges(3, -1, 1, 3, -3)
    @test_throws "is less than" validate_index_ranges(3, 3, -3, -1, 1)
    @test_throws "too small for this index type" validate_index_ranges(3, -1, -3, 3, -3)
    @test_throws "too large for this index type" validate_index_ranges(3, 3, 1, 3, -3)
    @test_throws "too small for this index type" validate_index_ranges(3, 3, -3, -1, -3)
    @test_throws "too large for this index type" validate_index_ranges(3, 3, -3, 3, 1)
    @test_throws "too large for ℓₘₐₓ" validate_index_ranges(2, 3, -2, 2, -2)
    @test_throws "too large for ℓₘₐₓ" validate_index_ranges(2, 2, -3, 2, -2)
    @test_throws "too large for ℓₘₐₓ" validate_index_ranges(2, 2, -2, 3, -2)
    @test_throws "too large for ℓₘₐₓ" validate_index_ranges(2, 2, -2, 2, -3)

    # The three-argument form, used where only the m′ range is constrained
    @test validate_index_ranges(3, 3, -3) === nothing
    @test_throws "must be non-negative" validate_index_ranges(-1, 0, 0)
    @test_throws "is less than" validate_index_ranges(3, -1, 1)
    @test_throws "too small for this index type" validate_index_ranges(3, -1, -3)
    @test_throws "too large for this index type" validate_index_ranges(3, 3, 1)
    @test_throws "too large for ℓₘₐₓ" validate_index_ranges(2, 3, -2)
    @test_throws "too large for ℓₘₐₓ" validate_index_ranges(2, 2, -3)

    # Half-odd-integer indices must bracket ±1/2, not 0: the recurrence seeds from both rows
    h(x) = HalfOddInteger(x)
    @test validate_index_ranges(h(5//2), h(5//2), h(-5//2)) === nothing
    @test_throws "too small for this index type" validate_index_ranges(h(5//2), h(-1//2), h(-5//2))
    @test_throws "too large for this index type" validate_index_ranges(h(5//2), h(5//2), h(1//2))

    # Every refusal is an `ArgumentError`, and the bracketing refusals state the rule
    @test_throws ArgumentError validate_index_ranges(-1, 0, 0, 0, 0)
    @test_throws ArgumentError validate_index_ranges(3, -1, 1, 3, -3)
    @test_throws ArgumentError validate_index_ranges(2, 3, -2, 2, -2)
    @test_throws ArgumentError validate_index_ranges(-1, 0, 0)
    @test_throws "the range of m′ must include 0, where the recurrence starts" validate_index_ranges(3, 3, 1)
    @test_throws "the range of m must include 0, where the recurrence starts" validate_index_ranges(3, 3, -3, -1, -3)
    @test_throws "the range of m′ must include both -1//2 and 1//2" validate_index_ranges(h(5//2), h(5//2), h(1//2))
    @test_throws "the range of m must include both -1//2 and 1//2" validate_index_ranges(h(5//2), h(5//2), h(-5//2), h(5//2), h(1//2))
end

@testitem "WignerMatrix: the container interface" begin
    import SphericalFunctions: WignerMatrix, ℓ, ℓₘᵢₙ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ
    import SphericalFunctions: ishalfinteger, HalfOddInteger

    for L ∈ (2, HalfOddInteger(3//2))
        n = Int(2L + 1)
        A = reshape(collect(1.0:n^2), n, n)
        w = WignerMatrix(copy(A), L)

        @test ℓ(w) == L
        @test parent(w) == A
        @test eltype(w) == Float64
        @test eltype(typeof(w)) == Float64
        @test ndims(w) == 2 && ndims(typeof(w)) == 2
        @test size(w) == (n, n)
        @test size(w, 1) == n && size(w, 3) == 1        # trailing dims behave as for arrays
        @test axes(w) == (m′ₘᵢₙ(w):m′ₘₐₓ(w), mₘᵢₙ(w):mₘₐₓ(w))
        @test axes(w, 1) == m′ₘᵢₙ(w):m′ₘₐₓ(w)
        @test axes(w, 3) == Base.OneTo(1)
        @test ishalfinteger(w) == !(L isa Integer)

        # Natural indexing reaches the 1-based storage in column-major order
        for (j, m) ∈ enumerate(mₘᵢₙ(w):mₘₐₓ(w)), (i, m′) ∈ enumerate(m′ₘᵢₙ(w):m′ₘₐₓ(w))
            @test w[m′, m] == A[i, j]
        end
        w[m′ₘᵢₙ(w), mₘₐₓ(w)] = -1.0
        @test w[m′ₘᵢₙ(w), mₘₐₓ(w)] == -1.0
        @test_throws BoundsError w[m′ₘₐₓ(w) + 1, mₘₐₓ(w)]
        @test_throws BoundsError w[m′ₘᵢₙ(w), mₘᵢₙ(w) - 1]
        @test_throws BoundsError w[m′ₘₐₓ(w) + 1, mₘₐₓ(w)] = 0.0

        # `Matrix`, `Array` and `collect` all give the same plain matrix
        M = Matrix(w)
        @test M isa Matrix{Float64}
        @test size(M) == (n, n)
        @test Array(w) == M
        @test collect(w) == M

        # Iteration is column-major over the block, matching `Matrix(w)`
        @test [x for x ∈ w] == vec(M)
        @test length(w) == n^2

        # `copy` is independent; `similar` shares only the axes
        c = copy(w)
        @test c == w
        c[m′ₘᵢₙ(c), mₘᵢₙ(c)] = 99.0
        @test w[m′ₘᵢₙ(w), mₘᵢₙ(w)] != 99.0
        s = similar(w)
        @test ℓ(s) == ℓ(w) && axes(s) == axes(w) && eltype(s) == Float64
        sc = similar(w, ComplexF64)
        @test eltype(sc) == ComplexF64 && axes(sc) == axes(w)

        # `==` compares ℓ, axes and values: the same zeros differ by ℓ alone, or by axes alone
        @test copy(w) == w
        c1 = copy(w)
        c1[m′ₘₐₓ(w), mₘₐₓ(w)] += 1
        @test c1 != w
        @test WignerMatrix(zeros(3, 3), 1) != WignerMatrix(zeros(3, 3), 2; m′ₘₐₓ=1, mₘₐₓ=1)
        @test WignerMatrix(zeros(3, 3), 2; m′ₘₐₓ=1, mₘₐₓ=1) !=
            WignerMatrix(zeros(3, 3), 2; m′ₘₐₓ=2, m′ₘᵢₙ=0, mₘₐₓ=1)

        # Display names the type and the ℓ, and the three-argument form prints the values
        summ = sprint(summary, w)
        @test occursin("WignerMatrix", summ)
        @test occursin("ℓ=$L", summ)
        @test sprint(show, w) == summ
        shown = sprint(show, MIME("text/plain"), w)
        @test occursin("WignerMatrix", shown)
        @test occursin(string(M[1, 1]), shown)
    end

    # A `Rational` ℓ is accepted and converted, and so are `Rational` indices
    w = WignerMatrix(zeros(4, 4), 3//2)
    @test ℓ(w) == HalfOddInteger(3//2)
    w[1//2, -1//2] = 7.0
    @test w[1//2, -1//2] == 7.0
    @test w[HalfOddInteger(1//2), HalfOddInteger(-1//2)] == 7.0

    # Storage too small for the requested ranges is refused, per dimension
    @test_throws "first dimension" WignerMatrix(zeros(3, 5), 2)
    @test_throws "second dimension" WignerMatrix(zeros(5, 3), 2)

    # The limits: each lower one defaults to minus the upper one, each has an ASCII spelling,
    # and where both spellings are given the Unicode one is used
    for (kw, axes′) in (
        ((; m′ₘₐₓ=1), (-1:1, -2:2)), ((; mp_max=1), (-1:1, -2:2)),
        ((; mₘₐₓ=1), (-2:2, -1:1)), ((; m_max=1, m_min=0), (-2:2, 0:1)),
        ((; mp_max=2, mp_min=-1, m_max=1, m_min=-2), (-1:2, -2:1)),
        ((; mp_max=0, m′ₘₐₓ=2, m′ₘᵢₙ=0), (0:2, -2:2)),
    )
        @test axes(WignerMatrix(zeros(5, 5), 2; kw...)) == axes′
    end
    h(x) = HalfOddInteger(x)
    # The limits are normalized against the kind of ℓ whatever the spelling of either: a
    # `Rational` limit with a `HalfOddInteger` ℓ, or the reverse, is converted
    for L in (3//2, h(3//2)), m in (1//2, h(1//2))
        wh = WignerMatrix(zeros(4, 4), L; m′ₘₐₓ=m, m_max=m)
        @test axes(wh) == (-1//2:1//2, -1//2:1//2) && ℓ(wh) === h(3//2)
        @test m′ₘₐₓ(wh) === h(1//2)
    end
    # A limit of the other kind, or of another integer type, is refused with the keyword named
    keyword = "must be an index of the same kind as the positional indices of the call"
    @test_throws ArgumentError WignerMatrix(zeros(4, 4), 3//2; m′ₘₐₓ=1)
    @test_throws "The keyword argument `m′ₘₐₓ` of `WignerMatrix`" WignerMatrix(zeros(4, 4), 3//2; m′ₘₐₓ=1)
    @test_throws keyword WignerMatrix(zeros(5, 5), 2; mₘᵢₙ=-1//2)
    @test_throws "`mp_max` of `WignerMatrix`" WignerMatrix(zeros(5, 5), 2; mp_max=1//2)
    @test_throws "`Int32` is narrower than `Int`" WignerMatrix(zeros(5, 5), 2; m_max=Int32(1))
    # ... and so is an ℓ that is not an `Int` or a half-odd-integer
    for (L, why) in (
        (Int32(2), "`Int32` is narrower than `Int`"), (UInt(2), "`UInt64` is unsigned"),
        (big(2), "`BigInt` is wider than `Int`"), (true, "A `Bool` is not an index"),
        (2//1, "2//1 is a whole number; write it as the integer 2"),
        (Int8(3)//Int8(2), "`Rational{Int8}` is not `Rational{Int}`"),
    )
        @test_throws ArgumentError WignerMatrix(zeros(5, 5), L)
        @test_throws why WignerMatrix(zeros(5, 5), L)
    end
    # The bracketing rule is stated when a range that misses the recurrence's seed is asked for
    @test_throws "the range of m′ must include 0" WignerMatrix(view(rand(5, 5), 4:5, :), 2; m′ₘₐₓ=2, m′ₘᵢₙ=1)
end

@testitem "WignerMatrixBatch: the container interface" begin
    import SphericalFunctions: WignerMatrixBatch, WignerMatrix, ℓ, Nᵣ
    import SphericalFunctions: m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, HalfOddInteger

    for L ∈ (2, HalfOddInteger(3//2))
        n, N = Int(2L + 1), 3
        A = reshape(collect(1.0:N*n^2), N, n, n)
        w = WignerMatrixBatch(copy(A), L)

        @test ℓ(w) == L && Nᵣ(w) == N
        @test parent(w) == A
        @test size(w) == (N, n, n)
        # `Array` is the storage, which holds exactly the block here, in `[iᵣ, m′, m]` order
        @test Array(w) == A

        for (k, m) ∈ enumerate(mₘᵢₙ(w):mₘₐₓ(w)), (j, m′) ∈ enumerate(m′ₘᵢₙ(w):m′ₘₐₓ(w)), iᵣ ∈ 1:N
            @test w[iᵣ, m′, m] == A[iᵣ, j, k]
        end
        w[2, m′ₘᵢₙ(w), mₘₐₓ(w)] = -5.0
        @test w[2, m′ₘᵢₙ(w), mₘₐₓ(w)] == -5.0
        @test_throws BoundsError w[0, m′ₘᵢₙ(w), mₘₐₓ(w)]
        @test_throws BoundsError w[N + 1, m′ₘᵢₙ(w), mₘₐₓ(w)]
        @test_throws BoundsError w[1, m′ₘₐₓ(w) + 1, mₘₐₓ(w)]
        @test_throws BoundsError w[1, m′ₘₐₓ(w) + 1, mₘₐₓ(w)] = 0.0

        # One rotor's matrix is a `WignerMatrix` view, so writing through it writes back
        w1 = w[1]
        @test w1 isa WignerMatrix
        @test ℓ(w1) == L && axes(w1) == (m′ₘᵢₙ(w):m′ₘₐₓ(w), mₘᵢₙ(w):mₘₐₓ(w))
        w1[m′ₘᵢₙ(w), mₘᵢₙ(w)] = 123.0
        @test w[1, m′ₘᵢₙ(w), mₘᵢₙ(w)] == 123.0
        @test_throws BoundsError w[0]
        @test_throws BoundsError w[N + 1]

        # `Array` follows the storage through writes; `Matrix` is refused as ambiguous
        @test Array(w) == parent(w)
        @test collect(w) == Array(w)
        @test_throws "3-dimensional" Matrix(w)

        # Iteration matches `Array(w)`, so the reducers work as they do on a plain array
        @test [x for x ∈ w] == vec(Array(w))
        @test sum(w) ≈ sum(Array(w))

        c = copy(w)
        @test Array(c) == Array(w)
        c[1, m′ₘᵢₙ(c), mₘᵢₙ(c)] = -77.0
        @test w[1, m′ₘᵢₙ(w), mₘᵢₙ(w)] != -77.0
        s = similar(w)
        @test size(s) == size(w) && Nᵣ(s) == N && ℓ(s) == L
        @test eltype(similar(w, ComplexF64)) == ComplexF64

        shown = sprint(show, MIME("text/plain"), w)
        @test occursin("WignerMatrixBatch", shown)
    end

    # Storage too small in either indexed dimension is refused
    @test_throws "second dimension" WignerMatrixBatch(zeros(2, 3, 5), 2)
    @test_throws "third dimension" WignerMatrixBatch(zeros(2, 5, 3), 2)

    # `Rational` indices reach a half-odd-integer batch
    w = WignerMatrixBatch(zeros(2, 4, 4), HalfOddInteger(3//2))
    w[1, 1//2, -1//2] = 4.0
    @test w[1, 1//2, -1//2] == 4.0
    @test w[1, HalfOddInteger(1//2), HalfOddInteger(-1//2)] == 4.0

    # A `Rational` ℓ and `Rational` keyword bounds are converted, as for `WignerMatrix`, and
    # each lower bound defaults to minus the upper one
    wr = WignerMatrixBatch(zeros(2, 4, 4), 3//2; mₘₐₓ=1//2)
    @test wr isa WignerMatrixBatch{typeof(HalfOddInteger(1//2))}
    @test ℓ(wr) == HalfOddInteger(3//2) && mₘₐₓ(wr) == HalfOddInteger(1//2)
    @test mₘᵢₙ(wr) == HalfOddInteger(-1//2)
    @test axes(WignerMatrixBatch(zeros(2, 4, 4), HalfOddInteger(3//2); m_max=1//2, mp_max=3//2)) ==
        (1:2, -3//2:3//2, -1//2:1//2)
    @test axes(WignerMatrixBatch(zeros(2, 5, 5), 2; mp_max=1, mp_min=0, m_min=-1)) ==
        (1:2, 0:1, -1:2)
    @test_throws "1//1 is a whole number; write it as the integer 1" WignerMatrixBatch(zeros(2, 3, 3), 1//1)
    @test_throws "`Int16` is narrower than `Int`" WignerMatrixBatch(zeros(2, 3, 3), Int16(1))
    @test_throws "the range of m must include 0" WignerMatrixBatch(zeros(2, 5, 5), 2; mₘₐₓ=2, mₘᵢₙ=1)
end

@testitem "DegreeBlock and DegreeBlockBatch: the container interface" begin
    import SphericalFunctions: DegreeBlock, DegreeBlockBatch, ℓ, Nᵣ, mₘᵢₙ, mₘₐₓ, HalfOddInteger

    for L ∈ (2, HalfOddInteger(3//2))
        n = Int(2L + 1)
        v = DegreeBlock(collect(1.0:n), L)

        @test ℓ(v) == L
        @test firstindex(v) == mₘᵢₙ(v) && lastindex(v) == mₘₐₓ(v)
        @test collect(keys(v)) == collect(mₘᵢₙ(v):mₘₐₓ(v))
        for (i, m) ∈ enumerate(mₘᵢₙ(v):mₘₐₓ(v))
            @test v[m] == i
        end
        v[mₘₐₓ(v)] = -3.0
        @test v[mₘₐₓ(v)] == -3.0
        @test_throws BoundsError v[mₘₐₓ(v) + 1]
        @test_throws BoundsError v[mₘᵢₙ(v) - 1]
        @test_throws BoundsError v[mₘₐₓ(v) + 1] = 0.0

        @test Vector(v) == parent(v)
        @test Array(v) == Vector(v)
        @test collect(v) == Vector(v)
        @test [x for x ∈ v] == Vector(v)

        c = copy(v)
        @test c == v
        c[mₘᵢₙ(c)] = 42.0
        @test v[mₘᵢₙ(v)] != 42.0
        @test ℓ(similar(v)) == L
        @test eltype(similar(v, ComplexF64)) == ComplexF64

        summ = sprint(summary, v)
        @test occursin("DegreeBlock", summ)
        @test sprint(show, v) == summ
        @test occursin("DegreeBlock", sprint(show, MIME("text/plain"), v))

        # The batched sibling, indexed `[iᵣ, m]`
        N = 3
        B = reshape(collect(1.0:N*n), N, n)
        b = DegreeBlockBatch(copy(B), L)
        @test ℓ(b) == L && Nᵣ(b) == N
        for (j, m) ∈ enumerate(mₘᵢₙ(b):mₘₐₓ(b)), iᵣ ∈ 1:N
            @test b[iᵣ, m] == B[iᵣ, j]
        end
        b[2, mₘₐₓ(b)] = -9.0
        @test b[2, mₘₐₓ(b)] == -9.0
        @test_throws BoundsError b[0, mₘₐₓ(b)]
        @test_throws BoundsError b[N + 1, mₘₐₓ(b)]
        @test_throws BoundsError b[1, mₘₐₓ(b) + 1]
        @test_throws BoundsError b[1, mₘₐₓ(b) + 1] = 0.0

        # One rotor's row is a `DegreeBlock` view
        b1 = b[1]
        @test b1 isa DegreeBlock
        @test ℓ(b1) == L
        b1[mₘᵢₙ(b)] = 55.0
        @test b[1, mₘᵢₙ(b)] == 55.0
        @test_throws BoundsError b[0]
        @test_throws BoundsError b[N + 1]

        @test Matrix(b) == [b[iᵣ, m] for iᵣ ∈ 1:N, m ∈ mₘᵢₙ(b):mₘₐₓ(b)]
        @test Array(b) == Matrix(b)
        @test collect(b) == Matrix(b)
        @test [x for x ∈ b] == vec(Matrix(b))

        cb = copy(b)
        @test Matrix(cb) == Matrix(b)
        cb[1, mₘᵢₙ(cb)] = -11.0
        @test b[1, mₘᵢₙ(b)] != -11.0
        @test Nᵣ(similar(b)) == N
        @test eltype(similar(b, ComplexF64)) == ComplexF64

        @test occursin("DegreeBlockBatch", sprint(summary, b))
        @test occursin("DegreeBlockBatch", sprint(show, MIME("text/plain"), b))
    end

    # Storage too small is refused
    @test_throws "length at least" DegreeBlock(zeros(4), 2)
    @test_throws "second dimension" DegreeBlockBatch(zeros(2, 4), 2)

    # A `Rational` ℓ, and `Rational` indexing
    v = DegreeBlock(zeros(4), 3//2)
    @test ℓ(v) == HalfOddInteger(3//2)
    v[1//2] = 6.0
    @test v[1//2] == 6.0
    @test v[HalfOddInteger(1//2)] == 6.0
    vb = DegreeBlockBatch(zeros(2, 4), 3//2; mₘᵢₙ=-1//2)
    @test ℓ(vb) == HalfOddInteger(3//2) && mₘᵢₙ(vb) == HalfOddInteger(-1//2)
    vb[2, 1//2] = 7.0
    @test vb[2, 1//2] == 7.0
    @test_throws "1//1 is a whole number; write it as the integer 1" DegreeBlockBatch(zeros(2, 3), 1//1)
    @test_throws "`Int32` is narrower than `Int`" DegreeBlock(zeros(3), Int32(1))

    # The limits of the m axis: each lower one defaults to minus the upper one, each has an
    # ASCII spelling, and they need not bracket zero, since no recurrence fills these blocks
    @test axes(DegreeBlock(zeros(5), 2; mₘₐₓ=1)) == (-1:1,)
    @test axes(DegreeBlock(zeros(5), 2; m_max=2, m_min=1)) == (1:2,)
    @test axes(DegreeBlock(zeros(5), 2; m_max=0, mₘₐₓ=-1, mₘᵢₙ=-2)) == (-2:-1,)
    @test axes(DegreeBlockBatch(zeros(3, 5), 2; m_min=0)) == (1:3, 0:2)
    @test axes(DegreeBlock(zeros(2), HalfOddInteger(3//2); mₘₐₓ=3//2, m_min=1//2)) == (1//2:3//2,)
    # ... but they must be in order and lie within ±ℓ, and ℓ must be non-negative
    for B in (DegreeBlock, DegreeBlockBatch)
        storage(n) = B === DegreeBlock ? zeros(n) : zeros(2, n)
        @test_throws ArgumentError B(storage(11), 1; mₘₐₓ=-1, mₘᵢₙ=1)
        @test_throws "mₘₐₓ=-1 is less than mₘᵢₙ=1" B(storage(11), 1; mₘₐₓ=-1, mₘᵢₙ=1)
        @test_throws "must lie within -ℓ:ℓ, which is -2:2" B(storage(11), 2; mₘₐₓ=5, mₘᵢₙ=-5)
        @test_throws "must lie within -ℓ:ℓ, which is -1:1" B(storage(11), 1; mₘₐₓ=1, mₘᵢₙ=-2)
        @test_throws "ℓ=-1 must be non-negative" B(storage(3), -1; mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws "ℓ=-1//2 must be non-negative" B(storage(3), -1//2; mₘₐₓ=1//2, mₘᵢₙ=-1//2)
        @test_throws "must lie within -ℓ:ℓ, which is -1//2:1//2" B(storage(3), 1//2; mₘₐₓ=3//2)
    end
end

@testitem "SpinMatrix and SpinMatrixBatch: the container interface" begin
    import SphericalFunctions: SpinMatrix, SpinMatrixBatch, DegreeBlock
    import SphericalFunctions: ℓ, Nᵣ, sₘᵢₙ, sₘₐₓ, mₘᵢₙ, mₘₐₓ, spins, HalfOddInteger

    for L ∈ (2, HalfOddInteger(3//2))
        n = Int(2L + 1)
        smin, smax = -L, L
        ns = Int(smax - smin) + 1
        A = reshape(collect(1.0:ns*n), ns, n)
        b = SpinMatrix(copy(A), L; sₘₐₓ=smax, sₘᵢₙ=smin)

        @test ℓ(b) == L
        @test sₘᵢₙ(b) == smin && sₘₐₓ(b) == smax
        @test collect(spins(b)) == collect(smin:smax)
        for (j, m) ∈ enumerate(mₘᵢₙ(b):mₘₐₓ(b)), (i, s) ∈ enumerate(smin:smax)
            @test b[s, m] == A[i, j]
        end
        b[smin, mₘₐₓ(b)] = -2.0
        @test b[smin, mₘₐₓ(b)] == -2.0
        @test_throws BoundsError b[smax + 1, mₘₐₓ(b)]
        @test_throws BoundsError b[smin, mₘᵢₙ(b) - 1]
        @test_throws BoundsError b[smax + 1, mₘₐₓ(b)] = 0.0

        # A whole spin row is a `DegreeBlock`, and it is a view
        row = b[smin, :]
        @test row isa DegreeBlock
        @test ℓ(row) == L
        row[mₘᵢₙ(b)] = 31.0
        @test b[smin, mₘᵢₙ(b)] == 31.0

        @test Matrix(b) == [b[s, m] for s ∈ smin:smax, m ∈ mₘᵢₙ(b):mₘₐₓ(b)]
        @test Array(b) == Matrix(b)
        @test collect(b) == Matrix(b)
        @test [x for x ∈ b] == vec(Matrix(b))

        c = copy(b)
        @test c == b
        c[smin, mₘᵢₙ(c)] = -13.0
        @test b[smin, mₘᵢₙ(b)] != -13.0
        @test ℓ(similar(b)) == L
        @test eltype(similar(b, ComplexF64)) == ComplexF64

        @test occursin("SpinMatrix", sprint(summary, b))
        @test occursin("SpinMatrix", sprint(show, MIME("text/plain"), b))

        # The batched sibling, indexed `[iᵣ, s, m]`
        N = 3
        C = reshape(collect(1.0:N*ns*n), N, ns, n)
        bb = SpinMatrixBatch(copy(C), L; sₘₐₓ=smax, sₘᵢₙ=smin)
        @test ℓ(bb) == L && Nᵣ(bb) == N
        for (k, m) ∈ enumerate(mₘᵢₙ(bb):mₘₐₓ(bb)), (j, s) ∈ enumerate(smin:smax), iᵣ ∈ 1:N
            @test bb[iᵣ, s, m] == C[iᵣ, j, k]
        end
        bb[2, smin, mₘₐₓ(bb)] = -6.0
        @test bb[2, smin, mₘₐₓ(bb)] == -6.0
        @test_throws BoundsError bb[0, smin, mₘₐₓ(bb)]
        @test_throws BoundsError bb[N + 1, smin, mₘₐₓ(bb)]
        @test_throws BoundsError bb[1, smax + 1, mₘₐₓ(bb)]
        @test_throws BoundsError bb[1, smax + 1, mₘₐₓ(bb)] = 0.0

        # One rotor's block is a `SpinMatrix` view
        bb1 = bb[1]
        @test bb1 isa SpinMatrix
        bb1[smin, mₘᵢₙ(bb)] = 88.0
        @test bb[1, smin, mₘᵢₙ(bb)] == 88.0
        @test_throws BoundsError bb[0]
        @test_throws BoundsError bb[N + 1]

        # A whole spin slice across all rotors
        sl = bb[:, smin, :]
        @test size(sl) == (N, n)

        @test Array(bb) == [bb[i, s, m] for i ∈ 1:N, s ∈ smin:smax, m ∈ mₘᵢₙ(bb):mₘₐₓ(bb)]
        @test collect(bb) == Array(bb)
        @test [x for x ∈ bb] == vec(Array(bb))

        cb = copy(bb)
        @test Array(cb) == Array(bb)
        cb[1, smin, mₘᵢₙ(cb)] = -21.0
        @test bb[1, smin, mₘᵢₙ(bb)] != -21.0
        @test Nᵣ(similar(bb)) == N
        @test eltype(similar(bb, ComplexF64)) == ComplexF64

        @test occursin("SpinMatrixBatch", sprint(summary, bb))
        @test occursin("SpinMatrixBatch", sprint(show, MIME("text/plain"), bb))
    end

    # Storage too small is refused, per dimension
    @test_throws "first dimension" SpinMatrix(zeros(3, 5), 2; sₘₐₓ=2, sₘᵢₙ=-2)
    @test_throws "second dimension" SpinMatrix(zeros(5, 3), 2; sₘₐₓ=2, sₘᵢₙ=-2)
    @test_throws "second dimension" SpinMatrixBatch(zeros(2, 3, 5), 2; sₘₐₓ=2, sₘᵢₙ=-2)
    @test_throws "third dimension" SpinMatrixBatch(zeros(2, 5, 3), 2; sₘₐₓ=2, sₘᵢₙ=-2)

    # A `Rational` ℓ is converted, along with the keyword spin bounds, which are indices like
    # any other and are normalized against the kind of ℓ.  These two containers are the only
    # ones that take spin bounds.
    br = SpinMatrix(zeros(4, 4), 3//2; sₘₐₓ=3//2, sₘᵢₙ=-3//2)
    @test br isa SpinMatrix
    @test ℓ(br) == HalfOddInteger(3//2)
    @test sₘₐₓ(br) == HalfOddInteger(3//2) && sₘᵢₙ(br) == HalfOddInteger(-3//2)
    bbr = SpinMatrixBatch(zeros(2, 4, 4), 3//2; sₘₐₓ=3//2, sₘᵢₙ=-3//2)
    @test bbr isa SpinMatrixBatch
    @test ℓ(bbr) == HalfOddInteger(3//2)
    @test sₘₐₓ(bbr) == HalfOddInteger(3//2) && sₘᵢₙ(bbr) == HalfOddInteger(-3//2)

    # ... and it agrees with the `HalfOddInteger` spelling, which is what they store
    b = SpinMatrix(zeros(4, 4), HalfOddInteger(3//2);
                   sₘₐₓ=HalfOddInteger(3//2), sₘᵢₙ=HalfOddInteger(-3//2))
    @test ℓ(b) == HalfOddInteger(3//2)
    b[1//2, -1//2] = 5.0
    @test b[1//2, -1//2] == 5.0
    bb = SpinMatrixBatch(zeros(2, 4, 4), HalfOddInteger(3//2);
                         sₘₐₓ=HalfOddInteger(3//2), sₘᵢₙ=HalfOddInteger(-3//2))
    @test ℓ(bb) == HalfOddInteger(3//2)
    bb[1, 1//2, -1//2] = 8.0
    @test bb[1, 1//2, -1//2] == 8.0
    # ... with a `Rational` bound on a `HalfOddInteger` ℓ, or the reverse
    @test spins(SpinMatrix(zeros(2, 4), HalfOddInteger(3//2); sₘₐₓ=1//2, sₘᵢₙ=-1//2)) ==
        spins(SpinMatrix(zeros(2, 4), 3//2; sₘₐₓ=HalfOddInteger(1//2), sₘᵢₙ=HalfOddInteger(-1//2)))

    # The spin bounds are required, in either spelling, and the m bounds are optional
    for B in (SpinMatrix, SpinMatrixBatch)
        storage(dims...) = B === SpinMatrix ? zeros(dims...) : zeros(2, dims...)
        @test spins(B(storage(3, 5), 2; s_max=1, s_min=-1)) == -1:1
        @test spins(B(storage(3, 5), 2; s_max=3, sₘₐₓ=1, s_min=-1)) == -1:1
        @test axes(B(storage(3, 5), 2; sₘₐₓ=1, sₘᵢₙ=-1, m_max=1))[end] == -1:1
        @test_throws UndefKeywordError B(storage(3, 5), 2; s_max=1)
        @test_throws UndefKeywordError(:sₘₐₓ) B(storage(3, 5), 2; sₘᵢₙ=-1)
        @test_throws UndefKeywordError(:sₘᵢₙ) B(storage(3, 5), 2; sₘₐₓ=1)
        # The spin axis may be any run of spin weights, beyond ±ℓ and on one side of zero, but
        # it must be in order; the m axis must be in order and within ±ℓ
        @test spins(B(storage(2, 5), 2; sₘₐₓ=4, sₘᵢₙ=3)) == 3:4
        @test_throws ArgumentError B(storage(3, 5), 2; sₘₐₓ=-1, sₘᵢₙ=1)
        @test_throws "sₘₐₓ=-1 is less than sₘᵢₙ=1" B(storage(3, 5), 2; sₘₐₓ=-1, sₘᵢₙ=1)
        @test_throws "must lie within -ℓ:ℓ, which is -1:1" B(storage(3, 5), 1; sₘₐₓ=1, sₘᵢₙ=-1, mₘₐₓ=2, mₘᵢₙ=0)
        @test_throws "mₘₐₓ=0 is less than mₘᵢₙ=1" B(storage(3, 5), 2; sₘₐₓ=1, sₘᵢₙ=-1, mₘₐₓ=0, mₘᵢₙ=1)
        # ... and every bound must be of the kind of ℓ
        @test_throws "The keyword argument `sₘᵢₙ`" B(storage(2, 4), 3//2; sₘₐₓ=1//2, sₘᵢₙ=0)
        @test_throws "`s_max` of `$(nameof(B))`" B(storage(3, 5), 2; s_max=1//2, s_min=-1)
    end
end

@testitem "Wigner containers: uninitialized BigFloat storage shows as #undef" begin
    import SphericalFunctions: WignerDMatrix, WignerMatrix, WignerMatrixBatch, DegreeBlock,
        DegreeBlockBatch, SpinMatrix, SpinMatrixBatch, HCalculator, HalfOddInteger

    # Storage of a non-bits type starts out unassigned, and reading such an element throws
    # an `UndefRefError`; displaying a container must not read it.
    for L ∈ (2, HalfOddInteger(3//2))
        n = Int(2L + 1)
        for x ∈ (
            WignerDMatrix(Complex{BigFloat}, L),
            WignerMatrix(Matrix{BigFloat}(undef, n, n), L),
            WignerMatrixBatch(Array{BigFloat}(undef, 2, n, n), L),
            DegreeBlock(Vector{BigFloat}(undef, n), L),
            DegreeBlockBatch(Matrix{BigFloat}(undef, 2, n), L),
            SpinMatrix(Matrix{BigFloat}(undef, n, n), L; sₘₐₓ=L, sₘᵢₙ=-L),
            SpinMatrixBatch(Array{BigFloat}(undef, 2, n, n), L; sₘₐₓ=L, sₘᵢₙ=-L),
        )
            @test occursin("#undef", sprint(show, MIME("text/plain"), x))
        end
    end
    @test occursin("#undef", sprint(show, MIME("text/plain"), HCalculator(big(0.3), 3).Hˡ))

    # Assigned elements print as usual, in the order of `Array(w)`
    w = WignerDMatrix(Complex{BigFloat}, 1)
    w[0, 1] = 7
    shown = sprint(show, MIME("text/plain"), w)
    @test occursin("#undef", shown) && occursin("7.0", shown)
    lines = split(shown, '\n')
    @test occursin("7.0", lines[3]) && endswith(rstrip(lines[3]), "im")
end

@testitem "WignerSeries: the blocks of every ℓ" begin
    import SphericalFunctions: WignerSeries, WignerMatrix, ℓ, HalfOddInteger

    blocks = [WignerMatrix(zeros(2ℓ + 1, 2ℓ + 1), ℓ) for ℓ ∈ 0:3]
    s = WignerSeries(blocks, 0, 3)

    @test length(s) == 4
    @test firstindex(s) == 0 && lastindex(s) == 3
    @test collect(keys(s)) == collect(0:3)
    @test axes(s) == (0:3,)
    @test axes(s, 1) == 0:3
    @test axes(s, 2) == Base.OneTo(1)
    @test ndims(s) == 1 && ndims(typeof(s)) == 1
    @test size(s) == (4,)
    @test size(s, 1) == 4 && size(s, 2) == 1
    @test eltype(s) == Pair{Int, eltype(blocks)}  # it iterates as ℓ => block, as a calculator does
    @test parent(s) === blocks && values(s) === blocks

    for ℓᵢ ∈ 0:3
        @test ℓ(s[ℓᵢ]) == ℓᵢ
    end
    @test_throws BoundsError s[4]
    @test_throws BoundsError s[-1]

    # Iteration gives ℓ => block pairs in order, as for a calculator; `first` and `last` are
    # blocks, as indexing is
    @test [(ℓᵢ, ℓ(b)) for (ℓᵢ, b) ∈ s] == [(ℓᵢ, ℓᵢ) for ℓᵢ ∈ 0:3]
    @test all(b === blocks[ℓᵢ + 1] for (ℓᵢ, b) ∈ s)
    @test collect(s) isa Vector{eltype(s)} && length(collect(s)) == 4
    @test first(s) === blocks[1] && last(s) === blocks[end]
    # ... and so do the forms that take a count, and `only`, which `Base` would otherwise build
    # from the pairs
    @test first(s, 2) == blocks[1:2] && last(s, 2) == blocks[3:4]
    @test first(s, 9) == blocks && isempty(last(s, 0))
    @test only(WignerSeries(blocks[1:1], 0, 0)) === blocks[1]
    @test_throws ArgumentError only(s)
    # Iteration covers every block over storage other than a `Vector` too, whose own
    # iteration state is not the integer position that the series counts
    sv = WignerSeries(view(blocks, 1:4), 0, 3)
    @test length(collect(sv)) == 4
    @test [ℓᵢ for (ℓᵢ, _) ∈ sv] == collect(0:3)
    @test all(b === blocks[ℓᵢ + 1] for (ℓᵢ, b) ∈ sv)

    c = copy(s)
    @test c == s
    @test length(similar(s)) == 4
    # `isequal`, `≈` and `hash` compare the range and the blocks, block by block, and agree
    @test isequal(c, s) && c ≈ s && hash(c) == hash(s) && length(Set([s, c])) == 1
    parent(c[2])[1] = NaN
    @test c != s && !isequal(c, s) && isequal(c, copy(c)) && hash(c) == hash(copy(c))
    @test WignerSeries(blocks[2:4], 1, 3) != s && !isequal(WignerSeries(blocks[2:4], 1, 3), s)
    rs = WignerSeries([WignerMatrix(1 .+ rand(2ℓ + 1, 2ℓ + 1), ℓ) for ℓ ∈ 0:3], 0, 3)
    near = copy(rs)
    parent(near[3])[1] += 1e-12
    @test near ≈ rs && !(near == rs) && !isapprox(near, rs; rtol=0, atol=1e-14)
    @test !SphericalFunctions.ishalfinteger(s)

    # An index of another kind or type is told what the series takes
    kind = "The indices of this `WignerSeries` are integers of type `Int`, like 3; got ℓ = "
    @test_throws ArgumentError s[2.0]
    @test_throws kind * "2.0::Float64.  `Float64` is not an index type." s[2.0]
    @test_throws "`Int8` is narrower than `Int`" s[Int8(2)]
    @test_throws "2//1 is a whole number" s[2//1]

    @test occursin("WignerSeries", sprint(show, s))
    @test occursin("WignerSeries", sprint(show, MIME("text/plain"), s))
    # Where the output is limited, as at the REPL, only the first two blocks and the last two
    # of more than four are shown
    long = WignerSeries([WignerMatrix(zeros(2ℓ + 1, 2ℓ + 1), ℓ) for ℓ ∈ 0:9], 0, 9)
    limited = sprint(
        show, MIME("text/plain"), long; context=(:limit => true, :displaysize => (24, 80))
    )
    full = sprint(show, MIME("text/plain"), long)
    @test all(ℓ -> occursin(" ℓ = $ℓ:", full), 0:9)
    @test all(ℓ -> occursin(" ℓ = $ℓ:", limited) == (ℓ ∈ (0, 1, 8, 9)), 0:9)
    @test occursin("⋮", limited) && !occursin("⋮", full)

    # The block count has to match the ℓ range it claims ...
    @test_throws DimensionMismatch WignerSeries(blocks, 0, 4)
    @test_throws "Got 4 blocks, but ℓ ∈ 0:4 needs 5." WignerSeries(blocks, 0, 4)
    @test_throws "Got 4 blocks, but ℓ ∈ 1:3 needs 3." WignerSeries(blocks, 1, 3)
    # ... each block must be the block of its position, of the index type of the bounds ...
    @test_throws ArgumentError WignerSeries(reverse(blocks), 0, 3)
    @test_throws "Block 1 is for ℓ=3, but in a series starting at ℓₘᵢₙ=0 block 1 must be for ℓ=0" WignerSeries(reverse(blocks), 0, 3)
    @test_throws "Block 1 is for ℓ=0, but in a series starting at ℓₘᵢₙ=5" WignerSeries(blocks, 5, 8)
    @test_throws "Block 2 is for ℓ=0" WignerSeries(blocks[[1, 1, 3, 4]], 0, 3)
    @test_throws "HalfOddInteger`, but the bounds of the series are of type `Int64`" WignerSeries([WignerMatrix(zeros(2, 2), 1//2)], 0, 0)
    @test_throws "block 1 is a `Matrix{Float64}`" WignerSeries([zeros(1, 1)], 0, 0)
    # ... and the bounds must be indices of one kind, `Rational`s being converted
    @test_throws "`Int32` is narrower than `Int`" WignerSeries(blocks, Int32(0), Int32(3))
    @test_throws "mixes integers (ℓₘₐₓ) with half-odd-integers (ℓₘᵢₙ)" WignerSeries(blocks, 1//2, 3)
    @test_throws MethodError WignerSeries(blocks, 0.0, 3.0)

    # Half-odd-integer series
    hblocks = [WignerMatrix(zeros(Int(2ℓ + 1), Int(2ℓ + 1)), ℓ)
               for ℓ ∈ HalfOddInteger(1//2):HalfOddInteger(5//2)]
    hs = WignerSeries(hblocks, HalfOddInteger(1//2), HalfOddInteger(5//2))
    @test length(hs) == 3
    @test ℓ(hs[HalfOddInteger(3//2)]) == HalfOddInteger(3//2)
    # ... whose bounds may be written as `Rational`s, which are converted
    hr = WignerSeries(hblocks, 1//2, 5//2)
    @test hr == hs && ℓ(hr[3//2]) == HalfOddInteger(3//2)
    @test SphericalFunctions.ℓₘᵢₙ(hr) === HalfOddInteger(1//2)
    @test SphericalFunctions.ishalfinteger(hs) && hash(hr) == hash(hs) && isequal(hr, hs)
    @test first(hs, 2) == hblocks[1:2] && last(hs) === hblocks[end]
    @test_throws "are half-odd-integers, each a `HalfOddInteger`" hs[1]
end

@testitem "WignerDMatrix and WignerdMatrix: the complex and real aliases" begin
    import SphericalFunctions: WignerDMatrix, WignerdMatrix, WignerMatrix
    import SphericalFunctions: ℓ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, HalfOddInteger

    # Both are aliases for `WignerMatrix`, so `show` names the underlying type
    D = WignerDMatrix(ComplexF64, 2)
    @test D isa WignerMatrix
    @test eltype(D) == ComplexF64
    @test ℓ(D) == 2
    @test size(D) == (5, 5)
    D[1, -2] = 3.0 + 0im
    @test D[1, -2] == 3.0 + 0im
    @test occursin("WignerMatrix", sprint(summary, D))

    d = WignerdMatrix(Float64, 2)
    @test d isa WignerMatrix
    @test eltype(d) == Float64
    @test ℓ(d) == 2
    d[1, -2] = 3.0
    @test d[1, -2] == 3.0

    # A restricted m′ range narrows the first axis only
    Dm = WignerDMatrix(ComplexF64, 3, 1)
    @test m′ₘᵢₙ(Dm) == -1 && m′ₘₐₓ(Dm) == 1
    @test mₘᵢₙ(Dm) == -3 && mₘₐₓ(Dm) == 3
    @test size(Dm) == (3, 7)
    # ... and the m range is given by keyword, in either spelling, with its lower limit
    # defaulting to minus the upper one; the storage is sized to the block
    Dmm = WignerDMatrix(ComplexF64, 3, 1; mₘₐₓ=2)
    @test axes(Dmm) == (-1:1, -2:2) && size(parent(Dmm)) == (3, 5)
    @test axes(WignerdMatrix(Float64, 3, 2; m_max=1, m_min=0)) == (-2:2, 0:1)
    @test axes(WignerdMatrix(Float64, 3//2, 1//2; mₘᵢₙ=-1//2)) == (-1//2:1//2, -1//2:3//2)
    # m′ₘₐₓ is positional in this form, so it is not a keyword; a symmetric m′ range is the
    # only one this form builds
    @test_throws MethodError WignerDMatrix(ComplexF64, 2; m′ₘₐₓ=1)
    @test_throws "m′ₘₐₓ=-1 is less than m′ₘᵢₙ=1" WignerDMatrix(ComplexF64, 2, -1)
    @test_throws "`Int8` is narrower than `Int`" WignerdMatrix(Float64, 2, Int8(1))

    # A `Rational` ℓ builds a half-odd-integer block
    Dh = WignerDMatrix(ComplexF64, 3//2)
    @test ℓ(Dh) == HalfOddInteger(3//2)
    @test size(Dh) == (4, 4)
    Dh[1//2, -1//2] = 1.0 + 2im
    @test Dh[1//2, -1//2] == 1.0 + 2im
    dh = WignerdMatrix(Float64, 3//2)
    @test ℓ(dh) == HalfOddInteger(3//2)
    @test size(dh) == (4, 4)

    # Wrapping existing storage, and the cross-type errors that catch the obvious mistake
    @test WignerDMatrix(zeros(ComplexF64, 5, 5), 2) isa WignerMatrix
    @test WignerdMatrix(zeros(5, 5), 2) isa WignerMatrix
    @test_throws "only supports complex types" WignerDMatrix(zeros(5, 5), 2)
    @test_throws "Perhaps you meant to use WignerdMatrix" WignerDMatrix(zeros(5, 5), 2)
    @test_throws "only supports real types" WignerdMatrix(zeros(ComplexF64, 5, 5), 2)
    @test_throws "Perhaps you meant to use WignerDMatrix" WignerdMatrix(zeros(ComplexF64, 5, 5), 2)

    # The storage-wrapping forms take the limits of `WignerMatrix`, in either spelling
    @test axes(WignerDMatrix(zeros(ComplexF64, 5, 5), 2; mp_max=1, m_max=0, m_min=-2)) ==
        (-1:1, -2:0)
    @test axes(WignerdMatrix(zeros(4, 4), 3//2; m′ₘₐₓ=1//2)) == (-1//2:1//2, -3//2:3//2)

    # An ℓ that is neither integer nor half-odd-integer is refused before the storage is
    # sized, rather than with a bare `InexactError`.  (The message is matched, not just
    # `Exception`, which the `InexactError` would satisfy too.)
    @test_throws ArgumentError WignerDMatrix(ComplexF64, 5//3)
    @test_throws "5//3 is neither an integer nor a half-odd-integer" WignerDMatrix(ComplexF64, 5//3)
    @test_throws "5//3 is neither an integer nor a half-odd-integer" WignerdMatrix(Float64, 5//3)
    @test_throws "only supports complex types" WignerDMatrix(zeros(5, 5), 2.5)
end

@testitem "SpinMatrix: the generic search functions agree with the dense matrix, or refuse" begin
    import SphericalFunctions: SpinMatrix, SpinMatrixBatch, sYlm, sYlmCalculator, recurrence!
    import SphericalFunctions: spins, sₘᵢₙ, sₘₐₓ, HalfOddInteger
    using Quaternionic: Rotor

    # A `SpinMatrix` iterates over every element, (s, m) in column-major order, so `pairs`,
    # and everything built on it — `findmax`, `argmax`, `findall`, `findfirst` — must cover
    # every element as well, or not be defined at all.  A `MethodError` is an acceptable
    # answer; one that looks at some of the elements, or reports a spin weight as the
    # position of an element, is not.
    function agrees_or_refuses(agrees, f)
        result = try
            f()
        catch e
            e isa MethodError && return true
            rethrow()
        end
        agrees(result)
    end
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    for b ∈ (
        SpinMatrix(reshape(collect(1.0:15.0), 3, 5), 2; sₘₐₓ=1, sₘᵢₙ=-1),
        sYlm(R, 3, -1:1)[2],
    )
        M = Matrix(b)
        @test agrees_or_refuses(r -> r == findmax(abs, M), () -> findmax(abs, b))
        @test agrees_or_refuses(r -> r == findmin(abs, M), () -> findmin(abs, b))
        @test agrees_or_refuses(r -> r == argmax(abs, M), () -> argmax(abs, b))
        @test agrees_or_refuses(r -> r == argmax(M), () -> argmax(b))
        large(x) = abs(x) > 0.1
        @test agrees_or_refuses(r -> r == findall(large, M), () -> findall(large, b))
        @test agrees_or_refuses(r -> r == findfirst(>(5) ∘ abs, M), () -> findfirst(>(5) ∘ abs, b))
        @test agrees_or_refuses(r -> length(r) == length(M), () -> collect(pairs(b)))
        # The reductions over the elements themselves are unaffected
        @test maximum(abs, b) == maximum(abs, M)
        @test sum(b) == sum(M)
    end

    # The spin axis is `spins(b)`, as it is for the calculator that produced the block
    b = SpinMatrix(zeros(3, 5), 2; sₘₐₓ=1, sₘᵢₙ=-1)
    @test collect(spins(b)) == [-1, 0, 1]
    @test collect(spins(SpinMatrixBatch(zeros(4, 3, 5), 2; sₘₐₓ=1, sₘᵢₙ=-1))) == [-1, 0, 1]
    @test collect(spins(SpinMatrix(zeros(2, 5), 2; sₘₐₓ=3, sₘᵢₙ=2))) == [2, 3]
    bₕ = SpinMatrix(zeros(2, 4), 3//2; sₘₐₓ=1//2, sₘᵢₙ=-1//2)
    @test collect(spins(bₕ)) == [HalfOddInteger(-1//2), HalfOddInteger(1//2)]
    @test first(spins(bₕ)) == sₘᵢₙ(bₕ) && last(spins(bₕ)) == sₘₐₓ(bₕ)
    calc = sYlmCalculator(R, 3, -1:1)
    @test collect(spins(recurrence!(calc, 2))) == collect(spins(calc))
end

@testitem "Wigner containers: the ASCII names of the accessors" begin
    import SphericalFunctions: ℓ, ℓₘᵢₙ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, sₘₐₓ, sₘᵢₙ
    import SphericalFunctions: ell, ell_min, mp_max, mp_min, m_max, m_min, s_max, s_min
    import SphericalFunctions: WignerMatrix, SpinMatrix

    # Each is the same function as its Unicode name, spelled as the keyword aliases are
    for (a, u) in ((ell, ℓ), (ell_min, ℓₘᵢₙ), (mp_max, m′ₘₐₓ), (mp_min, m′ₘᵢₙ),
                   (m_max, mₘₐₓ), (m_min, mₘᵢₙ), (s_max, sₘₐₓ), (s_min, sₘᵢₙ))
        @test a === u
    end
    w = WignerMatrix(zeros(5, 5), 2; mp_max=1, m_min=0)
    @test (ell(w), ell_min(w), mp_max(w), mp_min(w), m_max(w), m_min(w)) == (2, 0, 1, -1, 2, 0)
    b = SpinMatrix(zeros(3, 5), 2; s_max=1, s_min=-1)
    @test (s_max(b), s_min(b)) == (1, -1)
    # The spellings without the underscore are not defined
    for name in (:ellmin, :ellmax, :mpmax, :mpmin, :mmax, :mmin, :smax, :smin)
        @test !isdefined(SphericalFunctions, name)
    end
end

@testitem "Wigner containers: the labels in comparison, hashing, bounds and broadcasting" setup=[RefusalChecks] begin
    import SphericalFunctions: WignerMatrix, WignerMatrixBatch, DegreeBlock, DegreeBlockBatch,
        SpinMatrix, SpinMatrixBatch, isbatched, HalfOddInteger, m′ₘₐₓ, m′ₘᵢₙ, ℓₘᵢₙ, ℓ

    h = HalfOddInteger
    # For each block, the same storage under two sets of labels, of the same shape: another
    # ℓ, or another run of spin weights
    function twins(L)
        n = Int(2L + 1)
        A1, A2, A3 = collect(1.0:n), reshape(collect(1.0:3n), 3, n), reshape(collect(1.0:n^2), n, n)
        B2, B3 = reshape(collect(1.0:2n), 2, n), reshape(collect(1.0:2n^2), 2, n, n)
        B4 = reshape(collect(1.0:6n), 2, 3, n)
        L′ = L + 1
        [
            (WignerMatrix(copy(A3), L), WignerMatrix(copy(A3), L′; m′ₘₐₓ=L, mₘₐₓ=L)),
            (WignerMatrixBatch(copy(B3), L), WignerMatrixBatch(copy(B3), L′; m′ₘₐₓ=L, mₘₐₓ=L)),
            (DegreeBlock(copy(A1), L), DegreeBlock(copy(A1), L′; mₘₐₓ=L)),
            (DegreeBlockBatch(copy(B2), L), DegreeBlockBatch(copy(B2), L′; mₘₐₓ=L)),
            (SpinMatrix(copy(A2), L; sₘₐₓ=L, sₘᵢₙ=L-2), SpinMatrix(copy(A2), L; sₘₐₓ=L′, sₘᵢₙ=L-1)),
            (SpinMatrixBatch(copy(B4), L; sₘₐₓ=L, sₘᵢₙ=L-2),
             SpinMatrixBatch(copy(B4), L; sₘₐₓ=L′, sₘᵢₙ=L-1)),
        ]
    end
    labels = "must have the same labels"
    for L ∈ (2, h(3//2)), (w, other) ∈ twins(L)
        c = copy(w)
        # `==`, `isequal`, `≈` and `hash` agree, and count the labels as well as the numbers
        @test c == w && isequal(c, w) && c ≈ w && hash(c) == hash(w)
        @test haskey(Dict(w => 1), c) && length(unique([w, c, other])) == 2
        @test Array(other) == Array(w)
        @test other != w && !isequal(other, w) && !(other ≈ w)
        # The numbers compare as arrays do: NaN is `isequal` to itself but not `==`, the signed
        # zeros are `==` but not `isequal`, and a real and a complex block of the same values
        # are equal, with equal hashes
        wn = copy(w)
        parent(wn)[1] = NaN
        @test isequal(wn, copy(wn)) && wn != wn && hash(wn) == hash(copy(wn))
        w₊, w₋ = copy(w), copy(w)
        parent(w₊)[1], parent(w₋)[1] = 0.0, -0.0
        @test w₊ == w₋ && !isequal(w₊, w₋)
        wc = similar(w, ComplexF64)
        wc .= w
        @test wc == w && isequal(wc, w) && hash(wc) == hash(w)
        @test isapprox(w, similar(w) .= Array(w) .+ 1e-12; atol=1e-10)
        # `Array` and `collect` copy the block, whatever the storage
        @test Array(w) isa Array{Float64, ndims(w)} && collect(w) == Array(w)
        # `checkbounds` tests natural indices, as `getindex` does
        first_index, past_end = map(first, axes(w)), map(a -> last(a) + 1, axes(w))
        @test checkbounds(Bool, w, first_index...) && !checkbounds(Bool, w, past_end...)
        @test checkbounds(w, first_index...) === nothing
        @test_throws BoundsError checkbounds(w, past_end...)
        @test isbatched(w) == (ndims(w) == 3 || w isa DegreeBlockBatch)
        # Broadcasting gives a plain array, and pairs elements by position, so containers
        # combined in one broadcast, or written into one another, must have the same labels;
        # a plain array is paired as it is
        @test w .+ c isa Array{Float64, ndims(w)} && w .+ c == 2 .* Array(w)
        @test (w .== c) isa BitArray
        @test refuses(() -> w .- other, ArgumentError, labels)
        @test refuses(() -> copy(w) .= other, ArgumentError, labels)
        @test refuses(() -> copy(w) .= 2 .* other .+ w, ArgumentError, labels)
        @test refuses(() -> zeros(size(w)) .= w .+ other, ArgumentError, labels)
        @test (zeros(size(w)) .= w .+ Array(other)) == 2 .* Array(w)
        @test (copy(w) .= Array(other)) == w
        d = copy(w)
        d .= 2 .* d
        @test Array(d) == 2 .* Array(w)
    end
    # A block with an m′ axis and one with a spin axis are different blocks, even where their
    # axes are the same ranges
    A = rand(5, 5)
    w, b = WignerMatrix(copy(A), 2), SpinMatrix(copy(A), 2; sₘₐₓ=2, sₘᵢₙ=-2)
    @test axes(w) == axes(b) && Array(w) == Array(b)
    @test w != b && !isequal(w, b) && !(w ≈ b)
    @test refuses(() -> w .+ b, ArgumentError, labels)

    # The blocks with a single m axis have no m′, and say so
    for v ∈ (DegreeBlock(zeros(5), 2), DegreeBlockBatch(zeros(2, 5), 2),
             SpinMatrix(zeros(3, 5), 2; sₘₐₓ=1, sₘᵢₙ=-1),
             SpinMatrixBatch(zeros(2, 3, 5), 2; sₘₐₓ=1, sₘᵢₙ=-1))
        @test refuses(() -> m′ₘₐₓ(v), ArgumentError, "has no m′ axis")
        @test refuses(() -> m′ₘᵢₙ(v), ArgumentError, "has no m′ axis")
    end
    # `ℓₘᵢₙ` is defined for the index types, their values and the containers, and for nothing
    # else
    @test ℓₘᵢₙ(3) === 0 && ℓₘᵢₙ(h(5//2)) === h(1//2) && ℓₘᵢₙ(Int) === 0
    @test ℓₘᵢₙ(WignerMatrix(zeros(4, 4), 3//2)) === h(1//2)
    @test_throws MethodError ℓₘᵢₙ([1, 2])
    @test !applicable(ℓₘᵢₙ, "x") && !applicable(ℓₘᵢₙ, 1//2)

    # A half-integer index may be written as a `Rational{Int}` with denominator 2 or as a
    # `HalfOddInteger`, either way for each index, reading and writing; any other `Rational`
    # is refused with the reason
    wh = WignerMatrix(reshape(collect(1.0:16), 4, 4), 3//2)
    bh = WignerMatrixBatch(reshape(collect(1.0:32), 2, 4, 4), 3//2)
    sh = SpinMatrix(reshape(collect(1.0:8), 2, 4), 3//2; sₘₐₓ=1//2, sₘᵢₙ=-1//2)
    sbh = SpinMatrixBatch(reshape(collect(1.0:16), 2, 2, 4), 3//2; sₘₐₓ=1//2, sₘᵢₙ=-1//2)
    for (x, pre) ∈ ((wh, ()), (bh, (2,)), (sh, ()), (sbh, (2,)))
        expected = x[pre..., h(1//2), h(-3//2)]
        @test x[pre..., 1//2, -3//2] == x[pre..., 1//2, h(-3//2)] == x[pre..., h(1//2), -3//2] ==
            expected
        x[pre..., 1//2, h(-3//2)] = 100.0
        @test x[pre..., h(1//2), h(-3//2)] == 100.0
        x[pre..., h(1//2), -3//2] = 200.0
        @test x[pre..., h(1//2), h(-3//2)] == 200.0
        @test refuses(() -> x[pre..., 1//3, 1//2], ArgumentError, "1//3 is neither")
        @test refuses(
            () -> x[pre..., Int8(1)//Int8(2), 1//2], ArgumentError, "`Rational{Int8}` is not"
        )
        @test refuses(() -> x[pre..., 1//2, 1], ArgumentError, "are half-odd-integers")
        @test refuses(() -> (x[pre..., 1, 1//2] = 0.0), ArgumentError, " = 1::Int64")
    end
    @test sh[1//2, :] == sh[h(1//2), :] && sbh[:, 1//2, :] == sbh[:, h(1//2), :]
    @test refuses(() -> sh[1, :], ArgumentError, "are half-odd-integers")
    vh = DegreeBlock(collect(1.0:4), 3//2)
    @test vh[1//2] == vh[h(1//2)] && DegreeBlockBatch(ones(2, 4), 3//2)[1, -1//2] == 1
    @test refuses(() -> vh[1], ArgumentError, "indices of this `DegreeBlock` are half-odd")

    # The index of an integer container must be an `Int`: an index of the other kind, or of a
    # type that the index methods refuse, is refused with the reason, for reading and writing
    wi = WignerMatrix(reshape(collect(1.0:9), 3, 3), 1)
    bi = WignerMatrixBatch(reshape(collect(1.0:18), 2, 3, 3), 1)
    si = SpinMatrix(reshape(collect(1.0:6), 2, 3), 1; sₘₐₓ=1, sₘᵢₙ=0)
    sbi = SpinMatrixBatch(reshape(collect(1.0:12), 2, 2, 3), 1; sₘₐₓ=1, sₘᵢₙ=0)
    for (x, pre) ∈ ((wi, ()), (bi, (2,)), (si, ()), (sbi, (2,)))
        name = nameof(typeof(x))
        @test refuses(() -> x[pre..., Int32(1), 0], ArgumentError, "`Int32` is narrower")
        @test refuses(() -> x[pre..., Int32(1), 0], ArgumentError, "this `$name` are integers")
        @test refuses(() -> x[pre..., 1, 1//1], ArgumentError, "1//1 is a whole number")
        @test refuses(() -> x[pre..., 1//2, 0], ArgumentError, "are integers of type `Int`")
        @test refuses(() -> (x[pre..., 0, Int16(1)] = 1.0), ArgumentError, "`Int16` is narrower")
        @test x[pre..., 1, 0] == x[pre..., 1, 0]  # the natural form is untouched
    end
    @test refuses(() -> si[Int32(0), :], ArgumentError, "`Int32` is narrower")
    @test refuses(() -> sbi[:, Int32(0), :], ArgumentError, "`Int32` is narrower")
    vi = DegreeBlock(collect(1.0:3), 1)
    @test refuses(() -> vi[Int8(0)], ArgumentError, "`Int8` is narrower")
    @test refuses(() -> (vi[1//2] = 0.0), ArgumentError, "are integers of type `Int`")
    @test refuses(
        () -> DegreeBlockBatch(ones(2, 3), 1)[1, true], ArgumentError, "A `Bool` is not an index"
    )

    # The extent in a refusal is written with its signs, and an axis is shown as a unit range
    @test_throws "m′ₘₐₓ-m′ₘᵢₙ+1=2-(-2)+1=5; it is 2." WignerMatrix(zeros(2, 2), 2)
    @test repr(axes(WignerMatrix(zeros(5, 5), 2))) == "(-2:2, -2:2)"
    @test repr(axes(WignerMatrix(zeros(4, 4), 3//2))) == "(-3//2:3//2, -3//2:3//2)"
    @test repr(axes(WignerMatrixBatch(zeros(2, 5, 5), 2))) == "(1:2, -2:2, -2:2)"
end
