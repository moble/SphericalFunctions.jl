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
    import SphericalFunctions: WignerRange, half_integer, HalfOddInteger

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
    h = WignerRange(half_integer(-3//2):half_integer(5//2))
    @test first(h) == half_integer(-3//2) && last(h) == half_integer(5//2)
    @test step(h) == 1
    @test length(h) == 5
    @test firstindex(h) == 1 && lastindex(h) == 5
    @test h[1] == half_integer(-3//2)
    @test h[5] == half_integer(5//2)
    @test_throws BoundsError h[0]
    @test_throws BoundsError h[6]

    # An empty half-odd range has length 0 rather than a negative length
    @test length(WignerRange(half_integer(3//2):half_integer(1//2))) == 0

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
    @test half_integer(1//2) ∉ r
    @test 0.5 ∉ r
    @test half_integer(1//2) ∈ h && half_integer(-3//2) ∈ h
    @test half_integer(7//2) ∉ h
    @test 1 ∉ h
end

@testitem "Wigner containers: `validate_index_ranges` refuses every bad range" begin
    import SphericalFunctions: validate_index_ranges, half_integer

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
    h(x) = half_integer(x)
    @test validate_index_ranges(h(5//2), h(5//2), h(-5//2)) === nothing
    @test_throws "too small for this index type" validate_index_ranges(h(5//2), h(-1//2), h(-5//2))
    @test_throws "too large for this index type" validate_index_ranges(h(5//2), h(5//2), h(1//2))
end

@testitem "WignerMatrix: the container interface" begin
    import SphericalFunctions: WignerMatrix, ℓ, ℓₘᵢₙ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ
    import SphericalFunctions: ishalfinteger, half_integer

    for L ∈ (2, half_integer(3//2))
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

        # `==` compares ℓ, axes and values
        @test w != similar(w, ComplexF64)
        @test copy(w) == w

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
    @test ℓ(w) == half_integer(3//2)
    w[1//2, -1//2] = 7.0
    @test w[1//2, -1//2] == 7.0
    @test w[half_integer(1//2), half_integer(-1//2)] == 7.0

    # Storage too small for the requested ranges is refused, per dimension
    @test_throws "first dimension" WignerMatrix(zeros(3, 5), 2)
    @test_throws "second dimension" WignerMatrix(zeros(5, 3), 2)
end

@testitem "WignerMatrixBatch: the container interface" begin
    import SphericalFunctions: WignerMatrixBatch, WignerMatrix, ℓ, Nᵣ
    import SphericalFunctions: m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, half_integer

    for L ∈ (2, half_integer(3//2))
        n, N = Int(2L + 1), 3
        A = reshape(collect(1.0:N*n^2), N, n, n)
        w = WignerMatrixBatch(copy(A), L)

        @test ℓ(w) == L && Nᵣ(w) == N
        @test parent(w) == A
        @test size(w) == (N, n, n)

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

        # `Array` is the 3-d array in `[iᵣ, m′, m]` order; `Matrix` is refused as ambiguous
        @test Array(w) == permutedims(
            reshape([w[i, m′, m] for i ∈ 1:N, m′ ∈ m′ₘᵢₙ(w):m′ₘₐₓ(w), m ∈ mₘᵢₙ(w):mₘₐₓ(w)],
                    N, n, n),
            (1, 2, 3)
        )
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
    w = WignerMatrixBatch(zeros(2, 4, 4), half_integer(3//2))
    w[1, 1//2, -1//2] = 4.0
    @test w[1, 1//2, -1//2] == 4.0
    @test w[1, half_integer(1//2), half_integer(-1//2)] == 4.0
end

@testitem "DegreeBlock and DegreeBlockBatch: the container interface" begin
    import SphericalFunctions: DegreeBlock, DegreeBlockBatch, ℓ, Nᵣ, mₘᵢₙ, mₘₐₓ, half_integer

    for L ∈ (2, half_integer(3//2))
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
    @test ℓ(v) == half_integer(3//2)
    v[1//2] = 6.0
    @test v[1//2] == 6.0
    @test v[half_integer(1//2)] == 6.0
end

@testitem "SpinMatrix and SpinMatrixBatch: the container interface" begin
    import SphericalFunctions: SpinMatrix, SpinMatrixBatch, DegreeBlock
    import SphericalFunctions: ℓ, Nᵣ, sₘᵢₙ, sₘₐₓ, mₘᵢₙ, mₘₐₓ, half_integer

    for L ∈ (2, half_integer(3//2))
        n = Int(2L + 1)
        smin, smax = -L, L
        ns = Int(smax - smin) + 1
        A = reshape(collect(1.0:ns*n), ns, n)
        b = SpinMatrix(copy(A), L; sₘₐₓ=smax, sₘᵢₙ=smin)

        @test ℓ(b) == L
        @test sₘᵢₙ(b) == smin && sₘₐₓ(b) == smax
        @test collect(keys(b)) == collect(smin:smax)
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

    # A `Rational` ℓ is converted, along with the keyword spin bounds The spin bounds are
    # index keywords like any other, so `half_integer_kwargs` has to convert them: without
    # `:sₘₐₓ`/`:sₘᵢₙ` in `INDEX_KEYWORDS` the `Rational`-ℓ constructor hands an unconverted
    # `Rational` to a method typed on `IT` and throws a `TypeError`.  These two containers
    # are the only ones that take spin bounds.
    br = SpinMatrix(zeros(4, 4), 3//2; sₘₐₓ=3//2, sₘᵢₙ=-3//2)
    @test br isa SpinMatrix
    @test ℓ(br) == half_integer(3//2)
    @test sₘₐₓ(br) == half_integer(3//2) && sₘᵢₙ(br) == half_integer(-3//2)
    bbr = SpinMatrixBatch(zeros(2, 4, 4), 3//2; sₘₐₓ=3//2, sₘᵢₙ=-3//2)
    @test bbr isa SpinMatrixBatch
    @test ℓ(bbr) == half_integer(3//2)
    @test sₘₐₓ(bbr) == half_integer(3//2) && sₘᵢₙ(bbr) == half_integer(-3//2)

    # ... and it agrees with the `HalfOddInteger` spelling, which is what they store
    b = SpinMatrix(zeros(4, 4), half_integer(3//2);
                   sₘₐₓ=half_integer(3//2), sₘᵢₙ=half_integer(-3//2))
    @test ℓ(b) == half_integer(3//2)
    b[1//2, -1//2] = 5.0
    @test b[1//2, -1//2] == 5.0
    bb = SpinMatrixBatch(zeros(2, 4, 4), half_integer(3//2);
                         sₘₐₓ=half_integer(3//2), sₘᵢₙ=half_integer(-3//2))
    @test ℓ(bb) == half_integer(3//2)
    bb[1, 1//2, -1//2] = 8.0
    @test bb[1, 1//2, -1//2] == 8.0
end

@testitem "WignerSeries: the blocks of every ℓ" begin
    import SphericalFunctions: WignerSeries, WignerMatrix, ℓ, half_integer

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
    @test eltype(s) == eltype(blocks)
    @test parent(s) === blocks

    for ℓᵢ ∈ 0:3
        @test ℓ(s[ℓᵢ]) == ℓᵢ
    end
    @test_throws BoundsError s[4]
    @test_throws BoundsError s[-1]

    # Iteration gives the blocks in order
    @test [ℓ(b) for b ∈ s] == collect(0:3)

    c = copy(s)
    @test c == s
    @test length(similar(s)) == 4

    @test occursin("WignerSeries", sprint(show, s))
    @test occursin("WignerSeries", sprint(show, MIME("text/plain"), s))

    # The block count has to match the ℓ range it claims
    @test_throws "blocks" WignerSeries(blocks, 0, 4)
    @test_throws "blocks" WignerSeries(blocks, 1, 3)

    # Half-odd-integer series
    hblocks = [WignerMatrix(zeros(Int(2ℓ + 1), Int(2ℓ + 1)), ℓ)
               for ℓ ∈ half_integer(1//2):half_integer(5//2)]
    hs = WignerSeries(hblocks, half_integer(1//2), half_integer(5//2))
    @test length(hs) == 3
    @test ℓ(hs[half_integer(3//2)]) == half_integer(3//2)
end

@testitem "WignerDMatrix and WignerdMatrix: the complex and real aliases" begin
    import SphericalFunctions: WignerDMatrix, WignerdMatrix, WignerMatrix
    import SphericalFunctions: ℓ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, half_integer

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

    # A `Rational` ℓ builds a half-odd-integer block
    Dh = WignerDMatrix(ComplexF64, 3//2)
    @test ℓ(Dh) == half_integer(3//2)
    @test size(Dh) == (4, 4)
    Dh[1//2, -1//2] = 1.0 + 2im
    @test Dh[1//2, -1//2] == 1.0 + 2im
    dh = WignerdMatrix(Float64, 3//2)
    @test ℓ(dh) == half_integer(3//2)
    @test size(dh) == (4, 4)

    # Wrapping existing storage, and the cross-type errors that catch the obvious mistake
    @test WignerDMatrix(zeros(ComplexF64, 5, 5), 2) isa WignerMatrix
    @test WignerdMatrix(zeros(5, 5), 2) isa WignerMatrix
    @test_throws "only supports complex types" WignerDMatrix(zeros(5, 5), 2)
    @test_throws "Perhaps you meant to use WignerdMatrix" WignerDMatrix(zeros(5, 5), 2)
    @test_throws "only supports real types" WignerdMatrix(zeros(ComplexF64, 5, 5), 2)
    @test_throws "Perhaps you meant to use WignerDMatrix" WignerdMatrix(zeros(ComplexF64, 5, 5), 2)

    # An ℓ that is neither integer nor half-odd-integer is refused with the denominator rule
    # rather than a bare `InexactError` from sizing the storage
    @test_throws Exception WignerDMatrix(ComplexF64, 5//3)
    @test_throws Exception WignerdMatrix(Float64, 5//3)
end
