# Tests of `ModeWeights`, the wrapper for a vector of mode weights in the canonical ordering
#
#     [ f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ ],
#
# together with the spin weight `s` (see `src/mode_weights/mode_weights.jl`).  The oracles
# are the closed-form indexing functions `Ysize`, `Yindex` and `Yrange` (tested in
# `indexing.jl`), the operator matrices of `src/mode_weights/operators.jl` applied to the
# raw storage (tested against explicit differential operators in `test/operators.jl`), and —
# for the evaluation `w(R)` — `sYlm_closed_form(s, ℓ, m, θ, ϕ)` from the `Utilities` module,
# which is the explicit sum from `docs/src/30-conventions/01-summary.md` and shares no code
# with the package.
#
# `ℓₘᵢₙ` and `ℓₘₐₓ` are used as loop variables throughout, so the accessor functions of the
# same names are always called qualified, as `SphericalFunctions.ℓₘᵢₙ(w)`.

@testitem "ModeWeights construction" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yrange
    import DoubleFloats: Double64
    import OffsetArrays: OffsetArray
    import Random

    rng = Random.Xoshiro(20260910)

    # The four-argument form: the parameters come back unchanged, the storage is the very
    # same array (not a copy), and the array-like queries describe that storage.  The range
    # ℓₘₐₓ = ℓₘᵢₙ - 1 is the empty container.
    for s in -3:3, ℓₘᵢₙ in 0:3, ℓₘₐₓ in ℓₘᵢₙ-1:7
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
        @test w isa ModeWeights{ComplexF64}
        # This is an `AbstractModeContainer`, not an `AbstractVector`: `array_view` is the
        # route to the flat storage, and array semantics come with it.
        @test w isa SphericalFunctions.AbstractModeContainer{ComplexF64, Int}
        @test !(w isa AbstractVector)
        @test array_view(w) === parent(w)
        @test parent(w) === data
        @test spin(w) == s
        @test SphericalFunctions.ℓₘᵢₙ(w) == ℓₘᵢₙ
        @test SphericalFunctions.ℓₘₐₓ(w) == ℓₘₐₓ
        @test eltype(w) == ComplexF64
        @test length(w) == n
        @test size(w) == (n,)
        @test axes(w) == (Base.OneTo(n),)
        @test isempty(w) == (ℓₘₐₓ == ℓₘᵢₙ - 1)
        @test modes(w) == Yrange(ℓₘᵢₙ, ℓₘₐₓ)
        @test length(modes(w)) == n
    end

    # With only `ℓₘᵢₙ` given, `ℓₘₐₓ` is deduced from the length ...
    for ℓₘᵢₙ in 0:3, ℓₘₐₓ in ℓₘᵢₙ-1:8
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, n)
        for s in -3:3
            w = ModeWeights(data, s; ℓₘᵢₙ)
            @test parent(w) === data
            @test spin(w) == s
            @test SphericalFunctions.ℓₘᵢₙ(w) == ℓₘᵢₙ
            @test SphericalFunctions.ℓₘₐₓ(w) == ℓₘₐₓ
            @test length(w) == n
        end
        # ... `ℓₘᵢₙ` itself defaults to `abs(s)` ...
        for s in (-ℓₘᵢₙ, ℓₘᵢₙ)
            w = ModeWeights(data, s)
            @test spin(w) == s
            @test SphericalFunctions.ℓₘᵢₙ(w) == abs(s) == ℓₘᵢₙ
            @test SphericalFunctions.ℓₘₐₓ(w) == ℓₘₐₓ
            @test parent(w) === data
        end
    end
    # ... and `s` defaults to 0
    for ℓₘₐₓ in -1:8
        data = randn(rng, ComplexF64, (ℓₘₐₓ + 1)^2)
        w = ModeWeights(data)
        @test spin(w) == 0
        @test SphericalFunctions.ℓₘᵢₙ(w) == 0
        @test SphericalFunctions.ℓₘₐₓ(w) == ℓₘₐₓ
        @test parent(w) === data
    end
    # Every length that is not `Ysize(ℓₘᵢₙ, ℓₘₐₓ)` for some `ℓₘₐₓ ≥ ℓₘᵢₙ - 1` is rejected
    for ℓₘᵢₙ in 0:3
        valid = Set(Ysize(ℓₘᵢₙ, ℓₘₐₓ) for ℓₘₐₓ in ℓₘᵢₙ-1:12)
        for n in 0:Ysize(ℓₘᵢₙ, 12)
            if n ∈ valid
                @test ModeWeights(zeros(n), 0; ℓₘᵢₙ) isa ModeWeights
                @test ModeWeights(zeros(n), ℓₘᵢₙ) isa ModeWeights
            else
                @test_throws ArgumentError ModeWeights(zeros(n), 0; ℓₘᵢₙ)
                @test_throws ArgumentError ModeWeights(zeros(n), ℓₘᵢₙ)
            end
        end
    end
    @test_throws "length" ModeWeights(zeros(2), 0)
    @test_throws "length" ModeWeights(zeros(5), 1)

    # Uninitialized storage of a given element type, with and without `ℓₘᵢₙ`
    for T in (Float32, Float64, ComplexF64, Complex{Double64})
        for s in -2:2, ℓₘᵢₙ in 0:2, ℓₘₐₓ in ℓₘᵢₙ-1:5
            w = ModeWeights{T}(undef, s, ℓₘᵢₙ, ℓₘₐₓ)
            @test w isa ModeWeights{T}
            @test eltype(w) == T
            @test parent(w) isa Vector{T}
            @test length(w) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            @test spin(w) == s
            @test SphericalFunctions.ℓₘᵢₙ(w) == ℓₘᵢₙ
            @test SphericalFunctions.ℓₘₐₓ(w) == ℓₘₐₓ
            # The storage is writable
            if !isempty(w)
                w[1] = 1
                @test parent(w)[1] == 1
            end
        end
        for s in -2:2, ℓₘₐₓ in abs(s)-1:5
            w = ModeWeights{T}(undef, s, ℓₘₐₓ)
            @test w isa ModeWeights{T}
            @test length(w) == Ysize(abs(s), ℓₘₐₓ)
            @test spin(w) == s
            @test SphericalFunctions.ℓₘᵢₙ(w) == abs(s)
            @test SphericalFunctions.ℓₘₐₓ(w) == ℓₘₐₓ
        end
    end

    # Invalid parameters are `ArgumentError`s in every constructor form: a negative ℓₘᵢₙ ...
    for ℓₘᵢₙ in -3:-1
        @test_throws ArgumentError ModeWeights(zeros(3), 0, ℓₘᵢₙ, 1)
        @test_throws ArgumentError ModeWeights(zeros(3), 0; ℓₘᵢₙ)
        @test_throws ArgumentError ModeWeights{Float64}(undef, 0, ℓₘᵢₙ, 1)
    end
    @test_throws "ℓₘᵢₙ" ModeWeights(zeros(3), 0, -1, 1)
    @test_throws "ℓₘᵢₙ" ModeWeights(zeros(3), 0; ℓₘᵢₙ=-1)
    # ... ℓₘₐₓ < ℓₘᵢₙ - 1 (GitHub issue #52: `Ysize` would be negative) ...
    for ℓₘᵢₙ in 0:4, ℓₘₐₓ in -3:ℓₘᵢₙ-2
        @test_throws ArgumentError ModeWeights(Float64[], 0, ℓₘᵢₙ, ℓₘₐₓ)
        @test_throws ArgumentError ModeWeights{Float64}(undef, 0, ℓₘᵢₙ, ℓₘₐₓ)
    end
    @test_throws "ℓₘₐₓ" ModeWeights(zeros(3), 0, 3, 0)
    for s in -3:3, ℓₘₐₓ in -3:abs(s)-2
        @test_throws ArgumentError ModeWeights{Float64}(undef, s, ℓₘₐₓ)
    end
    # ... a length that does not match the given range ...
    for ℓₘᵢₙ in 0:3, ℓₘₐₓ in ℓₘᵢₙ-1:5
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        for bad in (n - 1, n + 1, 2n + 1)
            bad ≥ 0 || continue
            @test_throws ArgumentError ModeWeights(zeros(bad), 0, ℓₘᵢₙ, ℓₘₐₓ)
        end
    end
    @test_throws "length" ModeWeights(zeros(3), 0, 1, 2)
    # ... and storage that is not 1-based (the ordering is defined by 1-based positions)
    @test_throws ArgumentError ModeWeights(OffsetArray(zeros(3), 0:2), 1)
    @test_throws ArgumentError ModeWeights(OffsetArray(zeros(3), 0:2), 1, 1, 1)
    @test_throws ArgumentError ModeWeights(OffsetArray(zeros(3), -1:1), 0; ℓₘᵢₙ=1)

    # The empty container, ℓₘₐₓ = ℓₘᵢₙ - 1, in every form
    for s in -2:2, ℓₘᵢₙ in 0:3
        w = ModeWeights(ComplexF64[], s, ℓₘᵢₙ, ℓₘᵢₙ - 1)
        @test length(w) == 0
        @test isempty(w)
        @test isempty(modes(w))
        @test spin(w) == s
        @test SphericalFunctions.ℓₘᵢₙ(w) == ℓₘᵢₙ
        @test SphericalFunctions.ℓₘₐₓ(w) == ℓₘᵢₙ - 1
        @test sum(w) == 0
        w′ = ModeWeights(ComplexF64[], s; ℓₘᵢₙ)
        @test SphericalFunctions.ℓₘₐₓ(w′) == ℓₘᵢₙ - 1
        @test isempty(ModeWeights{Float32}(undef, s, ℓₘᵢₙ, ℓₘᵢₙ - 1))
    end
    @test SphericalFunctions.ℓₘₐₓ(ModeWeights(Float64[])) == -1
    @test SphericalFunctions.ℓₘₐₓ(ModeWeights(Float64[], 2)) == 1
    @test isempty(ModeWeights{Float64}(undef, 2, 1))

    # Integer indices are `Int`s.  Any other integer type is refused, in every form, with a
    # sentence saying why and how to write it: a narrow type overflows in `ℓ^2`, an unsigned
    # one wraps around in `-m`, and a wider one would reach code that stores `Int`s.
    data = zeros(Ysize(1, 3))
    @test ModeWeights(data, Int64(-1), Int64(1), Int64(3)) isa ModeWeights{Float64, Int}
    for (x, why) in (
        (Int8(1), "`Int8` is narrower than `Int`"), (Int32(1), "`Int32` is narrower than `Int`"),
        (UInt(1), "`UInt64` is unsigned"), (Int128(1), "`Int128` is wider than `Int`"),
        (big(1), "`BigInt` is wider than `Int`"), (true, "A `Bool` is not an index"),
    )
        @test_throws ArgumentError ModeWeights(data, -x, x, x)
        @test_throws why ModeWeights(data, -1, x, 3)
        @test_throws why ModeWeights(data, x)
        @test_throws why ModeWeights(data, -1; ℓₘᵢₙ=x)
        @test_throws why ModeWeights{Float64}(undef, x, 1, 3)
        @test_throws why ModeWeights{Float64}(undef, -1, x)
    end
    @test_throws "must all be integers of type `Int`" ModeWeights(data, Int8(-1), Int8(1), Int8(3))

    # The keyword `ℓₘᵢₙ` may also be spelled `ell_min`; where both are given, `ℓₘᵢₙ` is used
    for s in -2:2
        w = ModeWeights(zeros(Ysize(2, 4)), s; ell_min=2)
        @test SphericalFunctions.ℓₘᵢₙ(w) == 2 && SphericalFunctions.ℓₘₐₓ(w) == 4
        @test ModeWeights(zeros(Ysize(2, 4)), s; ell_min=0, ℓₘᵢₙ=2) == w
    end
    # A length that no ℓₘₐₓ can produce is refused rather than silently rounded
    @test_throws ArgumentError ModeWeights(randn(rng, ComplexF64, 7), 0; ℓₘᵢₙ=0)
    @test_throws "which is not Ysize(ℓₘᵢₙ=0, ℓₘₐₓ) for any ℓₘₐₓ" ModeWeights(zeros(7), 0; ℓₘᵢₙ=0)

    # A matrix holds one set of weights per column, and a `ModeWeights` labels one of them
    W = randn(rng, ComplexF64, Ysize(1, 3), 3)
    @test_throws ArgumentError ModeWeights(W, 1)
    @test_throws "ModeWeights(view(data, :, j), s)" ModeWeights(W, 1)
    @test_throws "15×3 matrix" ModeWeights(W, 1, 1, 3)
    @test parent(ModeWeights(view(W, :, 2), 1)) == W[:, 2]

    # Any 1-based AbstractVector serves as storage, without copying
    data = randn(rng, ComplexF64, 20)
    v = view(data, 3:11)  # Ysize(0, 2) = 9
    w = ModeWeights(v, 0)
    @test parent(w) === v
    @test SphericalFunctions.ℓₘₐₓ(w) == 2
    @test w[1] == data[3]
    w[1] = 0
    @test data[3] == 0
    r = ModeWeights(1:4, 0)
    @test eltype(r) == Int
    @test SphericalFunctions.ℓₘₐₓ(r) == 1
    @test r[1, 1] == 4
    @test collect(r) == 1:4
end


@testitem "ModeWeights indexing" setup=[RefusalChecks] begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yindex, Yrange, DegreeBlock
    import Random

    rng = Random.Xoshiro(20260911)

    for s in -2:2, ℓₘᵢₙ in unique((abs(s), 0)), ℓₘₐₓ in unique((ℓₘᵢₙ, 3, 6))
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(copy(data), s, ℓₘᵢₙ, ℓₘₐₓ)

        # Natural indexing reads the canonical position of the storage ...
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
            @test w[ℓ, m] === data[Yindex(ℓ, m, ℓₘᵢₙ)]
            @test w[ℓ, m] === parent(w)[Yindex(ℓ, m, ℓₘᵢₙ)]
        end
        # ... and linear indexing is the storage itself
        for i in 1:n
            @test w[i] === data[i]
        end
        @test firstindex(w) == 1
        @test lastindex(w) == n
        @test w[begin] === data[1]
        @test w[end] === data[n]
        @test eachindex(w) == 1:n
        @test collect(w) == data
        @test [w[i] for i in eachindex(w)] == data

        # `modes` lists the (ℓ, m) pairs in storage order, and pairs up with the values
        @test modes(w) == Yrange(ℓₘᵢₙ, ℓₘₐₓ)
        @test modes(w) == [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]
        @test length(modes(w)) == length(w)
        for (i, ((ℓ, m), v)) in enumerate(zip(modes(w), w))
            @test v === w[ℓ, m]
            @test v === w[i]
            @test Yindex(ℓ, m, ℓₘᵢₙ) == i
        end
        @test [w[ℓ, m] for (ℓ, m) in modes(w)] == data

        # setindex! by mode is visible linearly and in the storage ...
        for (i, (ℓ, m)) in enumerate(modes(w))
            v = ComplexF64(ℓ, m)  # a distinct value per mode
            w[ℓ, m] = v
            @test w[ℓ, m] === v
            @test w[i] === v
            @test parent(w)[i] === v
        end
        @test w == [ComplexF64(ℓ, m) for (ℓ, m) in modes(w)]
        # ... and setindex! linearly is visible by mode
        for (i, (ℓ, m)) in enumerate(modes(w))
            w[i] = data[i]
            @test w[ℓ, m] === data[i]
        end
        @test parent(w) == data
        # Assigning a value of another type converts it
        if n > 0
            w[ℓₘₐₓ, ℓₘₐₓ] = 3
            @test w[ℓₘₐₓ, ℓₘₐₓ] === ComplexF64(3)
            w[ℓₘₐₓ, ℓₘₐₓ] = data[end]
        end

        # `w[ℓ, :]` is a view over m ∈ -ℓ:ℓ
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ
            v = w[ℓ, :]
            @test axes(v) == (-ℓ:ℓ,)
            @test length(v) == 2ℓ + 1
            @test firstindex(v) == -ℓ
            @test lastindex(v) == ℓ
            for m in -ℓ:ℓ
                @test v[m] === w[ℓ, m]
            end
            @test collect(v) == data[Yindex(ℓ, -ℓ, ℓₘᵢₙ):Yindex(ℓ, ℓ, ℓₘᵢₙ)]
            @test_throws BoundsError v[ℓ + 1]
            @test_throws BoundsError v[-ℓ - 1]
            # Writing through the view changes `w` ...
            v[ℓ] = 7
            @test w[ℓ, ℓ] == 7
            @test parent(w)[Yindex(ℓ, ℓ, ℓₘᵢₙ)] == 7
            v .= ComplexF64(ℓ)
            @test all(w[ℓ, m] == ℓ for m in -ℓ:ℓ)
            # ... and writing to `w` is visible through the view
            w[ℓ, -ℓ] = 11
            @test v[-ℓ] == 11
            # Other ℓ blocks are untouched
            for ℓ′ in ℓₘᵢₙ:ℓₘₐₓ, m′ in -ℓ′:ℓ′
                if ℓ′ != ℓ
                    @test w[ℓ′, m′] === data[Yindex(ℓ′, m′, ℓₘᵢₙ)]
                end
            end
            # Restore
            for m in -ℓ:ℓ
                w[ℓ, m] = data[Yindex(ℓ, m, ℓₘᵢₙ)]
            end
        end
        @test parent(w) == data

        # Out-of-range modes are `BoundsError`s, reading and writing alike: ℓ outside
        # ℓₘᵢₙ:ℓₘₐₓ (with any m) ...
        for ℓ in unique((ℓₘᵢₙ - 1, ℓₘₐₓ + 1, -1, ℓₘₐₓ + 5)), m in (-1, 0, 1, ℓ, -ℓ)
            @test_throws BoundsError w[ℓ, m]
            @test_throws BoundsError w[ℓ, m] = 0
        end
        for ℓ in unique((ℓₘᵢₙ - 1, ℓₘₐₓ + 1, -1, ℓₘₐₓ + 5))
            @test_throws BoundsError w[ℓ, :]
        end
        # ... and m outside -ℓ:ℓ
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in (-ℓ - 1, ℓ + 1, -ℓ - 5, ℓ + 5)
            @test_throws BoundsError w[ℓ, m]
            @test_throws BoundsError w[ℓ, m] = 0
        end
        # Linear indices outside 1:n as well
        for i in (0, -1, n + 1, n + 5)
            @test_throws BoundsError w[i]
            @test_throws BoundsError w[i] = 0
        end
        # None of the failed writes touched the storage
        @test parent(w) == data
    end

    # The empty container has no valid mode at all
    w = ModeWeights(Float64[], 2)
    @test isempty(modes(w))
    for ℓ in -1:3, m in -3:3
        @test_throws BoundsError w[ℓ, m]
    end
    @test_throws BoundsError w[2, :]
    @test_throws BoundsError w[1]

    # The natural indices of `w[ℓ, m]` and `w[ℓ, :]` obey the rules of the constructors: an
    # integer container takes `Int`s, and an integer of another type is refused with the
    # reason, before any position is computed from it, since the index arithmetic is not
    # closed under it.  So are a `Rational`, even a whole one, and a half-odd-integer.  A
    # refused write leaves the storage alone.  Linear indices are positions in the storage,
    # and take any integer type, as the storage does.
    w = ModeWeights(collect(ComplexF64, 1:9), 0)
    data = copy(parent(w))
    narrow = "narrower than `Int`"
    for (index, why) ∈ (
        (Int32(2), narrow), (Int8(2), narrow), (UInt8(2), "unsigned"),
        (big(2), "wider than `Int`"), (Int128(2), "wider than `Int`"), (true, "`Bool`"),
        (2//1, "2//1 is a whole number; write it as the integer 2"),
        (1//2, "are integers of type `Int`, like 3; got ℓ = 1//2"),
    )
        @test refuses(() -> w[index, :], ArgumentError, why)
        @test refuses(() -> w[index, 0], ArgumentError, why)
        @test refuses(() -> w[index, 0] = 0, ArgumentError, why)
        @test refuses(() -> w[2, index], ArgumentError, "got m = ")
        @test refuses(() -> w[2, index] = 0, ArgumentError, "got m = ")
    end
    # A narrow index is refused for what it is, even where its value is out of range
    @test refuses(() -> w[Int32(3), :], ArgumentError, narrow)
    @test refuses(() -> w[big(300), :], ArgumentError, "wider than `Int`")
    @test parent(w) == data
    @test w[Int32(2)] === w[2] === parent(w)[2]
    w[UInt8(3)] = 30
    @test parent(w)[3] == 30
end


@testitem "ModeWeights array behavior" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, L², Lz, ð
    import LinearAlgebra: norm, dot, Diagonal, mul!
    import Random

    rng = Random.Xoshiro(20260912)

    same_range(a, b) = (
        spin(a) == spin(b)
        && SphericalFunctions.ℓₘᵢₙ(a) == SphericalFunctions.ℓₘᵢₙ(b)
        && SphericalFunctions.ℓₘₐₓ(a) == SphericalFunctions.ℓₘₐₓ(b)
    )

    for s in (-2, 0, 1), (ℓₘᵢₙ, ℓₘₐₓ) in ((abs(s), 5), (0, 4))
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
        ϵ = 100eps()

        # The broadcasts that are the arithmetic of mode weights — sums and differences of
        # weights with the same labels, products and quotients with numbers or with a plain
        # vector of factors, and a change of number type — return a new ModeWeights with the
        # same s and ℓ range, holding what the same broadcast on the storage gives
        for (result, expected) in (
            (2 .* w, 2 .* data),
            (w ./ 2, data ./ 2),
            (w .+ w, data .+ data),
            (w .- data, zeros(ComplexF64, n)),
            (data .* w, data .* data),
            (-w, -data),
            (ComplexF32.(w), ComplexF32.(data)),
            (2w, 2data), (w / 2, data / 2), (w + w, 2data), (w - w, zeros(ComplexF64, n)),
        )
            @test result isa ModeWeights
            @test same_range(result, w)
            @test parent(result) == expected
            @test eltype(result) == eltype(expected)
            @test parent(result) !== data
        end
        # Broadcasts that would label numbers which are not the mode weights of any function
        # are refused: the product of two sets of weights, a constant added to every weight,
        # and functions such as `conj`, `abs`, `real` and `==` applied elementwise
        for bad ∈ (() -> w .* conj.(w), () -> w .+ 1, () -> conj.(w), () -> abs.(w),
                   () -> real.(w), () -> w .== w, () -> w ./ w, () -> 1 ./ w)
            @test_throws ArgumentError bad()
        end
        if s != 0  # (for s = 0 the opposite spin weight is the same one)
            @test_throws "labels agree" w .+ ModeWeights(data, -s, ℓₘᵢₙ, ℓₘₐₓ)
            @test_throws "different labels" w + ModeWeights(data, -s, ℓₘᵢₙ, ℓₘₐₓ)
        end
        @test parent(w) === data
        @test parent(w) == data
        # In-place broadcast into a similar container, whose labels must match what is
        # computed
        w2 = similar(w)
        w2 .= 2 .* w .+ w
        @test w2 isa ModeWeights
        @test same_range(w2, w)
        @test parent(w2) == 2 .* data .+ data
        s != 0 && @test_throws "destination holds" ModeWeights(similar(data), -s, ℓₘᵢₙ, ℓₘₐₓ) .= w
        @test parent(w) == data
        w2 .= w
        @test w2 == w
        @test parent(w2) !== data

        # `similar(w)` and `similar(w, T)` keep the labels, while a size asks for plain
        # storage, even the size of `w`, so that the type does not depend on the size
        @test similar(w) isa ModeWeights{ComplexF64}
        @test same_range(similar(w), w)
        @test length(similar(w)) == n
        @test parent(similar(w)) !== data
        @test similar(w, Float32) isa ModeWeights{Float32}
        @test same_range(similar(w, Float32), w)
        @test similar(w, n) isa Vector{ComplexF64}
        @test similar(w, Float64, n) isa Vector{Float64}
        @test similar(w, (n,)) isa Vector{ComplexF64}
        @test similar(w, Float32, size(w)) isa Vector{Float32}
        @test similar(w, n + 1) isa Vector{ComplexF64}
        @test length(similar(w, n + 1)) == n + 1
        @test similar(w, Float64, n + 1) isa Vector{Float64}
        @test similar(w, (n, 2)) isa Matrix{ComplexF64}
        @test size(similar(w, (n, 2))) == (n, 2)
        # An operator *matrix* cannot check the labels of the weights, so its product with
        # them is refused; applied to the raw numbers it gives an unlabelled vector, and the
        # operator itself gives a correctly labelled one.  (`ℓₘᵢₙ` and `ℓₘₐₓ` are loop
        # variables here, not the accessor functions.)
        @test_throws "plain matrix cannot be checked" ð(s, ℓₘᵢₙ, ℓₘₐₓ) * w
        let wδ = ð(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w)
            @test wδ isa Vector
            @test !(wδ isa ModeWeights)
            @test wδ == parent(ð(w))
            @test spin(ð(w)) == s + 1
        end

        # `copy` is an independent ModeWeights with the same parameters
        c = copy(w)
        @test c isa ModeWeights{ComplexF64}
        @test same_range(c, w)
        @test c == w
        @test parent(c) !== data
        c[1] += 1
        c[ℓₘₐₓ, 0] = 0
        @test parent(w) == data
        @test w[1] == data[1]
        @test c != w
        # `map` applies an arbitrary function, so it returns plain numbers
        @test map(abs, w) isa Vector{Float64}
        @test map(abs, w) == abs.(data)

        # Reductions and equality agree with the storage
        @test sum(w) ≈ sum(data) atol=ϵ rtol=ϵ
        @test norm(w) ≈ norm(data) rtol=ϵ
        @test dot(w, w) ≈ dot(data, data) rtol=ϵ
        @test dot(w, data) ≈ dot(data, data) rtol=ϵ
        # `w'` is the adjoint of the storage, a plain row, so `w' * w` is a product with a
        # plain matrix and is refused; `dot` is the labelled inner product
        @test w' == data'
        @test_throws ArgumentError w' * w
        @test_throws "dot(w₁, w₂)" w' * w
        @test maximum(abs, w) == maximum(abs, data)
        @test w == data
        @test data == w
        @test isequal(w, data)
        @test collect(w) == data
        @test Vector(w) == data
        @test Vector(w) isa Vector{ComplexF64}
        @test w[1:min(3, n)] isa Vector{ComplexF64}

        # Products with plain matrices, of any shape, are refused on either side, and the
        # same products with the storage are ordinary products of arrays
        M = randn(rng, ComplexF64, n + 3, n)
        Msq = randn(rng, ComplexF64, n, n)
        for A in (M, Msq, Diagonal(ones(n)), L²(s, ℓₘᵢₙ, ℓₘₐₓ), Lz(s, ℓₘᵢₙ, ℓₘₐₓ))
            @test_throws "A * array_view(w)" A * w
        end
        @test_throws "A \\ array_view(w)" Msq \ w
        @test_throws "array_view(w) * A" w * transpose(M)
        for A in (M, Diagonal(ones(n)))
            y = zeros(ComplexF64, size(A, 1))
            @test_throws ArgumentError mul!(y, A, w)
            @test_throws "mul!(y, A, array_view(w))" mul!(y, A, w)
            @test_throws "mul!(y, A, array_view(w))" mul!(y, A, w, 2, 1)
            @test y == zeros(ComplexF64, size(A, 1))
        end
        @test mul!(zeros(ComplexF64, n + 3), M, array_view(w)) == M * data
        @test M * array_view(w) == M * data
        # The operator matrices applied by hand to the storage give the operators' numbers
        @test L²(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w) == parent(L²(w))
        @test Lz(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w) == parent(Lz(w))
        # Outer-product shapes cannot be a ModeWeights and must fall back to plain matrices
        @test w .+ transpose(w) isa Matrix{ComplexF64}
        @test_throws ArgumentError w * w'
        @test array_view(w) * w' isa Matrix{ComplexF64}
        @test parent(w) .+ w' isa Matrix{ComplexF64}
        @test size(parent(w) .+ w') == (n, n)

        # The two-argument `show` writes the labels, and heads the three-argument form,
        # which adds the weights
        @test repr(w) == "ModeWeights{ComplexF64} with s=$s, ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ"
        @test occursin("[$(repr(w))]", sprint(show, [w]))
        str = sprint(show, MIME("text/plain"), w)
        @test startswith(str, repr(w) * ":\n")
        @test occursin("ModeWeights{ComplexF64}", str)
        @test occursin("s=$s", str)
        @test occursin("ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ", str)
        @test count('\n', str) ≥ n
        limited = sprint(show, MIME("text/plain"), w; context=(:limit => true, :displaysize => (8, 80)))
        @test occursin("ModeWeights{ComplexF64}", limited)
        @test !isempty(sprint(show, w))
        @test !isempty(repr(w))
        @test !isempty(summary(w))
        @test !isempty(sprint(show, MIME("text/plain"), similar(w, Float32)))
    end

    # A vector combined with mode weights must hold one value per mode, so a length-1
    # ModeWeights is extended against a longer vector neither as a term nor as a factor: the
    # result of a broadcast over mode weights is always labelled, so its type does not
    # depend on the lengths, and one of another length is refused
    w1 = ModeWeights([3.0], 0)
    @test_throws ArgumentError w1 .+ [1.0, 2.0, 3.0]
    @test_throws DimensionMismatch w1 .* [1.0, 2.0, 3.0]
    @test_throws "must hold one value per mode" w1 .* [1.0, 2.0, 3.0]
    @test 2 .* w1 isa ModeWeights{Float64}
    @test (w1 .* [2.0]) isa ModeWeights{Float64}
    wc = ModeWeights(randn(rng, ComplexF64, 9), 0)
    @test @inferred(2 .* wc) isa ModeWeights
    @test @inferred(wc .+ wc) isa ModeWeights
    @test @inferred(wc ./ 2 .- wc) isa ModeWeights
    @test @inferred(similar(wc, 3)) isa Vector{ComplexF64}

    # The empty container reduces and prints like an empty vector
    w0 = ModeWeights(ComplexF64[], 1)
    @test sum(w0) == 0
    @test norm(w0) == 0
    @test isempty(2 .* w0)
    @test 2 .* w0 isa ModeWeights
    @test occursin("ℓ ∈ 1:0", sprint(show, MIME("text/plain"), w0))
end


@testitem "ModeWeights operators" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yindex
    import SphericalFunctions: L², Lz, L₊, L₋, R², Rz, R₊, R₋, ð, ð̄
    import DoubleFloats: Double64
    import Random

    rng = Random.Xoshiro(20260913)

    same_spin = (L², Lz, L₊, L₋, R², Rz)
    spin_changing = ((R₊, +1), (R₋, -1), (ð, +1), (ð̄, -1))

    for s in -2:2, ℓₘᵢₙ in unique((abs(s), 0)), ℓₘₐₓ in (3, 6)
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(copy(data), s, ℓₘᵢₙ, ℓₘₐₓ)

        # Every operator is the corresponding matrix applied to the storage, wrapped with
        # the same ℓ range; the input is untouched
        for O in same_spin
            Ow = O(w)
            @test Ow isa ModeWeights{ComplexF64}
            @test Ow == ModeWeights(O(s, ℓₘᵢₙ, ℓₘₐₓ) * data, s, ℓₘᵢₙ, ℓₘₐₓ)
            @test parent(Ow) == O(s, ℓₘᵢₙ, ℓₘₐₓ) * data
            @test spin(Ow) == s
            @test SphericalFunctions.ℓₘᵢₙ(Ow) == ℓₘᵢₙ
            @test SphericalFunctions.ℓₘₐₓ(Ow) == ℓₘₐₓ
            @test parent(Ow) !== parent(w)
            @test parent(w) == data
        end
        for (O, Δs) in spin_changing
            Ow = O(w)
            @test Ow isa ModeWeights{ComplexF64}
            @test Ow == ModeWeights(O(s, ℓₘᵢₙ, ℓₘₐₓ) * data, s + Δs, ℓₘᵢₙ, ℓₘₐₓ)
            @test parent(Ow) == O(s, ℓₘᵢₙ, ℓₘₐₓ) * data
            @test spin(Ow) == s + Δs
            @test SphericalFunctions.ℓₘᵢₙ(Ow) == ℓₘᵢₙ
            @test SphericalFunctions.ℓₘₐₓ(Ow) == ℓₘₐₓ
            @test parent(w) == data
            # Entries with ℓ below the input's |s| or the output's |s ± 1| vanish
            for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
                if ℓ < max(abs(s), abs(s + Δs))
                    @test Ow[ℓ, m] == 0
                end
            end
        end

        # Spin bookkeeping, and the relations ð = R₊, ð̄ = -R₋.  The latter two hold by
        # construction (one method of `*` in `src/mode_weights/mode_weights.jl` applies
        # every operator to a `ModeWeights`), so they pin the aliasing rather than the sign
        # convention; `test/operators.jl`'s Newman–Penrose item is what pins that.
        @test spin(ð(w)) == s + 1
        @test spin(R₊(w)) == s + 1
        @test spin(ð̄(w)) == s - 1
        @test spin(R₋(w)) == s - 1
        @test ð(w) == R₊(w)
        @test ð̄(w) == -R₋(w)
        @test spin(-R₋(w)) == spin(ð̄(w))
        @test spin(ð̄(ð(w))) == s
        @test spin(ð(ð(w))) == s + 2
        @test spin(ð̄(ð̄(w))) == s - 2

        # The documented actions on mode weights, read through natural indexing: for ℓ ≥ |s|
        #   {L² f}ₗₘ = ℓ(ℓ+1) fₗₘ,   {Lz f}ₗₘ = m fₗₘ,   {Rz f}ₗₘ = s fₗₘ,
        #   {L₊ f}ₗₘ = √((ℓ+m)(ℓ-m+1)) fₗ,ₘ₋₁,   {L₋ f}ₗₘ = √((ℓ-m)(ℓ+m+1)) fₗ,ₘ₊₁,
        #   {ð f}ₗₘ = √((ℓ-s)(ℓ+s+1)) fₗₘ,       {ð̄ f}ₗₘ = -√((ℓ+s)(ℓ-s+1)) fₗₘ,
        # and every operator maps the ℓ < |s| entries to zero
        L²w, Lzw, L₊w, L₋w, R²w, Rzw = L²(w), Lz(w), L₊(w), L₋(w), R²(w), Rz(w)
        R₊w, R₋w, ðw, ð̄w = R₊(w), R₋(w), ð(w), ð̄(w)
        ϵ = 10eps()
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
            if ℓ < abs(s)
                for Ow in (L²w, Lzw, L₊w, L₋w, R²w, Rzw, R₊w, R₋w, ðw, ð̄w)
                    @test Ow[ℓ, m] == 0
                end
            else
                @test L²w[ℓ, m] == ℓ * (ℓ + 1) * w[ℓ, m]
                @test R²w[ℓ, m] == ℓ * (ℓ + 1) * w[ℓ, m]
                @test Lzw[ℓ, m] == m * w[ℓ, m]
                @test Rzw[ℓ, m] == s * w[ℓ, m]
                if m == -ℓ
                    @test L₊w[ℓ, m] == 0
                else
                    @test L₊w[ℓ, m] ≈ √((ℓ + m) * (ℓ - m + 1)) * w[ℓ, m - 1] rtol=ϵ
                end
                if m == ℓ
                    @test L₋w[ℓ, m] == 0
                else
                    @test L₋w[ℓ, m] ≈ √((ℓ - m) * (ℓ + m + 1)) * w[ℓ, m + 1] rtol=ϵ
                end
                @test ðw[ℓ, m] ≈ √((ℓ - s) * (ℓ + s + 1)) * w[ℓ, m] rtol=ϵ
                @test R₊w[ℓ, m] ≈ √((ℓ - s) * (ℓ + s + 1)) * w[ℓ, m] rtol=ϵ
                @test ð̄w[ℓ, m] ≈ -√((ℓ + s) * (ℓ - s + 1)) * w[ℓ, m] rtol=ϵ
                @test R₋w[ℓ, m] ≈ √((ℓ + s) * (ℓ - s + 1)) * w[ℓ, m] rtol=ϵ
            end
        end
    end

    # The element type of the result follows the data: the operator matrix is built in the
    # real type of the data's element type, so nothing is promoted to Float64 or demoted
    for T in (Float32, Float64, Double64), CT in (T, Complex{T})
        for s in (-1, 2), ℓₘᵢₙ in unique((abs(s), 0))
            ℓₘₐₓ = 4
            data = randn(rng, CT, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
            w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
            for O in same_spin
                Ow = O(w)
                @test Ow isa ModeWeights{CT}
                @test eltype(Ow) == CT
                @test parent(Ow) == O(s, ℓₘᵢₙ, ℓₘₐₓ, T) * data
            end
            for (O, Δs) in spin_changing
                Ow = O(w)
                @test Ow isa ModeWeights{CT}
                @test eltype(Ow) == CT
                @test parent(Ow) == O(s, ℓₘᵢₙ, ℓₘₐₓ, T) * data
                @test spin(Ow) == s + Δs
            end
            # The coefficients are computed at the precision of T: the relative error of the
            # ladder coefficients is a few eps(T), which for `Double64` rules out coefficients
            # computed in `Float64`.  (For `Float32`, coefficients computed in `Float64` would
            # only be more accurate; the element-type checks above guard against promotion.)
            ϵ = 8eps(T)
            L₊w, ðw = L₊(w), ð(w)
            for ℓ in max(abs(s), ℓₘᵢₙ):ℓₘₐₓ, m in -ℓ+1:ℓ
                @test L₊w[ℓ, m] ≈ √(T((ℓ + m) * (ℓ - m + 1))) * w[ℓ, m - 1] rtol=ϵ
                @test ðw[ℓ, m] ≈ √(T((ℓ - s) * (ℓ + s + 1))) * w[ℓ, m] rtol=ϵ
            end
        end
    end
    # Integer data is promoted to Float64
    wi = ModeWeights(collect(1:Ysize(1, 2)), 1)
    @test L²(wi) isa ModeWeights{Float64}
    @test parent(L²(wi)) == L²(1, 1, 2) * parent(wi)
    @test ð(wi) isa ModeWeights{Float64}
    @test spin(ð(wi)) == 2
    @test parent(ð(wi)) == ð(1, 1, 2) * parent(wi)
end


@testitem "ModeWeights evaluation" setup=[Utilities] begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yindex, ð, ð̄, R₊, R₋
    import .Utilities: ℓmrange, sYlm_closed_form
    import Quaternionic: Rotor, from_spherical_coordinates, from_euler_angles
    import LinearAlgebra: norm
    import DoubleFloats: Double64
    import Random

    rng = Random.Xoshiro(20260914)

    # The independent reference for the harmonics at a rotor.  `sYlm_closed_form(s, ℓ, m, θ,
    # ϕ)` is the explicit sum of `docs/src/30-conventions/01-summary.md`, but it is only the
    # γ = 0 slice of the rotation group.  The same page gives the γ dependence as
    #
    #     ₛYₗₘ(𝐑 exp(γ𝐤/2)) = exp(-isγ) ₛYₗₘ(𝐑),
    #
    # so with 𝐑 = exp(α𝐤/2) exp(β𝐣/2) exp(γ𝐤/2) = `Quaternionic.from_euler_angles(α, β, γ)`,
    #
    #     ₛYₗₘ(𝐑) = exp(-isγ) ₛYₗₘ(θ=β, ϕ=α).
    #
    # Neither the closed form nor `from_euler_angles` runs any SphericalFunctions code, so
    # this is a truly independent oracle for `w(R)`.  Entries with ℓ < |s| are not harmonics
    # at all and are set to zero; the weights there must not contribute, which is checked
    # separately below.
    #
    # The rotors are *built* from Euler angles rather than decomposed into them: going the
    # other way, `to_euler_angles` loses about half the significant digits of β near the
    # poles (measured: the recovered β is off by 2.2e-16 at β = 0 and by 2.2e-30 at β =
    # 0.01, in Double64, the usual `acos` cancellation), which would cap this oracle at
    # √eps.
    function Yreference(s, ℓₘᵢₙ, ℓₘₐₓ, α::T, β::T, γ::T) where {T}
        [
            ℓ < abs(s) ? zero(Complex{T}) : cis(-s * γ) * sYlm_closed_form(s, ℓ, m, β, α)
            for (ℓ, m) in ℓmrange(ℓₘᵢₙ, ℓₘₐₓ)
        ]
    end

    for T in (Float64, Double64)
        # Measured below: `w(R)` differs from the closed-form sum by at most 16 eps(T) in
        # absolute value, for |w(R)| up to about 2.9, so 200 eps(T) leaves a factor ~13.
        ϵ = 200eps(T)
        euler_angles = [
            (T(1)/5,    T(3)/7,     T(9)/5),
            (T(23)/10,  T(11)/8,    T(1)/3),
            (T(41)/10,  2*T(π)/3,   T(57)/10),
            (-T(7)/4,   T(π)/7,     T(31)/10),
        ]
        rotors = [from_euler_angles(α, β, γ) for (α, β, γ) in euler_angles]
        @test all(R -> R isa Rotor{T}, rotors)

        for s in -2:2, ℓₘᵢₙ in unique((abs(s), 0)), ℓₘₐₓ in (3, 5)
            n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            data = randn(rng, Complex{T}, n)
            w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
            data2 = randn(rng, Complex{T}, n)
            w2 = ModeWeights(data2, s, ℓₘᵢₙ, ℓₘₐₓ)

            for (R, (α, β, γ)) in zip(rotors, euler_angles)
                # w(R) = Σᵢ wᵢ ₛYₗₘ(R)ᵢ over the canonical ordering, with the harmonics from
                # the closed form rather than from the package
                Y = Yreference(s, ℓₘᵢₙ, ℓₘₐₓ, α, β, γ)
                @test length(Y) == n
                f = w(R)
                @test f isa Complex{T}
                @test f ≈ sum(w[i] * Y[i] for i in 1:n) atol=ϵ rtol=ϵ
                @test f ≈ sum(w[ℓ, m] * Y[Yindex(ℓ, m, ℓₘᵢₙ)] for (ℓ, m) in modes(w)) atol=ϵ rtol=ϵ
                # (`Rotor(R)` renormalizes, so it can differ from `R` in the last bit.
                # Measured: no difference at all here.)
                @test w(Rotor(R)) ≈ f atol=ϵ rtol=ϵ

                # Linearity in the weights (measured: at most 0.6 eps(T)*scale, so 100 is a
                # factor ~170 above; this is metamorphic, but the absolute check above pins
                # the same quantity against the closed form)
                a, b = randn(rng, Complex{T}, 2)
                scale = abs(a) * norm(w) + abs(b) * norm(w2)
                @test (a .* w .+ b .* w2)(R) ≈ a * w(R) + b * w2(R) atol=100eps(T)*scale

                # With ℓₘᵢₙ = 0, the weights with ℓ < |s| belong to no harmonic and are not
                # read, so whatever they hold, a non-finite value included, the value is
                # exactly the same, and the ℓₘᵢₙ = |s| container gives the same value
                if ℓₘᵢₙ == 0 && s != 0
                    for junk in (randn(rng, Complex{T}), T(NaN), T(Inf))
                        wj = copy(w)
                        for ℓ in 0:abs(s)-1, m in -ℓ:ℓ
                            wj[ℓ, m] = junk
                        end
                        @test wj(R) == f
                    end
                    wt = ModeWeights(data[Yindex(abs(s), -abs(s), 0):end], s)
                    @test SphericalFunctions.ℓₘᵢₙ(wt) == abs(s)
                    @test wt(R) ≈ f atol=ϵ rtol=ϵ
                end

                # Evaluating ð(w) and ð̄(w) equals applying the explicit coefficients to the
                # weights pointwise, with the harmonics of the raised or lowered spin:
                #   (ð f)(R) = Σ √((ℓ-s)(ℓ+s+1)) fₗₘ ₛ₊₁Yₗₘ(R),
                #   (ð̄ f)(R) = -Σ √((ℓ+s)(ℓ-s+1)) fₗₘ ₛ₋₁Yₗₘ(R)
                Y₊ = Yreference(s + 1, ℓₘᵢₙ, ℓₘₐₓ, α, β, γ)
                Y₋ = Yreference(s - 1, ℓₘᵢₙ, ℓₘₐₓ, α, β, γ)
                ðf = sum(
                    √(T((ℓ - s) * (ℓ + s + 1))) * w[ℓ, m] * Y₊[Yindex(ℓ, m, ℓₘᵢₙ)]
                    for (ℓ, m) in modes(w) if ℓ ≥ abs(s)
                )
                ð̄f = -sum(
                    √(T((ℓ + s) * (ℓ - s + 1))) * w[ℓ, m] * Y₋[Yindex(ℓ, m, ℓₘᵢₙ)]
                    for (ℓ, m) in modes(w) if ℓ ≥ abs(s)
                )
                # Measured: at most 2.2 eps(T)*ℓₘₐₓ*norm(w), so 100 leaves a factor ~45
                tol = 100eps(T) * ℓₘₐₓ * norm(w)
                @test ð(w)(R) ≈ ðf atol=tol
                @test R₊(w)(R) ≈ ðf atol=tol
                @test ð̄(w)(R) ≈ ð̄f atol=tol
                @test R₋(w)(R) ≈ -ð̄f atol=tol
                @test ð(w)(R) isa Complex{T}
            end
        end

        # A single-mode container evaluates to the pointwise closed-form harmonic
        # sYlm_closed_form(s, ℓ, m, θ, ϕ) at R = from_spherical_coordinates(θ, ϕ) — which is
        # the γ = 0 slice, so the closed form applies with no rotation law at all —
        # including at the poles.  Measured: at most 4.1 eps(T) in absolute value, for
        # |ₛYₗₘ| up to 0.94.
        θs = T[0, T(π)/5, 2T(π)/3, T(π)]
        ϕs = T[0, 19//10, 51//10]
        for s in -2:2, ℓₘᵢₙ in unique((abs(s), 0))
            ℓₘₐₓ = 5
            n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            for θ in θs, ϕ in ϕs
                R = Rotor(from_spherical_coordinates(θ, ϕ))
                @test R isa Rotor{T}
                # `from_spherical_coordinates(θ, ϕ)` is the γ = 0 rotor with (α, β) = (ϕ, θ)
                Y = Yreference(s, ℓₘᵢₙ, ℓₘₐₓ, ϕ, θ, zero(T))
                for ℓ in abs(s):ℓₘₐₓ, m in -ℓ:ℓ
                    e = ModeWeights(zeros(Complex{T}, n), s, ℓₘᵢₙ, ℓₘₐₓ)
                    e[ℓ, m] = 1
                    @test e(R) ≈ sYlm_closed_form(s, ℓ, m, θ, ϕ) atol=ϵ rtol=ϵ
                    @test e(R) ≈ Y[Yindex(ℓ, m, ℓₘᵢₙ)] atol=ϵ rtol=ϵ
                end
            end
        end
    end

    # The result type follows the data and the rotor
    R = randn(rng, Rotor{Float32})
    @test ModeWeights(randn(rng, ComplexF32, Ysize(1, 3)), 1)(R) isa ComplexF32
    @test ModeWeights(randn(rng, Float32, Ysize(1, 3)), 1)(R) isa ComplexF32
end


### Half-integer indices
#
# The items below repeat the checks above for a `ModeWeights` whose spin weight and ℓ range
# are half-odd-integers.  The indices are written as `Rational`s with denominator 2, which
# is the public spelling, and the parameters that come back are compared against
# `HalfOddInteger`s, which is what the package stores.  Ranges such as `ℓₘᵢₙ-1:7//2` step
# through the half-odd-integers in the `Rational` spelling.

@testitem "ModeWeights half-integer construction" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yrange, HalfOddInteger
    import DoubleFloats: Double64
    import Random
    using Quaternionic: Rotor

    rng = Random.Xoshiro(20260916)
    h(x) = HalfOddInteger(x)

    # The four-argument form, with either spelling: the parameters come back as
    # `HalfOddInteger`s, the storage is the very same array, and ℓₘₐₓ = ℓₘᵢₙ - 1 is the
    # empty container
    for s in (-3//2, -1//2, 1//2, 3//2, 5//2), ℓₘᵢₙ in (1//2, 3//2, 5//2), ℓₘₐₓ in ℓₘᵢₙ-1:9//2
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        for w in (ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ), ModeWeights(data, h(s), h(ℓₘᵢₙ), h(ℓₘₐₓ)))
            @test w isa ModeWeights{ComplexF64, HalfOddInteger, Vector{ComplexF64}}
            @test w isa SphericalFunctions.AbstractModeContainer{ComplexF64, HalfOddInteger}
            @test !(w isa AbstractVector)
            @test parent(w) === data
            @test spin(w) === h(s)
            @test SphericalFunctions.ℓₘᵢₙ(w) === h(ℓₘᵢₙ)
            @test SphericalFunctions.ℓₘₐₓ(w) === h(ℓₘₐₓ)
            @test spin(w) == s
            @test length(w) == n
            @test size(w) == (n,)
            @test isempty(w) == (ℓₘₐₓ == ℓₘᵢₙ - 1)
            @test modes(w) == Yrange(ℓₘᵢₙ, ℓₘₐₓ)
            @test length(modes(w)) == n
        end
    end

    # With only `ℓₘᵢₙ` given, ℓₘₐₓ is deduced from the length, through the relation
    # (2ℓₘₐₓ+2)² = 4⋅length + (2ℓₘᵢₙ)² ...
    for ℓₘᵢₙ in (1//2, 3//2, 5//2), ℓₘₐₓ in ℓₘᵢₙ-1:11//2
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, n)
        for s in (-3//2, -1//2, 1//2, 3//2)
            w = ModeWeights(data, s; ℓₘᵢₙ)
            @test parent(w) === data
            @test spin(w) === h(s)
            @test SphericalFunctions.ℓₘᵢₙ(w) === h(ℓₘᵢₙ)
            @test SphericalFunctions.ℓₘₐₓ(w) === h(ℓₘₐₓ)
            @test length(w) == n
            w′ = ModeWeights(data, h(s); ℓₘᵢₙ=h(ℓₘᵢₙ))
            @test SphericalFunctions.ℓₘₐₓ(w′) === h(ℓₘₐₓ)
            @test parent(w′) === data
        end
        # ... and `ℓₘᵢₙ` itself defaults to `abs(s)`
        for s in (-ℓₘᵢₙ, ℓₘᵢₙ)
            w = ModeWeights(data, s)
            @test spin(w) === h(s)
            @test SphericalFunctions.ℓₘᵢₙ(w) === h(ℓₘᵢₙ)
            @test SphericalFunctions.ℓₘₐₓ(w) === h(ℓₘₐₓ)
        end
    end
    # Every length that is not `Ysize(ℓₘᵢₙ, ℓₘₐₓ)` for some half-odd ℓₘₐₓ ≥ ℓₘᵢₙ - 1 is
    # rejected.  That includes the lengths that *are* `Ysize` for an integer ℓₘₐₓ — the
    # perfect squares, for ℓₘᵢₙ = 1/2 — whose root has the wrong parity.
    for ℓₘᵢₙ in (1//2, 3//2, 5//2)
        valid = Set(Ysize(ℓₘᵢₙ, ℓₘₐₓ) for ℓₘₐₓ in ℓₘᵢₙ-1:25//2)
        for n in 0:Ysize(ℓₘᵢₙ, 25//2)
            if n ∈ valid
                @test ModeWeights(zeros(n), 1//2; ℓₘᵢₙ) isa ModeWeights{Float64, HalfOddInteger}
            else
                @test_throws ArgumentError ModeWeights(zeros(n), 1//2; ℓₘᵢₙ)
            end
        end
    end
    @test_throws "for any half-odd-integer ℓₘₐₓ" ModeWeights(zeros(4), 1//2)  # Ysize(0, 1)
    @test_throws "for any half-odd-integer ℓₘₐₓ" ModeWeights(zeros(1), 1//2)
    @test_throws "for any half-odd-integer ℓₘₐₓ" ModeWeights(zeros(3), 3//2)
    @test_throws "for any half-odd-integer ℓₘₐₓ" ModeWeights(zeros(9), 1//2)  # Ysize(0, 2)
    @test_throws "for any half-odd-integer ℓₘₐₓ" ModeWeights(zeros(16), 1//2)  # Ysize(0, 3)
    # The converse: a half-integer length is not an integer one
    @test_throws "for any ℓₘₐₓ" ModeWeights(zeros(2), 0)  # Ysize(1//2, 1//2)
    @test_throws "for any ℓₘₐₓ" ModeWeights(zeros(6), 1)  # Ysize(1//2, 3//2)

    # Uninitialized storage of a given element type, with and without `ℓₘᵢₙ`
    for T in (Float32, Float64, ComplexF64, Complex{Double64})
        for s in (-3//2, 1//2, 3//2), ℓₘᵢₙ in (1//2, 3//2), ℓₘₐₓ in ℓₘᵢₙ-1:7//2
            w = ModeWeights{T}(undef, s, ℓₘᵢₙ, ℓₘₐₓ)
            @test w isa ModeWeights{T, HalfOddInteger, Vector{T}}
            @test eltype(w) == T
            @test length(w) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            @test spin(w) === h(s)
            @test SphericalFunctions.ℓₘᵢₙ(w) === h(ℓₘᵢₙ)
            @test SphericalFunctions.ℓₘₐₓ(w) === h(ℓₘₐₓ)
            @test ModeWeights{T}(undef, h(s), h(ℓₘᵢₙ), h(ℓₘₐₓ)) isa ModeWeights{T, HalfOddInteger}
            if !isempty(w)
                w[1] = 1
                @test parent(w)[1] == 1
            end
        end
        for s in (-3//2, 1//2, 3//2), ℓₘₐₓ in abs(s)-1:7//2
            w = ModeWeights{T}(undef, s, ℓₘₐₓ)
            @test w isa ModeWeights{T, HalfOddInteger}
            @test length(w) == Ysize(abs(s), ℓₘₐₓ)
            @test spin(w) === h(s)
            @test SphericalFunctions.ℓₘᵢₙ(w) === h(abs(s))
            @test SphericalFunctions.ℓₘₐₓ(w) === h(ℓₘₐₓ)
        end
    end

    # Invalid parameters are `ArgumentError`s in every form: a negative ℓₘᵢₙ ...
    @test_throws "ℓₘᵢₙ" ModeWeights(zeros(2), 1//2, -1//2, 1//2)
    @test_throws "ℓₘᵢₙ" ModeWeights(zeros(2), 1//2; ℓₘᵢₙ=-1//2)
    @test_throws "ℓₘᵢₙ" ModeWeights{Float64}(undef, 1//2, -1//2, 1//2)
    # ... ℓₘₐₓ < ℓₘᵢₙ - 1 ...
    @test_throws "ℓₘₐₓ" ModeWeights(Float64[], 1//2, 3//2, -1//2)
    @test_throws "ℓₘₐₓ" ModeWeights{Float64}(undef, 1//2, 3//2, -1//2)
    @test_throws ArgumentError ModeWeights{Float64}(undef, 5//2, 1//2)
    # ... a length that does not match the given range ...
    for ℓₘᵢₙ in (1//2, 3//2), ℓₘₐₓ in ℓₘᵢₙ-1:5//2
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        for bad in (n - 1, n + 1, 2n + 1)
            bad ≥ 0 || continue
            @test_throws "length" ModeWeights(zeros(bad), 1//2, ℓₘᵢₙ, ℓₘₐₓ)
        end
    end
    # ... a `Rational` that is not a half-odd-integer or whose integer type is not `Int` ...
    @test_throws ArgumentError ModeWeights(zeros(2), 1//3)
    @test_throws "1//3 is neither an integer nor a half-odd-integer" ModeWeights(zeros(2), 1//3)
    @test_throws "1//1 is a whole number; write it as the integer 1" ModeWeights(zeros(4), 1//1)
    @test_throws "1//4 is neither an integer nor a half-odd-integer" ModeWeights(zeros(2), 1//2, 1//4, 1//2)
    @test_throws "3//1 is a whole number; write it as the integer 3" ModeWeights{Float64}(undef, 1//2, 3//1)
    @test_throws "`Rational{Int8}` is not `Rational{Int}`" ModeWeights(zeros(2), Int8(1)//Int8(2))
    @test_throws "`Rational{BigInt}` is not `Rational{Int}`" ModeWeights{Float64}(undef, big(1)//2, 7//2)
    # ... and, in every form, a mixture of integer and half-integer indices, which is
    # refused with a message naming both kinds and saying which argument is which
    mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
    for args in ((1//2, 0, 7//2), (1, 1//2, 7//2), (1//2, 1//2, 3), (h(1//2), 0, h(7//2)))
        @test_throws ArgumentError ModeWeights(zeros(20), args...)
        @test_throws mixed ModeWeights(zeros(20), args...)
        @test_throws mixed ModeWeights{Float64}(undef, args...)
    end
    @test_throws "and so mixes integers (ℓₘᵢₙ) with half-odd-integers (s, ℓₘₐₓ)" ModeWeights(zeros(20), 1//2, 0, 7//2)
    @test_throws mixed ModeWeights{Float64}(undef, 1//2, 3)
    @test_throws mixed ModeWeights{Float64}(undef, 1, 7//2)
    # A keyword `ℓₘᵢₙ` must be of the kind of the spin weight, in either spelling
    keyword = "The keyword argument `ℓₘᵢₙ` of `ModeWeights` must be an index of the same kind"
    @test_throws keyword ModeWeights(zeros(20), 1//2; ℓₘᵢₙ=0)
    @test_throws keyword ModeWeights(zeros(20), 1; ℓₘᵢₙ=1//2)
    @test_throws "`ell_min` of `ModeWeights`" ModeWeights(zeros(20), 1; ell_min=1//2)
    @test ModeWeights(zeros(Ysize(3//2, 7//2)), 1//2; ell_min=3//2) ==
        ModeWeights(zeros(Ysize(3//2, 7//2)), h(1//2); ℓₘᵢₙ=h(3//2))
end


@testitem "ModeWeights half-integer indexing" setup=[RefusalChecks] begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yindex, Yrange
    import SphericalFunctions: HalfOddInteger, DegreeBlock
    import Random

    rng = Random.Xoshiro(20260917)
    h(x) = HalfOddInteger(x)

    for s in (-1//2, 1//2, 3//2), ℓₘᵢₙ in unique((abs(s), 1//2)), ℓₘₐₓ in unique((ℓₘᵢₙ, 7//2, 11//2))
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(copy(data), s, ℓₘᵢₙ, ℓₘₐₓ)

        # `modes` lists the (ℓ, m) pairs in storage order, as `HalfOddInteger`s
        @test modes(w) == Yrange(ℓₘᵢₙ, ℓₘₐₓ)
        @test eltype(modes(w)) === Tuple{HalfOddInteger, HalfOddInteger}
        @test modes(w) == [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]
        @test length(modes(w)) == n
        for (i, ((ℓ, m), v)) in enumerate(zip(modes(w), w))
            @test ℓ isa HalfOddInteger && m isa HalfOddInteger
            # Natural indexing reads the canonical position of the storage, in either
            # spelling
            @test v === w[ℓ, m]
            @test v === w[Rational(ℓ), Rational(m)]
            @test v === w[i]
            @test v === data[Yindex(ℓ, m, ℓₘᵢₙ)]
            @test Yindex(ℓ, m, h(ℓₘᵢₙ)) == i
        end
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
            @test w[ℓ, m] === data[Yindex(ℓ, m, ℓₘᵢₙ)]
        end
        @test [w[ℓ, m] for (ℓ, m) in modes(w)] == data

        # setindex! by mode is visible linearly and in the storage, in either spelling ...
        for (i, (ℓ, m)) in enumerate(modes(w))
            v = ComplexF64(Rational(ℓ), Rational(m))  # a distinct value per mode
            w[ℓ, m] = v
            @test w[ℓ, m] === v
            @test w[i] === v
            w[Rational(ℓ), Rational(m)] = 2v
            @test w[ℓ, m] === 2v
            @test parent(w)[i] === 2v
        end
        # ... and setindex! linearly is visible by mode
        for (i, (ℓ, m)) in enumerate(modes(w))
            w[i] = data[i]
            @test w[ℓ, m] === data[i]
        end
        @test parent(w) == data

        # `w[ℓ, :]` is a `DegreeBlock` over m ∈ -ℓ:ℓ, indexed by half-odd-integers in either
        # spelling, and writing through it writes into `w`
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ
            v = w[ℓ, :]
            @test v isa DegreeBlock
            @test w[h(ℓ), :] isa DegreeBlock
            @test axes(v) == (-ℓ:ℓ,)
            @test length(v) == 2ℓ + 1
            @test firstindex(v) == -ℓ
            @test lastindex(v) == ℓ
            @test SphericalFunctions.ℓ(v) === h(ℓ)
            @test SphericalFunctions.mₘᵢₙ(v) === h(-ℓ)
            @test SphericalFunctions.mₘₐₓ(v) === h(ℓ)
            @test sprint(show, axes(v)) == "($(h(-ℓ)):$(h(ℓ)),)"  # the axis can be printed
            @test length(axes(v, 1)) == 2ℓ + 1
            for m in -ℓ:ℓ
                @test v[m] === w[ℓ, m]
                @test v[h(m)] === w[ℓ, m]
            end
            @test collect(v) == data[Yindex(ℓ, -ℓ, ℓₘᵢₙ):Yindex(ℓ, ℓ, ℓₘᵢₙ)]
            @test_throws BoundsError v[ℓ + 1]
            @test_throws BoundsError v[-ℓ - 1]
            v[ℓ] = 7
            @test w[ℓ, ℓ] == 7
            @test parent(w)[Yindex(ℓ, ℓ, ℓₘᵢₙ)] == 7
            w[ℓ, -ℓ] = 11
            @test v[-ℓ] == 11
            for ℓ′ in ℓₘᵢₙ:ℓₘₐₓ, m′ in -ℓ′:ℓ′
                if ℓ′ != ℓ
                    @test w[ℓ′, m′] === data[Yindex(ℓ′, m′, ℓₘᵢₙ)]
                end
            end
            for m in -ℓ:ℓ
                w[ℓ, m] = data[Yindex(ℓ, m, ℓₘᵢₙ)]
            end
        end
        @test parent(w) == data

        # Out-of-range modes are `BoundsError`s, reading and writing alike
        for ℓ in unique((ℓₘᵢₙ - 1, ℓₘₐₓ + 1, -1//2, ℓₘₐₓ + 5)), m in (-1//2, 1//2, ℓ, -ℓ)
            @test_throws BoundsError w[ℓ, m]
            @test_throws BoundsError w[ℓ, m] = 0
        end
        for ℓ in unique((ℓₘᵢₙ - 1, ℓₘₐₓ + 1, -1//2, ℓₘₐₓ + 5))
            @test_throws BoundsError w[ℓ, :]
        end
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in (-ℓ - 1, ℓ + 1, -ℓ - 5, ℓ + 5)
            @test_throws BoundsError w[ℓ, m]
            @test_throws BoundsError w[ℓ, m] = 0
        end
        @test parent(w) == data

        # Integer indices into a half-integer `w` are refused with an explanation that names
        # the index, even where only one of the two is an integer, as is a `Rational` that
        # is not a half-odd-integer or whose integer type is not `Int`; none of them touches
        # the storage
        kind = (
            "The indices of this `ModeWeights` are half-odd-integers, each a `HalfOddInteger` "
            * "or a `Rational{Int}` with denominator 2, like 7//2; got "
        )
        for (ℓ, m) in ((1, 0), (2, 1), (0, 0), (Int8(1), Int8(1)))
            @test refuses(() -> w[ℓ, m], ArgumentError, kind * "ℓ = $ℓ::$(typeof(ℓ))")
            @test refuses(() -> w[ℓ, m] = 0, ArgumentError, kind * "ℓ = ")
            @test refuses(() -> w[ℓ, :], ArgumentError, kind * "ℓ = ")
        end
        @test refuses(() -> w[3//2, 1], ArgumentError, kind * "m = 1::Int64")
        @test refuses(() -> w[1, 1//2], ArgumentError, kind * "ℓ = 1::Int64")
        @test refuses(() -> w[h(3//2), 1] = 0, ArgumentError, kind * "m = 1::Int64")
        @test refuses(() -> w[1//3, 1//3], ArgumentError, "1//3 is neither an integer nor a half")
        @test refuses(() -> w[1//1, :], ArgumentError, "1//1 is a whole number")
        @test refuses(
            () -> w[Int8(3)//Int8(2), 1//2], ArgumentError,
            "`Rational{Int8}` is not `Rational{Int}`"
        )
        @test parent(w) == data
    end

    # Conversely, an integer `w` refuses half-integer indices, and indices of another
    # integer type, and is otherwise as it was
    wi = ModeWeights(randn(rng, 9), 0)
    kind = "The indices of this `ModeWeights` are integers of type `Int`, like 3; got "
    @test refuses(() -> wi[3//2, 1//2], ArgumentError, kind * "ℓ = 3//2::Rational{Int64}")
    @test refuses(() -> wi[3//2, 1//2] = 0, ArgumentError, kind * "ℓ = 3//2")
    @test refuses(() -> wi[3//2, :], ArgumentError, kind * "ℓ = 3//2")
    @test refuses(() -> wi[1, 1//2], ArgumentError, kind * "m = 1//2")
    @test refuses(() -> wi[Int8(2), Int8(-1)], ArgumentError, "narrower than `Int`")
    @test wi[1, 0] === parent(wi)[Yindex(1, 0, 0)]
    # The integer path returns the same container as the half-integer one
    @test wi[1, :] isa DegreeBlock
    @test axes(wi[1, :]) == (-1:1,)

    # The empty container has no valid mode at all
    w = ModeWeights(Float64[], 3//2)
    @test isempty(modes(w))
    for ℓ in -1//2:5//2, m in -3//2:3//2
        @test_throws BoundsError w[ℓ, m]
    end
    @test_throws BoundsError w[3//2, :]
end


@testitem "ModeWeights half-integer array behavior" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, L², Lz, ð, HalfOddInteger
    import LinearAlgebra: norm, dot, Diagonal
    import Random

    rng = Random.Xoshiro(20260918)
    h(x) = HalfOddInteger(x)

    same_range(a, b) = (
        spin(a) === spin(b)
        && SphericalFunctions.ℓₘᵢₙ(a) === SphericalFunctions.ℓₘᵢₙ(b)
        && SphericalFunctions.ℓₘₐₓ(a) === SphericalFunctions.ℓₘₐₓ(b)
    )

    for s in (-1//2, 3//2), (ℓₘᵢₙ, ℓₘₐₓ) in ((abs(s), 9//2), (1//2, 7//2))
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
        ϵ = 100eps()

        # The arithmetic of mode weights keeps the wrapper, with the same `HalfOddInteger`
        # parameters
        for (result, expected) in (
            (2 .* w, 2 .* data),
            (w .+ w, data .+ data),
            (w .- data, zeros(ComplexF64, n)),
            (data .* w, data .* data),
            (-w, -data),
            (2w, 2data), (w + w, 2data),
        )
            @test result isa ModeWeights{eltype(expected), HalfOddInteger}
            @test same_range(result, w)
            @test parent(result) == expected
            @test parent(result) !== data
        end
        # Broadcasts that would label numbers which are not the mode weights of any function
        # are refused: the product of two sets of weights, a constant added to every weight,
        # and functions such as `conj`, `abs`, `real` and `==` applied elementwise
        for bad ∈ (() -> w .* conj.(w), () -> w .+ 1, () -> conj.(w), () -> abs.(w),
                   () -> real.(w), () -> w .== w, () -> w ./ w, () -> 1 ./ w)
            @test_throws ArgumentError bad()
        end
        @test_throws "labels agree" w .+ ModeWeights(data, -s, ℓₘᵢₙ, ℓₘₐₓ)
        @test parent(w) === data
        w2 = similar(w)
        w2 .= 2 .* w .+ w
        @test w2 isa ModeWeights{ComplexF64, HalfOddInteger}
        @test same_range(w2, w)
        @test parent(w2) == 2 .* data .+ data

        # `similar` and `copy` keep the wrapper and its parameters, while a size asks for
        # plain storage
        @test similar(w) isa ModeWeights{ComplexF64, HalfOddInteger}
        @test same_range(similar(w), w)
        @test similar(w, Float32) isa ModeWeights{Float32, HalfOddInteger}
        @test same_range(similar(w, Float32), w)
        @test similar(w, n) isa Vector{ComplexF64}
        @test similar(w, n + 1) isa Vector{ComplexF64}
        @test similar(w, (n, 2)) isa Matrix{ComplexF64}
        c = copy(w)
        @test c isa ModeWeights{ComplexF64, HalfOddInteger}
        @test same_range(c, w)
        @test c == w
        @test parent(c) !== data
        c[ℓₘₐₓ, 1//2] = 0
        @test parent(w) == data
        @test c != w
        @test map(abs, w) isa Vector{Float64}
        @test map(abs, w) == abs.(data)

        # Reductions, equality and conversion agree with the storage
        @test sum(w) ≈ sum(data) atol=ϵ rtol=ϵ
        @test norm(w) ≈ norm(data) rtol=ϵ
        @test dot(w, w) ≈ dot(data, data) rtol=ϵ
        @test w == data
        @test collect(w) == data
        @test Vector(w) isa Vector{ComplexF64}
        @test w[1:3] isa Vector{ComplexF64}

        # An operator *matrix* cannot check the labels, so it multiplies only the storage,
        # giving an unlabelled vector; the operator *function* gives a labelled
        # `ModeWeights`, with the same numbers
        @test L²(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w) isa Vector{ComplexF64}
        @test L²(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w) == parent(L²(w))
        @test Lz(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w) == parent(Lz(w))
        @test ð(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w) == parent(ð(w))
        @test spin(ð(w)) === h(s + 1)
        for A in (L²(s, ℓₘᵢₙ, ℓₘₐₓ), ð(s, ℓₘᵢₙ, ℓₘₐₓ), Diagonal(ones(n)))
            @test_throws "plain matrix cannot be checked" A * w
        end
        @test w .+ transpose(w) isa Matrix{ComplexF64}

        # `show` names the parameters in the `Rational` spelling
        str = sprint(show, MIME("text/plain"), w)
        @test occursin("ModeWeights{ComplexF64}", str)
        @test occursin("s=$s", str)
        @test occursin("ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ", str)
        @test count('\n', str) ≥ n
        @test !isempty(sprint(show, w))
        @test !isempty(summary(w))
    end

    # The empty container reduces and prints like an empty vector
    w0 = ModeWeights(ComplexF64[], 1//2)
    @test sum(w0) == 0
    @test isempty(2 .* w0)
    @test 2 .* w0 isa ModeWeights{ComplexF64, HalfOddInteger}
    @test occursin("ℓ ∈ 1//2:-1//2", sprint(show, MIME("text/plain"), w0))
end


@testitem "ModeWeights half-integer evaluation" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, sYlm, HalfOddInteger
    import Quaternionic: Rotor, Quaternion, components
    import LinearAlgebra: norm
    import DoubleFloats: Double64
    import Random

    rng = Random.Xoshiro(20260919)

    # The reference is the explicit sum Σ f_{ℓm} ₛYₗₘ(R), with the harmonics taken from the
    # documented definition ₛYₗₘ(R) = i^{2s} √((2ℓ+1)/4π) conj(𝔇ˡ_{m,-s}(R)) and the
    # factorization 𝔇ˡ_{m,-s} = e^{-imα} dˡ_{m,-s}(β) e^{isγ}, with Varshalovich's d from the
    # oracle, which shares no code with the package.  (`w(R)` itself is computed from `sYlm`,
    # so a reference built from `sYlm` would test only the summation.)  Everything is in
    # `BigFloat`, so that `Double64` is checked to its own precision.  The Euler angles
    # come from the half-angle phases of the rotor's components, so that the sign of a
    # half-integer 𝔇 is that of R and not of -R.  The weights are read through natural
    # indexing, so that the pairing of weights with harmonics is the one `modes` gives.
    function sYlm_oracle(R, ℓ, m, s)
        W, X, Y, Z = BigFloat.(components(Quaternion(R)))
        ϕₛ, ϕₐ = angle(Complex(W, Z)), angle(Complex(Y, X))
        α, β, γ = ϕₛ - ϕₐ, 2atan(abs(Complex(Y, X)), abs(Complex(W, Z))), ϕₛ + ϕₐ
        ℓ, m, s = Rational(ℓ), Rational(m), Rational(s)
        𝔇 = cis(-m * α) * HalfIntegerOracle.d_oracle(ℓ, m, -s, β) * cis(s * γ)
        cispi(BigFloat(s)) * √((2ℓ + 1) / (4BigFloat(π))) * conj(𝔇)
    end
    # Measured: at most 1.1 eps(T) times the norm of the weights, for both types; 10 is
    # asserted.  (Evaluated at -R instead, the reference misses by more than 10¹⁴ of those
    # units, so a sign error would be caught.)
    for T in (Float64, Double64)
        ϵ = 10eps(T)
        # (ℓₘₐₓ must be at least |s|: a container with ℓₘₐₓ < |s| holds no harmonics to sum.)
        for s in (-3//2, -1//2, 1//2, 3//2), ℓₘᵢₙ in unique((abs(s), 1//2)), ℓₘₐₓ in (max(abs(s), ℓₘᵢₙ), 7//2, 11//2)
            n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            data = randn(rng, Complex{T}, n)
            w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
            for _ in 1:3
                R = randn(rng, Rotor{T})
                expected = sum(
                    w[ℓ, m] * sYlm_oracle(R, ℓ, m, s) for (ℓ, m) in modes(w);
                    init=zero(Complex{BigFloat})
                )
                f = w(R)
                @test f isa Complex{T}
                @test f ≈ expected atol=ϵ*norm(data) rtol=ϵ
                # Linearity in the weights, with a second set of weights and complex
                # coefficients, so that the check is not exact in floating point
                w₂ = ModeWeights(randn(rng, Complex{T}, n), s, ℓₘᵢₙ, ℓₘₐₓ)
                a, b = randn(rng, Complex{T}, 2)
                @test (a .* w .+ b .* w₂)(R) ≈ a * f + b * w₂(R) atol=ϵ*norm(data) rtol=ϵ
            end
        end
    end

    # Weights below |s| belong to no harmonic, and are not read, so that whatever they hold,
    # a non-finite value included, the value of the function is exactly what it was
    for s in (3//2, -3//2), junk in (1e6, NaN, Inf)
        w = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), s, 1//2, 7//2)
        R = randn(rng, Rotor{Float64})
        f, F = w(R), w([R, R])
        @test isfinite(f) && F ≈ [f, f]
        for m in -1//2:1//2
            w[1//2, m] = junk
        end
        @test w(R) == f
        @test w([R, R]) == F
    end
    # A container with ℓₘₐₓ < |s| — the empty one included — holds no ℓ ≥ |s|, and so
    # describes the zero function, which evaluates to zero at every rotor, for either kind
    # of index; an empty vector of rotors gives an empty vector of values
    R = randn(rng, Rotor{Float64})
    for w in (ModeWeights(ComplexF64[], 1//2), ModeWeights(ComplexF64[], 1),
              ModeWeights(zeros(ComplexF64, 2), 3//2, 1//2, 1//2),
              ModeWeights(ones(ComplexF64, 2), 3//2, 1//2, 1//2))
        @test w(R) === zero(ComplexF64)
        @test w([R, R]) == zeros(ComplexF64, 2)
        @test w([R, R]) isa Vector{ComplexF64}
        @test w(Rotor{Float32}(R)) === zero(ComplexF64)
    end
    for w in (
        ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 5//2)), 1//2), ModeWeights(ComplexF32[], 0)
    )
        @test w(Rotor{Float64}[]) == ComplexF64[]
        @test w(Rotor{Float64}[]) isa Vector{ComplexF64}
    end
end


@testitem "ModeWeights half-integer operators" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, HalfOddInteger
    import SphericalFunctions: L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄
    import DoubleFloats: Double64
    import Random

    rng = Random.Xoshiro(20260920)
    h(x) = HalfOddInteger(x)

    same_spin = (L², Lz, L₊, L₋, Lx, Ly, R², Rz)
    spin_changing = ((R₊, +1), (R₋, -1), (ð, +1), (ð̄, -1))

    for T in (Float32, Float64, Double64)
        for s in (-3//2, -1//2, 1//2, 3//2), ℓₘᵢₙ in unique((abs(s), 1//2)), ℓₘₐₓ in (7//2, 11//2)
            n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            data = randn(rng, Complex{T}, n)
            w = ModeWeights(copy(data), s, ℓₘᵢₙ, ℓₘₐₓ)
            sh = h(s)

            # Every operator is the corresponding matrix, built in `T`, applied to the
            # storage, and wrapped with the same `HalfOddInteger` parameters; the input is
            # untouched
            for O in same_spin
                Ow = O(w)
                @test Ow isa ModeWeights{Complex{T}, HalfOddInteger}
                @test parent(Ow) == O(s, ℓₘᵢₙ, ℓₘₐₓ, T) * data
                @test spin(Ow) === sh
                @test SphericalFunctions.ℓₘᵢₙ(Ow) === h(ℓₘᵢₙ)
                @test SphericalFunctions.ℓₘₐₓ(Ow) === h(ℓₘₐₓ)
                @test parent(w) == data
            end
            for (O, Δs) in spin_changing
                Ow = O(w)
                @test Ow isa ModeWeights{Complex{T}, HalfOddInteger}
                @test parent(Ow) == O(s, ℓₘᵢₙ, ℓₘₐₓ, T) * data
                @test spin(Ow) === h(s + Δs)
                @test SphericalFunctions.ℓₘᵢₙ(Ow) === h(ℓₘᵢₙ)
                @test SphericalFunctions.ℓₘₐₓ(Ow) === h(ℓₘₐₓ)
                for (ℓ, m) in modes(w)
                    if ℓ < max(abs(sh), abs(sh + Δs))
                        @test Ow[ℓ, m] == 0
                    end
                end
            end
            @test ð(w) == R₊(w)
            @test ð̄(w) == -R₋(w)
            @test spin(ð̄(ð(w))) === sh
            @test spin(ð(ð(w))) === h(s + 2)
            @test spin(ð̄(ð̄(w))) === h(s - 2)

            # The documented actions on mode weights, read through natural indexing: for ℓ ≥
            # |s|
            #   {L² f}ₗₘ = ℓ(ℓ+1) fₗₘ,   {Lz f}ₗₘ = m fₗₘ,   {Rz f}ₗₘ = s fₗₘ,
            #   {L₊ f}ₗₘ = √((ℓ+m)(ℓ-m+1)) fₗ,ₘ₋₁,   {L₋ f}ₗₘ = √((ℓ-m)(ℓ+m+1)) fₗ,ₘ₊₁,
            #   {ð f}ₗₘ = √((ℓ-s)(ℓ+s+1)) fₗₘ,       {ð̄ f}ₗₘ = -√((ℓ+s)(ℓ-s+1)) fₗₘ.
            # The eigenvalues ℓ(ℓ+1), m and s are formed as `Rational`s here and converted,
            # independently of the package's numerator arithmetic; every one of them is
            # exactly representable, so the diagonal relations are exact.
            L²w, Lzw, L₊w, L₋w, R²w, Rzw = L²(w), Lz(w), L₊(w), L₋(w), R²(w), Rz(w)
            R₊w, R₋w, ðw, ð̄w = R₊(w), R₋(w), ð(w), ð̄(w)
            ϵ = 10eps(T)
            for (ℓ, m) in modes(w)
                ℓr, mr = Rational(ℓ), Rational(m)
                if ℓ < abs(sh)
                    for Ow in (L²w, Lzw, L₊w, L₋w, R²w, Rzw, R₊w, R₋w, ðw, ð̄w)
                        @test Ow[ℓ, m] == 0
                    end
                else
                    @test L²w[ℓ, m] == T(ℓr * (ℓr + 1)) * w[ℓ, m]
                    @test R²w[ℓ, m] == T(ℓr * (ℓr + 1)) * w[ℓ, m]
                    @test Lzw[ℓ, m] == T(mr) * w[ℓ, m]
                    @test Rzw[ℓ, m] == T(s) * w[ℓ, m]
                    if m == -ℓ
                        @test L₊w[ℓ, m] == 0
                    else
                        @test L₊w[ℓ, m] ≈ √(T((ℓ + m) * (ℓ - m + 1))) * w[ℓ, m - 1] rtol=ϵ
                    end
                    if m == ℓ
                        @test L₋w[ℓ, m] == 0
                    else
                        @test L₋w[ℓ, m] ≈ √(T((ℓ - m) * (ℓ + m + 1))) * w[ℓ, m + 1] rtol=ϵ
                    end
                    @test ðw[ℓ, m] ≈ √(T((ℓ - sh) * (ℓ + sh + 1))) * w[ℓ, m] rtol=ϵ
                    @test R₊w[ℓ, m] ≈ √(T((ℓ - sh) * (ℓ + sh + 1))) * w[ℓ, m] rtol=ϵ
                    @test ð̄w[ℓ, m] ≈ -√(T((ℓ + sh) * (ℓ - sh + 1))) * w[ℓ, m] rtol=ϵ
                    @test R₋w[ℓ, m] ≈ √(T((ℓ + sh) * (ℓ - sh + 1))) * w[ℓ, m] rtol=ϵ
                end
            end
        end
    end

    # Integer data is promoted to Float64
    wi = ModeWeights(collect(1:Ysize(1//2, 5//2)), 1//2)
    @test L²(wi) isa ModeWeights{Float64, HalfOddInteger}
    @test parent(L²(wi)) == L²(1//2, 1//2, 5//2) * parent(wi)
    @test ð(wi) isa ModeWeights{Float64, HalfOddInteger}
    @test spin(ð(wi)) === h(3//2)
end

# The items above cover construction, indexing, the operators and evaluation.  What is left
# is the small change of the array interface that generic code reaches for — the dimension
# queries, the conversions, the banded products and the unary operators.

@testitem "ModeWeights: the rest of the array interface" begin
    using LinearAlgebra: Bidiagonal, Tridiagonal, Diagonal, dot
    import SphericalFunctions: ℓₘᵢₙ, ℓₘₐₓ, spin, HalfOddInteger
    using Random

    rng = Random.Xoshiro(2026)
    s, lo, hi = -2, 2, 5
    n = Ysize(lo, hi)
    w = ModeWeights(randn(rng, ComplexF64, n), s, lo, hi)

    # Dimension queries: a `ModeWeights` is one-dimensional, and trailing dimensions behave
    # as they do for an ordinary vector
    @test ndims(w) == 1
    @test ndims(typeof(w)) == 1
    @test size(w) == (n,) && size(w, 1) == n && size(w, 2) == 1
    @test axes(w, 1) == Base.OneTo(n)
    @test axes(w, 2) == Base.OneTo(1)
    @test collect(keys(w)) == collect(1:n)
    @test eachindex(w) == Base.OneTo(n)
    @test firstindex(w) == 1 && lastindex(w) == n

    # The three conversions all give the same plain vector
    @test Array(w) == collect(w)
    @test Vector(w) == collect(w)
    @test Array(w) isa Vector{ComplexF64}
    @test w[2:5] == collect(w)[2:5]
    # A position may be a one-dimensional `CartesianIndex`, as broadcasting writes it
    @test w[CartesianIndex(3)] === w[3]
    let u = copy(w)
        u[CartesianIndex(3)] = 7
        @test u[3] == 7 && w[3] != 7
    end

    # Equality and `dot` both work with a plain vector on either side.  `dot` conjugates —
    # it is deliberately *not* what evaluating the weights does.
    v = collect(w)
    @test isequal(v, w)
    @test isequal(w, v)
    @test v == w && w == v
    @test dot(v, w) == dot(v, v)
    @test dot(w, v) == dot(v, v)
    # Between two ModeWeights the labels count, for `≈` and `dot` as for `==`: the same numbers
    # under the opposite spin weight are the weights of a different function
    wflip = ModeWeights(copy(v), -s, lo, hi)
    @test !(w ≈ wflip) && w != wflip
    @test_throws "different labels" dot(w, wflip)
    @test w ≈ copy(w) && dot(w, copy(w)) ≈ dot(v, v)

    # Unary plus is the identity, and unary minus negates
    @test +w === w
    @test collect(-w) == -v
    @test spin(-w) == s && ℓₘᵢₙ(-w) == lo && ℓₘₐₓ(-w) == hi

    # The banded matrices the operators produce multiply the storage of a `ModeWeights`, and a
    # product with the container itself is refused for each of them, on either side
    d = randn(rng, ComplexF64, n)
    B = Bidiagonal(d, randn(rng, ComplexF64, n-1), :U)
    T3 = Tridiagonal(randn(rng, ComplexF64, n-1), d, randn(rng, ComplexF64, n-1))
    for A in (Diagonal(d), B, T3)
        @test A * array_view(w) == A * v
        @test_throws ArgumentError A * w
        @test_throws ArgumentError A \ w
    end
    @test_throws ArgumentError w * reshape(v, 1, n)
    @test array_view(w) * reshape(v, 1, n) == v * reshape(v, 1, n)

    # The arithmetic of mode weights keeps the wrapper, and the labels with it
    b = w .+ v
    @test b isa ModeWeights
    @test spin(b) == s && ℓₘᵢₙ(b) == lo && ℓₘₐₓ(b) == hi
    @test collect(b) == 2v
    @test (2 .* w) isa ModeWeights
    @test (w .+ w) isa ModeWeights
    @test collect(w .* 2) == v .* 2
end

@testitem "ModeWeights broadcasting: a constant is refused however it is written" begin
    import SphericalFunctions: ModeWeights, spin, Ysize
    import Random

    rng = Random.Xoshiro(20260924)

    # Adding the same number to every mode weight gives the weights of no function, whether
    # the number is written as a literal, computed within the same broadcast, or given as a
    # vector or tuple of one element that broadcasting extends to every mode
    w = ModeWeights(randn(rng, ComplexF64, Ysize(0, 1)), 0)
    a, b = 2.0, 3.0
    @test_throws ArgumentError w .+ 1
    @test_throws ArgumentError w .+ a .* b
    @test_throws ArgumentError w .- sqrt.(4)
    @test_throws ArgumentError (1 .* 1) .- w
    @test_throws ArgumentError w .+ [1.0]
    @test_throws ArgumentError w .+ (1,)
    @test_throws ArgumentError w .+ fill(1.0)
    @test_throws ArgumentError w .+ Ref(1.0)

    # ... while sums of weights, and products and quotients with numbers or with a vector of
    # one factor per mode, keep the labels
    v = collect(1.0:4.0)
    for x ∈ (w .+ 2 .* w, (a .* w) .- w ./ b, v .* w .+ w, 2 .* w)
        @test x isa ModeWeights{ComplexF64}
        @test spin(x) == 0
    end
end


@testitem "ModeWeights broadcasting: complex.(a, b) compares the labels of both parts" begin
    import SphericalFunctions: ModeWeights, spin, Ysize
    import Random

    rng = Random.Xoshiro(20260925)

    # Real and imaginary parts combine into the weights of one function only if they are the
    # weights of the same spin weight and range of ℓ
    re = ModeWeights(randn(rng, Ysize(2, 4)), 2)
    im₋ = ModeWeights(randn(rng, Ysize(2, 4)), -2)
    @test_throws ArgumentError complex.(re, im₋)
    @test_throws ArgumentError Complex.(re, im₋)
    @test_throws ArgumentError ComplexF64.(re, im₋)
    # (ℓ ∈ 10:10 has as many modes as ℓ ∈ 2:4, so only the labels tell them apart)
    @test_throws ArgumentError complex.(re, ModeWeights(randn(rng, Ysize(10, 10)), 2, 10, 10))
    # A number as one of the parts adds a constant to every weight, which `+` refuses too
    @test_throws ArgumentError complex.(re, 1.0)
    @test_throws ArgumentError complex.(1.0, re)

    # Parts with the same labels combine, and the result has those labels
    im₊ = ModeWeights(randn(rng, Ysize(2, 4)), 2)
    for z ∈ (complex.(re, im₊), Complex.(re, im₊), ComplexF64.(re, im₊))
        @test z isa ModeWeights{ComplexF64}
        @test spin(z) == 2
        @test parent(z) == complex.(parent(re), parent(im₊))
    end
    # ... while a one-argument conversion keeps the labels of its argument
    @test complex.(re) isa ModeWeights{ComplexF64} && spin(complex.(re)) == 2
    @test float.(ModeWeights(collect(1:9), 0)) isa ModeWeights{Float64}
end


@testitem "ModeWeights: a rotor is not a factor, and mode weights are not quaternions" begin
    import SphericalFunctions: ModeWeights, D, spin, Ysize
    using Quaternionic: Rotor, Quaternion
    import Random

    rng = Random.Xoshiro(20260926)

    # A `Rotor` is a `Number`, but its product with mode weights is not the weights of any
    # function of the same spin weight; rotating the function is `D(R, ℓₘₐₓ) * w`
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    w = ModeWeights(randn(rng, ComplexF64, Ysize(0, 2)), 0)
    @test_throws ArgumentError R * w
    @test_throws ArgumentError w * R
    @test_throws ArgumentError R .* w
    @test_throws ArgumentError w .* R
    @test_throws ArgumentError Quaternion(1.0, 2.0, 3.0, 4.0) * w
    @test_throws ArgumentError ModeWeights(fill(Quaternion(1.0, 0.0, 0.0, 0.0), Ysize(0, 2)), 0)

    # Numbers scale the weights, and the rotation acts on them
    @test 2.0 * w isa ModeWeights{ComplexF64}
    @test (1 + 2im) * w isa ModeWeights{ComplexF64}
    @test w / 2 isa ModeWeights{ComplexF64}
    @test D(R, 2) * w isa ModeWeights{ComplexF64}
end


@testitem "ModeWeights: copying into another range of ℓ" begin
    import SphericalFunctions: ModeWeights, spin, Ysize, Yindex, HalfOddInteger, ð, ð̄, R₊, R₋
    import LinearAlgebra: mul!, dot
    import Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(20260927)
    h(x) = HalfOddInteger(x)
    lo(w) = SphericalFunctions.ℓₘᵢₙ(w)
    hi(w) = SphericalFunctions.ℓₘₐₓ(w)
    Rs = randn(rng, Rotor{Float64}, 4)

    # Every combination of a source range and a target range, for either kind of index: the
    # labels are those asked for, each weight the two ranges share is copied exactly, every
    # ℓ the source lacks is zero, and the ℓ < |s| of the target, which belong to no
    # harmonic, are zero whatever the source held there
    for (s, (a, b), targets) in (
        (2, (0, 6), ((0, 6), (2, 6), (3, 4), (0, 8), (5, 8), (7, 6))),
        (-1, (1, 4), ((0, 4), (1, 4), (2, 6), (1, 1))),
        (1//2, (1//2, 9//2), ((1//2, 9//2), (3//2, 7//2), (1//2, 11//2))),
        (-3//2, (1//2, 7//2), ((3//2, 7//2), (1//2, 9//2), (5//2, 3//2))),
    )
        w = ModeWeights(randn(rng, ComplexF64, Ysize(a, b)), s, a, b)
        for (c, d) in targets
            w′ = ModeWeights(w; ℓₘᵢₙ=c, ℓₘₐₓ=d)
            @test w′ isa ModeWeights{ComplexF64, typeof(spin(w))}
            @test spin(w′) === spin(w) && lo(w′) == c && hi(w′) == d
            @test parent(w′) !== parent(w)
            for ℓ in lo(w′):hi(w′), m in -ℓ:ℓ
                expected = ℓ < abs(spin(w)) || !(lo(w) ≤ ℓ ≤ hi(w)) ? zero(ComplexF64) : w[ℓ, m]
                @test w′[ℓ, m] === expected
            end
            # The same with the other spelling of the keywords, and with `HalfOddInteger`s
            @test ModeWeights(w; ell_min=c, ell_max=d) == w′
            c isa Rational && @test ModeWeights(w; ℓₘᵢₙ=h(c), ℓₘₐₓ=h(d)) == w′
            # The function is unchanged wherever the target range covers the source's
            # harmonics.  (Measured: exactly, since the weights added or dropped multiply
            # harmonics that are exactly zero; 4 eps allows for a change of summation
            # order.)
            if lo(w′) ≤ max(lo(w), abs(spin(w))) && hi(w′) ≥ hi(w)
                @test maximum(abs, w′(Rs) - w(Rs)) ≤ 4eps() * maximum(abs, w(Rs))
            end
        end
    end

    # The defaults give the range the transforms' analysis has, |s| ≤ ℓ ≤ ℓₘₐₓ(w), so a
    # container that already has it is copied, not relabelled
    w = ModeWeights(randn(rng, ComplexF64, Ysize(1, 6)), -1)
    c = ModeWeights(w)
    @test c == w && parent(c) !== parent(w)
    c[1, 0] = 0
    @test w[1, 0] != 0
    wz = ModeWeights(randn(rng, ComplexF64, Ysize(0, 6)), -2, 0, 6)
    @test lo(ModeWeights(wz)) == 2 && ModeWeights(wz)[2, 1] == wz[2, 1]
    @test ModeWeights(wz; ℓₘᵢₙ=5, ℓₘₐₓ=4) |> isempty

    # The spin-changing operators keep the range of their input; copying the result into the
    # range of the new spin weight makes it comparable with, and addable to, weights built
    # there, and it is the same function.  This is the case that needs truncation (ð on s ≥
    # 0) and the one that needs padding (ð̄ on s > 0).
    for (s, op, Δ) in ((1, ð, 1), (1, ð̄, -1), (-2, ð, 1), (-2, ð̄, -1), (0, R₊, 1), (1//2, ð, 1), (3//2, ð̄, -1))
        L = s isa Rational ? 11//2 : 6
        local w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), L)), s)
        ow = op * w
        @test lo(ow) == abs(s) && spin(ow) == s + Δ
        ow′ = ModeWeights(ow)
        @test lo(ow′) == abs(s + Δ) && hi(ow′) == L
        ref = ModeWeights(zeros(ComplexF64, Ysize(abs(s + Δ), L)), s + Δ)
        @test ow′ + ref == ow′
        @test ow′ ≈ ModeWeights(ow; ℓₘᵢₙ=abs(s + Δ))
        @test maximum(abs, ow′(Rs) - ow(Rs)) ≤ 4eps() * maximum(abs, ow(Rs))  # measured: 0
        # A destination of `mul!` must have the labels of the result, and the refusal says
        # how to allocate one and how to move the result into another range afterwards
        if abs(s + Δ) != abs(s)
            T = ComplexF64
            @test_throws ArgumentError mul!(ModeWeights{T}(undef, s + Δ, L), op, w)
            @test_throws "ModeWeights{ComplexF64}(undef, $(s + Δ), $(lo(w)), $L)" mul!(ModeWeights{T}(undef, s + Δ, L), op, w)
            @test_throws "ModeWeights(w′; ℓₘᵢₙ, ℓₘₐₓ)" mul!(ModeWeights{T}(undef, s + Δ, L), op, w)
            @test mul!(ModeWeights{T}(undef, s + Δ, lo(w), L), op, w) == ow
        end
    end

    # The keywords must be indices of the kind of `w`'s own, and the range must be valid
    w = ModeWeights(randn(rng, ComplexF64, Ysize(1, 4)), 1)
    wh = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), 1//2)
    kind = "must be an index of the kind of `w`'s own"
    @test_throws ArgumentError ModeWeights(w; ℓₘᵢₙ=1//2)
    @test_throws kind ModeWeights(w; ℓₘᵢₙ=1//2)
    @test_throws kind ModeWeights(wh; ℓₘᵢₙ=1)
    @test_throws kind ModeWeights(wh; ell_max=4)
    @test_throws "`Int32` is narrower than `Int`" ModeWeights(w; ℓₘᵢₙ=Int32(1))
    # ... and a refusal names the spelling of the keyword that the caller wrote
    @test_throws "The keyword argument `ell_min`" ModeWeights(w; ell_min=Int32(1))
    @test_throws "The keyword argument `ell_max`" ModeWeights(w; ell_max=Int32(3))
    @test_throws "The keyword argument `ℓₘₐₓ`" ModeWeights(w; ℓₘₐₓ=Int32(3))
    @test_throws "`Float64` is not an index type" ModeWeights(w; ℓₘₐₓ=4.0)
    @test_throws "1//3 is neither an integer nor a half-odd-integer" ModeWeights(wh; ℓₘᵢₙ=1//3)
    @test_throws "ℓₘᵢₙ=-1 must be non-negative" ModeWeights(w; ℓₘᵢₙ=-1)
    @test_throws "must be at least ℓₘᵢₙ-1" ModeWeights(w; ℓₘᵢₙ=4, ℓₘₐₓ=1)
    # ... and storage that has been resized since the labels were attached is refused
    v = randn(rng, ComplexF64, Ysize(1, 4))
    wv = ModeWeights(v, 1)
    resize!(v, 3)
    @test_throws DimensionMismatch ModeWeights(wv; ℓₘₐₓ=2)
end


@testitem "ModeWeights: products with plain matrices are refused" begin
    import SphericalFunctions: ModeWeights, spin, Ysize, sYlm, sYlm_matrix, relabel, ð
    import LinearAlgebra: Diagonal, Bidiagonal, Tridiagonal, Symmetric, UpperTriangular, dot
    import Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(20260928)
    Rs = randn(rng, Rotor{Float64}, 5)
    w = ModeWeights(randn(rng, ComplexF64, Ysize(1, 6)), 1)
    dw = ð * w                       # s = 2 over ℓ ∈ 1:6, as long as `w`
    n = length(w)

    # A synthesis matrix of one spin weight has the length of weights of another, so the product
    # would be silently wrong; it is refused, and the labelled evaluation checks the spin weight
    Y = sYlm_matrix(Rs, 6, 1)
    @test size(Y, 2) == length(dw)
    @test_throws ArgumentError Y * dw
    @test_throws "sYlm(R⃗, ℓₘₐₓ, s) * w" Y * dw
    @test_throws "spin weight s=2" sYlm(Rs, 6, 1) * dw
    @test Y * array_view(w) ≈ w(Rs) rtol=1e-13
    @test sYlm(Rs, 6, 1) * w ≈ w(Rs) rtol=1e-13

    # Every matrix type, and both sides of the product, and `\`
    for A in (
        randn(rng, n, n), Diagonal(ones(n)), Bidiagonal(ones(n), ones(n - 1), :L),
        Tridiagonal(ones(n - 1), ones(n), ones(n - 1)), Symmetric(randn(rng, n, n)),
        UpperTriangular(randn(rng, n, n)), ð(1, 1, 6), w', transpose(w),
    )
        @test_throws "A plain matrix cannot be checked against the labels" A * w
    end
    @test_throws "`w * A` is refused" w * randn(rng, n, 2)
    @test_throws "`A \\ w` is refused" randn(rng, n, n) \ w
    @test_throws "`A \\ w` is refused" Diagonal(ones(n)) \ w

    # `relabel` puts `w`'s own labels on the array, spin weight included, so the result of a
    # spin-changing operator matrix is labelled by the constructor with the new spin weight
    A = ð(1, 1, 6) * array_view(w)
    @test spin(relabel(w, A)) == 1
    @test ModeWeights(A, 2, 1, 6) == dw
end


@testitem "ModeWeights: equality, hashing and the linear algebra that keeps the labels" setup=[RefusalChecks] begin
    import SphericalFunctions: ModeWeights, spin, Ysize, relabel, array_view, D, HalfOddInteger
    import LinearAlgebra: I, rmul!, lmul!, axpy!, axpby!
    import Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(20260924)
    labels(x) = (spin(x), SphericalFunctions.ℓₘᵢₙ(x), SphericalFunctions.ℓₘₐₓ(x))
    h = HalfOddInteger
    for (s, lo, hi) ∈ ((-2, 2, 5), (0, 0, 3), (h(1//2), h(1//2), h(7//2)))
        n = Ysize(lo, hi)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(copy(data), s, lo, hi)
        wflip = ModeWeights(copy(data), -s, lo, hi)

        # `==`, `isequal` and `hash` agree.  Between two sets of weights they count the
        # labels; against a plain vector they compare the numbers alone, so the hash is that
        # of the numbers, and weights under other labels may share it
        @test isequal(w, copy(w)) && hash(w) == hash(copy(w)) == hash(data)
        @test length(Set([w, copy(w)])) == 1 && haskey(Dict(w => 1), copy(w))
        s != 0 && @test !isequal(w, wflip) && length(unique([w, copy(w), wflip])) == 2
        wn = copy(w)
        wn[1] = NaN
        @test isequal(wn, wn) && isequal(wn, copy(wn)) && wn != wn
        @test hash(wn) == hash(copy(wn))
        w₊, w₋ = copy(w), copy(w)
        w₊[1], w₋[1] = 0.0, -0.0
        @test w₊ == w₋ && !isequal(w₊, w₋)

        # Linear indexing by position, with any integer type, a colon or a vector of
        # positions; any index but one position gives plain numbers
        @test w[Int32(2)] === data[2] && w[UInt8(3)] === data[3]
        @test w[:] == data && w[:] isa Vector{ComplexF64} && w[[3, 1]] == data[[3, 1]]
        c = copy(w)
        c[1:2] = [10, 20]
        c[[n]] = [30]
        @test c[1] == 10 && c[2] == 20 && c[n] == 30 && labels(c) == labels(w)
        c[:] = data
        @test c == w

        # The linear-space operations keep the labels, and the in-place ones return their
        # destination
        z = zero(w)
        @test z isa ModeWeights && labels(z) == labels(w) && all(iszero, array_view(z))
        @test fill!(copy(w), 2) == ModeWeights(fill(2.0 + 0im, n), s, lo, hi)
        c = copy(w)
        @test rmul!(c, 2) === c && c == 2w
        @test lmul!(3, c) === c && c == 6w
        @test I * w == w && labels(I * w) == labels(w) && (2I) * w == 2w && w * (2I) == 2w
        c = copy(w)
        @test copyto!(c, 2 .* data) === c && array_view(c) == 2 .* data
        @test copyto!(c, w) === c && c == w
        @test copy!(c, 2w) === c && c == 2w
        y = copy(w)
        @test axpy!(2, w, y) === y && y ≈ 3w
        @test axpby!(2, w, 3, y) === y && y ≈ 11w
        # ... and refuse what would label numbers the labels do not describe
        @test refuses(() -> copyto!(copy(w), zeros(n + 1)), DimensionMismatch, "length $(n + 1)")
        if s != 0
            @test refuses(() -> copyto!(copy(w), wflip), ArgumentError, "different labels")
            @test refuses(() -> copy!(copy(w), wflip), ArgumentError, "different labels")
            @test refuses(() -> axpy!(2, wflip, copy(w)), ArgumentError, "different labels")
            @test refuses(() -> axpby!(2, wflip, 1, copy(w)), ArgumentError, "different labels")
        end
        # A sum is of two functions; `w .+ data` is the explicit form with a plain vector
        @test_throws MethodError w + data
        @test w .+ data == 2w

        # `relabel` wraps a vector of the same length under the labels of `w`, without copying
        # it, and refuses another length
        A = 2 .* data
        r = relabel(w, A)
        @test r isa ModeWeights && labels(r) == labels(w) && array_view(r) === A
        @test refuses(() -> relabel(w, zeros(n + 1)), DimensionMismatch, "length $(n + 1)")
    end

    # A rotation of weights that hold no ℓ at all needs no block of 𝔇, and gives weights
    # with the same labels and no entries
    R = randn(rng, Rotor{Float64})
    for (w, 𝔇) ∈ (
        (ModeWeights(ComplexF64[], 2), D(R, 0)), (ModeWeights(ComplexF64[], 3//2), D(R, 1//2))
    )
        rot = 𝔇 * w
        @test rot isa ModeWeights && isempty(array_view(rot)) && labels(rot) == labels(w)
    end
end
