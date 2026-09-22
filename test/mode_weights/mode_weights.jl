# Tests of `ModeWeights`, the wrapper for a vector of mode weights in the canonical ordering
#
#     [ f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ ],
#
# together with the spin weight `s` (see `src/mode_weights/mode_weights.jl`).  The oracles are
# the closed-form indexing functions `Ysize`, `Yindex` and `Yrange` (tested in `indexing.jl`),
# the operator matrices of `src/utilities/operators.jl` applied to the raw storage (tested
# against explicit differential operators in `test/operators.jl`), and — for the evaluation
# `w(R)` — the closed-form `sYlm(s, ℓ, m, θ, ϕ)` of the `Utilities` snippet, which is the
# explicit sum from `docs/src/30-conventions/01-summary.md` and shares no code with the package.
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
    # same array (not a copy), and the array-like queries describe that storage.  The
    # range ℓₘₐₓ = ℓₘᵢₙ - 1 is the empty container.
    for s in -3:3, ℓₘᵢₙ in 0:3, ℓₘₐₓ in ℓₘᵢₙ-1:7
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        data = randn(rng, ComplexF64, n)
        w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
        @test w isa ModeWeights{ComplexF64}
        # Since 3.0 this is an `AbstractModeContainer`, not an `AbstractVector`: `array_view` is
        # the route to the flat storage, and array semantics come with it.
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

    # Narrower integer types are kept, alone or promoted together
    for IT in (Int8, Int32, Int64)
        data = zeros(Ysize(1, 3))
        w = ModeWeights(data, IT(-1), IT(1), IT(3))
        @test w isa ModeWeights{Float64, IT}
        @test spin(w) === IT(-1)
        @test SphericalFunctions.ℓₘᵢₙ(w) === IT(1)
        @test SphericalFunctions.ℓₘₐₓ(w) === IT(3)
        w = ModeWeights(data, IT(-1))
        @test spin(w) == -1
        @test SphericalFunctions.ℓₘᵢₙ(w) == 1
        @test SphericalFunctions.ℓₘₐₓ(w) == 3
        @test typeof(spin(w)) == typeof(SphericalFunctions.ℓₘᵢₙ(w)) == typeof(SphericalFunctions.ℓₘₐₓ(w))
        @test length(ModeWeights{Float64}(undef, IT(2), IT(4))) == Ysize(2, 4)
    end

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


@testitem "ModeWeights indexing" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yindex, Yrange
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
end


@testitem "ModeWeights array behaviour" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, L², Lz, ð
    import LinearAlgebra: norm, dot, Diagonal
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

        # Shape-preserving broadcasts return a new ModeWeights with the same s and ℓ range,
        # holding what the same broadcast on the storage gives — including broadcasts that
        # change the element type or mix in a plain vector of the same length
        for (result, expected) in (
            (2 .* w, 2 .* data),
            (w .+ w, data .+ data),
            (w .* conj.(w), data .* conj.(data)),
            (w .+ 1, data .+ 1),
            (w .- data, zeros(ComplexF64, n)),
            (data .* w, data .* data),
            (-w, -data),
            (conj.(w), conj.(data)),
            (abs.(w), abs.(data)),
            (real.(w), real.(data)),
            (w .== w, trues(n)),
        )
            @test result isa ModeWeights
            @test same_range(result, w)
            @test parent(result) == expected
            @test eltype(result) == eltype(expected)
            @test parent(result) !== data
        end
        @test parent(w) === data
        @test parent(w) == data
        # In-place broadcast into a similar container
        w2 = similar(w)
        w2 .= 2 .* w .+ 1
        @test w2 isa ModeWeights
        @test same_range(w2, w)
        @test parent(w2) == 2 .* data .+ 1
        @test parent(w) == data
        w2 .= w
        @test w2 == w
        @test parent(w2) !== data

        # `similar` keeps the wrapper for the same length (with any element type), and falls
        # back to a plain array for any other shape
        @test similar(w) isa ModeWeights{ComplexF64}
        @test same_range(similar(w), w)
        @test length(similar(w)) == n
        @test parent(similar(w)) !== data
        @test similar(w, Float32) isa ModeWeights{Float32}
        @test same_range(similar(w, Float32), w)
        @test similar(w, n) isa ModeWeights{ComplexF64}
        @test similar(w, Float64, n) isa ModeWeights{Float64}
        @test similar(w, n + 1) isa Vector{ComplexF64}
        @test length(similar(w, n + 1)) == n + 1
        @test similar(w, Float64, n + 1) isa Vector{Float64}
        @test similar(w, (n, 2)) isa Matrix{ComplexF64}
        @test size(similar(w, (n, 2))) == (n, 2)
        # An operator *matrix* gives an unlabelled vector, because it cannot know what spin
        # weight its result has; the operator *function* gives a correctly labelled one.
        # (`ℓₘᵢₙ` and `ℓₘₐₓ` are loop variables here, not the accessor functions.)
        let wδ = ð(s, ℓₘᵢₙ, ℓₘₐₓ) * w
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
        # `map` keeps the wrapper, too
        @test map(abs, w) isa ModeWeights{Float64}
        @test parent(map(abs, w)) == abs.(data)

        # Reductions and equality agree with the storage
        @test sum(w) ≈ sum(data) atol=ϵ rtol=ϵ
        @test norm(w) ≈ norm(data) rtol=ϵ
        @test dot(w, w) ≈ dot(data, data) rtol=ϵ
        @test dot(w, data) ≈ dot(data, data) rtol=ϵ
        @test w' * w ≈ dot(data, data) rtol=ϵ
        @test maximum(abs, w) == maximum(abs, data)
        @test w == data
        @test data == w
        @test isequal(w, data)
        @test collect(w) == data
        @test Vector(w) == data
        @test Vector(w) isa Vector{ComplexF64}
        @test w[1:min(3, n)] isa Vector{ComplexF64}

        # Matrix products: a product that changes the length is a plain Vector; one that
        # keeps the length has the same values as the product with the storage
        M = randn(rng, ComplexF64, n + 3, n)
        @test M * w isa Vector{ComplexF64}
        @test M * w ≈ M * data rtol=ϵ
        Msq = randn(rng, ComplexF64, n, n)
        @test Msq * w isa AbstractVector{ComplexF64}
        @test length(Msq * w) == n
        @test Msq * w ≈ Msq * data rtol=ϵ
        # The operator matrices applied by hand give the operators' numbers
        @test L²(s, ℓₘᵢₙ, ℓₘₐₓ) * w == L²(w)
        @test Lz(s, ℓₘᵢₙ, ℓₘₐₓ) * w == Lz(w)
        @test ð(s, ℓₘᵢₙ, ℓₘₐₓ) * w == ð(w)
        @test Diagonal(ones(n)) * w == data
        # Outer-product shapes cannot be a ModeWeights and must fall back to plain matrices
        @test w .+ transpose(w) isa Matrix{ComplexF64}
        @test w * w' isa Matrix{ComplexF64}
        @test parent(w) .+ w' isa Matrix{ComplexF64}
        @test size(parent(w) .+ w') == (n, n)

        # `show` does not throw, and the three-argument form names the parameters
        str = sprint(show, MIME("text/plain"), w)
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

    # A length-1 ModeWeights broadcast against a longer vector changes shape, so the result
    # is a plain Vector
    w1 = ModeWeights([3.0], 0)
    @test w1 .+ [1.0, 2.0, 3.0] == [4.0, 5.0, 6.0]
    @test w1 .+ [1.0, 2.0, 3.0] isa Vector{Float64}
    @test 2 .* w1 isa ModeWeights{Float64}

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

        # Every operator is the corresponding matrix applied to the storage, wrapped with the
        # same ℓ range; the input is untouched
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
        # construction (all four `ModeWeights` methods come from one metaprogrammed loop in
        # `src/mode_weights/mode_weights.jl`), so they pin the aliasing rather than the sign
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
            # The coefficients are computed at the precision of T (no Float64 leaks): the
            # relative error of the ladder coefficients is a few eps(T)
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
    # (The `Utilities` snippet defines the closed-form `sYlm(s, ℓ, m, θ, ϕ)`, so `sYlm` here
    # is that reference function, never the package's own harmonics.)
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yindex, ð, ð̄, R₊, R₋
    import Quaternionic: Rotor, from_spherical_coordinates, from_euler_angles
    import LinearAlgebra: norm
    import DoubleFloats: Double64
    import Random

    rng = Random.Xoshiro(20260914)

    # The independent reference for the harmonics at a rotor.  `sYlm(s, ℓ, m, θ, ϕ)` from the
    # `Utilities` snippet is the explicit sum of `docs/src/30-conventions/01-summary.md`, but it is
    # only the γ = 0 slice of the rotation group.  The same page gives the γ dependence as
    #     ₛYₗₘ(𝐑 exp(γ𝐤/2)) = exp(-isγ) ₛYₗₘ(𝐑),
    # so with 𝐑 = exp(α𝐤/2) exp(β𝐣/2) exp(γ𝐤/2) = `Quaternionic.from_euler_angles(α, β, γ)`,
    #     ₛYₗₘ(𝐑) = exp(-isγ) ₛYₗₘ(θ=β, ϕ=α).
    # Neither the closed form nor `from_euler_angles` runs any SphericalFunctions code, so
    # this is a genuinely independent oracle for `w(R)`.  Entries with ℓ < |s| are not
    # harmonics at all and are set to zero; the weights there must not contribute, which is
    # checked separately below.
    #
    # The rotors are *built* from Euler angles rather than decomposed into them: going the
    # other way, `to_euler_angles` loses about half the significant digits of β near the
    # poles (measured: the recovered β is off by 2.2e-16 at β = 0 and by 2.2e-30 at β = 0.01,
    # in Double64, the usual `acos` cancellation), which would cap this oracle at √eps.
    function Yreference(s, ℓₘᵢₙ, ℓₘₐₓ, α::T, β::T, γ::T) where {T}
        [
            ℓ < abs(s) ? zero(Complex{T}) : cis(-s * γ) * sYlm(s, ℓ, m, β, α)
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

                # With ℓₘᵢₙ = 0, the weights with ℓ < |s| do not contribute, so junk there
                # is harmless and the ℓₘᵢₙ = |s| container gives the same value
                if ℓₘᵢₙ == 0 && s != 0
                    wj = copy(w)
                    for ℓ in 0:abs(s)-1, m in -ℓ:ℓ
                        wj[ℓ, m] = randn(rng, Complex{T})
                    end
                    @test wj(R) ≈ f atol=ϵ rtol=ϵ
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
        # sYlm(s, ℓ, m, θ, ϕ) at R = from_spherical_coordinates(θ, ϕ) — which is the γ = 0
        # slice, so the closed form applies with no rotation law at all — including at the
        # poles.  Measured: at most 4.1 eps(T) in absolute value, for |ₛYₗₘ| up to 0.94.
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
                    @test e(R) ≈ sYlm(s, ℓ, m, θ, ϕ) atol=ϵ rtol=ϵ
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
# The items below repeat the checks above for a `ModeWeights` whose spin weight and ℓ range are
# half-odd-integers.  The indices are written as `Rational`s with denominator 2, which is the
# public spelling, and the parameters that come back are compared against `HalfOddInteger`s,
# which is what the package stores.  Ranges such as `ℓₘᵢₙ-1:7//2` step through the
# half-odd-integers in the `Rational` spelling.

@testitem "ModeWeights half-integer construction" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yrange, HalfOddInteger
    import DoubleFloats: Double64
    import Random
    using Quaternionic: Rotor

    rng = Random.Xoshiro(20260916)
    h(x) = HalfOddInteger(x)

    # The four-argument form, with either spelling: the parameters come back as
    # `HalfOddInteger`s, the storage is the very same array, and ℓₘₐₓ = ℓₘᵢₙ - 1 is the empty
    # container
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
    # ... a `Rational` that is not a half-odd-integer ...
    @test_throws "denominator 2" ModeWeights(zeros(2), 1//3)
    @test_throws "denominator 2" ModeWeights(zeros(4), 1//1)
    @test_throws "denominator 2" ModeWeights(zeros(2), 1//2, 1//4, 1//2)
    @test_throws "denominator 2" ModeWeights{Float64}(undef, 1//2, 3//1)
    # ... and, in every form, a mixture of integer and half-integer indices, which is refused
    # with a message naming both spellings
    mixed = "all be integers, like 3, or all be half-odd-integers, like 7//2"
    for args in ((1//2, 0, 7//2), (1, 1//2, 7//2), (1//2, 1//2, 3), (h(1//2), 0, h(7//2)))
        @test_throws mixed ModeWeights(zeros(20), args...)
        @test_throws mixed ModeWeights{Float64}(undef, args...)
    end
    @test_throws mixed ModeWeights(zeros(20), 1//2; ℓₘᵢₙ=0)
    @test_throws mixed ModeWeights(zeros(20), 1; ℓₘᵢₙ=1//2)
    @test_throws mixed ModeWeights{Float64}(undef, 1//2, 3)
    @test_throws mixed ModeWeights{Float64}(undef, 1, 7//2)

    # The integer path is as it was: narrower integer types are kept when they agree, and
    # unified when they differ
    data = zeros(Ysize(1, 3))
    @test ModeWeights(data, Int8(-1), Int8(1), Int8(3)) isa ModeWeights{Float64, Int8}
    @test ModeWeights(data, Int8(-1)) isa ModeWeights{Float64, Int}
    @test ModeWeights(data, Int8(-1), 1, 3) isa ModeWeights{Float64, Int}
    @test ModeWeights(data, Int8(-1), 1, 3) == ModeWeights(data, -1, 1, 3)
    @test ModeWeights{Float64}(undef, Int8(1), 1, Int32(3)) isa ModeWeights{Float64, Int}
    # ... and evaluation, which builds a calculator with the narrow index type, gives the `Int`
    # result exactly
    R = randn(rng, Rotor{Float64})
    cdata = randn(rng, ComplexF64, Ysize(1, 3))
    @test ModeWeights(cdata, Int8(1), Int8(1), Int8(3))(R) == ModeWeights(cdata, 1, 1, 3)(R)
    @test ModeWeights(cdata, Int16(-1), Int16(1), Int16(3))(R) == ModeWeights(cdata, -1, 1, 3)(R)
end


@testitem "ModeWeights half-integer indexing" begin
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
            # Natural indexing reads the canonical position of the storage, in either spelling
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
            @test sprint(show, axes(v)) == "($(h(-ℓ)):1:$(h(ℓ)),)"  # the axis can be printed
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

        # Integer indices into a half-integer `w` are refused with an explanation, as is a
        # mixture of the two kinds among the indices themselves, or a `Rational` that is not
        # a half-odd-integer; none of them touches the storage
        kind = "indices of this `ModeWeights` are half-odd-integers, like 7//2"
        for (ℓ, m) in ((1, 0), (2, 1), (0, 0), (Int8(1), Int8(1)))
            @test_throws kind w[ℓ, m]
            @test_throws kind w[ℓ, m] = 0
            @test_throws kind w[ℓ, :]
        end
        mixed = "all be integers, like 3, or all be half-odd-integers, like 7//2"
        @test_throws mixed w[3//2, 1]
        @test_throws mixed w[1, 1//2]
        @test_throws mixed w[h(3//2), 1] = 0
        @test_throws "denominator 2" w[1//3, 1//3]
        @test_throws "denominator 2" w[1//1, :]
        @test parent(w) == data
    end

    # Conversely, an integer `w` refuses half-integer indices, and is otherwise as it was
    wi = ModeWeights(randn(rng, 9), 0)
    @test_throws "indices of this `ModeWeights` are integers, like 3" wi[3//2, 1//2]
    @test_throws "indices of this `ModeWeights` are integers, like 3" wi[3//2, 1//2] = 0
    @test_throws "indices of this `ModeWeights` are integers, like 3" wi[3//2, :]
    @test_throws "all be integers" wi[1, 1//2]
    @test wi[1, 0] === parent(wi)[Yindex(1, 0, 0)]
    @test wi[Int8(2), Int8(-1)] === parent(wi)[Yindex(2, -1, 0)]
    # Since 3.0 the integer path returns the same container as the half-integer one.
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


@testitem "ModeWeights half-integer array behaviour" begin
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

        # Shape-preserving broadcasts keep the wrapper, with the same `HalfOddInteger`
        # parameters
        for (result, expected) in (
            (2 .* w, 2 .* data),
            (w .+ w, data .+ data),
            (w .* conj.(w), data .* conj.(data)),
            (w .+ 1, data .+ 1),
            (w .- data, zeros(ComplexF64, n)),
            (data .* w, data .* data),
            (-w, -data),
            (conj.(w), conj.(data)),
            (abs.(w), abs.(data)),
            (real.(w), real.(data)),
            (w .== w, trues(n)),
        )
            @test result isa ModeWeights{eltype(expected), HalfOddInteger}
            @test same_range(result, w)
            @test parent(result) == expected
            @test parent(result) !== data
        end
        @test parent(w) === data
        w2 = similar(w)
        w2 .= 2 .* w .+ 1
        @test w2 isa ModeWeights{ComplexF64, HalfOddInteger}
        @test same_range(w2, w)
        @test parent(w2) == 2 .* data .+ 1

        # `similar`, `copy` and `map` keep the wrapper and its parameters
        @test similar(w) isa ModeWeights{ComplexF64, HalfOddInteger}
        @test same_range(similar(w), w)
        @test similar(w, Float32) isa ModeWeights{Float32, HalfOddInteger}
        @test same_range(similar(w, Float32), w)
        @test similar(w, n) isa ModeWeights{ComplexF64, HalfOddInteger}
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
        @test map(abs, w) isa ModeWeights{Float64, HalfOddInteger}
        @test parent(map(abs, w)) == abs.(data)

        # Reductions, equality and conversion agree with the storage
        @test sum(w) ≈ sum(data) atol=ϵ rtol=ϵ
        @test norm(w) ≈ norm(data) rtol=ϵ
        @test dot(w, w) ≈ dot(data, data) rtol=ϵ
        @test w == data
        @test collect(w) == data
        @test Vector(w) isa Vector{ComplexF64}
        @test w[1:3] isa Vector{ComplexF64}

        # An operator *matrix* gives an unlabelled vector; the operator *function* a labelled
        # `ModeWeights`, with the same numbers
        @test L²(s, ℓₘᵢₙ, ℓₘₐₓ) * w isa Vector{ComplexF64}
        @test L²(s, ℓₘᵢₙ, ℓₘₐₓ) * w == parent(L²(w))
        @test Lz(s, ℓₘᵢₙ, ℓₘₐₓ) * w == parent(Lz(w))
        @test ð(s, ℓₘᵢₙ, ℓₘₐₓ) * w == parent(ð(w))
        @test spin(ð(w)) === h(s + 1)
        @test Diagonal(ones(n)) * w == data
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


@testitem "ModeWeights half-integer evaluation" begin
    import SphericalFunctions: ModeWeights, modes, spin, Ysize, Yindex, sYlm, HalfOddInteger
    import Quaternionic: Rotor
    import LinearAlgebra: norm
    import DoubleFloats: Double64
    import Random

    rng = Random.Xoshiro(20260919)

    # The reference is the explicit sum Σ f_{ℓm} ₛYₗₘ(R) over the flat `sYlm`, which is checked
    # against the `sYlmCalculator` for half-integer indices in `test/sYlm/sYlm.jl`, read through
    # natural indexing so that the pairing of weights with harmonics is the one `modes` gives.
    for T in (Float64, Double64)
        ϵ = 100eps(T)
        # (ℓₘₐₓ must be at least |s|: a container with ℓₘₐₓ < |s| holds no harmonics to sum.)
        for s in (-3//2, -1//2, 1//2, 3//2), ℓₘᵢₙ in unique((abs(s), 1//2)), ℓₘₐₓ in (max(abs(s), ℓₘᵢₙ), 7//2, 11//2)
            n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            data = randn(rng, Complex{T}, n)
            w = ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
            for _ in 1:3
                R = randn(rng, Rotor{T})
                Y = array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ))
                expected = sum(
                    w[ℓ, m] * Y[Yindex(ℓ, m, ℓₘᵢₙ)] for (ℓ, m) in modes(w);
                    init=zero(Complex{T})
                )
                f = w(R)
                @test f isa Complex{T}
                @test f ≈ expected atol=ϵ*norm(data) rtol=ϵ
                # Linearity in the weights
                @test (2 .* w)(R) ≈ 2f atol=ϵ*norm(data) rtol=ϵ
            end
        end
    end

    # Weights below |s| do not contribute, because there are no such harmonics
    for s in (3//2, -3//2)
        w = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), s, 1//2, 7//2)
        R = randn(rng, Rotor{Float64})
        f = w(R)
        for m in -1//2:1//2
            w[1//2, m] = 1e6
        end
        @test w(R) ≈ f atol=1e-12 rtol=1e-12
    end
    # A container with ℓₘₐₓ < |s| — the empty one included — holds no harmonics and cannot be
    # evaluated; the refusal is the flat `sYlm`'s, and is the same for either kind of index
    @test_throws "exceeds ℓₘₐₓ" ModeWeights(ComplexF64[], 1//2)(randn(rng, Rotor{Float64}))
    @test_throws "exceeds ℓₘₐₓ" ModeWeights(ComplexF64[], 1)(randn(rng, Rotor{Float64}))
    @test_throws "exceeds ℓₘₐₓ" ModeWeights(zeros(ComplexF64, 2), 3//2, 1//2, 1//2)(randn(rng, Rotor{Float64}))
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

            # The documented actions on mode weights, read through natural indexing: for ℓ ≥ |s|
            #   {L² f}ₗₘ = ℓ(ℓ+1) fₗₘ,   {Lz f}ₗₘ = m fₗₘ,   {Rz f}ₗₘ = s fₗₘ,
            #   {L₊ f}ₗₘ = √((ℓ+m)(ℓ-m+1)) fₗ,ₘ₋₁,   {L₋ f}ₗₘ = √((ℓ-m)(ℓ+m+1)) fₗ,ₘ₊₁,
            #   {ð f}ₗₘ = √((ℓ-s)(ℓ+s+1)) fₗₘ,       {ð̄ f}ₗₘ = -√((ℓ+s)(ℓ-s+1)) fₗₘ.
            # The eigenvalues ℓ(ℓ+1), m and s are formed as `Rational`s here and converted,
            # independently of the package's numerator arithmetic; every one of them is exactly
            # representable, so the diagonal relations are exact.
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
# queries, the conversions, the banded products and the unary operators — together with the
# two-argument constructor that deduces `ℓₘₐₓ` from the length of the data.

@testitem "ModeWeights: the rest of the array interface" begin
    using LinearAlgebra: Bidiagonal, Tridiagonal, Diagonal, dot
    import SphericalFunctions: ℓₘᵢₙ, ℓₘₐₓ, spin, half_integer
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

    # Equality and `dot` both work with a plain vector on either side.  `dot` conjugates —
    # it is deliberately *not* what evaluating the weights does.
    v = collect(w)
    @test isequal(v, w)
    @test isequal(w, v)
    @test v == w && w == v
    @test dot(v, w) == dot(v, v)
    @test dot(w, v) == dot(v, v)

    # Unary plus is the identity, and unary minus negates
    @test +w === w
    @test collect(-w) == -v
    @test spin(-w) == s && ℓₘᵢₙ(-w) == lo && ℓₘₐₓ(-w) == hi

    # The banded matrices the operators produce all multiply a `ModeWeights`
    d = randn(rng, ComplexF64, n)
    @test Diagonal(d) * w == Diagonal(d) * v
    B = Bidiagonal(d, randn(rng, ComplexF64, n-1), :U)
    @test B * w == B * v
    T3 = Tridiagonal(randn(rng, ComplexF64, n-1), d, randn(rng, ComplexF64, n-1))
    @test T3 * w == T3 * v
    # ... and multiplying on the other side is the outer product
    @test w * reshape(v, 1, n) == v * reshape(v, 1, n)

    # Broadcasting keeps the wrapper when the shape is unchanged, and the labels with it
    b = w .+ 1
    @test b isa ModeWeights
    @test spin(b) == s && ℓₘᵢₙ(b) == lo && ℓₘₐₓ(b) == hi
    @test collect(b) == v .+ 1
    @test (2 .* w) isa ModeWeights
    @test (w .+ w) isa ModeWeights
    @test collect(w .* 2) == v .* 2
end

@testitem "ModeWeights: `ℓₘₐₓ` deduced from the length of the data" begin
    import SphericalFunctions: ℓₘᵢₙ, ℓₘₐₓ, spin, half_integer
    using Random

    rng = Random.Xoshiro(11)

    # Integer indices: (ℓₘₐₓ+1)² = length + ℓₘᵢₙ²
    for s ∈ (-2, 0, 3), lo ∈ (abs(s), abs(s) + 1), hi ∈ (abs(s) + 2, abs(s) + 4)
        data = randn(rng, ComplexF64, Ysize(lo, hi))
        w = ModeWeights(data, s; ℓₘᵢₙ=lo)
        @test spin(w) == s
        @test ℓₘᵢₙ(w) == lo
        @test ℓₘₐₓ(w) == hi                     # deduced, not given
        @test w == ModeWeights(data, s, lo, hi)
    end

    # The default `ℓₘᵢₙ` is `abs(s)`
    w = ModeWeights(randn(rng, ComplexF64, Ysize(2, 5)), -2)
    @test ℓₘᵢₙ(w) == 2 && ℓₘₐₓ(w) == 5

    # Half-odd-integer indices use the doubled form, and the root must square back exactly
    for twos ∈ (-3, 1), twolo ∈ (abs(twos), abs(twos) + 2), twohi ∈ (abs(twos) + 2, abs(twos) + 6)
        s, lo, hi = half_integer(twos//2), half_integer(twolo//2), half_integer(twohi//2)
        data = randn(rng, ComplexF64, Ysize(lo, hi))
        w = ModeWeights(data, s; ℓₘᵢₙ=lo)
        @test spin(w) == s && ℓₘᵢₙ(w) == lo && ℓₘₐₓ(w) == hi
    end

    # A `Rational` spelling reaches the same place
    wr = ModeWeights(randn(rng, ComplexF64, Ysize(half_integer(1//2), half_integer(7//2))), 1//2)
    @test ℓₘᵢₙ(wr) == half_integer(1//2) && ℓₘₐₓ(wr) == half_integer(7//2)

    # A length that no ℓₘₐₓ can produce is refused rather than silently rounded
    @test_throws Exception ModeWeights(randn(rng, ComplexF64, 7), 0; ℓₘᵢₙ=0)
end
