# Tests of `array_view` and `relabel` — the explicit route between the labelled containers
# and ordinary 1-based arrays, in `src/containers/array_view.jl`.
#
# The reason this route exists at all is the first test item below.  A block with integer
# indices could be an `OffsetArray`, but an `OffsetArray` with non-trivial offsets accepts
# `*` and `mul!` and returns *silently wrong* answers: a product of two blocks comes back as
# a 1-based `Matrix` of mostly zeros, and an adjoint product comes back holding
# uninitialized memory.  Refusing to be an `AbstractMatrix` turns that silence into a
# `MethodError`, and `array_view` is what a caller reaches for once they actually mean it.

@testitem "The composition law through array_view" begin
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(20260919)
    R₁ = randn(rng, Rotor{Float64})
    R₂ = randn(rng, Rotor{Float64})

    # 𝔇(R₁R₂) = 𝔇(R₁) 𝔇(R₂).  This is the property that `OffsetArray` blocks would break:
    # written as `𝔇₁[ℓ] * 𝔇₂[ℓ]` with such blocks, the product is wrong in the first
    # digit, with no error.
    for ℓₘₐₓ ∈ (4, 7//2)
        𝔇₁ = D(R₁, ℓₘₐₓ)
        𝔇₂ = D(R₂, ℓₘₐₓ)
        𝔇₁₂ = D(R₁ * R₂, ℓₘₐₓ)
        # Both ends from the series itself: a `HalfOddInteger` and a `Rational` deliberately
        # do not promote, so `ℓₘᵢₙ(𝔇₁):7//2` would be an error rather than a range.
        for ℓ ∈ SphericalFunctions.ℓₘᵢₙ(𝔇₁):SphericalFunctions.ℓₘₐₓ(𝔇₁)
            product = array_view(𝔇₁[ℓ]) * array_view(𝔇₂[ℓ])
            @test product ≈ array_view(𝔇₁₂[ℓ]) atol=100eps(Float64)
            # ... and the labelled form of the same answer
            relabelled = relabel(𝔇₁₂[ℓ], product)
            @test axes(relabelled) == axes(𝔇₁₂[ℓ])
            @test SphericalFunctions.ℓ(relabelled) == ℓ
            for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                @test relabelled[m′, m] ≈ 𝔇₁₂[ℓ][m′, m] atol=100eps(Float64)
            end
        end
    end
end

@testitem "Containers refuse linear algebra without array_view" begin
    using Quaternionic: Rotor
    using LinearAlgebra: LinearAlgebra, mul!, lu
    using Random

    rng = Random.Xoshiro(11)
    R = randn(rng, Rotor{Float64})
    𝔇 = D(R, 3)
    A = 𝔇[2]

    # The containers are deliberately not `AbstractArray`s, so every one of these is a
    # `MethodError` rather than a wrong answer.
    @test !(A isa AbstractArray)
    @test_throws MethodError A * A
    @test_throws MethodError A'
    @test_throws MethodError lu(A)
    @test_throws MethodError mul!(similar(A), A, A)
    # Nor is there linear indexing to get wrong
    @test_throws MethodError A[1]

    # Going through `array_view` is what makes them work
    @test array_view(A) * array_view(A) isa Matrix{ComplexF64}
    @test array_view(A)' isa AbstractMatrix
    @test lu(array_view(A)) isa LinearAlgebra.LU
end

@testitem "array_view aliases, Matrix copies" begin
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(22)
    R = randn(rng, Rotor{Float64})

    for ℓₘₐₓ ∈ (3, 5//2)
        𝔇 = D(R, ℓₘₐₓ)
        ℓ = ℓₘₐₓ
        w = 𝔇[ℓ]
        A = array_view(w)
        M = Matrix(w)
        @test A == M                      # same values ...
        @test A[1, 1] === w[-ℓ, -ℓ]       # ... and `array_view` is 1-based over the same block
        # `array_view` aliases: writing through it writes into the container
        A[1, 1] = 17
        @test w[-ℓ, -ℓ] == 17
        # `Matrix` does not: it was a copy taken before the write
        @test M[1, 1] != 17
    end
end

@testitem "array_view strides: BLAS eligibility" begin
    using Quaternionic: Rotor
    using LinearAlgebra: mul!
    using Random

    rng = Random.Xoshiro(33)
    ℓₘₐₓ = 5

    # An unbatched calculator's block is the leading entries of the calculator's buffer, so
    # its `array_view` is contiguous and linearly indexed, with the unit leading stride that
    # BLAS requires, whatever ℓ is.
    calc = DCalculator(randn(rng, Rotor{Float64}), ℓₘₐₓ)
    for ℓ ∈ 0:ℓₘₐₓ
        blk = recurrence!(calc, ℓ)
        A = array_view(blk)
        @test A isa StridedArray
        @test strides(A) == (1, size(A, 1))
        @test IndexStyle(A) === IndexLinear()
        @test pointer(A) == pointer(parent(blk))
    end

    # A batched block is contiguous as a whole, because the rotor axis leads, so it can be
    # reshaped to a matrix of `Nᵣ` times as many rows for BLAS without a copy ...
    N = 4
    rotors = randn(rng, Rotor{Float64}, N)
    batched = DCalculator(rotors, ℓₘₐₓ)
    blk = recurrence!(batched, ℓₘₐₓ)
    n = 2ℓₘₐₓ + 1
    @test strides(array_view(blk)) == (1, N, N * n)
    @test IndexStyle(array_view(blk)) === IndexLinear()
    x = randn(rng, ComplexF64, n, 3)
    y = mul!(zeros(ComplexF64, N * n, 3), reshape(array_view(blk), :, n), x)
    for iᵣ ∈ 1:N
        @test y[iᵣ:N:end, :] ≈ array_view(blk[iᵣ]) * x
    end
    # ... but a single rotor's slice out of it is strided by Nᵣ, so BLAS cannot take it.
    # `mul!` then falls back to the generic implementation: slower, never wrong.
    one_rotor = array_view(blk[2])
    @test strides(one_rotor) == (N, N * n)
    @test IndexStyle(one_rotor) === IndexLinear()
    single = DCalculator(rotors[2], ℓₘₐₓ)
    reference = array_view(recurrence!(single, ℓₘₐₓ))
    @test one_rotor == reference
    # The generic fallback gives the same answer BLAS would
    @test one_rotor * one_rotor ≈ reference * reference
end

@testitem "relabel round-trips every container" begin
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(44)
    ℓₘₐₓ, N = 4, 3
    R = randn(rng, Rotor{Float64})
    Rs = randn(rng, Rotor{Float64}, N)

    # One block of each shape the package hands out
    containers = Any[]
    push!(containers, D(R, ℓₘₐₓ)[3])                                   # WignerMatrix
    push!(containers, recurrence!(DCalculator(Rs, ℓₘₐₓ), 3))           # WignerMatrixBatch
    push!(containers, recurrence!(sYlmCalculator(R, ℓₘₐₓ, -2), 3))     # DegreeBlock
    push!(containers, recurrence!(sYlmCalculator(Rs, ℓₘₐₓ, -2), 3))    # DegreeBlockBatch
    push!(containers, recurrence!(sYlmCalculator(R, ℓₘₐₓ, -2:2), 3))   # SpinMatrix
    push!(containers, recurrence!(sYlmCalculator(Rs, ℓₘₐₓ, -2:2), 3))  # SpinMatrixBatch

    for w ∈ containers
        A = collect(array_view(w))          # an independent copy, so the round trip is visible
        r = relabel(w, A)
        @test typeof(r).name === typeof(w).name
        @test axes(r) == axes(w)
        @test size(r) == size(w)
        @test SphericalFunctions.ℓ(r) == SphericalFunctions.ℓ(w)
        @test array_view(r) == array_view(w)
        @test r == w
        # The storage of the result is `vec(A)`, which shares the memory of `A`
        @test parent(r) == vec(A)
        A[end] += 1
        @test array_view(r)[end] == A[end]
    end
end

@testitem "array_view of a ModeWeights is its flat storage" begin
    import OffsetArrays: OffsetArray
    using Random
    rng = Random.Xoshiro(55)

    for s ∈ (0, -2, 1//2)
        ℓₘₐₓ = s isa Rational ? 7//2 : 4
        w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
        @test array_view(w) === parent(w)
        @test array_view(w) isa Vector{ComplexF64}
        # The flat form is what a product with a synthesis matrix takes, so it must stay in
        # the canonical ordering
        for ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ
            @test array_view(w)[Yindex(ℓ, m, abs(s))] == w[ℓ, m]
        end
    end

    # An ordinary array is already in that form, which is what lets the transforms take
    # either a container or a plain array; an array with other axes is refused, because
    # every caller indexes the result from 1
    A = randn(rng, ComplexF64, 3, 4)
    @test array_view(A) === A
    @test array_view(view(A, :, 2:3)) == A[:, 2:3]
    for offset ∈ (OffsetArray(copy(A), 0:2, 1:4), OffsetArray(randn(rng, ComplexF64, 25), 0:24))
        @test_throws ArgumentError array_view(offset)
        @test_throws "offset arrays are not supported" array_view(offset)
    end
    # ... which the transforms, which take their input through `array_view`, inherit
    for method ∈ ("RS", "Minimal", "Matrix")
        𝒯 = SSHT(0, 4; method)
        @test_throws ArgumentError 𝒯 * OffsetArray(zeros(ComplexF64, 25), 0:24)
        @test_throws "offset arrays are not supported" 𝒯 * OffsetArray(zeros(ComplexF64, 25), 0:24)
        f = 𝒯 * zeros(ComplexF64, 25)
        @test_throws "offset arrays are not supported" 𝒯 \ OffsetArray(f, 0:length(f)-1)
    end
end

@testitem "Broadcast assignment writes through a container" begin
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(66)
    R = randn(rng, Rotor{Float64})
    c = sYlmCalculator(R, 4, -2)
    v = recurrence!(c, 3)

    v .= 5 + 0im
    @test all(v[m] == 5 for m ∈ -3:3)
    @test all(array_view(v) .== 5)

    # And a `ModeWeights` row view, which is the same container
    w = ModeWeights(zeros(ComplexF64, Ysize(0, 3)), 0)
    row = w[2, :]
    row .= 7 + 0im
    @test all(w[2, m] == 7 for m ∈ -2:2)
    @test w[1, 0] == 0    # other blocks untouched
end

@testitem "relabel refuses an array whose shape is not the block's" begin
    import SphericalFunctions: WignerMatrix, WignerMatrixBatch, DegreeBlock, DegreeBlockBatch
    import SphericalFunctions: SpinMatrix, SpinMatrixBatch, relabel, array_view
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(45)
    ℓₘₐₓ, N = 4, 3
    R = randn(rng, Rotor{Float64})
    Rs = randn(rng, Rotor{Float64}, N)

    # One block of each shape, built by hand and taken from a calculator, whose blocks are
    # the leading entries of a buffer sized for its largest ℓ
    containers = Any[
        WignerMatrix(zeros(ComplexF64, 3, 3), 1),
        WignerMatrixBatch(zeros(ComplexF64, 2, 3, 3), 1),
        DegreeBlock(zeros(ComplexF64, 3), 1),
        DegreeBlockBatch(zeros(ComplexF64, 2, 3), 1),
        SpinMatrix(zeros(ComplexF64, 3, 3), 1; sₘₐₓ=1, sₘᵢₙ=-1),
        SpinMatrixBatch(zeros(ComplexF64, 2, 3, 3), 1; sₘₐₓ=1, sₘᵢₙ=-1),
        recurrence!(DCalculator(R, ℓₘₐₓ), 2),
        recurrence!(DCalculator(Rs, ℓₘₐₓ), 2),
        recurrence!(sYlmCalculator(R, ℓₘₐₓ, -2), 2),
        recurrence!(sYlmCalculator(Rs, ℓₘₐₓ, -2), 2),
        recurrence!(sYlmCalculator(R, ℓₘₐₓ, -2:2), 2),
        recurrence!(sYlmCalculator(Rs, ℓₘₐₓ, -2:2), 2),
    ]

    # An array larger than the block, in any one dimension or in all of them, is refused
    # rather than covered in its first entries
    for w ∈ containers
        n = size(w)
        for d ∈ eachindex(n)
            larger = ntuple(i -> n[i] + (i == d), length(n))
            @test_throws DimensionMismatch relabel(w, zeros(ComplexF64, larger))
        end
        @test_throws DimensionMismatch relabel(w, zeros(ComplexF64, n .+ 2))
        # ... and so is a vector longer than the block, which its constructor would accept
        @test_throws DimensionMismatch relabel(w, zeros(ComplexF64, prod(n) + 5))

        # The block's own shape is accepted, and used as the storage
        A = randn(rng, ComplexF64, n)
        r = relabel(w, A)
        @test size(r) == n
        @test array_view(r) == A
    end
end

@testitem "Calculator blocks are contiguous views of one buffer" begin
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator, sλlmCalculator,
        recurrence!, array_view, WignerMatrix, WignerMatrixBatch, DegreeBlock,
        DegreeBlockBatch, SpinMatrix, SpinMatrixBatch, HalfOddInteger
    using Quaternionic: Rotor
    import Random

    # Every calculator keeps its block in the leading entries of one `Vector`, which is the
    # storage of the block it returns at every ℓ, so the block's type does not depend on ℓ
    # or on the limits, its `array_view` is contiguous and linearly indexed, and a sweep
    # allocates nothing.
    rng = Random.Xoshiro(20261001)
    R⃗ = randn(rng, Rotor{Float64}, 3)
    θ⃗ = [0.2, 1.1, 2.9]
    sweep(c) = (for (_, b) ∈ c; end; nothing)
    allocations(c) = (sweep(c); @allocated sweep(c))
    for (ℓmax, s, srange) ∈ ((5, -2, -2:2), (9//2, 1//2, -3//2:3//2))
        IT = ℓmax isa Integer ? Int : HalfOddInteger
        for (calc, Block, buffer) ∈ (
            (DCalculator(R⃗[1], ℓmax), WignerMatrix, :Wˡ),
            (DCalculator(R⃗, ℓmax), WignerMatrixBatch, :Wˡ),
            (dCalculator(θ⃗[1], ℓmax), WignerMatrix, :Wˡ),
            (dCalculator(θ⃗, ℓmax), WignerMatrixBatch, :Wˡ),
            (sYlmCalculator(R⃗[1], ℓmax, s), DegreeBlock, :Yˡ),
            (sYlmCalculator(R⃗, ℓmax, s), DegreeBlockBatch, :Yˡ),
            (sYlmCalculator(R⃗[1], ℓmax, srange), SpinMatrix, :Yˡ),
            (sYlmCalculator(R⃗, ℓmax, srange), SpinMatrixBatch, :Yˡ),
            (sλlmCalculator(θ⃗[1], ℓmax, s), DegreeBlock, :Yˡ),
            (sλlmCalculator(θ⃗, ℓmax, srange), SpinMatrixBatch, :Yˡ),
        )
            NT = eltype(getfield(calc, buffer))
            @test isconcretetype(eltype(calc))
            @test eltype(calc) === Pair{IT, Block{IT, NT, Vector{NT}}}
            for (ℓ, b) ∈ calc
                @test b isa Block{IT, NT, Vector{NT}}
                @test parent(b) === getfield(calc, buffer)
                A = array_view(b)
                @test A isa StridedArray && IndexStyle(A) === IndexLinear()
                @test strides(A) == Base.size_to_strides(1, size(A)...)
                @test pointer(A) == pointer(parent(b))
            end
            @test allocations(calc) == 0
        end
    end
end

@testitem "Slices of a block are views of its storage" begin
    import SphericalFunctions: DCalculator, sYlmCalculator, recurrence!, array_view,
        WignerMatrix, DegreeBlock, DegreeBlockBatch, SpinMatrix
    using Quaternionic: Rotor
    import Random

    # The slices along the leading axis — `w[iᵣ]` of a batch, `v[iᵣ]` and `b[iᵣ]`, and the
    # row `b[s, :]` of a `SpinMatrix` — take as storage the arithmetic progression of the
    # block's storage that holds their elements, as a strided view of the one-dimensional
    # storage.  The slice `b[:, s, :]` of a `SpinMatrixBatch` is not such a progression, and
    # its storage is the strided matrix of exactly its shape, which BLAS takes.
    rng = Random.Xoshiro(20261001)
    N = 3
    R⃗ = randn(rng, Rotor{Float64}, N)
    is_progression(p, b) =
        p isa SubArray && parent(p) === parent(b) && only(p.indices) isa StepRange
    for (ℓmax, s, srange) ∈ ((5, -2, -2:2), (9//2, 1//2, -3//2:3//2))
        ℓ = ℓmax - 1
        nₛ = length(srange)

        w = recurrence!(DCalculator(R⃗, ℓmax), ℓ)
        W = Array(w)
        for iᵣ ∈ 1:N
            @test w[iᵣ] isa WignerMatrix && is_progression(parent(w[iᵣ]), w)
            @test Array(w[iᵣ]) == W[iᵣ, :, :]
        end

        v = recurrence!(sYlmCalculator(R⃗, ℓmax, s), ℓ)
        V = Array(v)
        for iᵣ ∈ 1:N
            @test v[iᵣ] isa DegreeBlock && is_progression(parent(v[iᵣ]), v)
            @test Array(v[iᵣ]) == V[iᵣ, :]
        end

        b = recurrence!(sYlmCalculator(R⃗, ℓmax, srange), ℓ)
        B = Array(b)
        for iᵣ ∈ 1:N
            @test b[iᵣ] isa SpinMatrix && is_progression(parent(b[iᵣ]), b)
            @test Array(b[iᵣ]) == B[iᵣ, :, :]
        end
        for (j, sⱼ) ∈ enumerate(srange)
            slice = b[:, sⱼ, :]
            @test slice isa DegreeBlockBatch
            @test array_view(slice) isa StridedMatrix
            @test strides(array_view(slice)) == (1, N * nₛ)
            @test Array(slice) == B[:, j, :]
            for iᵣ ∈ 1:N
                # The two orders of slicing agree
                @test Array(b[iᵣ][sⱼ, :]) == Array(slice[iᵣ]) == B[iᵣ, j, :]
            end
        end

        b₁ = recurrence!(sYlmCalculator(R⃗[1], ℓmax, srange), ℓ)
        B₁ = Array(b₁)
        for (j, sⱼ) ∈ enumerate(srange)
            @test b₁[sⱼ, :] isa DegreeBlock && is_progression(parent(b₁[sⱼ, :]), b₁)
            @test Array(b₁[sⱼ, :]) == B₁[j, :]
        end
    end
end
