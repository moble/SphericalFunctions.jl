# Tests of `array_view` and `relabel` — the explicit route between the labelled containers
# and ordinary 1-based arrays, in `src/array_view.jl`.
#
# The reason this route exists at all is the first test item below.  Before version 3 the
# integer path returned `OffsetArray`s, and an `OffsetArray` with non-trivial offsets
# accepts `*` and `mul!` and returns *silently wrong* answers: a product of two blocks came
# back as a 1-based `Matrix` of mostly zeros, and an adjoint product came back holding
# uninitialized memory.  Refusing to be an `AbstractMatrix` turns that silence into a
# `MethodError`, and `array_view` is what a caller reaches for once they actually mean it.

@testitem "The composition law through array_view" begin
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(20260919)
    R₁ = randn(rng, Rotor{Float64})
    R₂ = randn(rng, Rotor{Float64})

    # 𝔇(R₁R₂) = 𝔇(R₁) 𝔇(R₂).  This is the property that the `OffsetArray` bug broke:
    # written as `𝔇₁[ℓ] * 𝔇₂[ℓ]` it used to give an answer wrong in the first digit, with
    # no error.
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
    using Random

    rng = Random.Xoshiro(33)
    ℓₘₐₓ = 5

    # An unbatched calculator's block has a unit leading stride, which is what BLAS
    # requires; contiguity is not required and is not present for ℓ < ℓₘₐₓ.
    calc = DCalculator(randn(rng, Rotor{Float64}), ℓₘₐₓ)
    for ℓ ∈ 0:ℓₘₐₓ
        A = array_view(recurrence!(calc, ℓ))
        @test A isa StridedArray
        @test stride(A, 1) == 1
    end

    # A batched block is also unit-strided as a whole, because the rotor axis leads ...
    N = 4
    rotors = randn(rng, Rotor{Float64}, N)
    batched = DCalculator(rotors, ℓₘₐₓ)
    blk = recurrence!(batched, ℓₘₐₓ)
    @test stride(array_view(blk), 1) == 1
    # ... but a single rotor's slice out of it is strided by Nᵣ, so BLAS cannot take it.
    # `mul!` then falls back to the generic implementation: slower, never wrong.
    one_rotor = array_view(blk[2])
    @test one_rotor isa StridedArray
    @test stride(one_rotor, 1) == N
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
    end
end

@testitem "array_view of a ModeWeights is its flat storage" begin
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
    # either a container or a plain array
    A = randn(rng, ComplexF64, 3, 4)
    @test array_view(A) === A
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
