# Tests of the places where an `@inbounds` access relies on an invariant established somewhere
# else — in a constructor, in a boundary method, or in a length check at the top of an
# operation.
#
# The kernels and the natural-index accessors of this package read and write their storage
# under `@inbounds`, and the accessors are `@propagate_inbounds`, so that a caller's own
# `@inbounds` removes their checks as well.  That is safe only if every value that reaches
# them has been validated where it entered: a block's limits must fit inside its storage, a
# rule's node count must give it a buffer to fill, and a vector that a container wraps without
# copying must still have the length its labels imply.  Each item below performs, as a user
# would, an operation whose input breaks one of these invariants, and expects the error that
# the boundary raises.  Without that error the access itself would be reached unchecked: an
# ordinary run would read or write outside the storage, which gives silently wrong numbers or
# corrupts memory, while a run with `--check-bounds=yes` would stop there with a
# `BoundsError`.  These items are meant to pass in both kinds of run, and they are tagged
# `:bounds` so that they can be run together, for example with
#
#     juliati . --check-bounds yes --filter ':bounds in tags'

@testitem "Bounds: the generic weight rules refuse too few nodes" tags=[:bounds] begin
    import SphericalFunctions: fejer2, clenshaw_curtis, salm2map, map2salm, Ysize
    import DoubleFloats: Double64

    # For a type that FFTW does not handle, Fejér's second rule fills a buffer of length n+1,
    # and the Clenshaw–Curtis rule one of length n-1, starting with the first element, so
    # n = -1 and n = 1 respectively would leave nothing to fill.
    @test_throws ArgumentError fejer2(-1, Double64)
    @test_throws ArgumentError clenshaw_curtis(1, Double64)

    # The same rule is what `salm2map` and `map2salm` build for a map with a single ring
    @test_throws ArgumentError salm2map(zeros(Complex{Double64}, Ysize(0, 3)), 0, 3, 7, 1)
    @test_throws ArgumentError map2salm(zeros(Complex{Double64}, 7, 1), 0, 3)
end

@testitem "Bounds: the machine-float weight rules refuse too few nodes" tags=[:bounds] begin
    import SphericalFunctions: fejer2, clenshaw_curtis, salm2map, Ysize

    # For the types that FFTW handles, the rules fill only the half spectrum that `irfft`
    # takes, of length (n+1)÷2+1 for Fejér's second rule and (n-1)÷2+1 for the Clenshaw–Curtis
    # rule.  Those lengths are zero for n ∈ {-3, -4} and n ∈ {-1, -2} respectively.
    for n ∈ (-3, -4)
        @test_throws ArgumentError fejer2(n)
        @test_throws ArgumentError fejer2(n, Float32)
    end
    for n ∈ (-1, -2)
        @test_throws ArgumentError clenshaw_curtis(n)
        @test_throws ArgumentError clenshaw_curtis(n, Float32)
    end

    # `salm2map` builds the Clenshaw–Curtis rule for its number of rings
    @test_throws ArgumentError salm2map(zeros(ComplexF64, Ysize(0, 3)), 0, 3, 7, -1)
    @test_throws ArgumentError salm2map(zeros(ComplexF64, Ysize(0, 3)), 0, 3, 7, -2)
end

@testitem "Bounds: the typed block constructors refuse storage smaller than the block" tags=[:bounds] begin
    import SphericalFunctions: WignerMatrix, WignerMatrixBatch, DegreeBlock, DegreeBlockBatch
    import SphericalFunctions: SpinMatrix, SpinMatrixBatch, HalfOddInteger

    # The natural-index accessors compare an index with the block's limits and then read the
    # storage under `@inbounds`, and iteration does the same for every element, so the limits
    # must fit inside the storage.  The typed constructors are the ones that `copy`, `similar`
    # and the views `w[iᵣ]`, `b[s, :]` and `b[iᵣ]` are built with.  Storage larger than the
    # block is legitimate, because a calculator's blocks sit in storage sized for its largest
    # ℓ; storage smaller than the block is not.
    @test_throws DimensionMismatch WignerMatrix{Int, Float64, Matrix{Float64}}(
        zeros(2, 2), 2, 2, -2, 2, -2
    )[2, 2]
    @test_throws DimensionMismatch WignerMatrixBatch{Int, Float64, Array{Float64, 3}}(
        zeros(1, 2, 2), 2, 2, -2, 2, -2, 3
    )[3, 2, 2]
    # ... including the extent of the rotor axis alone
    @test_throws DimensionMismatch WignerMatrixBatch{Int, Float64, Array{Float64, 3}}(
        zeros(1, 5, 5), 2, 2, -2, 2, -2, 2
    )[2, 0, 0]
    @test_throws DimensionMismatch DegreeBlock{Int, Float64, Vector{Float64}}(
        zeros(1), 3, 3, -3
    )[3]
    @test_throws DimensionMismatch DegreeBlockBatch{Int, Float64, Matrix{Float64}}(
        zeros(1, 1), 3, 3, -3, 2
    )[2, 3]
    @test_throws DimensionMismatch SpinMatrix{Int, Float64, Matrix{Float64}}(
        zeros(1, 1), 2, 1, -1, 2, -2
    )[1, 2]
    @test_throws DimensionMismatch SpinMatrixBatch{Int, Float64, Array{Float64, 3}}(
        zeros(1, 1, 1), 2, 1, -1, 2, -2, 2
    )[2, 1, 2]
    let h = HalfOddInteger
        @test_throws DimensionMismatch WignerMatrix{h, Float64, Matrix{Float64}}(
            zeros(1, 1), h(3//2), h(3//2), h(-3//2), h(3//2), h(-3//2)
        )[h(3//2), h(3//2)]
    end

    # Iteration, which is what `sum`, `maximum` and `collect` use, reads every element
    @test_throws DimensionMismatch sum(
        DegreeBlock{Int, Float64, Vector{Float64}}(zeros(1), 3, 3, -3)
    )

    # An axis of negative extent is refused, rather than compared with the storage, and so is
    # a negative number of rotors
    @test_throws ArgumentError WignerMatrix{Int, Float64, Matrix{Float64}}(
        zeros(5, 5), 2, -3, 2, 2, -2
    )
    @test_throws ArgumentError WignerMatrixBatch{Int, Float64, Array{Float64, 3}}(
        zeros(2, 5, 5), 2, 2, -2, 2, -2, -1
    )

    # Storage at least as large as the block is accepted, and the block keeps its own size
    @test size(WignerMatrix{Int, Float64, Matrix{Float64}}(zeros(7, 7), 2, 2, -2, 2, -2)) == (5, 5)
    @test size(
        WignerMatrixBatch{Int, Float64, Array{Float64, 3}}(zeros(3, 7, 7), 2, 2, -2, 2, -2, 2)
    ) == (2, 5, 5)
    @test size(DegreeBlock{Int, Float64, Vector{Float64}}(zeros(9), 3, 3, -3)) == (7,)
    @test size(DegreeBlockBatch{Int, Float64, Matrix{Float64}}(zeros(2, 9), 3, 3, -3, 2)) == (2, 7)
    @test size(SpinMatrix{Int, Float64, Matrix{Float64}}(zeros(4, 6), 2, 1, -1, 2, -2)) == (3, 5)
    @test size(
        SpinMatrixBatch{Int, Float64, Array{Float64, 3}}(zeros(2, 3, 5), 2, 1, -1, 2, -2, 2)
    ) == (2, 3, 5)
end

@testitem "Bounds: the dense reference recurrence refuses a restricted block" tags=[:bounds] begin
    import SphericalFunctions: WignerMatrix
    import SphericalFunctions:
        recurrence_step2!, recurrence_step3!, recurrence_step4!, recurrence_step5!,
        recurrence_step6!, convert_H_to_d!, convert_H_to_D!

    # These functions loop over every m of the block's ℓ, and take m′ₘᵢₙ to be -m′ₘₐₓ, under
    # `@inbounds`, so they need a block with the full range of m and a symmetric range of m′.
    # Each block here lies in the middle of a larger buffer, whose padding holds a value the
    # block does not, so that anything written outside the block can be seen.
    pad = 8
    function padded(::Type{NT}, ℓ; kwargs...) where {NT}
        n₁ = get(kwargs, :m′ₘₐₓ, ℓ) - get(kwargs, :m′ₘᵢₙ, -ℓ) + 1
        n₂ = get(kwargs, :mₘₐₓ, ℓ) - get(kwargs, :mₘᵢₙ, -ℓ) + 1
        buffer = fill(NT(7), n₁ + 2pad, n₂ + 2pad)
        block = view(buffer, pad+1:pad+n₁, pad+1:pad+n₂)
        fill!(block, one(NT))
        buffer, WignerMatrix(block, ℓ; kwargs...)
    end
    function padding_intact(buffer)
        outside = trues(size(buffer))
        outside[pad+1:end-pad, pad+1:end-pad] .= false
        all(==(7), buffer[outside])
    end
    sinβ, cosβ = sincos(0.3)

    # Steps 2 and 3 write one row of their first argument for every m ≥ 0 ...
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws ArgumentError recurrence_step2!(Hˡ, WignerMatrix(ones(5, 5), 2), sinβ, cosβ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws ArgumentError recurrence_step3!(Hˡ, WignerMatrix(ones(9, 9), 4), sinβ, cosβ)
        @test padding_intact(buffer)
    end
    # ... and read the m′ = 0 row of their second for every m ≥ 0 of its own ℓ
    let (buffer, Hˡ⁻¹) = padded(Float64, 2; mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws ArgumentError recurrence_step2!(WignerMatrix(zeros(7, 7), 3), Hˡ⁻¹, sinβ, cosβ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ⁺¹) = padded(Float64, 4; mₘₐₓ=2, mₘᵢₙ=-2)
        @test_throws ArgumentError recurrence_step3!(WignerMatrix(zeros(7, 7), 3), Hˡ⁺¹, sinβ, cosβ)
        @test padding_intact(buffer)
    end

    # Steps 4 and 5 run each row of m′ out to m = ℓ
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-2, mₘₐₓ=2, mₘᵢₙ=-2)
        @test_throws ArgumentError recurrence_step4!(Hˡ, sinβ, cosβ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=2, mₘᵢₙ=-2)
        @test_throws ArgumentError recurrence_step5!(Hˡ, sinβ, cosβ)
        @test padding_intact(buffer)
    end

    # Step 6 reflects through m → -m and m′ → -m′
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=3, mₘᵢₙ=-1)
        @test_throws ArgumentError recurrence_step6!(Hˡ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-1)
        @test_throws ArgumentError recurrence_step6!(Hˡ)
        @test padding_intact(buffer)
    end

    # The conversions multiply every element of the full block
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws ArgumentError convert_H_to_d!(Hˡ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-1)
        @test_throws ArgumentError convert_H_to_d!(Hˡ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(ComplexF64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws ArgumentError convert_H_to_D!(Hˡ, cis(0.2), cis(0.4))
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(ComplexF64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-1)
        @test_throws ArgumentError convert_H_to_D!(Hˡ, cis(0.2), cis(0.4))
        @test padding_intact(buffer)
    end
end

@testitem "Bounds: mul! with an operator refuses mode weights whose storage was resized" tags=[:bounds] begin
    import SphericalFunctions: ModeWeights, Δspin
    import SphericalFunctions: Lz, L₊, L₋, Lx
    using LinearAlgebra: mul!

    # A `ModeWeights` wraps its vector without copying it, so the vector can be resized after
    # the labels were checked against its length, while the operator kernels index up to the
    # length the labels imply, under `@inbounds`.  The operators here cover the four band
    # structures, whose kernels differ.
    for op ∈ (Lz, L₊, L₋, Lx)
        # An input shorter than its labels
        v = collect(1.0:16.0)
        w = ModeWeights(v, 0)
        resize!(v, 4)
        @test_throws DimensionMismatch mul!(
            ModeWeights(zeros(ComplexF64, 16), Δspin(op), 0, 3), op, w
        )
        @test_throws DimensionMismatch mul!(zeros(ComplexF64, 16), op, w)

        # A destination shorter than its labels
        u = zeros(ComplexF64, 16)
        w′ = ModeWeights(u, Δspin(op), 0, 3)
        resize!(u, 4)
        @test_throws DimensionMismatch mul!(w′, op, ModeWeights(collect(1.0:16.0), 0))
    end
end

@testitem "Bounds: an operator's product refuses mode weights whose storage was resized" tags=[:bounds] begin
    import SphericalFunctions: ModeWeights
    import SphericalFunctions: Lz, L₊, L₋, Lx

    # As for `mul!`, but here the result is allocated at the length of the input's storage and
    # then filled up to the length its labels imply, so the kernels would write past the end
    # of a fresh allocation.
    for op ∈ (Lz, L₊, L₋, Lx)
        v = collect(1.0:16.0)
        w = ModeWeights(v, 0)
        resize!(v, 4)
        @test_throws DimensionMismatch op * w
        @test_throws DimensionMismatch op(w)
    end
end

@testitem "Bounds: natural indexing refuses mode weights whose storage was resized" tags=[:bounds] begin
    import SphericalFunctions: ModeWeights, HalfOddInteger

    # `w[ℓ, m]` checks the mode against the labels and then reads the storage under
    # `@inbounds` at the position the labels give, so the storage must still have the length
    # the labels imply.
    v = collect(1.0:16.0)
    w = ModeWeights(v, 0)
    resize!(v, 4)
    @test_throws DimensionMismatch w[3, 0]
    @test_throws DimensionMismatch w[3, 3]
    @test_throws DimensionMismatch (w[3, 0] = 0.0)

    # ... for either kind of index
    vₕ = collect(1.0:20.0)
    wₕ = ModeWeights(vₕ, 1//2)
    resize!(vₕ, 4)
    @test_throws DimensionMismatch wₕ[7//2, 7//2]
    @test_throws DimensionMismatch wₕ[HalfOddInteger(7//2), HalfOddInteger(-1//2)]
    @test_throws DimensionMismatch (wₕ[7//2, 1//2] = 0.0)
end

@testitem "Bounds: HAxis refuses a largest ℓ below the smallest" tags=[:bounds] begin
    import SphericalFunctions: HAxis, HalfOddInteger

    # An axis starts at its smallest ℓ, and the natural-index accessors check an index only
    # against the current ℓ, so the storage must hold at least that one order.  With ℓₘₐₓ one
    # below the smallest ℓ it would hold nothing.
    @test_throws ArgumentError HAxis(Float64, 1, -1)[1, 0]
    @test_throws ArgumentError HAxis(Float64, 2, -1)[2, 0, 0]
    let h = HalfOddInteger
        @test_throws ArgumentError HAxis(Float64, 3, h(-1//2))[2, h(1//2)]
        @test_throws ArgumentError HAxis(Float64, 3, h(-1//2))[3, h(1//2), h(1//2)]
    end
end

@testitem "Bounds: a WignerSeries refuses to index blocks removed from it" tags=[:bounds] begin
    import SphericalFunctions: D
    using Quaternionic: Rotor

    # `values(s)` and `parent(s)` give the series' own vector of blocks, which can be resized,
    # while indexing checks ℓ against the labels of the series and then reads that vector
    # under `@inbounds` at the position the labels give.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    s = D(R, 3)
    pop!(values(s))
    @test_throws DimensionMismatch s[3]
    @test_throws DimensionMismatch last(s)

    # With the first block removed, every position would hold the block of the next ℓ
    s = D(R, 3)
    popfirst!(parent(s))
    @test_throws DimensionMismatch s[3]
    @test_throws DimensionMismatch s[0]
    @test_throws DimensionMismatch first(s)

    # Iteration pairs each ℓ with the block at the position its label gives, as indexing does
    @test_throws DimensionMismatch collect(s)
    @test_throws DimensionMismatch [ℓ for (ℓ, _) ∈ s]
end

@testitem "Bounds: a DegreeBlock refuses storage resized after its construction" tags=[:bounds] begin
    import SphericalFunctions: DegreeBlock, ModeWeights, relabel

    # A `DegreeBlock` may use a caller's vector as its storage, which can be resized after the
    # constructor has compared its length with the limits, while the accessors and iteration
    # read the storage under `@inbounds` at the positions the limits give.
    v = collect(1.0:5.0)
    b = DegreeBlock(v, 2)
    resize!(v, 1)
    @test_throws DimensionMismatch b[2]
    @test_throws DimensionMismatch (b[2] = 0.0)
    @test_throws DimensionMismatch sum(b)
    @test_throws DimensionMismatch collect(b)
    @test b[-2] == 1.0  # the one entry the storage still holds

    # ... including the block that `relabel` puts on a vector
    u = collect(1.0:7.0)
    b = relabel(ModeWeights(collect(1.0:16.0), 0)[3, :], u)
    resize!(u, 1)
    @test_throws DimensionMismatch b[3]
    @test_throws DimensionMismatch maximum(b)
end

@testitem "Bounds: the transforms refuse mode weights whose storage was resized" tags=[:bounds] begin
    import SphericalFunctions: SSHT, ModeWeights, nmodes, npixels
    using LinearAlgebra: mul!, ldiv!

    # A `ModeWeights` wraps its vector without copying it, so the transforms compare the
    # length of that vector with the labels before they read or fill it.
    for (method, inplace) ∈ (("RS", true), ("Minimal", true), ("Minimal", false),
                             ("Matrix", true), ("Matrix", false))
        𝒯 = method == "RS" ? SSHT(1, 4; method) : SSHT(1, 4; method, inplace)
        function resized_weights()
            v = zeros(ComplexF64, nmodes(𝒯))
            w = ModeWeights(v, 1, 1, 4)
            resize!(v, 3)
            w
        end
        f = zeros(ComplexF64, npixels(𝒯))
        @test_throws DimensionMismatch 𝒯 * resized_weights()
        @test_throws DimensionMismatch mul!(similar(f), 𝒯, resized_weights())
        @test_throws DimensionMismatch ldiv!(resized_weights(), 𝒯, f)
    end
end

@testitem "Bounds: HCalculator refuses rotor buffers shorter than its wedge" tags=[:bounds] begin
    import SphericalFunctions: HCalculator, HAxis, HWedge, recurrence!
    import SphericalFunctions: FixedSizeVector

    # The recurrence reads the phase e^{iβ} of every rotor of the wedge under `@inbounds`, so
    # the buffer of phases must be as long as the wedge has rotors.  The public constructors
    # size it from the rotor data, but the calculator's own constructor takes the buffers as
    # given.
    Hˡ = HWedge(Float64, 4, 2)
    h⃗ᵃ, h⃗ᵇ = HAxis(Float64, 4, 3), HAxis(Float64, 4, 3)
    h⃗ᵇ.ℓ = 1
    eⁱᵝ = FixedSizeVector{ComplexF64}(undef, 1)  # one phase for four rotors
    eⁱᵝ[1] = cis(0.3)
    none = FixedSizeVector{Float64}(undef, 0)  # the half angles, which integer indices do not use
    @test_throws DimensionMismatch recurrence!(
        HCalculator{Int, Float64, typeof(parent(Hˡ))}(
            h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, none, none, 2, 2, Ref(false), Ref(false)
        ),
        2
    )
end

@testitem "Bounds: the Wigner and harmonic calculators refuse buffers too small for them" tags=[:bounds] begin
    import SphericalFunctions: DCalculator, sYlmCalculator
    using Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(20260924)

    # `materialize!` writes the block and reads the power tables under `@inbounds`, for every
    # rotor, every (m′, m) or (s, m) the calculator serves, and every power up to 2ℓₘₐₓ, so the
    # calculators' own constructors compare the buffers they are given with all of those.
    # Each buffer is replaced here by one too small in a single dimension.
    R⃗ = randn(rng, Rotor{Float64}, 4)
    c = DCalculator(R⃗, 3)
    C = typeof(c)
    limits = (c.m′ₘₐₓ, c.m′ₘᵢₙ, c.mₘₐₓ, c.mₘᵢₙ)
    @test C(c.H, c.Wˡ, c.Z₊, c.Z₋, limits..., c.ℓ) isa C
    @test_throws DimensionMismatch C(c.H, zeros(ComplexF64, 1, 7, 7), c.Z₊, c.Z₋, limits..., c.ℓ)
    @test_throws DimensionMismatch C(c.H, zeros(ComplexF64, 4, 7, 6), c.Z₊, c.Z₋, limits..., c.ℓ)
    @test_throws DimensionMismatch C(c.H, c.Wˡ, zeros(ComplexF64, 2, 4), c.Z₋, limits..., c.ℓ)
    @test_throws DimensionMismatch C(c.H, c.Wˡ, c.Z₊, zeros(ComplexF64, 7, 1), limits..., c.ℓ)

    y = sYlmCalculator(R⃗, 3, -1:1)
    Y = typeof(y)
    state = (y.s, y.ℓ, y.phases)
    @test Y(y.H, y.Yˡ, y.Z₊, y.Z₋, state...) isa Y
    @test_throws DimensionMismatch Y(y.H, zeros(ComplexF64, 4, 1, 7), y.Z₊, y.Z₋, state...)
    @test_throws DimensionMismatch Y(y.H, zeros(ComplexF64, 4, 3, 5), y.Z₊, y.Z₋, state...)
    @test_throws DimensionMismatch Y(y.H, y.Yˡ, zeros(ComplexF64, 3, 4), y.Z₋, state...)
    @test_throws DimensionMismatch Y(y.H, y.Yˡ, y.Z₊, zeros(ComplexF64, 7, 2), state...)
end
