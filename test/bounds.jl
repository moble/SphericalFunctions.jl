# Tests of the places where an `@inbounds` access relies on an invariant established
# somewhere else — in a constructor, in a boundary method, or in a length check at the top
# of an operation.
#
# The kernels and the natural-index accessors of this package read and write their storage
# under `@inbounds`, and the accessors are `@propagate_inbounds`, so that a caller's own
# `@inbounds` removes their checks as well.  That is safe only if every value that reaches
# them has been validated where it entered: a block's limits must fit inside its storage, a
# rule's node count must give it a buffer to fill, and a vector that a container wraps
# without copying must still have the length its labels imply.  Each item below performs, as
# a user would, an operation whose input breaks one of these invariants, and expects the
# error that the boundary raises.  Without that error the access itself would be reached
# unchecked: an ordinary run would read or write outside the storage, which gives silently
# wrong numbers or corrupts memory, while a run with `--check-bounds=yes` would stop there
# with a `BoundsError`.  These items are meant to pass in both kinds of run, and they are
# tagged `:bounds` so that they can be run together, for example with
#
#     juliati . --check-bounds yes --filter ':bounds in tags'

@testitem "Bounds: the generic weight rules refuse too few nodes" tags=[:bounds] begin
    import SphericalFunctions: fejer2, clenshaw_curtis, salm2map, map2salm, Ysize
    import DoubleFloats: Double64

    # For a type that FFTW does not handle, the transform is GenericFFT's, but the rules
    # fill the same half spectrum as for the other types (see the next item), whose length is
    # zero for these values of n; they are refused before anything is filled.
    for n ∈ (-3, -4)
        @test_throws ArgumentError fejer2(n, Double64)
    end
    for n ∈ (-1, -2)
        @test_throws ArgumentError clenshaw_curtis(n, Double64)
    end

    # `salm2map` and `map2salm` refuse a map with fewer than two rings before they build the
    # Clenshaw–Curtis rule for it
    @test_throws "needs at least two rings" salm2map(zeros(Complex{Double64}, Ysize(0, 3)), 0, 3, 7, 1)
    @test_throws "needs at least two rings" map2salm(zeros(Complex{Double64}, 7, 1), 0, 3)
end

@testitem "Bounds: the machine-float weight rules refuse too few nodes" tags=[:bounds] begin
    import SphericalFunctions: fejer2, clenshaw_curtis, salm2map, Ysize

    # The rules fill only the half spectrum that `irfft` takes, of length (n+1)÷2+1 for
    # Fejér's second rule and (n-1)÷2+1 for the Clenshaw–Curtis rule.  Those lengths are
    # zero for n ∈ {-3, -4} and n ∈ {-1, -2} respectively.
    for n ∈ (-3, -4)
        @test_throws ArgumentError fejer2(n)
        @test_throws ArgumentError fejer2(n, Float32)
    end
    for n ∈ (-1, -2)
        @test_throws ArgumentError clenshaw_curtis(n)
        @test_throws ArgumentError clenshaw_curtis(n, Float32)
    end

    # `salm2map` refuses these numbers of rings before it builds the Clenshaw–Curtis rule
    @test_throws "needs at least two rings" salm2map(zeros(ComplexF64, Ysize(0, 3)), 0, 3, 7, -1)
    @test_throws "needs at least two rings" salm2map(zeros(ComplexF64, Ysize(0, 3)), 0, 3, 7, -2)
end

@testitem "Bounds: the typed block constructors refuse storage smaller than the block" tags=[:bounds] begin
    import SphericalFunctions: WignerMatrix, WignerMatrixBatch, DegreeBlock, DegreeBlockBatch
    import SphericalFunctions: SpinMatrix, SpinMatrixBatch, HalfOddInteger

    # The natural-index accessors compare an index with the block's limits and then read the
    # storage under `@inbounds`, and iteration does the same for every element, so the
    # limits must fit inside the storage.  The typed constructors are the ones that `copy`,
    # `similar` and the views `w[iᵣ]`, `b[s, :]` and `b[iᵣ]` are built with.  Storage larger
    # than the block is legitimate, because a calculator's blocks sit in storage sized for
    # its largest ℓ; storage smaller than the block is not.
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

@testitem "Bounds: the dense reference recurrence refuses a restricted block" setup=[DenseRecurrence] tags=[:bounds] begin
    import SphericalFunctions: WignerMatrix
    import .DenseRecurrence:
        recurrence_step2!, recurrence_step3!, recurrence_step4!, recurrence_step5!,
        recurrence_step6!, convert_H_to_d!, convert_H_to_D!

    # This item applies the same kind of check to the dense implementation of the recurrence
    # in the test module `DenseRecurrence` (in `test/wigner/recurrence.jl`), against which
    # the suite checks the engine.  No user calls these functions, but their `@inbounds`
    # loops rely on the block they are given just as the package's kernels do.
    #
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
        @test_throws "needs a block with the full range" recurrence_step2!(Hˡ, WignerMatrix(ones(5, 5), 2), sinβ, cosβ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws "needs a block with the full range" recurrence_step3!(Hˡ, WignerMatrix(ones(9, 9), 4), sinβ, cosβ)
        @test padding_intact(buffer)
    end
    # ... and read the m′ = 0 row of their second for every m ≥ 0 of its own ℓ
    let (buffer, Hˡ⁻¹) = padded(Float64, 2; mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws "needs a block with the full range" recurrence_step2!(WignerMatrix(zeros(7, 7), 3), Hˡ⁻¹, sinβ, cosβ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ⁺¹) = padded(Float64, 4; mₘₐₓ=2, mₘᵢₙ=-2)
        @test_throws "needs a block with the full range" recurrence_step3!(WignerMatrix(zeros(7, 7), 3), Hˡ⁺¹, sinβ, cosβ)
        @test padding_intact(buffer)
    end

    # Steps 4 and 5 run each row of m′ out to m = ℓ
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-2, mₘₐₓ=2, mₘᵢₙ=-2)
        @test_throws "needs a block with the full range" recurrence_step4!(Hˡ, sinβ, cosβ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=2, mₘᵢₙ=-2)
        @test_throws "needs a block with the full range" recurrence_step5!(Hˡ, sinβ, cosβ)
        @test padding_intact(buffer)
    end

    # Step 6 reflects through m → -m and m′ → -m′
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=3, mₘᵢₙ=-1)
        @test_throws "needs a block with the full range" recurrence_step6!(Hˡ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-1)
        @test_throws "needs a block with the full range" recurrence_step6!(Hˡ)
        @test padding_intact(buffer)
    end

    # The conversions multiply every element of the full block
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws "needs a block with the full range" convert_H_to_d!(Hˡ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(Float64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-1)
        @test_throws "needs a block with the full range" convert_H_to_d!(Hˡ)
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(ComplexF64, 3; m′ₘₐₓ=1, m′ₘᵢₙ=-1, mₘₐₓ=1, mₘᵢₙ=-1)
        @test_throws "needs a block with the full range" convert_H_to_D!(Hˡ, cis(0.2), cis(0.4))
        @test padding_intact(buffer)
    end
    let (buffer, Hˡ) = padded(ComplexF64, 3; m′ₘₐₓ=2, m′ₘᵢₙ=-1)
        @test_throws "needs a block with the full range" convert_H_to_D!(Hˡ, cis(0.2), cis(0.4))
        @test padding_intact(buffer)
    end
end

@testitem "Bounds: mul! with an operator refuses mode weights whose storage was resized" tags=[:bounds] begin
    import SphericalFunctions: ModeWeights, Δspin
    import SphericalFunctions: Lz, L₊, L₋, Lx
    using LinearAlgebra: mul!

    # A `ModeWeights` wraps its vector without copying it, so the vector can be resized
    # after the labels were checked against its length, while the operator kernels index up
    # to the length the labels imply, under `@inbounds`.  The operators here cover the four
    # band structures, whose kernels differ.
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

    # As for `mul!`, but here the result is allocated at the length of the input's storage
    # and then filled up to the length its labels imply, so the kernels would write past the
    # end of a fresh allocation.
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
    import SphericalFunctions: HAxis

    # An axis starts at order 0, whose elements the recurrence writes under `@inbounds`, so
    # the storage must hold at least that one order.  With ℓₘₐₓ one below 0 it would hold
    # nothing.
    @test_throws ArgumentError HAxis(Float64, 1, -1)
    @test_throws ArgumentError HAxis(Float64, 2, -1)
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
    d̄ₗ = FixedSizeVector{Float64}(undef, 2)  # the table of √δ², one entry per step of ℓ
    @test_throws DimensionMismatch recurrence!(
        HCalculator{Int, Float64, typeof(parent(Hˡ))}(
            h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, none, none, d̄ₗ, 2, 2, Ref(false), Ref(false)
        ),
        2
    )
end

@testitem "Bounds: the Wigner and harmonic calculators refuse buffers too small for them" tags=[:bounds] begin
    import SphericalFunctions: DCalculator, sYlmCalculator, SphericalFunctionsEngine
    using Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(20260924)

    # `materialize!` writes the block and reads the power tables under `@inbounds`, for
    # every rotor, every (m′, m) or (s, m) the calculator serves, and every power up to
    # 2ℓₘₐₓ, so the calculators' own constructors compare the buffers they are given with
    # all of those.  Each buffer is replaced here by one too small in a single dimension:
    # the power tables of the engine are laid out [iᵣ, k+1], with a row for each of the 4
    # rotors and 7 columns.  The calculator's copy of its rotors must likewise hold one for
    # each.
    R⃗ = randn(rng, Rotor{Float64}, 4)
    c = DCalculator(R⃗, 3)
    C = typeof(c)
    e = c.engine
    with_tables(e, Z₊, Z₋) = SphericalFunctionsEngine(e.H, Z₊, Z₋)
    limits = (c.m′ₘₐₓ, c.m′ₘᵢₙ, c.mₘₐₓ, c.mₘᵢₙ, c.m′ₘₐₓˢ, c.m′ₘᵢₙˢ, c.mₘₐₓˢ, c.mₘᵢₙˢ)
    @test C(e, c.Wˡ, c.rotors, limits..., c.ℓ, c.lift) isa C
    @test_throws DimensionMismatch C(e, zeros(ComplexF64, 1, 7, 7), c.rotors, limits..., c.ℓ, c.lift)
    @test_throws DimensionMismatch C(e, zeros(ComplexF64, 4, 7, 6), c.rotors, limits..., c.ℓ, c.lift)
    @test_throws DimensionMismatch C(with_tables(e, zeros(ComplexF64, 3, 7), e.Z₋), c.Wˡ, c.rotors, limits..., c.ℓ, c.lift)
    @test_throws DimensionMismatch C(with_tables(e, e.Z₊, zeros(ComplexF64, 4, 6)), c.Wˡ, c.rotors, limits..., c.ℓ, c.lift)
    @test_throws DimensionMismatch C(e, c.Wˡ, c.rotors[1:3], limits..., c.ℓ, c.lift)

    y = sYlmCalculator(R⃗, 3, -1:1)
    Y = typeof(y)
    e = y.engine
    state = (y.s, y.ℓ, y.phases, y.lift)
    @test Y(e, y.Yˡ, y.rotors, state...) isa Y
    @test_throws DimensionMismatch Y(e, zeros(ComplexF64, 4, 1, 7), y.rotors, state...)
    @test_throws DimensionMismatch Y(e, zeros(ComplexF64, 4, 3, 5), y.rotors, state...)
    @test_throws DimensionMismatch Y(with_tables(e, zeros(ComplexF64, 3, 7), e.Z₋), y.Yˡ, y.rotors, state...)
    @test_throws DimensionMismatch Y(with_tables(e, e.Z₊, zeros(ComplexF64, 4, 6)), y.Yˡ, y.rotors, state...)
    @test_throws DimensionMismatch Y(e, y.Yˡ, y.rotors[1:3], state...)
end


# The items above concern invariants of the storage; those below concern the type of an
# integer index, which must be `Int`.  The index arithmetic is done in the type of the
# indices: the flat layout puts the modes of ℓ from position ℓ² - ℓₘᵢₙ² + 1 onward, a block
# finds m at position m - mₘᵢₙ + 1 of its storage, and nearly every range of m runs from -ℓ to
# ℓ.  That arithmetic is closed under `Int` but not under the other integer types.  In `Int8`,
# ℓₘᵢₙ² overflows once ℓₘᵢₙ ≥ 12, and m - mₘᵢₙ once ℓ ≥ 64; in `Int16`, ℓₘᵢₙ² overflows once
# ℓₘᵢₙ ≥ 182; in an unsigned type, -ℓ and ℓₘᵢₙ - 1 wrap around to huge values; and a type
# wider than `Int` reaches code that stores `Int`s.  The boundary methods therefore refuse an
# integer index of any other type, a `Bool` included, with an `ArgumentError` that names the
# argument and says how to write it.  Each item below performs, as a user would, operations
# with indices of one of these types, and expects that refusal, which it recognizes by the
# sentence that explains it.  Without the refusal, some of these operations would compute a
# position outside the storage they index, so that an ordinary run would read outside it and
# a run with `--check-bounds=yes` would stop with a `BoundsError`; the others would return
# sizes or values that are silently wrong, or fail far from the cause.  These items are tagged
# `:narrow_integers` as well as `:bounds`.

@testsnippet IndexTypeRefusals begin
    # Each of these matches the `ArgumentError` with which a boundary method refuses an index
    # of one kind of wrong type, by the sentence that says what is wrong with it.  A function
    # given to `@test_throws` is applied to the error as `showerror` prints it, which begins
    # with the name of the exception's type.
    refusal(sentence) =
        message -> startswith(message, "ArgumentError") && occursin(sentence, message)
    narrow_refusal = refusal("narrower than `Int`")
    unsigned_refusal = refusal("is unsigned")
    wide_refusal = refusal("wider than `Int`")
    bool_refusal = refusal("A `Bool` is not an index")
end

@testitem "Bounds: products of harmonic values and mode weights refuse narrow labels" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: sYlm, ModeWeights, Ysize
    using Quaternionic: Rotor

    # The product reads the entries of the weights against a range of the storage of the
    # harmonic values, which starts at the position of the weights' smallest ℓ and is as long
    # as the weights.  In `Int8` with ℓₘᵢₙ = 12, the overflow of ℓₘᵢₙ² would put that position
    # at 257 and make both containers 308 entries long, rather than 1 and 52, so that the range
    # would run 256 entries past the end of the storage.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws narrow_refusal sYlm(R, Int8(13), Int8(12)) * ModeWeights(
        zeros(ComplexF64, Ysize(Int8(12), Int8(13))), Int8(12), Int8(12), Int8(13)
    )
    @test_throws narrow_refusal (
        sYlm(R, Int16(201), Int16(200))
        * ModeWeights{ComplexF64}(undef, Int16(200), Int16(200), Int16(201))
    )

    # ... for several rotors, or a range of spin weights, as well
    @test_throws narrow_refusal (
        sYlm([R, R], Int8(15), Int8(2); ℓₘᵢₙ=Int8(14))
        * ModeWeights{ComplexF64}(undef, Int8(2), Int8(14), Int8(15))
    )
    @test_throws narrow_refusal (
        sYlm(R, Int8(15), Int8(1):Int8(2); ℓₘᵢₙ=Int8(14))
        * ModeWeights{ComplexF64}(undef, Int8(2), Int8(14), Int8(15))
    )

    # With only the weights labelled in `Int8`, the range would be as long as their 308
    # entries, in harmonic values that hold 52
    @test_throws narrow_refusal (
        sYlm(R, 13, 12) * ModeWeights{ComplexF64}(undef, Int8(12), Int8(12), Int8(13))
    )
end

@testitem "Bounds: Yindex refuses narrow indices, which would point past the end of the modes" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: Yindex, Ysize

    # A vector that holds the modes of ℓ ∈ 12:13 has `Ysize(12, 13)` = 52 entries.  In `Int8`
    # the overflow of ℓₘᵢₙ² would put the first of them at position 257 and the last at 308,
    # and in `Int16`, with ℓ ∈ 200:201, the last of 804 at position 66340.
    f̃ = zeros(ComplexF64, Ysize(12, 13))
    @test_throws narrow_refusal f̃[Yindex(Int8(12), Int8(-12), Int8(12))]
    @test_throws narrow_refusal f̃[Yindex(Int8(13), Int8(13), Int8(12))]
    g̃ = zeros(ComplexF64, Ysize(200, 201))
    @test_throws narrow_refusal g̃[Yindex(Int16(201), Int16(201), Int16(200))]
end

@testitem "Bounds: ModeWeights refuses narrow labels whose layout would overflow" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: ModeWeights, Ysize

    # The constructors compare the length of the data with the `Ysize` of the labels, and the
    # `undef` form allocates that many entries.  With ℓₘᵢₙ² overflowing, data of the right
    # length would be refused and data of the wrong length accepted, and the form that deduces
    # ℓₘₐₓ from the length of the data would take the square root of a negative number.
    @test_throws narrow_refusal ModeWeights(zeros(Ysize(12, 13)), Int8(0), Int8(12), Int8(13))
    @test_throws narrow_refusal ModeWeights(
        zeros(Ysize(200, 201)), Int16(0), Int16(200), Int16(201)
    )
    @test_throws narrow_refusal ModeWeights{ComplexF64}(undef, Int8(0), Int8(12), Int8(20))
    @test_throws narrow_refusal ModeWeights(zeros(Ysize(12, 13)), Int8(0); ℓₘᵢₙ=Int8(12))
end

@testitem "Bounds: sYlm refuses a narrow ℓₘᵢₙ whose square would overflow" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: sYlm, sYlm_matrix
    using Quaternionic: Rotor

    # The harmonic values would be allocated with the overflowed `Ysize`, 316 entries for
    # ℓ ∈ 14:15 in `Int8` rather than 60, and each would be written at the overflowed
    # `Yindex`, which leaves the first 256 entries unwritten.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws narrow_refusal sYlm(R, Int8(15), Int8(2); ℓₘᵢₙ=Int8(14))
    @test_throws narrow_refusal sYlm(R, Int16(201), Int16(2); ℓₘᵢₙ=Int16(200))
    @test_throws narrow_refusal sYlm_matrix([R], Int8(15), Int8(2); ℓₘᵢₙ=Int8(14))
end

@testitem "Bounds: the harmonics refuse Int8 indices at ℓ ≥ 64" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: sYlm, sYlmCalculator, recurrence!
    using Quaternionic: Rotor

    # A block of ℓ finds m at position m - mₘᵢₙ + 1 of its storage, and reads or writes it
    # there under `@inbounds`.  In `Int8`, m - mₘᵢₙ exceeds 127 for some m once ℓ ≥ 64, and
    # wraps around to a negative position; the extent 2ℓ + 1 of the block wraps in the same
    # way.  The flat values would be wrong as well once ℓₘₐₓ + |s| ≥ 128, because the exponents
    # m ± s of their phases would wrap.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws narrow_refusal sYlm(R, Int8(100), Int8(0))[Int8(100)][Int8(50)]
    @test_throws narrow_refusal (sYlm(R, Int8(70), Int8(0))[Int8(70)][Int8(70)] = 0)
    @test_throws narrow_refusal (
        recurrence!(sYlmCalculator(R, Int8(127), Int8(-3)), Int8(64))[Int8(64)]
    )
    @test_throws narrow_refusal sYlm(R, Int8(126), Int8(2))
    @test_throws narrow_refusal sYlm(R, Int8(127), Int8(-3))
end

@testitem "Bounds: D and d refuse Int8 indices at ℓₘₐₓ ≥ 64" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: D, d, DCalculator, dCalculator, recurrence!
    using Quaternionic: Rotor

    # The calculators allocate m′ₘₐₓ - m′ₘᵢₙ + 1 rows and mₘₐₓ - mₘᵢₙ + 1 columns for their
    # blocks, and find (m′, m) at the offsets m′ - m′ₘᵢₙ and m - mₘᵢₙ.  In `Int8` each of these
    # wraps around once ℓₘₐₓ ≥ 64.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws narrow_refusal D(R, Int8(64))
    @test_throws narrow_refusal d(0.3, Int8(100))
    @test_throws narrow_refusal (
        recurrence!(DCalculator(R, Int8(100)), Int8(64))[Int8(64), Int8(-64)]
    )
    @test_throws narrow_refusal dCalculator(0.3, Int8(127))
end

@testitem "Bounds: blocks and mode weights refuse Int8 labels at ℓ ≥ 64" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: DegreeBlock, DegreeBlockBatch, WignerMatrix, ModeWeights, Ysize

    # As for the harmonics, the offset m - mₘᵢₙ of an element and the extent of an axis would
    # wrap around in `Int8` once ℓ ≥ 64.
    @test_throws narrow_refusal DegreeBlock(zeros(ComplexF64, 129), Int8(64))[Int8(64)]
    @test_throws narrow_refusal (DegreeBlock(zeros(ComplexF64, 129), Int8(64))[Int8(64)] = 1)
    @test_throws narrow_refusal sum(DegreeBlockBatch(zeros(ComplexF64, 2, 129), Int8(64)))
    @test_throws narrow_refusal (
        WignerMatrix(zeros(ComplexF64, 129, 129), Int8(64))[Int8(64), Int8(-64)]
    )
    @test_throws narrow_refusal let w = ModeWeights(
            zeros(ComplexF64, Ysize(0, 64)), Int8(0), Int8(0), Int8(64)
        )
        w[Int8(64), :][Int8(64)]
    end
end

@testitem "Bounds: the operators refuse narrow labels" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: ModeWeights, Ysize, L₊, L₋, Lx

    # The raising and lowering operators weight the mode (ℓ, m) with √((ℓ ∓ m)(ℓ ± m + 1)),
    # and in `Int8` the factor ℓ ∓ m wraps around to a negative number for some m once ℓ ≥ 64.
    # The diagonal of each band is as long as the `Ysize` of the labels, which overflows in
    # `Int8` once ℓₘᵢₙ ≥ 12, while the off-diagonals hold one entry for each pair of
    # neighbouring modes, so that their lengths would disagree.
    for op ∈ (L₊, L₋, Lx)
        @test_throws narrow_refusal op * ModeWeights(
            zeros(ComplexF64, Ysize(0, 100)), Int8(0), Int8(0), Int8(100)
        )
        @test_throws narrow_refusal op(Int8(0), Int8(0), Int8(100))
        @test_throws narrow_refusal op(Int8(0), Int8(12), Int8(13))
    end
end

@testitem "Bounds: Ysize and Yrange refuse unsigned indices" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: Ysize, Yrange

    # In an unsigned type the check ℓₘₐₓ ≥ ℓₘᵢₙ - 1 would wrap around at ℓₘᵢₙ = 0, and so
    # refuse every ℓₘₐₓ, and the range -ℓ:ℓ of m would be empty for ℓ ≥ 1, because -ℓ wraps
    # around to a huge value; `Yrange` would then fill only the first of the entries that it
    # allocates, and return the others unwritten.
    @test_throws unsigned_refusal Ysize(UInt(0), UInt(3))
    @test_throws unsigned_refusal Ysize(UInt(3))
    @test_throws unsigned_refusal Yrange(UInt(0), UInt(2))
    @test_throws unsigned_refusal Yrange(UInt(2))
end

@testitem "Bounds: ModeWeights refuses unsigned labels" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: ModeWeights, L²

    # With unsigned labels the range -ℓ:ℓ of m would be empty, so that every mode would be
    # reported out of bounds, every block would hold nothing, and the product with an operator
    # would give the wrong weights; and ℓₘᵢₙ = 0 would be refused by the wrapped check
    # ℓₘₐₓ ≥ ℓₘᵢₙ - 1.
    @test_throws unsigned_refusal ModeWeights(
        zeros(ComplexF64, 24), UInt(1), UInt(1), UInt(4)
    )[UInt(2), UInt(1)]
    @test_throws unsigned_refusal ModeWeights(
        zeros(ComplexF64, 24), UInt(1), UInt(1), UInt(4)
    )[UInt(2), :][UInt(1)]
    @test_throws unsigned_refusal L² * ModeWeights(
        collect(ComplexF64, 1:24), UInt(1), UInt(1), UInt(4)
    )
    @test_throws unsigned_refusal ModeWeights{Float64}(undef, 0, UInt(4))
    @test_throws unsigned_refusal ModeWeights(zeros(16), UInt(0))
end

@testitem "Bounds: the operators refuse unsigned indices" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: L², Lz, L₊, ð

    # The operators build their matrices mode by mode over the range -ℓ:ℓ of m, which would
    # be empty in an unsigned type, so that the matrix for ℓ ∈ 1:4 would have no rows at all.
    # Promoted to a common type with an unsigned ℓₘₐₓ, s = 0 and ℓₘᵢₙ = 0 would be unsigned as
    # well, and would be refused by the wrapped check ℓₘₐₓ ≥ ℓₘᵢₙ - 1.
    @test_throws unsigned_refusal L²(UInt(1), UInt(1), UInt(4))
    @test_throws unsigned_refusal ð(UInt(1), UInt(1), UInt(4))
    @test_throws unsigned_refusal L₊(UInt(1), UInt(1), UInt(4))
    @test_throws unsigned_refusal Lz(0, 0, UInt(2))
    @test_throws unsigned_refusal L²(0, UInt(2))
end

@testitem "Bounds: the Wigner and harmonic functions refuse unsigned indices" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: D, d, DCalculator, HCalculator, sYlm
    using Quaternionic: Rotor

    # The limits of m′ and m default to -ℓₘₐₓ, which wraps around in an unsigned type, and the
    # recurrences are written for signed indices, so these calls would fail inside the
    # validation of their ranges.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws unsigned_refusal D(R, UInt(3))
    @test_throws unsigned_refusal d(0.3, UInt(3))
    @test_throws unsigned_refusal DCalculator(R, UInt(3))
    @test_throws unsigned_refusal HCalculator(0.3, UInt(6))
    @test_throws unsigned_refusal sYlm(R, UInt(3), UInt(1))
    @test_throws unsigned_refusal sYlm(R, UInt(3), 1)
end

@testitem "Bounds: the pixelizations refuse narrow and unsigned indices" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
    import SphericalFunctions: leja_pixels, leja_rotors

    # A pixelization has `Ysize(|s|, ℓₘₐₓ)` points, one for each mode, so the overflow of s²
    # in `Int8` for |s| ≥ 12 would give 308 points for the 52 modes of s = 12 and ℓₘₐₓ = 13,
    # and the wrapped check of that count would refuse every unsigned ℓₘₐₓ at s = 0.
    @test_throws narrow_refusal golden_ratio_spiral_pixels(Int8(12), Int8(13))
    @test_throws narrow_refusal golden_ratio_spiral_rotors(Int8(12), Int8(13))
    @test_throws narrow_refusal leja_pixels(Int8(12), Int8(13))
    @test_throws narrow_refusal leja_rotors(Int8(12), Int8(13))
    @test_throws unsigned_refusal golden_ratio_spiral_pixels(UInt(0), UInt(3))
    @test_throws unsigned_refusal leja_pixels(UInt(0), UInt(3))
end

@testitem "Bounds: Bool indices are refused" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: ModeWeights, Ysize, L², D, sYlm
    using Quaternionic: Rotor

    # A `Bool` is an `Integer`, but -true is the `Int` -1, which a `Bool` label cannot hold, so
    # the block of ℓ = true could not be built, and the Wigner and harmonic functions would
    # fail inside the validation of their ranges.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws bool_refusal ModeWeights(zeros(4), false, false, true)[true, :]
    @test_throws bool_refusal L² * ModeWeights(zeros(4), false, false, true)
    @test_throws bool_refusal L²(false, false, true)
    @test_throws bool_refusal Ysize(false, true)
    @test_throws bool_refusal D(R, true)
    @test_throws bool_refusal sYlm(R, true, false)
end

@testitem "Bounds: calculators and harmonic values refuse a narrow index type" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: DCalculator, dCalculator, HCalculator, sYlmCalculator, sYlm, ℓ
    using Quaternionic: Rotor

    # Iteration passes ℓ + 1 on as the next state, which is an `Int` whatever the integer
    # type of ℓ, so a calculator or harmonic values of a narrower index type would stop with
    # a `MethodError` at the second step, and `ℓ` of an `HCalculator` would not be inferred to
    # a single type.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws narrow_refusal collect(DCalculator(R, Int8(3)))
    @test_throws narrow_refusal for (ℓ′, 𝔡ˡ) ∈ dCalculator(0.5, Int16(3)) end
    @test_throws narrow_refusal collect(dCalculator(0.5, Int32(3)))
    @test_throws narrow_refusal collect(sYlmCalculator(R, Int16(3), Int16(1)))
    @test_throws narrow_refusal [ℓ′ for (ℓ′, Yˡ) ∈ sYlm(R, Int32(3), Int32(1))]
    @test_throws narrow_refusal @inferred(ℓ(HCalculator(0.3, Int32(4))))
end

@testitem "Bounds: blocks refuse a narrow index type, which Int indices could not read" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: D, d, DCalculator, recurrence!, DegreeBlock, sYlm, ModeWeights
    using Quaternionic: Rotor

    # The accessors of a block take indices of its own index type only, while a caller writes
    # literals, which are `Int`s, and the iteration of a Wigner matrix or a spin matrix, which
    # every reduction uses, forms its indices as a lower limit plus an `Int`; with any other
    # integer index type, neither would find a method.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws narrow_refusal sum(D(R, Int32(2))[2])
    @test_throws narrow_refusal D(R, Int32(2))[2][1, 0]
    @test_throws narrow_refusal maximum(abs, d(0.3, Int32(2))[2])
    @test_throws narrow_refusal all(isfinite, recurrence!(DCalculator(R, Int32(3)), 2))
    @test_throws narrow_refusal DegreeBlock(zeros(5), Int32(2))[1]
    @test_throws narrow_refusal sYlm(R, Int32(3), Int32(0))[2][1]
    @test_throws narrow_refusal sum(sYlm(R, Int32(3), Int32(-1):Int32(1))[2])
    @test_throws narrow_refusal ModeWeights(
        collect(ComplexF64, 1:15), Int8(1), Int8(1), Int8(3)
    )[2, :][1]
end

@testitem "Bounds: indices wider than Int are refused" tags=[:bounds, :narrow_integers] setup=[IndexTypeRefusals] begin
    import SphericalFunctions: D, d, HCalculator, recurrence!, sYlm
    using Quaternionic: Rotor

    # The H recurrence computes the offsets of its wedge in the type of the indices, and the
    # wedge is indexed with `Int`s only, so an `Int128` or a `BigInt` index would reach it and
    # find no method.
    R = Rotor(1.0, 2.0, 3.0, 4.0)
    @test_throws wide_refusal D(R, Int128(3))
    @test_throws wide_refusal D(R, big(3))
    @test_throws wide_refusal d(0.3, Int128(10))
    @test_throws wide_refusal sYlm(R, big(3), big(1))
    @test_throws wide_refusal recurrence!(HCalculator(0.3, big(6)), big(6))
end
