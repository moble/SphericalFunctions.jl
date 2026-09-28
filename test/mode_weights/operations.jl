# Tests of the container products in `src/mode_weights/operations.jl`:
#
#     𝔇 * w   rotates mode weights
#     Y * w   evaluates the function at the rotor(s)
#
# Both are bilinear — neither conjugates — which is the single most consequential detail
# here, and is what the "does not conjugate" item below pins against a future "fix".

@testitem "Rotating mode weights: the defining property" begin
    using Quaternionic: Rotor, RotorF64
    import DoubleFloats: Double64
    import SphericalFunctions: HalfOddInteger
    using Random

    rng = Random.Xoshiro(20260920)
    # `w(R)` is already pinned against the closed-form harmonics elsewhere, so it serves as
    # an independent oracle here: no new reference is needed.  The half-integer cases matter
    # because the group structure checked in the next item holds for the conjugate of 𝔇 as
    # well, so only this property distinguishes 𝔇 from it; for the same reason each case is
    # also compared with the rotation the other way, `w(R * Q)`, which must differ.
    h = HalfOddInteger
    cases = [
        [(s, ℓₘᵢₙ, 4) for s ∈ (-2, 0, 1) for ℓₘᵢₙ ∈ unique((abs(s), 0))];
        [(h(1//2), h(1//2), h(7//2)), (h(-3//2), h(3//2), h(9//2)), (h(3//2), h(1//2), h(3//2))]
    ]
    for T ∈ (Float32, Float64, Double64), (s, ℓₘᵢₙ, ℓₘₐₓ) ∈ cases
        w = ModeWeights(randn(rng, Complex{T}, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ)
        R = randn(rng, Rotor{T})
        # Measured: at most 15 eps(T) over thirty seeds of this sweep, against 150 to 450
        # eps(T) allowed, while the rotation the other way differs by at least 260ϵ (in
        # Float32, and by more than 10¹⁰ϵ in the other types)
        ϵ = 100ℓₘₐₓ * eps(T)
        rot = D(R, ℓₘₐₓ) * w
        wrong_way = zero(T)
        for Q ∈ randn(rng, Rotor{T}, 3)
            @test isapprox(rot(Q), w(inv(R) * Q); atol=ϵ, rtol=ϵ)
            wrong_way = max(wrong_way, abs(rot(Q) - w(R * Q)))
        end
        @test wrong_way > 100ϵ
        # the rotation changes neither the spin weight nor the range of ℓ, and computes in
        # the precision of the rotor and the weights
        @test rot isa ModeWeights{Complex{T}}
        @test spin(rot) === s
        @test SphericalFunctions.ℓₘᵢₙ(rot) === ℓₘᵢₙ && SphericalFunctions.ℓₘₐₓ(rot) === ℓₘₐₓ
    end
end

@testitem "Rotating mode weights: group structure" begin
    using Quaternionic: Rotor, RotorF64
    using LinearAlgebra: norm
    using Random

    rng = Random.Xoshiro(77)
    for ℓₘₐₓ ∈ (4, 7//2)
        # Measured: at most 17 eps over thirty seeds, against 350 or 400 eps allowed
        ϵ = 100ℓₘₐₓ * eps()
        s = ℓₘₐₓ isa Rational ? 1//2 : -2
        ℓₘᵢₙ = abs(s)
        w = ModeWeights(randn(rng, ComplexF64, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ)
        R₁ = randn(rng, RotorF64); R₂ = randn(rng, RotorF64)

        # Left multiplication by a representation composes with no transpose or inverse
        @test isapprox(array_view(D(R₁,ℓₘₐₓ) * (D(R₂,ℓₘₐₓ) * w)),
                       array_view(D(R₁*R₂, ℓₘₐₓ) * w); atol=ϵ, rtol=ϵ)
        @test isapprox(array_view(D(one(RotorF64), ℓₘₐₓ) * w), array_view(w); atol=ϵ, rtol=ϵ)
        @test isapprox(array_view(D(inv(R₁), ℓₘₐₓ) * (D(R₁, ℓₘₐₓ) * w)),
                       array_view(w); atol=ϵ, rtol=ϵ)
        # 𝔇 is unitary, so each ℓ block keeps its norm
        rot = D(R₁, ℓₘₐₓ) * w
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
            @test isapprox(norm(array_view(rot[ℓ, :])), norm(array_view(w[ℓ, :])); atol=ϵ, rtol=ϵ)
        end
        # The double cover: 𝔇(-R) = (-1)^{2ℓ} 𝔇(R), exactly
        sign = ℓₘₐₓ isa Rational ? -1 : 1
        @test array_view(D(-R₁, ℓₘₐₓ) * w) == sign .* array_view(D(R₁, ℓₘₐₓ) * w)
        # A rotation is block-diagonal in ℓ and touches nothing but m, so it commutes with ð
        @test isapprox(array_view(ð(D(R₁,ℓₘₐₓ) * w)), array_view(D(R₁,ℓₘₐₓ) * ð(w)); atol=ϵ, rtol=ϵ)
    end
end

@testitem "Rotating mode weights: refusals and the calculator" begin
    import SphericalFunctions: spin, Nᵣ
    using Quaternionic: Rotor, RotorF64
    using LinearAlgebra: mul!
    using Random

    rng = Random.Xoshiro(88)
    ℓₘₐₓ = 4
    w = ModeWeights(randn(rng, ComplexF64, Ysize(2, ℓₘₐₓ)), -2, 2, ℓₘₐₓ)
    R = randn(rng, RotorF64)

    # `D` starts at ℓ=0 while `w` starts at |s|, so containment is what is required
    @test D(R, ℓₘₐₓ) * w isa ModeWeights
    @test array_view(D(R, ℓₘₐₓ + 3) * w) == array_view(D(R, ℓₘₐₓ) * w)
    @test_throws "ℓ range of" D(R, 2) * w
    # A restricted block cannot rotate: every m mixes into every m′
    @test_throws "needs the whole" D(R, ℓₘₐₓ; m′ₘₐₓ=2) * w
    # Kinds may not be mixed
    wh = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), 1//2)
    @test_throws "must be of one kind" D(R, ℓₘₐₓ) * wh
    # A batched calculator rotates by one rotor, not many, and a batch of one is still a
    # batch, whose blocks have a rotor axis, so it is refused in the same way, with the
    # reason, rather than failing inside the product
    for calc ∈ (DCalculator(randn(rng, RotorF64, 3), ℓₘₐₓ), DCalculator([R], ℓₘₐₓ),
                dCalculator([0.3], ℓₘₐₓ))
        @test_throws ArgumentError calc * w
        @test_throws "built for a batch of $(Nᵣ(calc)) rotor" calc * w
        @test_throws "even of length one" mul!(similar(w), calc, w)
    end
    # A destination with other labels is refused, with the way to allocate one and the way
    # to copy weights into another range of ℓ
    @test_throws ArgumentError mul!(ModeWeights(zeros(ComplexF64, Ysize(0, ℓₘₐₓ)), -2, 0, ℓₘₐₓ), D(R, ℓₘₐₓ), w)
    @test_throws "`ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` copies" mul!(ModeWeights(zeros(ComplexF64, Ysize(0, ℓₘₐₓ)), -2, 0, ℓₘₐₓ), D(R, ℓₘₐₓ), w)

    # The calculator streams, and agrees with the series exactly: `D` copies the same blocks
    @test array_view(DCalculator(R, ℓₘₐₓ) * w) == array_view(D(R, ℓₘₐₓ) * w)
    # `mul!` writes into a correctly labelled destination ...
    dst = similar(w)
    @test array_view(mul!(dst, D(R, ℓₘₐₓ), w)) == array_view(D(R, ℓₘₐₓ) * w)
    @test array_view(mul!(similar(w), DCalculator(R, ℓₘₐₓ), w)) == array_view(dst)
    # ... but not into a mislabelled one, and not in place
    @test_throws "changes neither" mul!(ModeWeights(zeros(ComplexF64, Ysize(2, ℓₘₐₓ)), 1, 2, ℓₘₐₓ),
                                        D(R, ℓₘₐₓ), w)
    @test_throws "aliases the input" mul!(w, D(R, ℓₘₐₓ), w)
    # A bare vector, at least as long as the result, is accepted too, and the result comes
    # back labelled, as a ModeWeights over its first entries; a shorter one is refused
    for calc ∈ (D(R, ℓₘₐₓ), DCalculator(R, ℓₘₐₓ))
        out = zeros(ComplexF64, length(array_view(w)) + 2)
        w′ = mul!(out, calc, w)
        @test w′ isa ModeWeights && parent(array_view(w′)) === out
        @test (spin(w′), SphericalFunctions.ℓₘᵢₙ(w′), SphericalFunctions.ℓₘₐₓ(w′)) ==
            (spin(w), SphericalFunctions.ℓₘᵢₙ(w), SphericalFunctions.ℓₘₐₓ(w))
        @test array_view(w′) == array_view(D(R, ℓₘₐₓ) * w) && all(iszero, out[end-1:end])
        @test_throws "at least" mul!(zeros(ComplexF64, length(array_view(w)) - 1), calc, w)
    end
    @test_throws "aliases the input" mul!(parent(w), D(R, ℓₘₐₓ), w)
end

@testitem "Evaluating mode weights: the four shapes" begin
    using Quaternionic: Rotor, RotorF64
    import DoubleFloats: Double64
    using Random

    rng = Random.Xoshiro(123)
    s, ℓₘₐₓ = -2, 4
    w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
    R = randn(rng, RotorF64); R⃗ = randn(rng, RotorF64, 3)
    # Measured: at most 7 eps over thirty seeds, against 400 eps allowed
    ϵ = 100ℓₘₐₓ * eps()

    # One rotor, one spin: bit-exact against `w(R)`, which uses the same reduction
    @test sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w == w(R)
    # Many rotors
    @test isapprox(sYlm(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w, [w(r) for r ∈ R⃗]; atol=ϵ, rtol=ϵ)
    # A spin range selects `w`'s own row, in both the single and batched shapes
    @test isapprox(sYlm(R, ℓₘₐₓ, -2:2) * w, w(R); atol=ϵ, rtol=ϵ)
    @test isapprox(sYlm(R⃗, ℓₘₐₓ, -2:2) * w, [w(r) for r ∈ R⃗]; atol=ϵ, rtol=ϵ)
    # Harmonics wider in ℓ than the weights give the same answer
    @test isapprox(sYlm(R, ℓₘₐₓ + 2, s; ℓₘᵢₙ=abs(s)) * w, w(R); atol=ϵ, rtol=ϵ)
    # The calculator streams; its per-ℓ partial sums associate differently, hence ≈
    @test isapprox(sYlmCalculator(R, ℓₘₐₓ, s) * w, w(R); atol=ϵ, rtol=ϵ)
    @test isapprox(sYlmCalculator(R⃗, ℓₘₐₓ, s) * w, [w(r) for r ∈ R⃗]; atol=ϵ, rtol=ϵ)
    # And it is the labelled form of the product `sYlm_matrix` documents as
    # `f = Y * array_view(f̃)`, which is the same arithmetic (measured: exactly)
    @test sYlm_matrix(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * array_view(w) == sYlm(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w
    # ... while the product of the plain matrix with the labelled weights is refused, since it
    # could not check them
    @test_throws "sYlm(R⃗, ℓₘₐₓ, s) * w" sYlm_matrix(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w

    # Refusals
    # same ℓ range, so only the spin weight differs and only that check can fire
    w1 = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), 1, abs(s), ℓₘₐₓ)
    @test_throws ArgumentError sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w1
    @test_throws "spin weight s=1, but these harmonic values serve only s=-2." sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w1
    @test_throws "but these harmonic values serve only s ∈ -1:0." sYlm(R, ℓₘₐₓ, -1:0) * w1
    @test_throws "ℓ range of" sYlm(R, 3, s; ℓₘᵢₙ=abs(s)) * w

    # Weights stored from ℓ = 0 < |s| need harmonics only from |s| up, since the entries below
    # belong to no harmonic; those entries are not read by any of the shapes, so a
    # non-finite value there leaves every result exactly as it was
    w0 = ModeWeights(randn(rng, ComplexF64, Ysize(0, ℓₘₐₓ)), s, 0, ℓₘₐₓ)
    products = (
        () -> sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w0, () -> sYlm(R⃗, ℓₘₐₓ, s) * w0,
        () -> sYlm(R, ℓₘₐₓ, -2:2) * w0, () -> sYlm(R⃗, ℓₘₐₓ, -2:2) * w0,
        () -> sYlmCalculator(R, ℓₘₐₓ, s) * w0, () -> sYlmCalculator(R⃗, ℓₘₐₓ, s) * w0,
        () -> w0(R), () -> w0(R⃗),
    )
    before = map(f -> f(), products)
    @test isapprox(before[1], w0(R); atol=ϵ, rtol=ϵ)
    @test isapprox(before[2], [w0(r) for r ∈ R⃗]; atol=ϵ, rtol=ϵ)
    for ℓ ∈ 0:1, m ∈ -ℓ:ℓ
        w0[ℓ, m] = NaN
    end
    @test map(f -> f(), products) == before
    @test_throws "of which ℓ ∈ 2:4 are needed" sYlm(R, 3, s; ℓₘᵢₙ=abs(s)) * w0

    # The same shapes in other precisions, where the harmonics and the result are of the
    # precision of the rotor
    for T ∈ (Float32, Double64)
        wT = ModeWeights(randn(rng, Complex{T}, Ysize(abs(s), ℓₘₐₓ)), s)
        RT, R⃗T = randn(rng, Rotor{T}), randn(rng, Rotor{T}, 3)
        # Measured: at most 5.2 eps(T) over thirty seeds, against 400 eps(T) allowed
        ϵT = 100ℓₘₐₓ * eps(T)
        @test wT(RT) isa Complex{T}
        @test sYlm(RT, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * wT == wT(RT)
        @test isapprox(sYlm(R⃗T, ℓₘₐₓ, s) * wT, [wT(r) for r ∈ R⃗T]; atol=ϵT, rtol=ϵT)
        @test isapprox(sYlm(RT, ℓₘₐₓ, -2:2) * wT, wT(RT); atol=ϵT, rtol=ϵT)
        @test isapprox(sYlm(R⃗T, ℓₘₐₓ, -2:2) * wT, [wT(r) for r ∈ R⃗T]; atol=ϵT, rtol=ϵT)
        @test isapprox(sYlmCalculator(RT, ℓₘₐₓ, s) * wT, wT(RT); atol=ϵT, rtol=ϵT)
        @test isapprox(sYlmCalculator(R⃗T, ℓₘₐₓ, s) * wT, [wT(r) for r ∈ R⃗T]; atol=ϵT, rtol=ϵT)
    end
end

@testitem "Evaluating mode weights does not conjugate" begin
    using Quaternionic: Rotor, RotorF64
    using LinearAlgebra: dot
    using Random

    rng = Random.Xoshiro(4321)
    s, ℓₘₐₓ = -2, 3
    w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
    R = randn(rng, RotorF64)
    Y = sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

    # The bilinear product is the right one ...
    @test Y * w == w(R)
    # ... and the conjugating one is demonstrably different, by far more than any tolerance
    @test abs(dot(array_view(Y), array_view(w)) - w(R)) > 1e-6
    # ... so `dot` on these types errors rather than quietly answering with the wrong phase.
    # This is the tombstone: it stops a future maintainer "fixing" a MethodError by adding
    # the harmful non-conjugating method.  The message points to the labelled product.
    @test_throws "conjugates its first argument" dot(Y, w)
    @test_throws "Use `Y * w`, the bilinear product" dot(Y, w)
    @test_throws "conjugates its first argument" dot(w, Y)
    @test_throws "conjugates its first argument" dot(sYlmCalculator(R, ℓₘₐₓ, s), w)
    # The conjugating inner product of two sets of weights is still available and unchanged
    @test dot(w, w) ≈ sum(abs2, array_view(w))
end


@testitem "Mode-weight operations allocate only their result" begin
    import SphericalFunctions: ModeWeights, DCalculator, sYlmCalculator, D, sYlm, ð, array_view
    import LinearAlgebra: mul!
    import Quaternionic: Rotor
    import Random

    # Measured inside functions, as a user's inner loop sees it; at top level the boxing of
    # a dynamically dispatched call would be counted and would prove nothing.
    inplace(dst, A, src) = (mul!(dst, A, src); @allocated mul!(dst, A, src))
    product(A, src) = (A * src; @allocated A * src)

    rng = Random.Xoshiro(1234)
    R = randn(rng, Rotor{Float64})
    s, ℓₘₐₓ = -2, 8
    w = ModeWeights(randn(rng, ComplexF64, SphericalFunctions.Ysize(abs(s), ℓₘₐₓ)), s)

    # Rotation: `mul!` into a correctly labelled destination writes only into it.  The hint
    # of the message in `check_ℓ_covers` names `ℓₘₐₓ(w)`, so it is built only when the check
    # fails; built on every call, it would be some 350 bytes.
    w′ = similar(w)
    @test inplace(w′, D(R, ℓₘₐₓ), w) == 0
    @test inplace(w′, DCalculator(R, ℓₘₐₓ), w) == 0

    # Evaluation returns a scalar, so it has nothing to allocate at all
    @test product(sYlm(R, ℓₘₐₓ, s), w) == 0
    @test product(sYlmCalculator(R, ℓₘₐₓ, s), w) == 0

    # An operator applied in place, where the destination takes the *output* spin weight
    out = ModeWeights{ComplexF64}(undef, s + 1, abs(s), ℓₘₐₓ)
    @test inplace(out, ð, w) == 0
    @test out == ð * w
end

# Evaluation at many rotors at once, and the refusals on the evaluation side — the rotation
# side is covered above, but `check_evaluation` has its own ℓ-range and spin-range checks,
# and the kind mismatch is reported with the name of the kind the container actually holds.

@testitem "Evaluating mode weights: many rotors, and the refusals" begin
    using Quaternionic: Rotor, RotorF64
    import SphericalFunctions: HalfOddInteger
    using Random

    rng = Random.Xoshiro(2026)
    ℓₘₐₓ = 4
    s, ℓₘᵢₙ = -2, 2
    w = ModeWeights(randn(rng, ComplexF64, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ)

    # A vector of rotors evaluates at each, and agrees with evaluating one at a time
    R⃗ = randn(rng, RotorF64, 5)
    vals = w(R⃗)
    @test length(vals) == length(R⃗)
    for (i, R) ∈ enumerate(R⃗)
        @test vals[i] ≈ w(R)
    end
    # A one-element vector still gives a vector, not a scalar
    @test length(w(R⃗[1:1])) == 1

    R = R⃗[1]

    # A calculator whose ℓ range does not cover the weights says so, and names the fix
    @test_throws "ℓ range of" sYlmCalculator(R, 2, s) * w
    @test_throws "Build it with ℓₘₐₓ=4" sYlmCalculator(R, 2, s) * w

    # A calculator that does not serve this spin weight says which ones it does
    @test_throws ArgumentError sYlmCalculator(R, ℓₘₐₓ, 1) * w
    @test_throws "spin weight s=-2, but this calculator serves only s=1." sYlmCalculator(R, ℓₘₐₓ, 1) * w
    @test_throws "this calculator serves only s ∈ -1:1." sYlmCalculator(R, ℓₘₐₓ, -1:1) * w

    # Mixing the two kinds of index is refused, and the message names the kind the container
    # actually holds, which is what `index_kind_name` is for
    wh = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), 1//2)
    @test_throws "must be of one kind" sYlmCalculator(R, ℓₘₐₓ, s) * wh
    @test_throws "integers" sYlmCalculator(R, ℓₘₐₓ, s) * wh
    @test_throws "half-odd-integers" sYlmCalculator(R, HalfOddInteger(7//2), HalfOddInteger(1//2)) * w

    # The real harmonics are refused, whether flat or as a calculator, and for either kind
    # of index: they omit the phase i^{2s}, and depend on θ alone.  (For a half-odd spin
    # weight they would give f(θ, 0) times a constant ±i.)
    import SphericalFunctions: sλlm, sλlmCalculator
    @test_throws "real harmonics" sλlm(0.9, ℓₘₐₓ, s; ℓₘᵢₙ) * w
    @test_throws "real harmonics" sλlmCalculator(0.9, ℓₘₐₓ, s) * w
    @test_throws "real harmonics" sλlm(0.9, 7//2, 1//2) * wh
    @test_throws "real harmonics" sλlmCalculator(0.9, 7//2, 1//2) * wh
end

# The units composed as a user composes them: analysis, a differential operator, a change of
# the range of ℓ, and synthesis; the same transform applied to several functions at once,
# one column at a time; and a rotation followed by synthesis, analysis and evaluation.  Each
# result is compared with an evaluation of the weights at the transform's own rotors, which
# goes through none of the transform's code.

@testitem "Pipelines: analysis, operators and synthesis compose" begin
    import SphericalFunctions: ModeWeights, SSHT, D, ð, ð̄, R₊, spin, Ysize, rotors
    import LinearAlgebra: mul!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(20260930)

    # Analysis → operator → the range of ℓ of the new spin weight → synthesis, for the cases
    # that truncate the range (ð for s ≥ 0, ð̄ for s ≤ 0), that pad it with zeros (ð̄ for s
    # > 0, ð for s < 0), and that leave it alone (ð for s = -1/2).  Measured: at most 37 eps
    # relative to the largest value, against 100ℓₘₐₓ eps allowed.
    for method ∈ ("RS", "Matrix"), (s, op, Δ) ∈ (
        (1, ð, 1), (1, ð̄, -1), (0, R₊, 1), (-2, ð, 1), (-2, ð̄, -1),
        (1//2, ð, 1), (3//2, ð̄, -1), (-1//2, ð, 1),
    )
        L = s isa Rational ? 11//2 : 6
        ϵ = 100L * eps()
        𝒯, 𝒯′ = SSHT(s, L; method), SSHT(s + Δ, L; method)
        w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), L)), s)
        w₁ = 𝒯 \ (𝒯 * w)
        @test maximum(abs, array_view(w₁) - array_view(w)) ≤ ϵ * maximum(abs, array_view(w))
        ow = op * w₁
        @test spin(ow) == s + Δ && SphericalFunctions.ℓₘᵢₙ(ow) == abs(s)
        w₂ = ModeWeights(ow)
        @test SphericalFunctions.ℓₘᵢₙ(w₂) == abs(s + Δ)
        g = 𝒯′ * w₂
        @test maximum(abs, g - (op * w)(rotors(𝒯′))) ≤ ϵ * maximum(abs, g)
        # ... and the copy of the result into the range of the new spin weight has the
        # labels of the new transform's analysis, so the two can be compared and subtracted
        g̃ = 𝒯′ \ g
        @test g̃ ≈ w₂ rtol=ϵ
        @test (g̃ - w₂) isa ModeWeights
    end
end

@testitem "Pipelines: several functions at once, one column at a time" begin
    import SphericalFunctions: ModeWeights, SSHT, D, ð, Ysize
    import LinearAlgebra: mul!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(20260931)
    for method ∈ ("RS", "Matrix"), s ∈ (1, 3//2)
        L = s isa Rational ? 9//2 : 5
        𝒯 = SSHT(s, L; method)
        f = 𝒯 * ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), L)), s)
        # A multi-column analysis is a plain array, whose columns are labelled one at a time
        # without copying, so an operator or a rotation applies to each with no allocation,
        # and exactly as it applies to a copy of the column
        W = 𝒯 \ hcat(f, 2f, 3f)
        @test W isa Matrix{ComplexF64}
        @test_throws "ModeWeights(view(data, :, j), s)" ModeWeights(W, s)
        R = randn(rng, Rotor{Float64})
        W′, W″ = similar(W), similar(W)
        for j ∈ axes(W, 2)
            mul!(view(W′, :, j), ð, ModeWeights(view(W, :, j), s))
            mul!(view(W″, :, j), D(R, L), ModeWeights(view(W, :, j), s))
        end
        for j ∈ axes(W, 2)
            @test W′[:, j] == array_view(ð * ModeWeights(W[:, j], s))
            @test W″[:, j] == array_view(D(R, L) * ModeWeights(W[:, j], s))
        end
    end
end

@testitem "Pipelines: rotation, synthesis, analysis and evaluation" begin
    import SphericalFunctions: ModeWeights, SSHT, D, Ysize, rotors
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(20260932)
    # Measured: at most 34 eps against 100ℓₘₐₓ eps allowed
    for method ∈ ("RS", "Matrix"), s ∈ (1, -2, 1//2, -3//2)
        L = s isa Rational ? 9//2 : 5
        ϵ = 100L * eps()
        𝒯 = SSHT(s, L; method)
        w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), L)), s)
        R = randn(rng, Rotor{Float64})
        w′ = D(R, L) * w
        g = 𝒯 * w′
        @test maximum(abs, g - w′(rotors(𝒯))) ≤ ϵ * maximum(abs, g)
        w″ = 𝒯 \ g
        @test maximum(abs, array_view(w″) - array_view(w′)) ≤ ϵ * maximum(abs, array_view(w′))
        Q = randn(rng, Rotor{Float64}, 4)
        @test maximum(abs, w″(Q) - w(inv(R) .* Q)) ≤ ϵ * maximum(abs, w(inv(R) .* Q))
    end
end
