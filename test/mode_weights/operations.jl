# Tests of the container products in `src/mode_weights/operations.jl`:
#
#     𝔇 * w   rotates mode weights
#     Y * w   evaluates the function at the rotor(s)
#
# Both are bilinear — neither conjugates — which is the single most consequential detail here,
# and is what the "does not conjugate" item below pins against a future "fix".

@testitem "Rotating mode weights: the defining property" begin
    using Quaternionic: Rotor, RotorF64
    using Random

    rng = Random.Xoshiro(20260920)
    # `w(R)` is already pinned against the closed-form harmonics elsewhere, so it serves as an
    # independent oracle here: no new reference is needed.
    for T ∈ (Float64,), s ∈ (-2, 0, 1), ℓₘᵢₙ ∈ (abs(s), 0)
        ℓₘₐₓ = 4
        w = ModeWeights(randn(rng, Complex{T}, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ)
        R = randn(rng, Rotor{T}); Q = randn(rng, Rotor{T})
        # measured worst case 3e-15 over this sweep; 1e-10 is a wide safety margin
        ϵ = 1e-10
        @test isapprox((D(R, ℓₘₐₓ) * w)(Q), w(inv(R) * Q); atol=ϵ, rtol=ϵ)
        # the rotation changes neither the spin weight nor the range of ℓ
        rot = D(R, ℓₘₐₓ) * w
        @test rot isa ModeWeights
        @test spin(rot) == s
        @test SphericalFunctions.ℓₘᵢₙ(rot) === ℓₘᵢₙ && SphericalFunctions.ℓₘₐₓ(rot) === ℓₘₐₓ
    end
end

@testitem "Rotating mode weights: group structure" begin
    using Quaternionic: Rotor, RotorF64
    using LinearAlgebra: norm
    using Random

    rng = Random.Xoshiro(77)
    ϵ = 1e-10
    for ℓₘₐₓ ∈ (4, 7//2)
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
    # A batched calculator rotates by one rotor, not many
    @test_throws "Nᵣ=" DCalculator(randn(rng, RotorF64, 3), ℓₘₐₓ) * w

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
end

@testitem "Evaluating mode weights: the four shapes" begin
    using Quaternionic: Rotor, RotorF64
    using Random

    rng = Random.Xoshiro(123)
    s, ℓₘₐₓ = -2, 4
    w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
    R = randn(rng, RotorF64); R⃗ = randn(rng, RotorF64, 3)
    ϵ = 1e-10

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
    # And it is the same product `sYlm_matrix` documents as `f = Y * f̃`
    @test isapprox(sYlm_matrix(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * array_view(w),
                   sYlm(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w; atol=ϵ, rtol=ϵ)

    # Refusals
    # same ℓ range, so only the spin weight differs and only that check can fire
    w1 = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), 1, abs(s), ℓₘₐₓ)
    @test_throws "spin weight" sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s)) * w1
    @test_throws "ℓ range of" sYlm(R, 3, s; ℓₘᵢₙ=abs(s)) * w
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
    # This is the tombstone: it stops a future maintainer "fixing" a MethodError by adding the
    # harmful non-conjugating method.
    @test_throws "conjugates its first argument" dot(Y, w)
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

    # Measured inside functions, as a user's inner loop sees it; at top level the boxing of a
    # dynamically dispatched call would be counted and would prove nothing.
    inplace(dst, A, src) = (mul!(dst, A, src); @allocated mul!(dst, A, src))
    product(A, src) = (A * src; @allocated A * src)

    rng = Random.Xoshiro(1234)
    R = randn(rng, Rotor{Float64})
    s, ℓₘₐₓ = -2, 8
    w = ModeWeights(randn(rng, ComplexF64, SphericalFunctions.Ysize(abs(s), ℓₘₐₓ)), s)

    # Rotation: `mul!` into a correctly labelled destination writes only into it.  These were
    # ~350 bytes each until the hint message in `check_ℓ_covers` was made lazy — it named
    # `ℓₘₐₓ(w)`, so it was built on every call whether or not the check failed.
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
    import SphericalFunctions: half_integer
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
    @test_throws "spin weight" sYlmCalculator(R, ℓₘₐₓ, 1) * w

    # Mixing the two kinds of index is refused, and the message names the kind the container
    # actually holds, which is what `index_kind_name` is for
    wh = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), 1//2)
    @test_throws "must be of one kind" sYlmCalculator(R, ℓₘₐₓ, s) * wh
    @test_throws "integers" sYlmCalculator(R, ℓₘₐₓ, s) * wh
    @test_throws "half-odd-integers" sYlmCalculator(R, half_integer(7//2), half_integer(1//2)) * w
end
