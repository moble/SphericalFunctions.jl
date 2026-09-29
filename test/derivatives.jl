# Tests of the rules for automatic differentiation that the extensions attach to `D_array`
# and `sYlm_array` (see `src/derivatives.jl`).  The references are the explicit polynomial of
# the `ExplicitWignerMatrices` module, which shares no code with the package, evaluated in
# `BigFloat` and differentiated by ForwardDiff.  The polynomial is smooth everywhere, so its
# derivatives are good references at the poles too.

@testsnippet DerivativeTools begin
    import ForwardDiff
    using Quaternionic: Quaternion, Rotor, from_euler_angles
    import SphericalFunctions
    import SphericalFunctions: D, sYlm, array_view, pole_radius

    # A rotor with the components of `v`, not normalized, so that derivatives are taken in all
    # four directions of the space of quaternions, including the one that changes only the
    # norm.  Normalizing with `Rotor(v...)` would be differentiable too, but it goes through
    # Quaternionic's `abs`, which some backends cannot follow; see the note on the backends in
    # "Derivatives: first order, with every backend".
    rotorof(v) = Rotor{eltype(v)}(v[1], v[2], v[3], v[4])

    # Complex values as one real vector, of their real parts followed by their imaginary parts
    realified(z) = [real.(vec(z)); imag.(vec(z))]

    # The values of `D` and of `sYlm` as real vectors, from the package and from the
    # polynomial, in the order in which the package stores them.  These are callable structs
    # rather than closures, so that Enzyme compiles each kind once for all of the points at
    # which it is tested.
    struct DValues{T, L}
        n::T
        limits::L  # (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
    DValues(n) = DValues(n, (n, -n, n, -n))
    function (f::DValues)(v)
        m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ = f.limits
        S = D(rotorof(v), f.n; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        realified(reduce(vcat, [vec(parent(b)) for b ∈ values(S)]))
    end
    function polynomial(f::DValues, v)
        m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ = f.limits
        realified([
            ExplicitWignerMatrices.D_polynomial(ℓ, m′, m, v)[1]
            for ℓ ∈ (f.n isa Integer ? 0 : 1//2):f.n
            for m ∈ max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ) for m′ ∈ max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)
        ])
    end
    struct YValues{T, S}
        n::T
        s::S
        ℓₘᵢₙ::T
    end
    (f::YValues)(v) = realified(array_view(sYlm(rotorof(v), f.n, f.s; ℓₘᵢₙ=f.ℓₘᵢₙ)))
    function polynomial(f::YValues, v)
        realified([
            abs(s) ≤ ℓ ? ExplicitWignerMatrices.sYlm_polynomial(ℓ, m, s, v)[1] : zero(Complex{eltype(v)})
            for ℓ ∈ f.ℓₘᵢₙ:f.n for m ∈ -ℓ:ℓ for s ∈ (f.s isa AbstractRange ? f.s : (f.s,))
        ])
    end
    Base.show(io::IO, f::DValues) = print(io, "D(R, ", f.n, "; limits=", f.limits, ")")
    Base.show(io::IO, f::YValues) = print(io, "sYlm(R, ", f.n, ", ", f.s, "; ℓₘᵢₙ=", f.ℓₘᵢₙ, ")")

    reference_jacobian(f, v) = Float64.(ForwardDiff.jacobian(w -> polynomial(f, w), big.(v)))

    # The functions tested, with integer and half-integer indices, whole and restricted blocks
    # of 𝔇, and one spin weight or a range of them.
    functions = [
        DValues(3),
        DValues(4, (2, -1, 3, -2)),
        DValues(5//2),
        DValues(7//2, (3//2, -1//2, 7//2, -5//2)),
        YValues(4, -1, 1),
        YValues(4, -2:1, 0),
        YValues(7//2, 1//2, 1//2),
    ]
    largest_ℓ(f) = Float64(f.n)

    # Rotors at the poles, near them, and away from both, the last not normalized.  One is
    # just outside the radius within which the calculators use the expansion about a pole in
    # place of the recurrence.  There, for the functions with integer indices, the
    # recurrence's first derivatives are wrong by 70 to 180 times `eps()` relative to the
    # largest, several times the tolerance used below, so that a backend that bypassed the
    # rules would fail.  (With half-integer indices the recurrence happens to be accurate
    # there.)
    function test_points(ℓₘₐₓ)
        components(R) = [R[1], R[2], R[3], R[4]]
        near(r) = components(from_euler_angles(0.4, 2asin(r), -1.3))
        [
            "the identity" => [1.0, 0.0, 0.0, 0.0],
            "a rotation about z" => [cos(0.35), 0.0, 0.0, sin(0.35)],
            "a rotor at β = π" => components(
                Quaternion(0.0, 0.0, 1.0, 0.0) * Quaternion(cos(0.35), 0.0, 0.0, sin(0.35))
            ),
            "a rotor near β = 0" => near(1e-7),
            "a rotor just outside the expansion" => near(1.05 * pole_radius(Float64, ℓₘₐₓ)),
            "a generic rotor" => 1.7 .* [0.3, -0.5, 0.7, 0.2],
        ]
    end

    # Agreement to within `rtol` times the largest element of the reference, in the maximum
    # norm.  DifferentiationInterfaceTest calls this with its own keywords.
    function within(a, b; atol=0, rtol)
        maximum(abs, a - b; init=zero(real(eltype(b)))) ≤
            atol + rtol * max(1, maximum(abs, b; init=zero(real(eltype(b)))))
    end
end


@testitem "Derivatives: the generators and their adjoints" setup=[DerivativeTools] begin
    import SphericalFunctions: D_array, sYlm_array, D_array_widened, D_is_widened, D_narrowed,
        D_pushforward, D_pullback, sYlm_pushforward, sYlm_pullback, rotor_generator,
        rotor_cotangent, IntegerHalf, HalfOddInteger
    using Random: Xoshiro

    # The derivative through the recurrence, which is accurate away from the poles, is found
    # by calling the methods of `D_array` and `sYlm_array` for a general rotor, which the
    # rules do not intercept.
    recurrence_D(R, limits...) = invoke(D_array, Tuple{Rotor, ntuple(_ -> typeof(limits[1]), 5)...}, R, limits...)
    recurrence_Y(R, ℓₘₐₓ, s, ℓₘᵢₙ) = invoke(sYlm_array, Tuple{Rotor, typeof(ℓₘₐₓ), Any, typeof(ℓₘᵢₙ)}, R, ℓₘₐₓ, s, ℓₘᵢₙ)
    h(x) = HalfOddInteger(x)

    rng = Xoshiro(1234)
    q = 1.7 .* [0.3, -0.5, 0.7, 0.2]
    for q̇ ∈ ([0.1, 0.4, -0.3, 0.9], q, randn(rng, 4))  # the second changes only the norm
        path(t) = rotorof(q .+ t .* q̇)
        v = rotor_generator(q, q̇)
        for limits ∈ ((4, 4, -4, 4, -4), (4, 2, -1, 3, -2), h.((7//2, 7//2, -7//2, 7//2, -7//2)), h.((7//2, 5//2, -1//2, 7//2, -3//2)))
            Aʷ = D_array_widened(rotorof(q), limits...)
            A = D_array(rotorof(q), limits...)
            # The block restricted in m′ is read from the widened values bit for bit.
            @test D_narrowed(Aʷ, limits...) == A
            @test D_is_widened(limits[1:3]...) == (limits[2] != limits[1] || limits[3] != -limits[1])
            Ȧ = D_pushforward(Aʷ, v, limits...)
            @test within(Ȧ, ForwardDiff.derivative(t -> recurrence_D(path(t), limits...), 0.0); rtol=16eps())
            # The pullback is the adjoint of the pushforward: Σ Re(conj(Ā) Ȧ) = R̄ ⋅ Ṙ.
            Ā = randn(rng, ComplexF64, length(A))
            R̄ = rotor_cotangent(q, D_pullback(Aʷ, Ā, limits...))
            @test sum(real(conj.(Ā) .* Ȧ)) ≈ sum(R̄ .* q̇) atol=16eps() * (sum(abs, conj.(Ā) .* Ȧ) + sum(abs, R̄ .* q̇))
            # The cotangent of a function of R/‖R‖ is orthogonal to R.
            @test abs(sum(R̄ .* q)) ≤ 16eps() * sum(abs, R̄ .* q)
        end
        for (ℓₘₐₓ, s, ℓₘᵢₙ) ∈ ((4, -1, 1), (4, -2:1, 0), (4, 2, 0), (h(7//2), h(1//2), h(1//2)), (h(7//2), h(-3//2):h(1//2), h(1//2)))
            Y = sYlm_array(rotorof(q), ℓₘₐₓ, s, ℓₘᵢₙ)
            Ẏ = sYlm_pushforward(Y, v, ℓₘᵢₙ, ℓₘₐₓ)
            @test within(Ẏ, ForwardDiff.derivative(t -> recurrence_Y(path(t), ℓₘₐₓ, s, ℓₘᵢₙ), 0.0); rtol=16eps())
            Ȳ = randn(rng, ComplexF64, size(Y))
            R̄ = rotor_cotangent(q, sYlm_pullback(Y, Ȳ, ℓₘᵢₙ, ℓₘₐₓ))
            @test sum(real(conj.(Ȳ) .* Ẏ)) ≈ sum(R̄ .* q̇) atol=16eps() * (sum(abs, conj.(Ȳ) .* Ẏ) + sum(abs, R̄ .* q̇))
            @test abs(sum(R̄ .* q)) ≤ 16eps() * sum(abs, R̄ .* q)
        end
    end

    # `rotor_cotangent` is the adjoint of `rotor_generator`, for rotors of any norm.
    for _ ∈ 1:10
        R, Ṙ, g = randn(rng, 4), randn(rng, 4), randn(rng, 3)
        @test sum(g .* rotor_generator(R, Ṙ)) ≈ sum(rotor_cotangent(R, g) .* Ṙ)
    end
end


@testitem "Derivatives: D and sYlm label the same values as the calculators" setup=[DerivativeTools] begin
    import SphericalFunctions: DCalculator, sYlmCalculator, recurrence!, sYlm!, HarmonicValues
    # `D` of a rotor is computed through `D_array`, which the rules attach to, and labelled
    # by `D_series`.  Its blocks are those a calculator returns, bit for bit, and of their
    # types.
    R = rotorof(1.7 .* [0.3, -0.5, 0.7, 0.2])
    for (n, limits) ∈ ((4, (;)), (4, (m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=3, mₘᵢₙ=-2)), (7//2, (;)), (7//2, (m′ₘₐₓ=3//2, m′ₘᵢₙ=-1//2)))
        S = D(R, n; limits...)
        c = DCalculator(R, n; limits...)
        for (ℓ, b) ∈ c
            @test S[ℓ] == b
            @test typeof(S[ℓ]) == typeof(copy(b))
        end
    end
end


@testitem "Derivatives: ChainRules rules" setup=[DerivativeTools] begin
    import ChainRulesCore
    using ChainRulesCore: NoTangent, ZeroTangent, Tangent, @thunk
    import SphericalFunctions: D_array, sYlm_array, rotor_generator, D_pushforward,
        D_array_widened, sYlm_pushforward
    using StaticArrays: SVector

    q = 1.7 .* [0.3, -0.5, 0.7, 0.2]
    q̇ = [0.1, 0.4, -0.3, 0.9]
    R = rotorof(q)
    Y = sYlm_array(R, 4, -1, 1)
    Ẏ = sYlm_pushforward(Y, rotor_generator(q, q̇), 1, 4)
    # A tangent of the rotor is read by its components, whatever its type.
    tangents = (
        Quaternion(q̇...),
        rotorof(q̇),
        q̇,
        Tangent{typeof(R)}(components=SVector{4}(q̇...)),
        Tangent{typeof(R)}(components=Tangent{SVector{4, Float64}}(data=Tuple(q̇))),
        @thunk(Quaternion(q̇...)),
    )
    for Ṙ ∈ tangents
        Ω, Ω̇ = ChainRulesCore.frule((NoTangent(), Ṙ, NoTangent(), NoTangent(), NoTangent()), sYlm_array, R, 4, -1, 1)
        @test Ω == Y
        @test Ω̇ == Ẏ
    end
    @test ChainRulesCore.frule((NoTangent(), ZeroTangent(), NoTangent(), NoTangent(), NoTangent()), sYlm_array, R, 4, -1, 1)[2] isa ZeroTangent
    @test ChainRulesCore.frule((NoTangent(), Tangent{typeof(R)}(), NoTangent(), NoTangent(), NoTangent()), sYlm_array, R, 4, -1, 1)[2] isa ZeroTangent

    limits = (4, 2, -1, 3, -2)
    A = D_array(R, limits...)
    Ȧ = D_pushforward(D_array_widened(R, limits...), rotor_generator(q, q̇), limits...)
    Ω, Ω̇ = ChainRulesCore.frule((NoTangent(), Quaternion(q̇...), ntuple(_ -> NoTangent(), 5)...), D_array, R, limits...)
    @test Ω == A
    @test Ω̇ == Ȧ

    # The rules check the limits they are given, as `D` does, although only the widened
    # limits reach a calculator.
    for bad ∈ ((4, -1, -2, 4, -4), (4, 5, -4, 4, -4))
        @test_throws ArgumentError D_array(R, bad...)
        @test_throws ArgumentError ChainRulesCore.rrule(D_array, R, bad...)
        @test_throws ArgumentError ChainRulesCore.frule((NoTangent(), Quaternion(q̇...), ntuple(_ -> NoTangent(), 5)...), D_array, R, bad...)
        @test_throws ArgumentError ForwardDiff.derivative(t -> D_array(rotorof(q .+ t .* q̇), bad...), 0.0)
    end

    # The cotangent is a `Quaternion`, the adjoint of the pushforward, and orthogonal to R.
    for (Ω, back, Ω̇) ∈ (
        (ChainRulesCore.rrule(sYlm_array, R, 4, -1, 1)..., Ẏ),
        (ChainRulesCore.rrule(D_array, R, limits...)..., Ȧ),
    )
        Ω̄ = randn(ComplexF64, size(Ω))
        cotangents = back(Ω̄)
        R̄ = cotangents[2]
        @test R̄ isa Quaternion{Float64}
        @test all(c -> c isa NoTangent, cotangents[[1; 3:end]])
        @test sum(real(conj.(Ω̄) .* Ω̇)) ≈ sum(R̄[i] * q̇[i] for i ∈ 1:4)
        @test abs(sum(R̄[i] * q[i] for i ∈ 1:4)) ≤ 64eps()
        @test back(ZeroTangent())[2] isa ZeroTangent
        @test back(@thunk(Ω̄))[2] == R̄
    end
end


@testitem "Derivatives: first order, with every backend" setup=[ExplicitWignerMatrices, DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoReverseDiff, AutoEnzyme,
        AutoMooncake, AutoMooncakeForward, AutoZygote
    using DifferentiationInterfaceTest: Scenario, test_differentiation
    import Enzyme, Mooncake, ReverseDiff, Zygote

    # Every backend that DifferentiationInterface offers for the packages with rules.  Each
    # reaches the rules its own way: ForwardDiff, ReverseDiff, Enzyme, and Mooncake through
    # their own, and Zygote through ChainRules'.  The rotor is built without normalizing it
    # (see `rotorof`) because the normalization of `Rotor(v...)` fails in two of them, for
    # reasons that have nothing to do with the rules: when ReverseDiff replays its tape, as
    # it does for a Jacobian, the broadcast in Quaternionic's `rotor` writes into an
    # `SVector`, and Enzyme's own rule for `hypot`, which Quaternionic's `abs` calls, fails
    # in batched reverse mode.
    backends = [
        AutoForwardDiff(),
        AutoReverseDiff(),
        AutoReverseDiff(compile=true),
        AutoEnzyme(mode=Enzyme.Forward),
        AutoEnzyme(mode=Enzyme.Reverse),
        AutoMooncake(config=nothing),
        AutoMooncakeForward(config=nothing),
        AutoZygote(),
    ]
    scenarios = [
        Scenario{:jacobian, :out}(
            f, v; res1=reference_jacobian(f, v), prep_args=(; x=copy(v), contexts=()),
            name="$f at $point"
        )
        for f ∈ functions for (point, v) ∈ test_points(largest_ℓ(f))
    ]
    # The largest error measured is about 3eps() times the largest element.
    test_differentiation(
        backends, scenarios; correctness=true, isapprox=within, atol=0, rtol=16eps(),
        detailed=true, testset_name="Jacobians"
    )
end


@testitem "Derivatives: second order" setup=[ExplicitWignerMatrices, DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoReverseDiff, AutoZygote, SecondOrder
    using DifferentiationInterfaceTest: Scenario, test_differentiation
    import ReverseDiff, Zygote
    using Random: Xoshiro

    # Hessians of a real combination of the values.  Nested ForwardDiff applies the rule to
    # dual numbers of dual numbers, and forward-over-reverse applies the forward rule to the
    # reverse rule's own arithmetic.  The combinations that are missing here fail in the
    # tools themselves rather than in the rules.  Mooncake refuses dual numbers, and cannot
    # yet differentiate its own reverse pass in forward mode; reverse-over-forward, with
    # either ReverseDiff or Zygote outside, fails on a conversion within the tools; and
    # Enzyme over Enzyme, although it gives the right Hessian of the harmonics when runtime
    # activity is enabled in the outer pass, aborts the whole process for 𝔇 with
    # half-integer indices, on an assertion in Enzyme's handling of Julia's calling
    # convention (as of Enzyme 0.13.205 on Julia 1.13).
    backends = [
        AutoForwardDiff(),
        SecondOrder(AutoForwardDiff(), AutoReverseDiff()),
        SecondOrder(AutoForwardDiff(), AutoZygote()),
    ]
    struct Combination{F}
        f::F
        c::Vector{Float64}
    end
    (h::Combination)(v) = sum(h.c .* h.f(v))
    rng = Xoshiro(5678)
    scenarios = map([functions[1], functions[3], functions[5], functions[7]]) do f
        c = randn(rng, length(f(test_points(largest_ℓ(f))[end][2])))
        h = Combination(f, c)
        [
            Scenario{:hessian, :out}(
                h, v;
                res1=Float64.(ForwardDiff.gradient(w -> sum(c .* polynomial(f, w)), big.(v))),
                res2=Float64.(ForwardDiff.hessian(w -> sum(c .* polynomial(f, w)), big.(v))),
                prep_args=(; x=copy(v), contexts=()), name="Hessian of $f at $point"
            )
            for (point, v) ∈ test_points(largest_ℓ(f))
        ]
    end
    test_differentiation(
        backends, reduce(vcat, scenarios); correctness=true, isapprox=within, atol=0,
        rtol=64eps(), detailed=true, testset_name="Hessians"
    )
end
