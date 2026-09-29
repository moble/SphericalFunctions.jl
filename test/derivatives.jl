# Tests of the rules for automatic differentiation (see `src/derivatives.jl`): the rules for
# the calculators' steps, which ForwardDiff, Enzyme, and Mooncake use, and through which
# they differentiate `D`, `sYlm`, and `sYlm_matrix` too; and the rules for those functions'
# arrays, which ChainRules and ReverseDiff use.  The references are the explicit polynomial
# of the `ExplicitWignerMatrices` module, which shares no code with the package, evaluated
# in `BigFloat` and differentiated by ForwardDiff.  The polynomial is smooth everywhere, so
# its derivatives are good references at the poles too.

@testsnippet DerivativeTools begin
    import ForwardDiff
    using Quaternionic: Quaternion, Rotor, from_euler_angles
    import SphericalFunctions
    import SphericalFunctions: D, sYlm, sYlm_matrix, DCalculator, sYlmCalculator, array_view,
        Ysize

    # The rotation data of rotor `k` from the components `v[4k-3:4k]`, as a `Quaternion` or
    # as a `Rotor` whose norm need not be 1, so that derivatives are taken in all four
    # directions of the space of quaternions, including the one that changes only the norm.
    # Normalizing with `Rotor(v...)` would be differentiable too, but it goes through
    # Quaternionic's `abs`, which some backends cannot follow; see the note on the backends
    # in "Derivatives: first order, with every backend".
    quaternionof(v, k=1) = Quaternion{eltype(v)}(v[4k - 3], v[4k - 2], v[4k - 1], v[4k])
    rotorof(v, k=1) = Rotor{eltype(v)}(v[4k - 3], v[4k - 2], v[4k - 1], v[4k])

    # Complex values as one real vector, of their real parts followed by their imaginary parts
    realified(z) = [real.(vec(z)); imag.(vec(z))]

    # The functions tested: each is a callable struct that gives some values as a real
    # vector, with a method of `polynomial` that gives the same values from the polynomial.
    # They are callable structs rather than closures so that Enzyme compiles each kind once,
    # for all of the points at which it is tested.
    ℓs(n) = (n isa Integer ? 0 : 1//2):n
    rows(ℓ, (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)) = max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)
    cols(ℓ, (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)) = max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ)
    full(n) = (n, -n, n, -n)
    spins(s) = s isa AbstractRange ? s : (s,)
    Y_polynomial(ℓ, m, s, v) =
        abs(s) ≤ ℓ ? ExplicitWignerMatrices.sYlm_polynomial(ℓ, m, s, v)[1] : zero(Complex{eltype(v)})
    components(v, k) = v[(4k - 3):4k]

    # All the blocks of `D` of one quaternion
    struct DValues{T, L}
        n::T
        limits::L  # (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
    DValues(n) = DValues(n, full(n))
    function (f::DValues)(v)
        m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ = f.limits
        S = D(quaternionof(v), f.n; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        realified(reduce(vcat, [vec(parent(b)) for b ∈ values(S)]))
    end
    polynomial(f::DValues, v) = realified([
        ExplicitWignerMatrices.D_polynomial(ℓ, m′, m, v)[1]
        for ℓ ∈ ℓs(f.n) for m ∈ cols(ℓ, f.limits) for m′ ∈ rows(ℓ, f.limits)
    ])

    # `sYlm` of one rotor
    struct YValues{T, S}
        n::T
        s::S
        ℓₘᵢₙ::T
    end
    (f::YValues)(v) = realified(array_view(sYlm(rotorof(v), f.n, f.s; ℓₘᵢₙ=f.ℓₘᵢₙ)))
    polynomial(f::YValues, v) = realified([
        Y_polynomial(ℓ, m, s, v) for ℓ ∈ f.ℓₘᵢₙ:f.n for m ∈ -ℓ:ℓ for s ∈ spins(f.s)
    ])

    # `sYlm_matrix` of `nᵣ` quaternions
    struct YMatrixValues{T, S}
        n::T
        s::S
        ℓₘᵢₙ::T
        nᵣ::Int
    end
    function (f::YMatrixValues)(v)
        realified(sYlm_matrix([quaternionof(v, k) for k ∈ 1:f.nᵣ], f.n, f.s; ℓₘᵢₙ=f.ℓₘᵢₙ))
    end
    polynomial(f::YMatrixValues, v) = realified([
        Y_polynomial(ℓ, m, s, components(v, k))
        for ℓ ∈ f.ℓₘᵢₙ:f.n for m ∈ -ℓ:ℓ for s ∈ spins(f.s) for k ∈ 1:f.nᵣ
    ])

    # An element of a block, for the rotor iᵣ, whether or not the block has a rotor axis and
    # a spin axis.  The choice is by dispatch, as a user's code would make it, so that the
    # loops below are type-stable, which Enzyme needs.
    entry(b::SphericalFunctions.WignerMatrix, iᵣ, m′, m) = b[m′, m]
    entry(b::SphericalFunctions.WignerMatrixBatch, iᵣ, m′, m) = b[iᵣ, m′, m]
    entry(b::SphericalFunctions.DegreeBlock, iᵣ, s, m) = b[m]
    entry(b::SphericalFunctions.DegreeBlockBatch, iᵣ, s, m) = b[iᵣ, m]
    entry(b::SphericalFunctions.SpinMatrix, iᵣ, s, m) = b[s, m]
    entry(b::SphericalFunctions.SpinMatrixBatch, iᵣ, s, m) = b[iᵣ, s, m]

    # Every block of a `DCalculator` of `nᵣ` quaternions, in turn, as a loop over the
    # calculator reads them; batched when `nᵣ > 1`.  The number of rotors `N` is a type
    # parameter, so that whether the calculator is batched is settled at compile time, as it
    # is for a user's own code.
    struct DCalcValues{T, L, N}
        n::T
        limits::L
        length::Int  # of the values, before they are made real
    end
    DCalcValues(n, limits, N) = DCalcValues{typeof(n), typeof(limits), N}(
        n, limits, sum(length(rows(ℓ, limits)) * length(cols(ℓ, limits)) for ℓ ∈ ℓs(n)) * N
    )
    rotor_data(v, ::Val{1}, make) = make(v)
    rotor_data(v, ::Val{N}, make) where {N} = [make(v, k) for k ∈ 1:N]
    function (f::DCalcValues{T, L, N})(v) where {T, L, N}
        m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ = f.limits
        c = DCalculator(rotor_data(v, Val(N), quaternionof), f.n; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        z = Vector{eltype(eltype(values(c)))}(undef, f.length)
        k = 0
        # The calculator's ℓ is of its own index type, a `HalfOddInteger` for half-integer ℓ,
        # which does not mix with the `Rational` limits, so they are converted to it.
        for (ℓ, b) ∈ c
            limits = map(x -> typeof(ℓ)(x), f.limits)
            for m ∈ cols(ℓ, limits), m′ ∈ rows(ℓ, limits), iᵣ ∈ 1:N
                z[k += 1] = entry(b, iᵣ, m′, m)
            end
        end
        realified(z)
    end
    polynomial_values(f::DCalcValues{T, L, N}, v) where {T, L, N} = [
        ExplicitWignerMatrices.D_polynomial(ℓ, m′, m, components(v, k))[1]
        for ℓ ∈ ℓs(f.n) for m ∈ cols(ℓ, f.limits) for m′ ∈ rows(ℓ, f.limits) for k ∈ 1:N
    ]
    polynomial(f::DCalcValues, v) = realified(polynomial_values(f, v))

    # Every block of an `sYlmCalculator` of `N` rotors, in turn; batched when `N > 1`
    struct YCalcValues{T, S, N}
        n::T
        s::S
        length::Int  # of the values, before they are made real
    end
    YCalcValues(n, s, N) =
        YCalcValues{typeof(n), typeof(s), N}(n, s, Ysize(minimum(abs, spins(s)), n) * length(spins(s)) * N)
    function (f::YCalcValues{T, S, N})(v) where {T, S, N}
        c = sYlmCalculator(rotor_data(v, Val(N), rotorof), f.n, f.s)
        z = Vector{eltype(eltype(values(c)))}(undef, f.length)
        k = 0
        for (ℓ, b) ∈ c, m ∈ -ℓ:ℓ, s ∈ spins(f.s), iᵣ ∈ 1:N
            z[k += 1] = entry(b, iᵣ, s, m)
        end
        realified(z)
    end
    polynomial_values(f::YCalcValues{T, S, N}, v) where {T, S, N} = [
        Y_polynomial(ℓ, m, s, components(v, k))
        for ℓ ∈ minimum(abs, spins(f.s)):f.n for m ∈ -ℓ:ℓ for s ∈ spins(f.s) for k ∈ 1:N
    ]
    polynomial(f::YCalcValues, v) = realified(polynomial_values(f, v))

    for F ∈ (DValues, YValues, YMatrixValues, DCalcValues, YCalcValues)
        @eval Base.show(io::IO, f::$F) = print(io, $(string(F)), "(", join(map(string, Tuple(getfield(f, k) for k ∈ fieldnames($F))), ", "), ")")
    end
    nrotors(f) = hasfield(typeof(f), :nᵣ) ? f.nᵣ : 1
    nrotors(::DCalcValues{T, L, N}) where {T, L, N} = N
    nrotors(::YCalcValues{T, S, N}) where {T, S, N} = N
    largest_ℓ(f) = Float64(f.n)
    reference_jacobian(f, v) = Float64.(ForwardDiff.jacobian(w -> polynomial(f, w), big.(v)))

    # The functions of whole arrays, with integer and half-integer indices, whole and
    # restricted blocks of 𝔇, one spin weight or a range of them, and one rotor or several;
    # and the loops over calculators, of one rotor or a batch.
    array_functions = [
        DValues(3),
        DValues(4, (2, -1, 4, -4)),
        DValues(4, (4, -4, 2, -1)),
        DValues(7//2, (3//2, -1//2, 7//2, -5//2)),
        YValues(4, -1, 1),
        YValues(4, -2:1, 0),
        YValues(7//2, 1//2, 1//2),
        YMatrixValues(3, -1, 1, 2),
        YMatrixValues(5//2, -1//2:1//2, 1//2, 2),
    ]
    calculator_functions = [
        DCalcValues(3, full(3), 1),
        DCalcValues(4, (2, -1, 3, -2), 2),
        DCalcValues(5//2, full(5//2), 1),
        YCalcValues(4, -2:1, 2),
        YCalcValues(7//2, 1//2, 1),
    ]

    # Rotors at the poles, near them, and away from both, the last not normalized.  At a
    # distance r from a pole, the first derivatives of the recurrence itself are wrong by
    # about `eps()`/r relative to their size, and at a pole they are NaN, so that a backend
    # that bypassed the rules would fail at these points.  A function of several rotors is
    # tested with each of these as its first rotor, and a generic rotor for the others.
    function test_points(ℓₘₐₓ, nᵣ=1)
        components(R) = [R[1], R[2], R[3], R[4]]
        near(r, β₀=0) = components(from_euler_angles(0.4, β₀ + 2asin(r), -1.3))
        others = reduce(vcat, [[0.2, 0.9 - 0.1k, -0.3, 0.1k] for k ∈ 1:(nᵣ - 1)]; init=Float64[])
        [
            name => vcat(v, others) for (name, v) ∈ (
                "the identity" => [1.0, 0.0, 0.0, 0.0],
                "a rotation about z" => [cos(0.35), 0.0, 0.0, sin(0.35)],
                "a rotor at β = π" => components(
                    Quaternion(0.0, 0.0, 1.0, 0.0) * Quaternion(cos(0.35), 0.0, 0.0, sin(0.35))
                ),
                "a rotor near β = 0" => near(1e-7),
                "a rotor near β = π" => near(-1e-3, π),
                "a generic rotor" => 1.7 .* [0.3, -0.5, 0.7, 0.2],
            )
        ]
    end
    points(f) = test_points(largest_ℓ(f), nrotors(f))

    # Agreement to within `rtol` times the largest element of the reference, in the maximum
    # norm.  DifferentiationInterfaceTest calls this with its own keywords.
    function within(a, b; atol=0, rtol)
        maximum(abs, a - b; init=zero(real(eltype(b)))) ≤
            atol + rtol * max(1, maximum(abs, b; init=zero(real(eltype(b)))))
    end
end


@testitem "Derivatives: the kernels are adjoint to themselves" setup=[DerivativeTools] begin
    import SphericalFunctions: rotor_generator, rotor_cotangent, wigner_block_pushforward!,
        wigner_block_pullback!, harmonic_block_pushforward!, harmonic_block_pullback!,
        HalfOddInteger
    using Random: Xoshiro

    # Each pullback is the adjoint of its pushforward, Σ Re(conj(z̄) ż) = Σ Ḡ ⋅ G, whatever
    # the values, which here are random, as are the generators and the cotangents; the
    # comparisons with the polynomial below test that the pushforwards are the right ones.
    rng = Xoshiro(1234)
    h(x) = x isa Integer ? x : HalfOddInteger(x)
    for Nᵣ ∈ (1, 3), (ℓ, rows, cols, outrows, outcols) ∈ (
        (4, -4:4, -4:4, -4:4, -4:4),
        (4, -2:3, -4:4, -2:3, -4:4),
        (4, -4:4, -2:3, -4:4, -2:3),
        (4, -3:4, -2:3, -2:3, -2:3),
        (4, -2:3, -3:4, -2:3, -2:3),
        (h(7//2), h(-5//2):h(5//2), h(-7//2):h(7//2), h(-5//2):h(5//2), h(-7//2):h(7//2)),
    ), left ∈ (true, false)
        # A block is differentiated from the side whose neighbors it holds.
        can = left ? (first(outrows) == -ℓ || first(outrows) > first(rows)) && (last(outrows) == ℓ || last(outrows) < last(rows)) :
            (first(outcols) == -ℓ || first(outcols) > first(cols)) && (last(outcols) == ℓ || last(outcols) < last(cols))
        can || continue
        A = randn(rng, ComplexF64, Nᵣ, length(rows), length(cols))
        G = randn(rng, 3, Nᵣ)
        Ȧ = similar(A, Nᵣ, length(outrows), length(outcols))
        wigner_block_pushforward!((x, ẋ) -> only(ẋ), Ȧ, A, ℓ, rows, cols, outrows, outcols, left, G, Val(1))
        Ā = randn(rng, ComplexF64, size(Ȧ))
        Ḡ = wigner_block_pullback!(zeros(3, Nᵣ), A, Ā, ℓ, rows, cols, outrows, outcols, left)
        @test sum(real(conj.(Ā) .* Ȧ)) ≈ sum(Ḡ .* G)
    end
    for Nᵣ ∈ (1, 3), (ℓ, nₛ) ∈ ((4, 1), (4, 3), (h(5//2), 2))
        n = Int(2ℓ) + 1
        A = randn(rng, ComplexF64, Nᵣ, nₛ, n)
        G = randn(rng, 3, Nᵣ)
        Ȧ = similar(A)
        harmonic_block_pushforward!((x, ẋ) -> only(ẋ), Ȧ, A, ℓ, 1:nₛ, G, Val(1))
        Ā = randn(rng, ComplexF64, size(Ȧ))
        Ḡ = harmonic_block_pullback!(zeros(3, Nᵣ), A, Ā, ℓ, 1:nₛ)
        @test sum(real(conj.(Ā) .* Ȧ)) ≈ sum(Ḡ .* G)
    end

    # `rotor_cotangent` is the adjoint of `rotor_generator` from either side, for rotors of
    # any norm, and the cotangent of a function of R/‖R‖ is orthogonal to R.
    for _ ∈ 1:10, left ∈ (true, false)
        R, Ṙ, g = randn(rng, 4), randn(rng, 4), randn(rng, 3)
        R̄ = rotor_cotangent(left, R, g)
        @test sum(g .* rotor_generator(left, R, Ṙ)) ≈ sum(R̄ .* Ṙ)
        @test abs(sum(R̄ .* R)) ≤ 16eps() * sum(abs, R̄ .* R)
    end
end


@testitem "Derivatives: D and sYlm label the same values as the calculators" setup=[DerivativeTools] begin
    import SphericalFunctions: recurrence!
    # `D` of a rotor is computed through `D_array`, which the rules for ChainRules and
    # ReverseDiff attach to, and labelled by `D_series`.  Its blocks are those a calculator
    # returns, bit for bit, and of their types, whether or not the calculator stores more rows
    # or columns than the block.
    R = rotorof(1.7 .* [0.3, -0.5, 0.7, 0.2])
    for (n, limits) ∈ (
        (4, (;)), (4, (m′ₘₐₓ=2, m′ₘᵢₙ=-1)), (4, (mₘₐₓ=2, mₘᵢₙ=-1)),
        (4, (m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=3, mₘᵢₙ=-2)), (4, (m′ₘₐₓ=3, m′ₘᵢₙ=-2, mₘₐₓ=2, mₘᵢₙ=-1)),
        (7//2, (;)), (7//2, (m′ₘₐₓ=3//2, m′ₘᵢₙ=-1//2))
    )
        S = D(R, n; limits...)
        c = DCalculator(R, n; limits...)
        for (ℓ, b) ∈ c
            @test S[ℓ] == b
            @test typeof(S[ℓ]) == typeof(copy(b))
        end
    end
end


@testitem "Derivatives: calculators of dual numbers allocate nothing when stepped" setup=[DerivativeTools] begin
    import SphericalFunctions: set_R!
    # A calculator whose rotors are dual numbers runs the recurrence on their values, and
    # lifts each block, into storage allocated with it, so a sweep over its blocks allocates
    # nothing, and neither does setting new rotors.  Nor does a calculator of floats, whose
    # sweep is exactly as it was.
    T = ForwardDiff.Tag{Nothing, Float64}
    dual(x, k) = ForwardDiff.Dual{T}(x, ntuple(i -> Float64(i == k), 4)...)
    q = 1.7 .* [0.3, -0.5, 0.7, 0.2]
    Rd = Rotor{ForwardDiff.Dual{T, Float64, 4}}((dual(q[k], k) for k ∈ 1:4)...)
    function sweep(c)
        s = 0.0
        for b ∈ values(c)
            s += ForwardDiff.value(real(first(b)))
        end
        s
    end
    reset!(c, R) = (set_R!(c, R); nothing)
    # Measured behind a function barrier, since `c` has a different type in each iteration
    # of the loop below, and Julia 1.10 boxes the result of a dynamic call.
    sweep_allocations(c) = @allocated sweep(c)
    for c ∈ (
        DCalculator(Rd, 8), DCalculator(Rd, 8; m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=3, mₘᵢₙ=-2),
        DCalculator([Rd, Rd], 8), sYlmCalculator(Rd, 8, -1), sYlmCalculator([Rd, Rd], 8, -2:2),
        DCalculator(rotorof(q), 8), sYlmCalculator(rotorof(q), 8, -1),
    )
        R = eltype(c.rotors) <: Quaternion{Float64} ? rotorof(q) : Rd
        R = SphericalFunctions.isbatched(c) ? [R, R] : R
        sweep(c); reset!(c, R)
        @test sweep_allocations(c) == 0
        @test @allocated(reset!(c, R)) == 0
    end
end


@testitem "Derivatives: ChainRules rules" setup=[DerivativeTools] begin
    import ChainRulesCore
    using ChainRulesCore: NoTangent, ZeroTangent, Tangent, @thunk
    import SphericalFunctions: D_array, sYlm_array, sYlm_matrix_array, D_array_with_stored,
        D_array_pushforward, derivatives_from_left, rotor_generator, rotor_generators,
        harmonic_array_pushforward
    using StaticArrays: SVector

    q = 1.7 .* [0.3, -0.5, 0.7, 0.2]
    q̇ = [0.1, 0.4, -0.3, 0.9]
    R = rotorof(q)
    Y = sYlm_array(R, 4, -1, 1)
    Ẏ = harmonic_array_pushforward(
        (y, ẏ) -> only(ẏ), Y, false, 1, 4, rotor_generators(true, [R], [q̇]), Val(1)
    )
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
    _, stored, calc = D_array_with_stored(R, limits...)
    Ȧ = D_array_pushforward(calc, stored, rotor_generator(derivatives_from_left(calc), R, q̇))
    Ω, Ω̇ = ChainRulesCore.frule((NoTangent(), Quaternion(q̇...), ntuple(_ -> NoTangent(), 5)...), D_array, R, limits...)
    @test Ω == A
    @test Ω̇ == Ȧ

    # The rules check the limits they are given, as `D` does.
    for bad ∈ ((4, -1, -2, 4, -4), (4, 5, -4, 4, -4))
        @test_throws ArgumentError D_array(R, bad...)
        @test_throws ArgumentError ChainRulesCore.rrule(D_array, R, bad...)
        @test_throws ArgumentError ChainRulesCore.frule((NoTangent(), Quaternion(q̇...), ntuple(_ -> NoTangent(), 5)...), D_array, R, bad...)
        @test_throws ArgumentError ForwardDiff.derivative(t -> D_array(rotorof(q .+ t .* q̇), bad...), 0.0)
    end

    # The cotangent is a `Quaternion`, the adjoint of the pushforward, and orthogonal to R;
    # that of a vector of rotors is a vector of them.
    flat(x::AbstractVector{<:AbstractArray}) = reduce(vcat, vec.(x))
    flat(x::AbstractArray) = vec(x)
    for (Ω, back, Ω̇) ∈ (
        (ChainRulesCore.rrule(sYlm_array, R, 4, -1, 1)..., Ẏ),
        (ChainRulesCore.rrule(D_array, R, limits...)..., Ȧ),
    )
        Ω̄ = Ω isa AbstractVector{<:AbstractArray} ? [randn(ComplexF64, size(b)) for b ∈ Ω] : randn(ComplexF64, size(Ω))
        cotangents = back(Ω̄)
        R̄ = cotangents[2]
        @test R̄ isa Quaternion{Float64}
        @test all(c -> c isa NoTangent, cotangents[[1; 3:end]])
        @test sum(real(conj.(flat(Ω̄)) .* flat(Ω̇))) ≈ sum(R̄[i] * q̇[i] for i ∈ 1:4)
        @test abs(sum(R̄[i] * q[i] for i ∈ 1:4)) ≤ 64eps()
        @test back(ZeroTangent())[2] isa ZeroTangent
        @test back(@thunk(Ω̄))[2] == R̄
    end
    R⃗ = [R, rotorof(q̇)]
    Ω, back = ChainRulesCore.rrule(sYlm_matrix_array, R⃗, 3, -1:1, 1)
    R̄⃗ = back(randn(ComplexF64, size(Ω)))[2]
    @test R̄⃗ isa Vector{Quaternion{Float64}}
    @test length(R̄⃗) == 2
end


@testitem "Derivatives: first order, with every backend" setup=[ExplicitWignerMatrices, DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoReverseDiff, AutoEnzyme,
        AutoMooncake, AutoMooncakeForward, AutoZygote
    using DifferentiationInterfaceTest: Scenario, test_differentiation
    import Enzyme, Mooncake, ReverseDiff, Zygote

    # Every backend that DifferentiationInterface offers for the packages with rules.  Each
    # reaches the rules its own way: ForwardDiff, Enzyme, and Mooncake through those for the
    # calculators, and ReverseDiff and Zygote through those for the arrays of `D`, `sYlm`, and
    # `sYlm_matrix`.  The rotors are built without normalizing them (see `quaternionof`),
    # because the normalization of `Rotor(v...)` fails in two of them, for reasons that have
    # nothing to do with the rules: when ReverseDiff replays its tape, as it does for a
    # Jacobian, the broadcast in Quaternionic's `rotor` writes into an `SVector`, and Enzyme's
    # own rule for `hypot`, which Quaternionic's `abs` calls, fails in batched reverse mode.
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
        for f ∈ array_functions for (point, v) ∈ points(f)
    ]
    # The largest error measured is about 3eps() times the largest element.
    test_differentiation(
        backends, scenarios; correctness=true, isapprox=within, atol=0, rtol=16eps(),
        detailed=true, testset_name="Jacobians"
    )
end


@testitem "Derivatives: first order, through calculators" setup=[ExplicitWignerMatrices, DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoReverseDiff, AutoEnzyme,
        AutoMooncake, AutoMooncakeForward
    using DifferentiationInterfaceTest: Scenario, test_differentiation
    import Enzyme, Mooncake, ReverseDiff

    # Loops over calculators, as the documentation recommends for many ℓ or many rotors,
    # which every backend here differentiates by the rules for the calculators' steps.
    # Zygote cannot follow the mutation of a calculator's buffers at all, and is not tested
    # here.  Nor is Enzyme's reverse mode on Julia 1.10, where its type analysis cannot
    # deduce the type of an integer in the loop over the blocks when bounds are checked, as
    # `Pkg.test` checks them.
    backends = [
        AutoForwardDiff(),
        AutoReverseDiff(),
        AutoEnzyme(mode=Enzyme.Forward),
        AutoMooncake(config=nothing),
        AutoMooncakeForward(config=nothing),
    ]
    VERSION ≥ v"1.11" && insert!(backends, 4, AutoEnzyme(mode=Enzyme.Reverse))
    scenarios = [
        Scenario{:jacobian, :out}(
            f, v; res1=reference_jacobian(f, v), prep_args=(; x=copy(v), contexts=()),
            name="$f at $point"
        )
        for f ∈ calculator_functions for (point, v) ∈ points(f)
    ]
    # The rules are as accurate here as for the arrays.
    test_differentiation(
        backends, scenarios; correctness=true, isapprox=within, atol=0, rtol=16eps(),
        detailed=true, testset_name="Jacobians through calculators"
    )
end


@testitem "Derivatives: second order" setup=[ExplicitWignerMatrices, DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoReverseDiff, AutoZygote, SecondOrder
    using DifferentiationInterfaceTest: Scenario, test_differentiation
    import ReverseDiff, Zygote
    using Random: Xoshiro

    # Hessians of a real combination of the values.  Nested ForwardDiff lifts the blocks of a
    # calculator of dual numbers, whose values are themselves lifted from a calculator of
    # floats, and forward-over-reverse applies the forward rules to the reverse rules' own
    # arithmetic.  The combinations that are missing here fail in the tools themselves rather
    # than in the rules.  Mooncake refuses dual numbers, and cannot yet differentiate its own
    # reverse pass in forward mode; reverse-over-forward, with either ReverseDiff or Zygote
    # outside, fails on a conversion within the tools; and Enzyme over Enzyme, although it
    # gives the right Hessian of the harmonics when runtime activity is enabled in the outer
    # pass, aborts the whole process for 𝔇 with half-integer indices, on an assertion in
    # Enzyme's handling of Julia's calling convention (as of Enzyme 0.13.205 on Julia 1.13).
    struct Combination{F}
        f::F
        c::Vector{Float64}
    end
    (h::Combination)(v) = sum(h.c .* h.f(v))
    rng = Xoshiro(5678)
    function scenarios(functions)
        reduce(vcat, map(functions) do f
            v₀ = points(f)[end][2]
            c = randn(rng, length(f(v₀)))
            h = Combination(f, c)
            [
                Scenario{:hessian, :out}(
                    h, v;
                    res1=Float64.(ForwardDiff.gradient(w -> sum(c .* polynomial(f, w)), big.(v))),
                    res2=Float64.(ForwardDiff.hessian(w -> sum(c .* polynomial(f, w)), big.(v))),
                    prep_args=(; x=copy(v), contexts=()), name="Hessian of $f at $point"
                )
                for (point, v) ∈ points(f)
            ]
        end)
    end
    test_differentiation(
        [AutoForwardDiff(), SecondOrder(AutoForwardDiff(), AutoReverseDiff()), SecondOrder(AutoForwardDiff(), AutoZygote())],
        scenarios(array_functions[[1, 4, 5, 7, 8]]); correctness=true, isapprox=within,
        atol=0, rtol=64eps(), detailed=true, testset_name="Hessians"
    )
    test_differentiation(
        [AutoForwardDiff()], scenarios(calculator_functions[[1, 3, 5]]); correctness=true,
        isapprox=within, atol=0, rtol=64eps(), detailed=true,
        testset_name="Hessians through calculators"
    )
end
