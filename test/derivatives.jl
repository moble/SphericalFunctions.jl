# Tests of the rules for automatic differentiation (see `src/derivatives/kernels.jl`): the
# rules for the calculators' steps, of 𝔇, `d`, ₛYₗₘ, and ₛλₗₘ, which ForwardDiff, Enzyme,
# Mooncake, and ReverseDiff use, and through which the first three differentiate `D`, `d`,
# `sYlm`, and `sYlm_matrix` too; and the rules for those functions' arrays, which ChainRules
# and ReverseDiff use.  The references are the explicit polynomial of the
# `ExplicitWignerMatrices` module, which shares no code with the package, evaluated in
# `BigFloat` and differentiated by ForwardDiff.  The polynomial is smooth everywhere, so its
# derivatives are good references at the poles too.

@testsnippet DerivativeTools begin
    import ForwardDiff
    using Quaternionic: Quaternion, Rotor, from_euler_angles
    import SphericalFunctions
    import SphericalFunctions: D, d, sYlm, sYlm_matrix, DCalculator, sYlmCalculator,
        dCalculator, sλlmCalculator, array_view, Ysize

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

    # All the blocks of `D` of one quaternion.  Each block is copied as an `Array`, because
    # Enzyme cannot prove the types in a vector of the views of one vector that are the
    # blocks' storage.
    struct DValues{T, L}
        n::T
        limits::L  # (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
    DValues(n) = DValues(n, full(n))
    function (f::DValues)(v)
        m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ = f.limits
        S = D(quaternionof(v), f.n; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        realified(reduce(vcat, [vec(Array(b)) for b ∈ values(S)]))
    end
    polynomial(f::DValues, v) = realified([
        ExplicitWignerMatrices.D_polynomial(ℓ, m′, m, v)[1]
        for ℓ ∈ ℓs(f.n) for m ∈ cols(ℓ, f.limits) for m′ ∈ rows(ℓ, f.limits)
    ])

    # All the blocks of `d` of the one angle `x[1]`, which are those of `D` of the rotor
    # (cos(β/2), 0, sin(β/2), 0), copied as those of `D` are
    struct dValues{T, L}
        n::T
        limits::L  # (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
    dValues(n) = dValues(n, full(n))
    function (f::dValues)(x)
        m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ = f.limits
        S = d(x[1], f.n; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        reduce(vcat, [vec(Array(b)) for b ∈ values(S)])
    end
    function polynomial(f::dValues, x)
        c, s, z = cos(x[1] / 2), sin(x[1] / 2), zero(x[1])
        [
            real(ExplicitWignerMatrices.D_polynomial(ℓ, m′, m, [c, z, s, z])[1])
            for ℓ ∈ ℓs(f.n) for m ∈ cols(ℓ, f.limits) for m′ ∈ rows(ℓ, f.limits)
        ]
    end

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

    # All the blocks of `d` of the phase e^{ix[1]}, which are those of `d` of the angle x[1]
    struct dPhaseValues{T}
        n::T
    end
    (f::dPhaseValues)(x) = reduce(vcat, [vec(Array(b)) for b ∈ values(d(cis(x[1]), f.n))])
    polynomial(f::dPhaseValues, x) = polynomial(dValues(f.n), x)

    # All the blocks of `d` of a quaternion, which are those of `d` of the β of its Euler
    # decomposition
    struct dRotorValues{T}
        n::T
    end
    (f::dRotorValues)(v) = reduce(vcat, [vec(Array(b)) for b ∈ values(d(quaternionof(v), f.n))])
    polynomial(f::dRotorValues, v) =
        polynomial(dValues(f.n), [2atan(sqrt(v[2]^2 + v[3]^2), sqrt(v[1]^2 + v[4]^2))])

    # The angles of a calculator of `d` or of ₛλₗₘ, one or a vector of `N` of them
    angle_data(x, ::Val{1}) = x[1]
    angle_data(x, ::Val{N}) where {N} = x[1:N]
    # The rotor of a rotation through β about y
    y_rotor(β) = [cos(β / 2), zero(β), sin(β / 2), zero(β)]

    # Every block of a `dCalculator` of `N` angles, in turn; batched when `N > 1`
    struct dCalcValues{T, L, N}
        n::T
        limits::L
        length::Int
    end
    dCalcValues(n, limits, N) = dCalcValues{typeof(n), typeof(limits), N}(
        n, limits, sum(length(rows(ℓ, limits)) * length(cols(ℓ, limits)) for ℓ ∈ ℓs(n)) * N
    )
    function (f::dCalcValues{T, L, N})(x) where {T, L, N}
        m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ = f.limits
        c = dCalculator(angle_data(x, Val(N)), f.n; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        z = Vector{eltype(eltype(values(c)))}(undef, f.length)
        k = 0
        for (ℓ, b) ∈ c
            limits = map(x -> typeof(ℓ)(x), f.limits)
            for m ∈ cols(ℓ, limits), m′ ∈ rows(ℓ, limits), iᵣ ∈ 1:N
                z[k += 1] = entry(b, iᵣ, m′, m)
            end
        end
        z
    end
    polynomial(f::dCalcValues{T, L, N}, x) where {T, L, N} = [
        real(ExplicitWignerMatrices.D_polynomial(ℓ, m′, m, y_rotor(x[k]))[1])
        for ℓ ∈ ℓs(f.n) for m ∈ cols(ℓ, f.limits) for m′ ∈ rows(ℓ, f.limits) for k ∈ 1:N
    ]

    # Every block of an `sλlmCalculator` of `N` angles, in turn; batched when `N > 1`.  The
    # real harmonics are the harmonics at (θ, 0), divided by i^{2s} for a half-odd s.
    struct λCalcValues{T, S, N}
        n::T
        s::S
        length::Int
    end
    λCalcValues(n, s, N) =
        λCalcValues{typeof(n), typeof(s), N}(n, s, Ysize(minimum(abs, spins(s)), n) * length(spins(s)) * N)
    function (f::λCalcValues{T, S, N})(x) where {T, S, N}
        c = sλlmCalculator(angle_data(x, Val(N)), f.n, f.s)
        z = Vector{eltype(eltype(values(c)))}(undef, f.length)
        k = 0
        for (ℓ, b) ∈ c, m ∈ -ℓ:ℓ, s ∈ spins(f.s), iᵣ ∈ 1:N
            z[k += 1] = entry(b, iᵣ, s, m)
        end
        z
    end
    function λ_polynomial(ℓ, m, s, θ)
        Y = Y_polynomial(ℓ, m, s, y_rotor(θ))
        real(isinteger(s) ? Y : Y / (1, im, -1, -im)[mod(Int(2s), 4) + 1])
    end
    polynomial(f::λCalcValues{T, S, N}, x) where {T, S, N} = [
        λ_polynomial(ℓ, m, s, x[k])
        for ℓ ∈ minimum(abs, spins(f.s)):f.n for m ∈ -ℓ:ℓ for s ∈ spins(f.s) for k ∈ 1:N
    ]

    for F ∈ (DValues, dValues, YValues, YMatrixValues, DCalcValues, YCalcValues, dPhaseValues,
             dRotorValues, dCalcValues, λCalcValues)
        @eval Base.show(io::IO, f::$F) = print(io, $(string(F)), "(", join(map(string, Tuple(getfield(f, k) for k ∈ fieldnames($F))), ", "), ")")
    end
    nrotors(f) = hasfield(typeof(f), :nᵣ) ? f.nᵣ : 1
    nrotors(::DCalcValues{T, L, N}) where {T, L, N} = N
    nrotors(::YCalcValues{T, S, N}) where {T, S, N} = N
    nrotors(::dCalcValues{T, L, N}) where {T, L, N} = N
    nrotors(::λCalcValues{T, S, N}) where {T, S, N} = N
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
        DValues(4, (2, -1, 3, -2)),
    ]
    angle_array_functions = [
        dValues(4), dValues(4, (2, -1, 3, -2)), dValues(7//2), dPhaseValues(4),
        dRotorValues(7//2),
    ]
    calculator_functions = [
        DCalcValues(3, full(3), 1),
        DCalcValues(4, (2, -1, 3, -2), 2),
        DCalcValues(5//2, full(5//2), 1),
        YCalcValues(4, -2:1, 2),
        YCalcValues(7//2, 1//2, 1),
    ]
    angle_calculator_functions = [
        dCalcValues(4, full(4), 1), dCalcValues(4, (2, -1, 3, -2), 3),
        dCalcValues(7//2, full(7//2), 2), λCalcValues(4, -2:1, 3), λCalcValues(7//2, 1//2, 1),
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
    # Angles at the poles, near them, and away from both, and for a calculator of several
    # angles, triples of them in which each of those comes first.
    angle_points = [
        "β = $β" => [β] for β ∈ (0.0, 1e-7, π / 2, π - 1e-3, Float64(π), 2.2)
    ]
    points(f::Union{dValues, dPhaseValues}) = angle_points
    angle_triples = [
        "(0, π, 2.2)" => [0.0, π, 2.2], "(π, 0.4, 0)" => [π, 0.4, 0.0],
        "(1e-7, π - 1e-3, 1.3)" => [1e-7, π - 1e-3, 1.3],
        "(π - 1e-3, 1e-7, 2.2)" => [π - 1e-3, 1e-7, 2.2], "(0.4, 1.3, 2.2)" => [0.4, 1.3, 2.2],
    ]
    points(f::Union{dCalcValues, λCalcValues}) =
        [name => x[1:nrotors(f)] for (name, x) ∈ angle_triples]
    # A quaternion's β has no derivative at the poles, and neither do some elements of its
    # `d` (see "Derivatives: d of a rotor at the poles"), so its points are near them.
    points(f::dRotorValues) = test_points(largest_ℓ(f))[4:6]

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
        HalfOddInteger, angle_generators, angle_cotangent, zero_cotangents, AngleCotangents,
        rotation_angle, rotation_angle_gradient, rotation_angle_cotangent
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
        (0, 0:0, 0:0, 0:0, 0:0),
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
        # The same for a real block of `d`, along the generators of its angles
        a = randn(rng, Nᵣ, length(rows), length(cols))
        β̇ = randn(rng, Nᵣ)
        ȧ = similar(a, Nᵣ, length(outrows), length(outcols))
        wigner_block_pushforward!((x, ẋ) -> only(ẋ), ȧ, a, ℓ, rows, cols, outrows, outcols, left, angle_generators(β̇), Val(1))
        ā = randn(rng, size(ȧ))
        Ḡ = wigner_block_pullback!(zero_cotangents(Float64, Nᵣ), a, ā, ℓ, rows, cols, outrows, outcols, left)
        @test sum(ā .* ȧ) ≈ sum(angle_cotangent(Ḡ, i) * β̇[i] for i ∈ 1:Nᵣ) rtol=1e-13 atol=1e-14
    end
    for Nᵣ ∈ (1, 3), (ℓ, nₛ) ∈ ((4, 1), (4, 3), (h(5//2), 2))
        n = Int(2ℓ) + 1
        A = randn(rng, ComplexF64, Nᵣ, nₛ, n)
        G = randn(rng, 3, Nᵣ)
        Ȧ = similar(A)
        harmonic_block_pushforward!((x, ẋ) -> only(ẋ), Ȧ, A, ℓ, G, Val(1))
        Ā = randn(rng, ComplexF64, size(Ȧ))
        Ḡ = harmonic_block_pullback!(zeros(3, Nᵣ), A, Ā, ℓ)
        @test sum(real(conj.(Ā) .* Ȧ)) ≈ sum(Ḡ .* G)
        a = randn(rng, Nᵣ, nₛ, n)
        β̇ = randn(rng, Nᵣ)
        ȧ = similar(a)
        harmonic_block_pushforward!((x, ẋ) -> only(ẋ), ȧ, a, ℓ, angle_generators(β̇), Val(1))
        ā = randn(rng, size(ȧ))
        Ḡ = harmonic_block_pullback!(zero_cotangents(Float64, Nᵣ), a, ā, ℓ)
        @test sum(ā .* ȧ) ≈ sum(angle_cotangent(Ḡ, i) * β̇[i] for i ∈ 1:Nᵣ)
    end

    # The gradient of the angle of a phase or of a rotor is that of `rotation_angle`, and is
    # not finite at the poles, where a zero cotangent of the angle still gives zero.
    z = 1.3cis(0.7)
    gz = rotation_angle_gradient(z)
    @test [real(gz), imag(gz)] ≈ ForwardDiff.gradient(v -> rotation_angle(Complex(v...)), [real(z), imag(z)])
    q = 1.7 .* [0.3, -0.5, 0.7, 0.2]
    @test collect(rotation_angle_gradient(quaternionof(q))) ≈ ForwardDiff.gradient(v -> rotation_angle(quaternionof(v)), q)
    for R ∈ (Quaternion(1.0, 0.0, 0.0, 0.0), Quaternion(0.0, 0.0, 1.0, 0.0))
        @test !all(isfinite, rotation_angle_gradient(R))
        @test all(iszero, rotation_angle_cotangent(R, 0.0))
    end
    @test_throws DimensionMismatch AngleCotangents(zeros(1, 2), zeros(1, 3))

    # `rotor_cotangent` is the adjoint of `rotor_generator` from either side, for rotors of
    # any norm, and the cotangent of a function of R/‖R‖ is orthogonal to R.
    for _ ∈ 1:10, left ∈ (true, false)
        R, Ṙ, g = randn(rng, 4), randn(rng, 4), randn(rng, 3)
        R̄ = rotor_cotangent(left, R, g)
        @test sum(g .* rotor_generator(left, R, Ṙ)) ≈ sum(R̄ .* Ṙ)
        @test abs(sum(R̄ .* R)) ≤ 16eps() * sum(abs, R̄ .* R)
    end
end


@testitem "Derivatives: the kernels allocate nothing" setup=[DerivativeTools] begin
    import SphericalFunctions: wigner_block_pushforward!, harmonic_block_pushforward!,
        angle_generators
    # The kernels are specialized on the function that combines a value with its
    # derivatives, as the rules of Enzyme, Mooncake, and ChainRules pass it, so that it is
    # not called dynamically.  Before Julia 1.12, a real block with the generators of an
    # angle allocates 16 bytes per call, however many rotations there are.
    combine = (x, ẋ) -> only(ẋ)
    wigner(f, Ȧ, A, G) = @allocated wigner_block_pushforward!(f, Ȧ, A, 4, -4:4, -4:4, -4:4, -4:4, true, G, Val(1))
    harmonic(f, Ȧ, A, G) = @allocated harmonic_block_pushforward!(f, Ȧ, A, 4, G, Val(1))
    for Nᵣ ∈ (1, 3)
        for (A, G) ∈ ((randn(ComplexF64, Nᵣ, 9, 9), randn(3, Nᵣ)), (randn(Nᵣ, 9, 9), angle_generators(randn(Nᵣ))))
            Ȧ = similar(A)
            wigner(combine, Ȧ, A, G)
            @test wigner(combine, Ȧ, A, G) == 0 skip=(VERSION < v"1.12" && eltype(A) <: Real)
        end
        for (A, G) ∈ ((randn(ComplexF64, Nᵣ, 2, 9), randn(3, Nᵣ)), (randn(Nᵣ, 2, 9), angle_generators(randn(Nᵣ))))
            Ȧ = similar(A)
            harmonic(combine, Ȧ, A, G)
            @test harmonic(combine, Ȧ, A, G) == 0 skip=(VERSION < v"1.12" && eltype(A) <: Real)
        end
    end
end


@testitem "Derivatives: D and sYlm label the same values as the calculators" setup=[DerivativeTools] begin
    import SphericalFunctions: recurrence!
    # `D` of a rotor is computed through `D_array`, which the rules for ChainRules and
    # ReverseDiff attach to, and labelled by `D_series`.  Its blocks are those a calculator
    # returns, bit for bit, over views of one vector, and are of one type, whether or not
    # the derivatives of the block need rows or columns beyond it.
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
            @test typeof(S[ℓ]) == typeof(D(R, n)[ℓ])
            @test parent(S[ℓ]) isa SubArray{ComplexF64, 1, Vector{ComplexF64}}
        end
        s = d(0.7, n; limits...)
        for (ℓ, b) ∈ dCalculator(0.7, n; limits...)
            @test s[ℓ] == b
            @test typeof(s[ℓ]) == typeof(d(0.7, n)[ℓ])
        end
    end
end


@testitem "Derivatives: calculators of dual numbers allocate nothing when stepped" setup=[DerivativeTools] begin
    import SphericalFunctions: set_R!, set_β!, set_θ!
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
            s += float(ForwardDiff.value(ForwardDiff.value(real(first(b)))))
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
    # The same for calculators of `d` and of ₛλₗₘ, of angles that are dual numbers, single
    # and nested
    βd = dual(0.7, 1)
    βdd = ForwardDiff.Dual{ForwardDiff.Tag{Int, typeof(βd)}}(βd, βd)
    set_angles!(c::SphericalFunctions.dCalculator, β) = (set_β!(c, β); nothing)
    set_angles!(c, θ) = (set_θ!(c, θ); nothing)
    for β ∈ (βd, βdd), c ∈ (
        dCalculator(β, 8), dCalculator([β, β], 8; m′ₘₐₓ=2, mₘₐₓ=3), sλlmCalculator(β, 8, -2),
        sλlmCalculator([β, β], 8, -2:2),
    )
        x = SphericalFunctions.isbatched(c) ? [β, β] : β
        sweep(c); set_angles!(c, x)
        @test sweep_allocations(c) == 0
        @test @allocated(set_angles!(c, x)) == 0
    end
end


@testitem "Derivatives: ReverseDiff's angles give one calculator type" setup=[DerivativeTools] begin
    import ReverseDiff
    import SphericalFunctions: floattype, set_θ!, set_β!
    # ReverseDiff records in the type of each of its numbers where the number came from, but
    # a calculator stores its rotor data in one type, that of an element of a tracked
    # `Vector`.  So an angle that is an input of the tape and one computed on it, each alone
    # or in a vector, give a calculator of that one type, and a calculator built from one of
    # them accepts any other through `set_θ!`.
    x = ReverseDiff.track([0.3, 1.1])
    types = (
        floattype(x[1]), floattype(2x[1]), floattype([x[1]]), floattype([2x[1]]),
        floattype(typeof(2x[1])), floattype(sYlmCalculator(2x[1], 4, -1)),
        floattype(dCalculator(2x[1], 4)), floattype(sλlmCalculator(2x[1], 4, -1)),
    )
    @test all(==(first(types)), types)
    total(c) = sum(Y -> sum(z -> real(z) + 2imag(z), array_view(Y)), values(c))
    function f(x)
        c = sYlmCalculator(2x[1], 4, -1)
        s = total(c)
        for θ ∈ (x[2], 3x[2])
            set_θ!(c, θ)
            s += total(c)
        end
        s
    end
    @test ReverseDiff.gradient(f, [0.3, 1.1]) ≈ ForwardDiff.gradient(f, [0.3, 1.1])
    function g(x)
        c = dCalculator(2x[1], 4)
        s = total(c)
        for β ∈ (x[2], 3x[2])
            set_β!(c, β)
            s += total(c)
        end
        s
    end
    @test ReverseDiff.gradient(g, [0.3, 1.1]) ≈ ForwardDiff.gradient(g, [0.3, 1.1])
end


@testitem "Derivatives: ChainRules rules" setup=[DerivativeTools] begin
    import ChainRulesCore
    using ChainRulesCore: NoTangent, ZeroTangent, Tangent, @thunk
    import SphericalFunctions: D_array, d_array, sYlm_array, sYlm_matrix_array,
        wigner_arrays_with_derivative_values, wigner_arrays_pushforward, derivatives_from_left,
        rotor_generator, rotor_generators, harmonic_array_pushforward
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
    _, derivative_vals, calc = wigner_arrays_with_derivative_values(ComplexF64, R, limits...)
    Ȧ = wigner_arrays_pushforward(
        calc, derivative_vals,
        reshape(collect(rotor_generator(derivatives_from_left(calc), R, q̇)), 3, 1)
    )
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

    # `d_array` of an angle, a phase, and a rotor: the tangent is ForwardDiff's, the
    # cotangent is adjoint to it, and that of a phase is tangent to the unit circle, as that
    # of a rotor is orthogonal to the rotor.
    dlimits = (4, 2, -1, 3, -2)
    for (x, ẋ, path) ∈ (
        (0.7, 0.3, t -> 0.7 + 0.3t),
        (cis(0.7), 0.4 + 0.9im, t -> cis(0.7) + (0.4 + 0.9im) * t),
        (Quaternion(q...), Quaternion(q̇...), t -> Quaternion((q .+ t .* q̇)...)),
    )
        local Ω, Ω̇, back
        Ω, Ω̇ = ChainRulesCore.frule((NoTangent(), ẋ, ntuple(_ -> NoTangent(), 5)...), d_array, x, dlimits...)
        @test Ω == d_array(x, dlimits...)
        forward = ForwardDiff.derivative(t -> reduce(vcat, vec.(d_array(path(t), dlimits...))), 0.0)
        @test flat(Ω̇) ≈ forward
        Ω, back = ChainRulesCore.rrule(d_array, x, dlimits...)
        Ω̄ = [randn(size(b)) for b ∈ Ω]
        x̄ = back(Ω̄)[2]
        @test sum(flat(Ω̄) .* flat(Ω̇)) ≈ (x isa Real ? x̄ * ẋ : x isa Complex ? real(conj(x̄) * ẋ) : sum(x̄[i] * ẋ[i] for i ∈ 1:4))
        if x isa Complex
            @test abs(real(conj(x̄) * x)) ≤ 16eps() * abs(x̄)
        elseif x isa Quaternion
            @test abs(sum(x̄[i] * x[i] for i ∈ 1:4)) ≤ 64eps() * sqrt(sum(abs2, x̄[i] for i ∈ 1:4))
        end
        @test back(ZeroTangent())[2] isa ZeroTangent
        @test back(@thunk(Ω̄))[2] == x̄
    end
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


@testitem "Derivatives: d and ₛλₗₘ, with every backend" setup=[ExplicitWignerMatrices, DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoReverseDiff, AutoEnzyme,
        AutoMooncake, AutoMooncakeForward, AutoZygote
    using DifferentiationInterfaceTest: Scenario, test_differentiation
    import Enzyme, Mooncake, ReverseDiff, Zygote

    # `d` of angles, phases, and rotors, and loops over the calculators of `d` and of ₛλₗₘ,
    # whole and restricted, with integer and half-integer indices, of one angle and of
    # several, at and near the poles.  Every backend reaches the rules, as for 𝔇 and ₛYₗₘ;
    # Zygote, which cannot follow a calculator, is tested on the arrays alone.  Enzyme's
    # reverse mode is tested on the calculators alone, and not on Julia 1.10, as for the
    # calculators of 𝔇 and ₛYₗₘ: on the arrays, its batched Jacobians of `d` intermittently
    # crash Enzyme 0.13.209 in the garbage collector.
    backends = [
        AutoForwardDiff(),
        AutoReverseDiff(),
        AutoReverseDiff(compile=true),
        AutoEnzyme(mode=Enzyme.Forward),
        AutoMooncake(config=nothing),
        AutoMooncakeForward(config=nothing),
    ]
    calculator_backends = copy(backends)
    VERSION ≥ v"1.11" && insert!(calculator_backends, 5, AutoEnzyme(mode=Enzyme.Reverse))
    scenarios(functions) = [
        Scenario{:jacobian, :out}(
            f, x; res1=reference_jacobian(f, x), prep_args=(; x=copy(x), contexts=()),
            name="$f at $point"
        )
        for f ∈ functions for (point, x) ∈ points(f)
    ]
    test_differentiation(
        [backends; AutoZygote()], scenarios(angle_array_functions); correctness=true,
        isapprox=within, atol=0, rtol=16eps(), detailed=true, testset_name="Jacobians of d"
    )
    test_differentiation(
        calculator_backends, scenarios(angle_calculator_functions); correctness=true,
        isapprox=within, atol=0, rtol=16eps(), detailed=true,
        testset_name="Jacobians through calculators of d and ₛλₗₘ"
    )
end


@testitem "Derivatives: ReverseDiff's recorded tapes" setup=[ExplicitWignerMatrices, DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoReverseDiff, prepare_gradient, gradient
    import ReverseDiff
    import SphericalFunctions: recurrence!, set_θ!
    using Random: Xoshiro

    # A tape that ReverseDiff records at one point and compiles gives the gradient at
    # others, for every function of whole arrays and every loop over a calculator, as well
    # as for loops that compute blocks out of order, run twice, stop early, interleave two
    # calculators, or set new rotor data partway through.
    struct Twice end
    (::Twice)(x) = (c = dCalculator(x, 4); [[sum(b) for (_, b) ∈ c]; [2sum(b) for (_, b) ∈ c]])
    struct Interleaved end
    function (::Interleaved)(x)
        c, y = DCalculator(quaternionof(x), 3), sλlmCalculator(x[2], 3, 0)
        [real(a[0, 0]) + b[0] for ((_, a), (_, b)) ∈ zip(c, y)]
    end
    struct EarlyExit end
    function (::EarlyExit)(x)
        z = eltype(x)[]
        for (ℓ, b) ∈ sYlmCalculator(rotorof(x), 6, -1)
            ℓ > 2 && break
            push!(z, real(b[0]), imag(b[ℓ]))
        end
        z
    end
    struct OutOfOrder end
    function (::OutOfOrder)(x)
        c = dCalculator(x[1:2], 5)
        [sum(recurrence!(c, 4)), sum(recurrence!(c, 1)), sum(recurrence!(c, 5)), sum(recurrence!(c, 5))]
    end
    struct Reset end
    function (::Reset)(x)
        c = sλlmCalculator(x[1], 4, -1)
        z = [b[0] for (_, b) ∈ c]
        set_θ!(c, 2x[2])
        [z; [b[0] for (_, b) ∈ c]]
    end
    rng = Xoshiro(2468)
    backend = AutoReverseDiff(compile=true)
    for f ∈ [
        DValues(4), YValues(4, -1, 1), YMatrixValues(3, -1, 1, 2), dValues(4), dPhaseValues(4),
        dRotorValues(7//2), DCalcValues(3, full(3), 1), DCalcValues(4, (2, -1, 3, -2), 2),
        YCalcValues(4, -2:1, 2), YCalcValues(7//2, 1//2, 1), angle_calculator_functions...
    ]
        x = [last(p) for p ∈ points(f)[end-2:end]]
        w = randn(rng, length(f(x[1])))
        h(v) = sum(w .* f(v))
        prep = prepare_gradient(h, backend, x[1])
        for v ∈ x
            @test within(gradient(h, prep, backend, v), ForwardDiff.gradient(h, v); rtol=1e-13)
        end
    end
    for f ∈ (Twice(), Interleaved(), EarlyExit(), OutOfOrder(), Reset())
        x = f isa Union{Interleaved, EarlyExit} ?
            [[0.3, -0.5, 0.7, 0.2], [0.9, 0.1, -0.2, 0.4], [-0.2, 0.6, 0.6, -0.1]] :
            [[0.4, 1.3], [0.9, 2.0], [2.5, 0.1]]
        w = randn(rng, length(f(x[1])))
        h(v) = sum(w .* f(v))
        prep = prepare_gradient(h, backend, x[1])
        for v ∈ x
            @test within(gradient(h, prep, backend, v), ForwardDiff.gradient(h, v); rtol=1e-13)
        end
    end
end


@testitem "Derivatives: d of a rotor at the poles" setup=[DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoReverseDiff, AutoMooncake,
        AutoZygote, jacobian
    import Mooncake, ReverseDiff, Zygote

    # At a pole, β of a rotor has no derivative, and neither have the elements of `d` with
    # m′ - m = ±1 at β = 0, or m′ + m = ±1 at β = π, which are `NaN`; the others have the
    # derivative zero there.
    labels = [(m′, m) for ℓ ∈ 0:4 for m ∈ -ℓ:ℓ for m′ ∈ -ℓ:ℓ]
    f = dRotorValues(4)
    for (v, β) ∈ (([1.0, 0, 0, 0], 0), ([0.6, 0, 0, 0.8], 0), ([0.0, 0, 1, 0], π), ([0.0, 0.6, 0.8, 0], π))
        singular = [abs(β == 0 ? m′ - m : m′ + m) == 1 for (m′, m) ∈ labels]
        @test count(singular) == 40
        for backend ∈ (AutoForwardDiff(), AutoReverseDiff(), AutoZygote(), AutoMooncake(config=nothing))
            J = jacobian(f, backend, v)
            @test all(i -> any(isnan, J[i, :]), findall(singular))
            @test all(iszero, J[.!singular, :])
        end
    end
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


@testitem "Derivatives: array_view in a loop, under Enzyme's batched reverse mode" setup=[DerivativeTools] tags=[:derivatives] begin
    using DifferentiationInterface: AutoForwardDiff, AutoEnzyme, jacobian
    import Enzyme

    # Each element of each block of a calculator, read through `array_view(b)[i]` in the
    # loop over its blocks, under Enzyme's batched reverse mode, which
    # DifferentiationInterface's `jacobian` uses.  This is the pattern of
    # EnzymeAD/Enzyme#3316, which crashes Enzyme as it compiles the loop when a branch in
    # the loop's body leads Enzyme's optimizer to peel the first iteration; the extension
    # for Enzyme keeps the branches of `array_view` out of the loop (see the note on
    # `check_storage` there).
    struct ViewReads{T, S}
        n::T
        s::S
    end
    function (f::ViewReads)(v)
        c = sYlmCalculator(rotorof(v), f.n, f.s)
        z = zeros(eltype(eltype(values(c))), Ysize(f.n))
        k = 0
        for (ℓ, b) ∈ c, m ∈ -ℓ:ℓ
            z[k += 1] = array_view(b)[Int(m + ℓ) + 1]
        end
        realified(z)
    end
    v = [0.3, 0.4, -0.5, 0.6]
    for f ∈ (ViewReads(4, 1), ViewReads(7//2, 1//2))
        J = jacobian(f, AutoEnzyme(mode=Enzyme.Reverse), v)
        @test within(J, jacobian(f, AutoForwardDiff(), v); rtol=16eps())
    end
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
        scenarios([array_functions[[1, 4, 5, 7, 8]]; dValues(4)]); correctness=true, isapprox=within,
        atol=0, rtol=64eps(), detailed=true, testset_name="Hessians"
    )
    test_differentiation(
        [AutoForwardDiff()],
        scenarios([calculator_functions[[1, 3, 5]]; dCalcValues(4, full(4), 1); λCalcValues(4, -2:1, 3)]);
        correctness=true,
        isapprox=within, atol=0, rtol=64eps(), detailed=true,
        testset_name="Hessians through calculators"
    )
end
