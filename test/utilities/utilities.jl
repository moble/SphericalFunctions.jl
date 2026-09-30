# Helpers shared by many test items and by the literature-comparison pages: ranges of mode
# indices, samples of angles, directions and rotors, and a closed-form ₛYₗₘ that owes
# nothing to the package.  They are defined in a `@testmodule`, which is evaluated once in
# each test process rather than once in every item that uses it.  The module exports
# nothing, so that no name it defines can shadow one defined in an item; each item imports
# the helpers it uses by name, which also shows where they come from.
#
# Every sampling helper takes the random-number generator as its first argument.  Each item
# creates its own, as in `rng = Random.Xoshiro(1234)`, so that its samples are the same
# whichever items have already run in the same process.
@testmodule Utilities begin
    import Random: AbstractRNG
    import Quaternionic
    import Quaternionic: QuatVec, Rotor, 𝐢, 𝐣, 𝐤

    ℓmrange(ℓₘᵢₙ, ℓₘₐₓ) = [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]
    ℓmrange(ℓₘₐₓ) = ℓmrange(0, ℓₘₐₓ)
    function sℓmrange(ℓₘₐₓ, sₘₐₓ)
        sₘₐₓ = min(abs(sₘₐₓ), ℓₘₐₓ)
        [
            (s, ℓ, m)
            for s in -sₘₐₓ:sₘₐₓ
            for ℓ in abs(s):ℓₘₐₓ
            for m in -ℓ:ℓ
        ]
    end
    function ℓm′mrange(ℓₘₐₓ)
        [
            (ℓ, m′, m)
            for ℓ in 0:ℓₘₐₓ
            for m′ in -ℓ:ℓ
            for m in -ℓ:ℓ
        ]
    end

    # `n` samples drawn uniformly from the interval [a, b].  They are computed in `T`
    # itself, so that any type `rand` supports can be sampled, `Double64` and `BigFloat`
    # included, and they are clamped to the interval, so that rounding cannot put one
    # outside it.
    uniform_samples(rng::AbstractRNG, ::Type{T}, a, b, n) where {T} =
        clamp.(T(a) .+ (T(b) - T(a)) .* rand(rng, T, n), T(a), T(b))

    # Samples of an azimuthal angle in [0, 2π]: the ends, 0, π and 2π together with their
    # neighbors, and `n÷2` random values in each half of the range.
    αrange(rng::AbstractRNG, ::Type{T}=Float64, n=15) where {T} = T[
        0; nextfloat(T(0)); uniform_samples(rng, T, 0, T(π), n÷2); prevfloat(T(π)); T(π);
        nextfloat(T(π)); uniform_samples(rng, T, T(π), 2T(π), n÷2); prevfloat(2T(π)); 2T(π)
    ]
    # Samples of a polar angle in [0, π]: both poles, their neighbors, and `n` random values
    # between them.  `avoid_poles` keeps every sample, random ones included, at least that
    # far from either pole; the formulas of several reference pages are singular there.  The
    # edges of the band are computed in `T`, so that the sample just inside the upper edge
    # is distinct from the edge itself in every precision.
    βrange(rng::AbstractRNG, ::Type{T}=Float64, n=15; avoid_poles=0) where {T} = T[
        T(avoid_poles); nextfloat(T(avoid_poles));
        uniform_samples(rng, T, T(avoid_poles), T(π)-T(avoid_poles), n);
        prevfloat(T(π)-T(avoid_poles)); T(π)-T(avoid_poles)
    ]
    γrange(rng::AbstractRNG, ::Type{T}=Float64, n=15) where {T} = αrange(rng, T, n)
    αβγrange(rng::AbstractRNG, ::Type{T}=Float64, n=15; avoid_poles=0) where {T} =
        vec(collect(Iterators.product(
            αrange(rng, T, n), βrange(rng, T, n; avoid_poles), γrange(rng, T, n)
        )))

    # The spherical coordinates (θ, ϕ) are sampled like the Euler angles (β, α)
    const θrange = βrange
    θϕrange(rng::AbstractRNG, ::Type{T}=Float64, n=15; avoid_poles=0) where {T} = vec(collect(
        Iterators.product(θrange(rng, T, n; avoid_poles), αrange(rng, T, n))
    ))

    # Unit vectors: the axes, their negatives, and `n` random directions
    v̂range(rng::AbstractRNG, ::Type{T}=Float64, n=15) where {T} = QuatVec{T}[
        𝐢; 𝐣; 𝐤;
        -𝐢; -𝐣; -𝐤;
        Quaternionic.normalize.(randn(rng, QuatVec{T}, n))
    ]
    # Rotors: the identity, and the rotations by π and by ±π/2 about each axis, each with both
    # signs, followed by `n` random rotors
    function Rrange(rng::AbstractRNG, ::Type{T}=Float64, n=15) where {T}
        invsqrt2 = inv(√T(2))
        [
            [
                # `sign*R` promotes to `Quaternion`; these are rotations, so say so.
                Rotor{T}(sign*R)
                for R in [
                    Rotor{T}(1);
                    [Rotor{T}(𝐯) for 𝐯 in (𝐢,𝐣,𝐤)];
                    [Rotor{T}(invsqrt2 + invsqrt2*𝐯) for 𝐯 in (𝐢,𝐣,𝐤)];
                    [Rotor{T}(invsqrt2 - invsqrt2*𝐯) for 𝐯 in (𝐢,𝐣,𝐤)]
            ]
            for sign in (1,-1)
            ];
            randn(rng, Rotor{T}, n)
        ]
    end

    """
        array_equal(a1, a2, equal_nan=false)

    Ensure that arrays have same types and shapes, and all elements are the same.  If
    `equal_nan` is `true`, NaNs in the same place in each array will be considered to be
    equal.

    Note that this is slightly stricter than the numpy version of this function, because
    arrays of different type will not be considered equal.

    """
    function array_equal(a1::T1, a2::T2, equal_nan=false) where {T1, T2}
        if T1 !== T2 || size(a1) != size(a2)
            return false
        end
        all(e->e[1]==e[2] || (equal_nan && isnan(e[1]) && isnan(e[2])), zip(a1, a2))
    end

    """
        sYlm_closed_form(s, ℓ, m, θ, ϕ)

    The spin-weighted spherical harmonic ₛYₗₘ(θ, ϕ), evaluated from the explicit sum over
    factorials given on the conventions pages.  It shares no code with the package, which is
    what makes it a reference for the package's harmonics.
    """
    function sYlm_closed_form(s::Int, ell::Int, m::Int, theta::T, phi::T) where {T<:Real}
        # Eqs. (II.7) and (II.8) of https://arxiv.org/abs/0709.0093v3 [AjithEtAl_2011](@cite)
        # Note their weird definition w.r.t. `-s`
        k_min = max(0, m + s)
        k_max = min(ell + m, ell + s)
        sin_half_theta, cos_half_theta = sincos(theta / 2)
        return (-1)^(-s) * sqrt((2 * ell + 1) / (4 * T(π))) *
            T(sum(
                (-1) ^ (k)
                * sqrt(factorial(big(ell + m)) * factorial(big(ell - m)) * factorial(big(ell - s)) * factorial(big(ell + s)))
                * (cos_half_theta ^ (2 * ell + m + s - 2 * k))
                * (sin_half_theta ^ (2 * k - s - m))
                / (factorial(big(ell + m - k)) * factorial(big(ell + s - k)) * factorial(big(k)) * factorial(big(k - s - m)))
                for k in k_min:k_max
            )) *
            cis(m * phi)
    end

    """
        sYlm_closed_form_pixels(s, ℓ, m, pixels)

    The same closed form as `sYlm_closed_form`, evaluated on a whole list of `(θ, ϕ)`
    pixels, with the pixel-independent factorials hoisted out of the loop.  That is about 40
    times faster, which is what makes pixel-by-pixel comparisons against a transform
    affordable; it is the same transcription of the conventions-page formula, so it owes
    nothing to the package.

    Checked against `sYlm_closed_form` itself in the "SSHT synthesis" test item.
    """
    function sYlm_closed_form_pixels(
        s::Int, ℓ::Int, m::Int, p::AbstractVector{<:AbstractVector{T}}
    ) where {T}
        kmin, kmax = max(0, m + s), min(ℓ + m, ℓ + s)
        𝒩 = sqrt(
            factorial(big(ℓ + m)) * factorial(big(ℓ - m))
            * factorial(big(ℓ - s)) * factorial(big(ℓ + s))
        )
        c = T[
            (-1)^k * 𝒩 / (
                factorial(big(ℓ + m - k)) * factorial(big(ℓ + s - k))
                * factorial(big(k)) * factorial(big(k - s - m))
            )
            for k in kmin:kmax
        ]
        pre = T(-1)^(-s) * sqrt((2ℓ + 1) / (4 * T(π)))
        map(p) do θϕ
            sθ, cθ = sincos(θϕ[1] / 2)
            pre * sum(
                c[k-kmin+1] * cθ^(2ℓ + m + s - 2k) * sθ^(2k - s - m)
                for k in kmin:kmax
            ) * cis(m * θϕ[2])
        end
    end

    # The Levi-Civita symbol
    ε(j,k,l) = ifelse(
        (j,k,l)∈((1,2,3),(2,3,1),(3,1,2)),
        1,
        ifelse(
            (j,k,l)∈((2,1,3),(1,3,2),(3,2,1)),
            -1,
            0
        )
    )

end  # @testmodule Utilities


# Two checks of what the compiler makes of a call, which work alike on every Julia version
# the package supports.  `Base.infer_return_type` exists only from Julia 1.11 on; on 1.10
# the same answer comes from `Core.Compiler.return_type`, which is what `Base.promote_op`
# uses.  The printed optimized code marks a dynamically dispatched call with the word
# "dynamic" only from Julia 1.12 on, so `dynamic_calls` reads the code itself: a call that
# inference has resolved is an `:invoke`, or a `:call` of a builtin function, and any other
# `:call` is dispatched at run time.
@testmodule InferenceChecks begin
    inferred_type(f, types) = @static if isdefined(Base, :infer_return_type)
        Base.infer_return_type(f, types)
    else
        Core.Compiler.return_type(f, Tuple{types...})
    end

    builtin(f::GlobalRef) =
        isdefined(f.mod, f.name) && getglobal(f.mod, f.name) isa Core.Builtin
    builtin(f) = f isa Core.Builtin
    # Calls that only construct an exception, on a branch that throws it.  Julia 1.10, under
    # `Pkg.test` with bounds checking forced on, leaves the `BoundsError(A, i)` of each
    # bounds check as a `:call`; these say nothing about how the code runs when it does not
    # throw.
    exception_type(f::GlobalRef) = isdefined(f.mod, f.name) && exception_type(getglobal(f.mod, f.name))
    exception_type(f) = f isa Type && f <: Exception
    # The statements of the optimized code that are dispatched at run time, and their number
    dynamic_call_list(f, types) = filter(
        ex -> Meta.isexpr(ex, :call) && !builtin(ex.args[1]) && !exception_type(ex.args[1]),
        only(Base.code_typed(f, types; optimize=true)).first.code
    )
    dynamic_calls(f, types) = length(dynamic_call_list(f, types))
end

@testitem "Utilities: the sampling helpers and array_equal" setup=[Utilities] begin
    import .Utilities: αrange, βrange, γrange, αβγrange, θϕrange, v̂range, Rrange, array_equal
    using Random
    using DoubleFloats: Double64
    using Quaternionic: QuatVec, Rotor

    rng = Random.Xoshiro(3)

    # `avoid_poles` excludes a band around each pole from the random samples as well as from
    # the fixed endpoints.  Excluding it from the endpoints alone would put about one sample
    # in 300 calls inside the band.
    for T ∈ (Float64, Float32, Float16), _ ∈ 1:300
        β = βrange(rng, T, 7; avoid_poles=1e-3)
        @test all(T(1e-3) ≤ b ≤ T(π) - T(1e-3) for b ∈ β)
    end
    # The edges of the band are computed in `T`, so the sample just inside the upper edge is
    # distinct from the edge itself in every precision
    for T ∈ (Float64, Float32, Float16)
        β = βrange(rng, T, 3; avoid_poles=1e-3)
        @test β[end-1] == prevfloat(β[end]) < β[end] == T(π) - T(1e-3)
        @test β[1] == T(1e-3) && β[2] == nextfloat(β[1])
    end

    # The helpers draw from the generator they are given and from no other, so the same seed
    # gives the same samples, and the global generator is left as it was; an item's samples
    # therefore do not depend on which items ran before it in the same process
    global_state = copy(Random.default_rng())
    for f ∈ (αrange, βrange, γrange, αβγrange, θϕrange, v̂range, Rrange)
        @test f(Random.Xoshiro(5), Float64, 6) == f(Random.Xoshiro(5), Float64, 6)
        @test f(Random.Xoshiro(5), Float64, 6) != f(Random.Xoshiro(6), Float64, 6)
    end
    @test copy(Random.default_rng()) == global_state

    # A wider type is sampled in that type, not converted from `Float64` samples
    for T ∈ (Double64, BigFloat)
        β = βrange(rng, T, 20; avoid_poles=1e-3)
        @test eltype(β) === T
        @test length(β) == 24
        @test all(T(1e-3) ≤ b ≤ T(π) - T(1e-3) for b ∈ β)
        @test allunique(β)
        @test all(b -> T(Float64(b)) != b, β[3:end-2])
        α = αrange(rng, T, 20)
        @test eltype(α) === T
        @test all(0 ≤ x ≤ 2T(π) for x ∈ α)
        @test count(x -> 0 < x < T(π), α) ≥ 10 && count(x -> T(π) < x < 2T(π), α) ≥ 10
        @test eltype(αβγrange(rng, T, 3)) === NTuple{3, T}
        @test eltype(θϕrange(rng, T, 3)) === NTuple{2, T}
        @test eltype(v̂range(rng, T, 3)) === QuatVec{T}
        @test eltype(Rrange(rng, T, 3)) === Rotor{T}
    end

    # Both edges below π and 2π are sampled, once each
    for T ∈ (Float64, Float32)
        α = αrange(rng, T, 4)
        @test count(==(prevfloat(T(π))), α) == 1
        @test count(==(prevfloat(2T(π))), α) == 1
        @test extrema(α) == (0, 2T(π))
    end

    # `array_equal` compares type, shape and elements, with NaNs equal only on request
    @test array_equal([1.0, NaN], [1.0, NaN], true)
    @test !array_equal([1.0, NaN], [1.0, NaN])
    @test !array_equal([1.0, NaN], [2.0, NaN], true)
    @test !array_equal([1.0], [1.0f0])
    @test !array_equal([1.0], [1.0, 1.0])
end


@testitem "Utilities: InferenceChecks" setup=[InferenceChecks] begin
    import .InferenceChecks: inferred_type, dynamic_calls

    # `inferred_type` agrees with inference on every supported Julia version, including a
    # call that cannot return
    @test inferred_type(x -> 2x + 1, (Int,)) === Int
    @test inferred_type(x -> x / 2, (Int,)) === Float64
    @test inferred_type(x -> throw(ArgumentError("no")), (Int,)) === Union{}
    @test inferred_type(r -> r[], (Base.RefValue{Any},)) === Any

    # `dynamic_calls` finds a call that is dispatched at run time, and none in the same code
    # once the argument's type is known, so that a count of zero means something
    f(x) = x + 1
    g(r) = f(r[])
    @test dynamic_calls(g, (Base.RefValue{Any},)) ≥ 1
    @test dynamic_calls(g, (Base.RefValue{Int},)) == 0
    @test dynamic_calls(x -> sum(x), (Vector{Float64},)) == 0
end
