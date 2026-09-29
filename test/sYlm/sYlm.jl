# Tests of the spin-weighted spherical harmonics layer: `sYlmCalculator`, `sYlm`, `sYlm!`,
# and `sYlm_matrix`, against the closed-form expression for ₛYₗₘ and the defining relation
# to the Wigner 𝔇 matrices.

@testitem "sYlm vs the closed form, for arbitrary rotors" setup=[Utilities] begin
    import SphericalFunctions
    import SphericalFunctions: sYlm, Yindex, Ysize
    import .Utilities: Rrange, sYlm_closed_form
    using Quaternionic: Quaternion, Rotor, components, 𝐢, 𝐣, 𝐤
    using DoubleFloats: Double64
    using Random

    rng = Random.Xoshiro(17)

    # Reference, category 1 (a closed-form formula, from the conventions pages).
    # `sYlm_closed_form(s, ℓ, m, θ, ϕ)` from the `Utilities` module is the explicit sum of
    # Ajith et al., Eqs. (II.7)-(II.8), evaluated on the sphere — that is, at the rotor
    # `from_spherical_coordinates(θ, ϕ)`, whose Euler angles are (ϕ, θ, 0).  A general rotor
    # has a third angle, and ₛYₗₘ depends on it through the spin weight alone:
    #
    #     ₛYₗₘ(R_{αβγ}) = (-1)^s √((2ℓ+1)/4π) conj(𝔇ˡ_{m,-s}) = ₛYₗₘ(θ=β, ϕ=α) e^{-i s γ}.
    #
    # The angles are taken from the quaternion components (`atan(|Rₐ|, |Rₛ|)` is scale-free
    # and needs no normalization, unlike the `acos` inside `to_euler_angles`), and the
    # closed form is evaluated at four times the working precision and rounded to `T`, so
    # the reference contributes at most half an ulp of its own error.
    function euler_angles(R)
        w, x, y, z = BigFloat.(components(Quaternion(R)))
        ϕₛ, ϕₐ = angle(Complex(w, z)), angle(Complex(y, x))
        (ϕₛ - ϕₐ, 2 * atan(abs(Complex(y, x)), abs(Complex(w, z))), ϕₛ + ϕₐ)
    end
    function Yref(::Type{T}, R, s, ℓₘₐₓ) where {T}
        v = setprecision(BigFloat, 4 * precision(T) + 64) do
            α, β, γ = euler_angles(R)
            [
                sYlm_closed_form(s, ℓ, m, β, α) * cis(-s * γ)
                for ℓ in abs(s):ℓₘₐₓ for m in -ℓ:ℓ
            ]
        end
        Complex{T}.(v)  # in `Yindex(ℓ, m, abs(s))` order
    end

    for T ∈ (Float64, Double64)
        # Measured worst error over everything below: 5.2 eps (Float64), 2.9 eps (Double64).
        ϵ = 20 * eps(T)
        for ℓₘₐₓ ∈ (0, 1, 2, 5, 9)
            for s ∈ -min(3, ℓₘₐₓ):min(3, ℓₘₐₓ)
                for R ∈ Rrange(rng, T, 8)
                    Y = array_view(sYlm(R, ℓₘₐₓ, s))
                    @test eltype(Y) === Complex{T}
                    @test length(Y) == Ysize(abs(s), ℓₘₐₓ)
                    # Accumulated and asserted once per rotor: an engine that is wrong
                    # everywhere would otherwise print one failure per (ℓ, m).
                    errY = maximum(abs, Y .- Yref(T, R, s, ℓₘₐₓ))
                    @test errY ≤ ϵ
                    Y₀ = array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=0))
                    @test length(Y₀) == Ysize(0, ℓₘₐₓ)
                    @test all(iszero, Y₀[1:s^2])
                    @test Y₀[s^2+1:end] == Y
                end
            end
        end
    end
end

@testitem "sYlm vs closed form and definition" setup=[Utilities] begin
    import SphericalFunctions
    import SphericalFunctions: sYlm, Yindex, D
    import .Utilities: θϕrange, sYlm_closed_form
    using Quaternionic
    using Random
    rng = Random.Xoshiro(42)
    ℓₘₐₓ = 5
    for (θ, ϕ) ∈ θϕrange(rng, Float64, 6)
        R = Rotor(from_spherical_coordinates(θ, ϕ))
        for s ∈ -2:2
            Y = array_view(sYlm(R, ℓₘₐₓ, s))
            for ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ
                # Closed form on the sphere
                @test Y[Yindex(ℓ, m, abs(s))] ≈ sYlm_closed_form(s, ℓ, m, θ, ϕ) atol=1e-13 rtol=1e-13
            end
        end
    end
    # The definition ₛYₗₘ = (-1)^s √((2ℓ+1)/4π) conj(𝔇ˡₘ,₋ₛ) for random rotors
    for R ∈ randn(rng, Rotor{Float64}, 5)
        𝔇 = D(R, ℓₘₐₓ)
        for s ∈ -2:2
            Y = array_view(sYlm(R, ℓₘₐₓ, s))
            for ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ
                @test Y[Yindex(ℓ, m, abs(s))] ≈ (-1)^s * √((2ℓ+1)/(4π)) * conj(𝔇[ℓ][m, -s]) atol=1e-14
            end
        end
    end
end

@testitem "sYlm_matrix" setup=[Utilities, RefusalChecks] begin
    import SphericalFunctions
    import SphericalFunctions: sYlm, sYlm_matrix, Ysize
    import .Utilities: sYlm_closed_form
    using Quaternionic: Rotor, Quaternion, components
    using DoubleFloats: Double64
    using Random
    rng = Random.Xoshiro(7)

    # Reference, category 1 (a closed-form formula, from the conventions pages): the same
    # closed form and the same Euler-angle bookkeeping as the "sYlm vs the closed form" item
    # above, evaluated at four times the working precision and rounded.
    function euler_angles(R)
        w, x, y, z = BigFloat.(components(Quaternion(R)))
        ϕₛ, ϕₐ = angle(Complex(w, z)), angle(Complex(y, x))
        (ϕₛ - ϕₐ, 2 * atan(abs(Complex(y, x)), abs(Complex(w, z))), ϕₛ + ϕₐ)
    end
    function Yref(::Type{T}, Rs, s, ℓₘₐₓ) where {T}
        rows = setprecision(BigFloat, 4 * precision(T) + 64) do
            map(Rs) do R
                α, β, γ = euler_angles(R)
                [sYlm_closed_form(s, ℓ, m, β, α) * cis(-s * γ) for ℓ in abs(s):ℓₘₐₓ for m in -ℓ:ℓ]
            end
        end
        Complex{T}[rows[i][j] for i in eachindex(rows), j in eachindex(first(rows))]
    end

    for T ∈ (Float64, Double64)
        # Measured worst error over the grid below: 3.0 eps (Float64), 1.2 eps (Double64).
        ϵ = 20 * eps(T)
        Rs = randn(rng, Rotor{T}, 7)
        for ℓₘₐₓ ∈ (2, 6), s ∈ -2:2
            M = sYlm_matrix(Rs, ℓₘₐₓ, s)
            @test M isa Matrix{Complex{T}}
            @test size(M) == (7, Ysize(abs(s), ℓₘₐₓ))
            for (i, R) ∈ enumerate(Rs)
                @test M[i, :] == array_view(sYlm(R, ℓₘₐₓ, s))
            end
            errM = maximum(abs, M .- Yref(T, Rs, s, ℓₘₐₓ))
            @test errM ≤ ϵ
            M₀ = sYlm_matrix(Rs, ℓₘₐₓ, s; ℓₘᵢₙ=0)
            @test size(M₀) == (7, Ysize(0, ℓₘₐₓ))
            @test all(iszero, M₀[:, 1:s^2])
            @test M₀[:, s^2+1:end] == M
        end
        # Converting the rotors is the way to raise the computation type, and the result is
        # then accurate to the *raised* type.  The components of these rotors are exactly
        # representable in Float64, so one reference serves both computations: measured
        # 5.2e-32 = 1.1 eps(Double64) against the closed form, where the Float64 computation
        # of the very same rotors is off by ~1e-16.
        R64 = Rotor{Float64}.(Rs)
        MD = sYlm_matrix(Rotor{Double64}.(R64), 4, 1)
        @test eltype(MD) === Complex{Double64}
        errMD = maximum(abs, MD .- Yref(Double64, R64, 1, 4))
        @test errMD ≤ 8 * eps(Double64)
    end
    # Errors
    Rs = randn(rng, Rotor{Float64}, 3)
    @test refuses(() -> sYlm_matrix(Rs, 2, 3), ArgumentError, "exceeds ℓₘₐₓ")
    @test refuses(() -> sYlm(Rs[1], 2, 3), ArgumentError, "exceeds ℓₘₐₓ")
    # The keyword may also be spelled `ell_min`
    @test sYlm_matrix(Rs, 4, 1; ell_min=0) == sYlm_matrix(Rs, 4, 1; ℓₘᵢₙ=0)
end

@testitem "sYlmCalculator all spins and batches" begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, sYlm, recurrence!, Yindex
    using Quaternionic: Rotor
    import SphericalFunctions: DegreeBlockBatch, array_view
    using Random
    rng = Random.Xoshiro(11)
    ℓₘₐₓ, sₘₐₓ, Nᵣ = 7, 3, 5
    Rs = randn(rng, Rotor{Float64}, Nᵣ)
    calc = sYlmCalculator(Rs, ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
    @test SphericalFunctions.Nᵣ(calc) == Nᵣ
    @test SphericalFunctions.ℓₘₐₓ(calc) == ℓₘₐₓ
    @test SphericalFunctions.spins(calc) == -sₘₐₓ:sₘₐₓ
    @test_throws MethodError SphericalFunctions.spin(calc)
    # Every spin weight from one calculator, batched, equals the single-rotor results exactly
    singles = Dict((i, s) => array_view(sYlm(Rs[i], ℓₘₐₓ, s; ℓₘᵢₙ=0)) for i ∈ 1:Nᵣ for s ∈ -sₘₐₓ:sₘₐₓ)
    for ℓ ∈ 0:ℓₘₐₓ
        block = recurrence!(calc, ℓ)
        @test SphericalFunctions.ℓ(calc) == ℓ
        for s ∈ -sₘₐₓ:sₘₐₓ
            blk = block[:, s, :]
            @test blk isa DegreeBlockBatch
            @test axes(blk) == (1:Nᵣ, -ℓ:ℓ)
            for i ∈ 1:Nᵣ, m ∈ -ℓ:ℓ
                if ℓ < abs(s)
                    @test iszero(blk[i, m])
                else
                    @test blk[i, m] == singles[(i, s)][Yindex(ℓ, m)]
                end
            end
        end
    end
    # Arbitrary ℓ order gives the same results as the sequential order
    for ℓ ∈ (0, 3, 1, 7, 7, 4, 0, 2)
        block = recurrence!(calc, ℓ)
        for s ∈ (-2, 0, 3)
            blk = block[:, s, :]
            for i ∈ 1:Nᵣ, m ∈ -ℓ:ℓ
                @test blk[i, m] == (ℓ < abs(s) ? 0 : singles[(i, s)][Yindex(ℓ, m)])
            end
        end
    end
    # copy keeps the natural axes and survives the next recurrence!; collect is 1-based
    v = recurrence!(calc, 4)[:, 1, :]
    c = copy(v)
    a = collect(v)
    @test axes(c) == axes(v)
    @test a isa Matrix{ComplexF64} && size(a) == (Nᵣ, 9)
    recurrence!(calc, 5)
    @test array_view(c) == [singles[(i, 1)][Yindex(4, m)] for i ∈ 1:Nᵣ, m ∈ -4:4]
    # Nᵣ == 1: a single rotor and a vector view
    c1 = sYlmCalculator(Rs[2], ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
    row = recurrence!(c1, 3)[-1, :]
    @test axes(row) == (-3:3,)
    @test collect(row) == singles[(2, -1)][Yindex(3, -3):Yindex(3, 3)]
    # Angle input, given at construction, evaluates at (θ, ϕ=0)
    θs = [0.0, 0.7, 1.9, π]
    cθ = sYlmCalculator(θs, 4, -2:2)
    blkθ = recurrence!(cθ, 2)
    using Quaternionic: from_spherical_coordinates
    for (i, θ) ∈ enumerate(θs), s ∈ -2:2, m ∈ -2:2
        Yref = array_view(sYlm(Rotor(from_spherical_coordinates(θ, 0.0)), 4, s; ℓₘᵢₙ=0))[Yindex(2, m)]
        @test blkθ[i, s, m] ≈ Yref atol=1e-15
        @test imag(blkθ[i, s, m]) == 0
    end
    # similar and show
    c2 = similar(calc)
    @test SphericalFunctions.Nᵣ(c2) == Nᵣ && SphericalFunctions.ℓₘₐₓ(c2) == ℓₘₐₓ && SphericalFunctions.spins(c2) == -sₘₐₓ:sₘₐₓ
    @test c2.Yˡ !== calc.Yˡ
    str = sprint(show, calc)
    @test occursin("sYlmCalculator", str) && occursin("ℓₘₐₓ=$ℓₘₐₓ", str) && occursin("Nᵣ=$Nᵣ", str)
    @test occursin("batched", str)
    str2 = sprint(show, MIME("text/plain"), c2)
    @test occursin("sYlmCalculator", str2) && occursin("ℓₘₐₓ=$ℓₘₐₓ", str2)
    @test occursin("nothing computed yet", str2)
end

@testitem "sYlmCalculator errors" setup=[RefusalChecks] begin
    import SphericalFunctions: sYlmCalculator, sYlm, sYlm!, sYlm_matrix, recurrence!
    using Quaternionic: Rotor
    using Random
    rng = Random.Xoshiro(3)
    R = randn(rng, Rotor{Float64})
    @test refuses(() -> sYlmCalculator(R, 3, 4), ArgumentError, "|s|=4, which exceeds ℓₘₐₓ=3")
    @test refuses(() -> sYlmCalculator(R, 3, -4:4), ArgumentError, "exceeds ℓₘₐₓ")
    # A negative spin weight is an ordinary one; what is refused is a range that is not a
    # consecutive run from low to high, with a note saying how to write one
    @test SphericalFunctions.spin(sYlmCalculator(R, 3, -1)) == -1
    unit = "A range of indices must be a unit range, running upward in steps of 1"
    @test refuses(() -> sYlmCalculator(R, 3, -2:2:2), ArgumentError, unit)
    @test refuses(() -> sYlmCalculator(R, 3, 2:-1:-2), ArgumentError, unit)
    @test refuses(() -> sYlmCalculator(R, 7//2, 3//2:-1:-3//2), ArgumentError, unit)
    @test refuses(() -> sYlmCalculator(R, 3, 0:1:2), ArgumentError, "may be written 0:2")
    @test refuses(() -> sYlmCalculator(R, 3, Base.OneTo(2)), ArgumentError, "as 1:2")
    # A unit range written downward is empty, and that is refused as well, by the flat
    # functions too
    for f ∈ (
        () -> sYlmCalculator(R, 3, 2:-2), () -> sYlmCalculator(R, 7//2, 3//2:-3//2),
        () -> sYlm(R, 3, 2:1), () -> sYlm_matrix([R], 3, 2:1),
        () -> sYlm!(zeros(ComplexF64, 2, 16), R, 3, 2:1),
    )
        @test refuses(f, ArgumentError, "runs downward or is empty")
    end
    calc = sYlmCalculator(R, 4, -2:2)
    # A spin weight the calculator does not serve is out of bounds of the block it returns
    @test_throws BoundsError recurrence!(calc, 2)[3, :]
    @test refuses(() -> recurrence!(calc, 5), ArgumentError, "out of bounds")
    @test refuses(() -> recurrence!(calc, R, -1), ArgumentError, "out of bounds")
    @test refuses(() -> recurrence!(calc, 2.0), ArgumentError, "so ℓ must be one too")
    # A complex "phase" is not a valid rotor for an sYlmCalculator
    @test refuses(() -> recurrence!(calc, cis(0.3), 2), ArgumentError, "rotors")
    @test refuses(() -> recurrence!(calc, [cis(0.3)], 2), ArgumentError, "rotors")
    @test refuses(() -> recurrence!(calc, [R, R], 2), DimensionMismatch, "Expected 1 rotors")
    batched = sYlmCalculator([R, R, R], 4, -2:2)
    @test refuses(() -> recurrence!(batched, R, 2), DimensionMismatch, "expects Nᵣ=3")
    @test refuses(() -> recurrence!(batched, [R, R], 2), DimensionMismatch, "Expected 3 rotors")
    @test refuses(() -> sYlm(R, 2, 3), ArgumentError, "|s|=3 exceeds ℓₘₐₓ=2")
    # A negative ℓₘₐₓ is named as such, rather than as a spin weight too large for it
    @test refuses(() -> sYlm(R, -1, 0), ArgumentError, "ℓₘₐₓ=-1 must be at least 0")
    # The message about ℓₘᵢₙ names the floor of the integer kind, and the bound ℓₘₐₓ
    @test refuses(
        () -> sYlm(R, 2, 0; ℓₘᵢₙ=-1), ArgumentError,
        "ℓₘᵢₙ=-1 must satisfy 0 ≤ ℓₘᵢₙ ≤ ℓₘₐₓ=2."
    )
    @test refuses(() -> sYlm(R, 2, 1; ℓₘᵢₙ=3), ArgumentError, "0 ≤ ℓₘᵢₙ ≤ ℓₘₐₓ=2")
    Y = zeros(ComplexF64, 5)
    @test refuses(() -> sYlm!(Y, R, 3, 0), DimensionMismatch, "Output vector has length")
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 25), sYlmCalculator(R, 4, 1), R, 2), ArgumentError,
        "not among them"
    )
    # A calculator built from a vector of rotors, even of one, has blocks with a rotor index,
    # so it cannot fill the vector of one rotor's values
    for c ∈ (batched, sYlmCalculator([R], 4, -2:2))
        @test refuses(
            () -> sYlm!(zeros(ComplexF64, 25), c, R, 1), ArgumentError,
            "`sYlm!` needs a calculator built for a single rotor"
        )
    end
    # The output's element type must be the calculator's own; it does not decide the type
    @test refuses(
        () -> sYlm!(zeros(ComplexF32, 25), calc, R, 1; ℓₘᵢₙ=0), ArgumentError,
        "element type must be Complex{Float64}"
    )
    @test refuses(
        () -> sYlm!(zeros(ComplexF32, 25), R, 4, 1; ℓₘᵢₙ=0), ArgumentError,
        "element type must be Complex{Float64}"
    )
    # ... and a rotor of another float type cannot be pushed through a calculator.
    # (`check_rotor_type` owns this message.)
    @test_throws "given data would give Float32" sYlm!(zeros(ComplexF64, 25), calc, Rotor{Float32}(R), 1; ℓₘᵢₙ=0)
    # An output whose shape does not suit the spin argument is refused, saying what is needed:
    # a vector for one spin weight, and a matrix for a range of them
    shape = "fills a vector for one spin weight, and a matrix with length(s) rows"
    @test refuses(() -> sYlm!(zeros(ComplexF64, 2, 25), R, 4, 0), ArgumentError, shape)
    @test refuses(() -> sYlm!(zeros(ComplexF64, 125), R, 4, -2:2), ArgumentError, shape)
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 2, 25), sYlmCalculator(R, 4, 0), R), ArgumentError, shape
    )
    @test refuses(() -> sYlm!(zeros(ComplexF64, 125), calc, R), ArgumentError, shape)
    @test refuses(() -> sYlm!(zeros(ComplexF64, 5, 25), calc, R, 1), ArgumentError, shape)
end

@testitem "sYlm! reuses a calculator" begin
    import SphericalFunctions: sYlmCalculator, sYlm, sYlm!, Ysize
    using Quaternionic: Rotor
    using Random
    rng = Random.Xoshiro(5)
    ℓₘₐₓ = 6
    Rs = randn(rng, Rotor{Float64}, 4)
    calc = sYlmCalculator(Rs[1], ℓₘₐₓ, -2:2)
    Y = Vector{ComplexF64}(undef, Ysize(0, ℓₘₐₓ))
    for R ∈ Rs, s ∈ -2:2
        # The values come back labelled, as a `HarmonicValues` over (a view of) `Y` itself
        @test parent(array_view(sYlm!(Y, calc, R, s; ℓₘᵢₙ=0))) === Y
        @test Y == array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=0))
        n = Ysize(abs(s), ℓₘₐₓ)
        sYlm!(Y, calc, R, s)
        @test Y[1:n] == array_view(sYlm(R, ℓₘₐₓ, s))
    end
    # Allocation-free after warm-up, whether or not the labelled result is used.  This is
    # measured inside a function, as a loop reusing the calculator runs: a call from global
    # scope must box the `HarmonicValues` it returns, which in compiled code is never built.
    R = randn(rng, Rotor{Float64})
    ignore_result!(Y, calc, R) = (sYlm!(Y, calc, R, 1; ℓₘᵢₙ=0); nothing)
    use_result!(Y, calc, R) = sum(abs2, array_view(sYlm!(Y, calc, R, 1; ℓₘᵢₙ=0)))
    ignore_result!(Y, calc, R); use_result!(Y, calc, R)
    @test @allocated(ignore_result!(Y, calc, R)) == 0
    # Before Julia 1.12 the labelled result is not elided once it is used, and costs one small
    # allocation (measured 16 bytes on 1.10 and 1.11); from 1.12 on it costs nothing.
    @test @allocated(use_result!(Y, calc, R)) ≤ (VERSION ≥ v"1.12" ? 0 : 16)
end

@testitem "sYlm generic types" begin
    import SphericalFunctions: sYlmCalculator, sYlm, recurrence!, Yindex, Ysize
    using Quaternionic: Rotor, from_spherical_coordinates
    using DoubleFloats: Double64
    import ForwardDiff
    import MathChecker: checked, unchecked
    using Random
    rng = Random.Xoshiro(13)
    R64 = randn(rng, Rotor{Float64})
    # Float32, against Float64 at the same (Float32-rounded) rotor, so that the difference
    # is the Float32 arithmetic alone.  Measured at most 11.7 eps(Float32) over thirty seeds.
    R32 = Rotor{Float32}(R64)
    Y32 = array_view(sYlm(R32, 20, -1))
    Y64 = array_view(sYlm(Rotor{Float64}(R32), 20, -1))
    @test eltype(Y32) === ComplexF32
    @test all(isfinite, Y32)
    @test maximum(abs, Y32 .- Y64) < 50eps(Float32)
    # BigFloat vs Double64
    YB = array_view(sYlm(Rotor{BigFloat}(R64), 4, 2))
    YD = array_view(sYlm(Rotor{Double64}(R64), 4, 2))
    @test maximum(abs, YB .- YD) < 1e-30
    # ForwardDiff through the rotor.  ϕ = 0 is included deliberately: there the spinor phase
    # `z₊` is exactly 1, and a `sqrt` of an exact zero inside `complex_powers!` used to make
    # every derivative NaN.  It is also the case the ring-based transforms use.
    θ₀ = 0.8
    for ϕ ∈ (0.0, 0.3)
        f(θ) = real(array_view(sYlm(Rotor(from_spherical_coordinates(θ, ϕ)), 3, 1))[Yindex(3, 2, 1)])
        dual = ForwardDiff.derivative(f, θ₀)
        h = 1e-6
        fd = (f(θ₀ + h) - f(θ₀ - h)) / 2h
        @test isfinite(dual)
        # The central difference is accurate to about eps/h ≈ 2e-10; measured 4e-11
        @test abs(dual - fd) < 1e-8
    end
    # ... and directly through `complex_powers!` at the exact phase 1
    let dz3 = ForwardDiff.derivative(
            x -> (Z = Vector{Complex{typeof(x)}}(undef, 6);
                  SphericalFunctions.complex_powers!(Z, Complex(one(x), zero(x) * x));
                  real(Z[3])),
            0.0
        )
        @test isfinite(dz3)
    end
    # No uninitialized memory is read: signaling NaNs everywhere, then a full sweep.
    # `MathChecker.Checked` with only the NaN check enabled throws as soon as a NaN takes
    # part in an operation, which is what turns an unwritten entry into a test failure.
    NC = checked(Float64; precision=false, nan=true, inf=false)
    # The half-integer cases reach the half-angle buffers and the seed of the rows m′ =
    # ±1/2, which only they use.
    for (ℓₘₐₓ, sₘₐₓ, Nᵣ) ∈ ((0, 0, 1), (2, 1, 1), (5, 2, 3), (9, 3, 2), (7//2, 3//2, 2), (1//2, 1//2, 1))
        Rs = randn(rng, Rotor{Float64}, Nᵣ)
        # A calculator works in the float type of its rotors, so the checked type is applied
        # to them; `ref` keeps the plain-Float64 rotors and is the value to compare against.
        calc = sYlmCalculator(Rotor{NC}.(Rs), ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
        fill!(calc, NaN)  # after construction, which stores the rotor data `fill!` preserves
        ref = sYlmCalculator(Rs, ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
        lo = ℓₘₐₓ isa Integer ? 0 : 1//2
        for ℓ ∈ [collect(lo:1:ℓₘₐₓ); lo + (ℓₘₐₓ - lo) ÷ 2]
            blk = recurrence!(calc, ℓ)
            refblk = recurrence!(ref, ℓ)
            for s ∈ -sₘₐₓ:sₘₐₓ
                # Every entry must have been written: an untouched one still holds the
                # sentinel NaN, which turns `err` into NaN and fails the comparison.  The
                # values are only approximately equal to the plain-`Float64` run because the
                # `muladd` in `complex_powers!` is fused for `Float64` but not for a wrapper
                # type, so the two round differently.
                err = 0.0
                for i ∈ 1:Nᵣ, m ∈ -ℓ:ℓ
                    # A vector of rotors, even of one, gives batched blocks
                    z = blk[i, s, m]
                    zref = refblk[i, s, m]
                    err = max(err, abs(unchecked(real(z)) - real(zref)), abs(unchecked(imag(z)) - imag(zref)))
                end
                @test err < 1e-13
            end
        end
        # `fill!` keeps the stored rotor data, as its docstring promises, so the recurrence
        # can be re-run without re-supplying it — and must still write every element.
        fill!(calc, NaN)
        blk = recurrence!(calc, ℓₘₐₓ)
        refblk = recurrence!(ref, ℓₘₐₓ)
        for s ∈ -sₘₐₓ:sₘₐₓ
            err = 0.0
            for i ∈ 1:Nᵣ, m ∈ -ℓₘₐₓ:ℓₘₐₓ
                z = blk[i, s, m]
                zref = refblk[i, s, m]
                err = max(err, abs(unchecked(real(z)) - real(zref)), abs(unchecked(imag(z)) - imag(zref)))
            end
            @test err < 1e-13
        end
    end
end


### Half-integer indices in the flat functions.
#
# The first item checks the half-integer values themselves, against their definition
# evaluated with the independent Wigner 𝔇 oracle of `test/wigner/half_integer_oracle.jl`.
# The rest take `sYlmCalculator` as their reference, and check that the flat functions lay
# its values out in the canonical ordering, accept the `Rational` spelling, refuse a mixture
# of the two kinds of index, and keep the two properties the transforms rest on —
# orthonormality on the sphere and antiperiodicity in ϕ.  Angles are fixed wherever a value
# is asserted exactly, so that a failure is reproducible.

@testitem "sYlm half-integer against the Wigner 𝔇 oracle" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: sYlm, sYlmCalculator, sλlm, recurrence!
    using Quaternionic: to_spherical_coordinates, from_spherical_coordinates

    # The documented definition, ₛYₗₘ(R) = i^{2s} √((2ℓ+1)/4π) conj(𝔇ˡ_{m,-s}(R)),
    # evaluated with the oracle's 𝔇 — a verbatim copy of Boyle (2016), cross-checked
    # against Varshalovich et al. — rather than with this package's `D`.  Nothing else in
    # the suite compares half-integer harmonics with an independent reference:
    # orthonormality, antiperiodicity, round trips and comparisons between the package's own
    # functions are all blind to a sign or phase error, such as a flipped sign for s = ±3/2.
    # Measured worst case 4.4 eps, and 2.5 eps for the real flavor below; 20 eps is
    # asserted.  i^{2s} is the principal branch e^{iπs} the documentation settles on
    oracle(R, ℓ, m, s) = cispi(Float64(s)) * √((2ℓ + 1) / 4π) * conj(HalfIntegerOracle.D_oracle(R, ℓ, m, -s))
    ℓₘₐₓ = 9//2
    for R ∈ HalfIntegerOracle.rotors(), s ∈ (-5//2, -3//2, -1//2, 1//2, 3//2, 5//2)
        Y = sYlm(R, ℓₘₐₓ, s)
        calc = sYlmCalculator(R, ℓₘₐₓ, s)
        for ℓ ∈ abs(s):ℓₘₐₓ
            block = recurrence!(calc, ℓ)
            for m ∈ -ℓ:ℓ
                expected = oracle(R, ℓ, m, s)
                @test Y[ℓ][m] ≈ expected atol=20eps()
                @test block[m] ≈ expected atol=20eps()
            end
        end
    end

    # The real flavor at θ is ₛY(θ, 0) / i^{2s}, so the same oracle checks it at the rotor
    # of (θ, 0)
    for θ ∈ (0.3, 1.1, 2.9), s ∈ (-3//2, 1//2, 5//2)
        R = from_spherical_coordinates(θ, 0.0)
        Λ = sλlm(θ, ℓₘₐₓ, s)
        for ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ
            @test Λ[ℓ][m] ≈ oracle(R, ℓ, m, s) / cispi(Float64(s)) atol=20eps()
        end
    end
end

@testitem "sYlm half-integer vs sYlmCalculator blocks" begin
    import SphericalFunctions: sYlm, sYlmCalculator, recurrence!, Ysize, Yindex, HalfOddInteger
    using Quaternionic: Rotor, from_euler_angles

    ℓₘₐₓ = 9//2
    Rs = [Rotor(from_euler_angles(α, β, γ)) for (α, β, γ) ∈ ((0.0, 0.0, 0.0), (0.7, 1.1, 2.3), (2.9, 0.4, 5.1), (4.0, 2.2, 0.3))]
    for R ∈ Rs
        calc = sYlmCalculator(R, ℓₘₐₓ, -3//2:3//2)
        for s ∈ (-3//2, -1//2, 1//2, 3//2), ℓₘᵢₙ ∈ (abs(s), 1//2)
            Y = array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ))
            @test eltype(Y) === ComplexF64
            @test length(Y) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            # The flat function and the calculator are the same engine, so the values agree
            # exactly, whatever the calculator's own sₘₐₓ.
            for ℓ ∈ ℓₘᵢₙ:1:ℓₘₐₓ
                blk = recurrence!(calc, ℓ)
                for m ∈ -ℓ:ℓ
                    @test Y[Yindex(ℓ, m, ℓₘᵢₙ)] == blk[s, m]
                end
            end
        end
        # The default ℓₘᵢₙ is |s|, for the half-integer kind as for the integer one.
        @test array_view(sYlm(R, ℓₘₐₓ, 3//2)) == array_view(sYlm(R, ℓₘₐₓ, 3//2; ℓₘᵢₙ=3//2))
        @test length(array_view(sYlm(R, ℓₘₐₓ, 3//2))) == Ysize(3//2, ℓₘₐₓ)
    end
    # For half-integer s the values include the phase i^{2s} = ±i: the ϕ = γ = 0 values,
    # which are real for integer s, are here purely imaginary.
    Rθ = Rotor(from_euler_angles(0.0, 1.1, 0.0))
    for s ∈ (-1//2, 1//2, 3//2)
        Yθ = array_view(sYlm(Rθ, ℓₘₐₓ, s))
        @test maximum(abs ∘ real, Yθ) == 0
        @test maximum(abs ∘ imag, Yθ) > 0.1
    end
end

@testitem "sYlm half-integer ℓₘᵢₙ = 1//2 gives zeros below |s|" setup=[RefusalChecks] begin
    import SphericalFunctions: sYlm, sYlm_matrix, Ysize, Yindex
    using Quaternionic: Rotor, from_spherical_coordinates

    ℓₘₐₓ = 9//2
    R = from_spherical_coordinates(0.7, 1.2)
    for s ∈ (-3//2, 3//2, 5//2)
        Y = array_view(sYlm(R, ℓₘₐₓ, s))
        Y₀ = array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=1//2))
        @test length(Y₀) == Ysize(1//2, ℓₘₐₓ)
        # The entries for ℓ < |s| are the first Ysize(1//2, |s| - 1) of them, and all zero.
        n₀ = Ysize(1//2, abs(s) - 1)
        @test n₀ == Yindex(abs(s), -abs(s), 1//2) - 1
        @test all(iszero, Y₀[1:n₀])
        @test Y₀[n₀+1:end] == Y
        M₀ = sYlm_matrix([R, -R], ℓₘₐₓ, s; ℓₘᵢₙ=1//2)
        @test all(iszero, M₀[:, 1:n₀])
        @test M₀[:, n₀+1:end] == sYlm_matrix([R, -R], ℓₘₐₓ, s)
    end
    # The floor of ℓₘᵢₙ is 1/2; anything below it, or above ℓₘₐₓ, is refused, with a message
    # that names the floor of this kind of index rather than the integer one.
    @test refuses(() -> sYlm(R, ℓₘₐₓ, 1//2; ℓₘᵢₙ=-1//2), ArgumentError, "must satisfy")
    @test refuses(() -> sYlm(R, ℓₘₐₓ, 1//2; ℓₘᵢₙ=11//2), ArgumentError, "must satisfy")
    @test refuses(
        () -> sYlm(R, ℓₘₐₓ, 1//2; ℓₘᵢₙ=-1//2), ArgumentError,
        "ℓₘᵢₙ=-1//2 must satisfy 1//2 ≤ ℓₘᵢₙ ≤ ℓₘₐₓ=9//2."
    )
    @test refuses(
        () -> sYlm_matrix([R, -R], ℓₘₐₓ, 3//2; ℓₘᵢₙ=11//2), ArgumentError, "1//2 ≤ ℓₘᵢₙ"
    )
    @test refuses(() -> sYlm(R, 1//2, 3//2), ArgumentError, "exceeds ℓₘₐₓ")
    @test refuses(() -> sYlm(R, -1//2, 1//2), ArgumentError, "ℓₘₐₓ=-1//2 must be at least 1//2")
end

@testitem "sYlm! half-integer, both forms, equals sYlm" setup=[RefusalChecks] begin
    import SphericalFunctions: sYlm, sYlm!, sYlmCalculator, Ysize
    using Quaternionic: Rotor, from_spherical_coordinates

    ℓₘₐₓ = 9//2
    Rs = [from_spherical_coordinates(θ, ϕ) for (θ, ϕ) ∈ ((0.0, 0.0), (0.7, 1.2), (2.2, 4.0), (π, 0.3))]
    calc = sYlmCalculator(Rs[1], ℓₘₐₓ, -3//2:3//2)
    Y = Vector{ComplexF64}(undef, Ysize(1//2, ℓₘₐₓ))
    for R ∈ Rs, s ∈ (-3//2, -1//2, 1//2, 3//2)
        @test parent(array_view(sYlm!(Y, calc, R, s; ℓₘᵢₙ=1//2))) === Y
        @test Y == array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=1//2))
        n = Ysize(abs(s), ℓₘₐₓ)
        sYlm!(Y, calc, R, s)
        @test Y[1:n] == array_view(sYlm(R, ℓₘₐₓ, s))
        # The allocating form, with the indices spelled as `Rational`s
        Y′ = Vector{ComplexF64}(undef, n)
        @test parent(array_view(sYlm!(Y′, R, ℓₘₐₓ, s))) === Y′
        @test Y′ == array_view(sYlm(R, ℓₘₐₓ, s))
    end
    # Allocation-free after warm-up, as for the integer kind (and measured inside a function
    # for the same reason)
    R = Rs[2]
    reuse!(Y, calc, R) = sum(abs2, array_view(sYlm!(Y, calc, R, 1//2; ℓₘᵢₙ=1//2)))
    reuse!(Y, calc, R)
    @test @allocated(reuse!(Y, calc, R)) ≤ (VERSION ≥ v"1.12" ? 0 : 16)
    # Errors: a spin weight the calculator does not serve, the output length, and a spin
    # weight or ℓₘᵢₙ of the wrong kind, which is refused with a message naming the kind of
    # the calculator, or that of the spin weight, rather than with a bare conversion error.
    @test refuses(() -> sYlm!(Y, calc, R, 5//2), ArgumentError, "not among them")
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 3), calc, R, 1//2), DimensionMismatch,
        "Output vector has length"
    )
    @test refuses(
        () -> sYlm!(Y, calc, R, 1), ArgumentError,
        "indices are half-odd-integers, like 7//2, so the spin weight s must be one too"
    )
    @test refuses(
        () -> sYlm!(Y, calc, R, 1//2; ℓₘᵢₙ=0), ArgumentError, "keyword argument `ℓₘᵢₙ`"
    )
    # Without a spin weight the calculator alone fixes the kind of ℓₘᵢₙ
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 4, Ysize(1//2, ℓₘₐₓ)), calc, R; ℓₘᵢₙ=0), ArgumentError,
        "The indices of this `sYlmCalculator` are half-odd-integers"
    )
    icalc = sYlmCalculator(R, 4, 1)
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 25), icalc, R, 1//2), ArgumentError,
        "indices are integers, like 3, so the spin weight s must be one too"
    )
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 25), icalc, R, 1; ℓₘᵢₙ=1//2), ArgumentError,
        "keyword argument `ℓₘᵢₙ`"
    )
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 25), icalc, R; ell_min=1//2), ArgumentError,
        "The indices of this `sYlmCalculator` are integers of type `Int`, like 3; got ell_min"
    )
end

@testitem "sYlm_matrix half-integer rows equal sYlm" setup=[RefusalChecks] begin
    import SphericalFunctions: sYlm, sYlm_matrix, Ysize
    using Quaternionic: Rotor, from_spherical_coordinates

    ℓₘₐₓ = 7//2
    Rs = [from_spherical_coordinates(θ, ϕ) for θ ∈ (0.0, 0.7, 2.2, π) for ϕ ∈ (0.0, 1.2, 4.0)]
    for s ∈ (-3//2, -1//2, 1//2, 3//2), ℓₘᵢₙ ∈ (abs(s), 1//2)
        M = sYlm_matrix(Rs, ℓₘₐₓ, s; ℓₘᵢₙ)
        @test M isa Matrix{ComplexF64}
        @test size(M) == (length(Rs), Ysize(ℓₘᵢₙ, ℓₘₐₓ))
        for (i, R) ∈ enumerate(Rs)
            @test M[i, :] == array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ))
        end
    end
    @test refuses(() -> sYlm_matrix(Rs, 1//2, 3//2), ArgumentError, "exceeds ℓₘₐₓ")
end

@testitem "sYlm half-integer spellings and mixed kinds" setup=[RefusalChecks] begin
    import SphericalFunctions: sYlm, sYlm!, sYlm_matrix, sYlmCalculator, Ysize, HalfOddInteger
    using Quaternionic: Rotor, from_spherical_coordinates

    R = from_spherical_coordinates(0.7, 1.2)
    Rs = [R, -R]
    ℓₘₐₓ, s, ℓₘᵢₙ = HalfOddInteger(7//2), HalfOddInteger(1//2), HalfOddInteger(1//2)
    # The `Rational` spelling is normalized at the boundary and gives exactly what the
    # `HalfOddInteger` spelling gives, keyword included.
    @test array_view(sYlm(R, 7//2, 1//2)) == array_view(sYlm(R, ℓₘₐₓ, s))
    @test array_view(sYlm(R, 7//2, 3//2; ℓₘᵢₙ=1//2)) == array_view(sYlm(R, ℓₘₐₓ, HalfOddInteger(3//2); ℓₘᵢₙ))
    @test array_view(sYlm(R, 7//2, HalfOddInteger(1//2))) == array_view(sYlm(R, ℓₘₐₓ, s))
    @test sYlm_matrix(Rs, 7//2, 1//2) == sYlm_matrix(Rs, ℓₘₐₓ, s)
    @test sYlm_matrix(Rs, 7//2, 1//2; ℓₘᵢₙ=1//2) == sYlm_matrix(Rs, ℓₘₐₓ, s; ℓₘᵢₙ)
    Yr = array_view(sYlm!(Vector{ComplexF64}(undef, Ysize(1//2, 7//2)), R, 7//2, 1//2))
    Yh = array_view(sYlm!(Vector{ComplexF64}(undef, Ysize(1//2, 7//2)), R, ℓₘₐₓ, s))
    @test Yr == Yh
    calc = sYlmCalculator(R, 7//2, 1//2)
    @test array_view(sYlm!(similar(Yr), calc, R, 1//2; ℓₘᵢₙ=1//2)) ==
        array_view(sYlm!(similar(Yr), calc, R, s; ℓₘᵢₙ))
    # ... as are the keyword's ASCII spelling and the calculator form's
    @test array_view(sYlm(R, 7//2, 3//2; ell_min=1//2)) ==
        array_view(sYlm(R, ℓₘₐₓ, HalfOddInteger(3//2); ℓₘᵢₙ))
    @test array_view(sYlm!(similar(Yr), calc, R, 1//2; ell_min=1//2)) ==
        array_view(sYlm!(similar(Yr), calc, R, s; ℓₘᵢₙ))

    # A mixture of the two kinds of positional index is refused with a message that lists
    # the indices, names both spellings of a half-odd-integer, and says which kind each
    # index is
    mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
    for f ∈ (
        () -> sYlm(R, 7//2, 1), () -> sYlm(R, 4, 1//2), () -> sYlm(R, 4, HalfOddInteger(1//2)),
        () -> sYlm_matrix(Rs, 7//2, 1), () -> sYlm_matrix(Rs, 4, 1//2),
        () -> sYlm!(similar(Yr), R, 7//2, 1), () -> sYlm!(similar(Yr), R, 4, 1//2),
    )
        @test refuses(f, ArgumentError, mixed)
    end
    @test refuses(
        () -> sYlm(R, 7//2, 1), ArgumentError,
        "and so mixes integers (s) with half-odd-integers (ℓₘₐₓ)"
    )
    # A keyword of the other kind is refused by name, in either spelling
    @test refuses(() -> sYlm(R, 7//2, 1//2; ℓₘᵢₙ=0), ArgumentError, "keyword argument `ℓₘᵢₙ`")
    @test refuses(() -> sYlm(R, 4, 1; ℓₘᵢₙ=1//2), ArgumentError, "keyword argument `ℓₘᵢₙ`")
    @test refuses(() -> sYlm(R, 4, 1; ell_min=1//2), ArgumentError, "keyword argument `ell_min`")
    # A `Rational` that is not a half-odd-integer is refused as such
    @test refuses(() -> sYlm(R, 7//3, 1//3), ArgumentError, "7//3 is neither an integer nor")
    @test refuses(() -> sYlm(R, 4//1, 1//1), ArgumentError, "4//1 is a whole number")

    # An integer index of any type but `Int` is refused, because the index arithmetic is not
    # closed under it: ℓ² overflows a narrow type, and -m wraps around in an unsigned one.
    for (IT, sentence) ∈ (
        (Int8, "narrower than `Int`"), (Int16, "narrower than `Int`"),
        (Int32, "narrower than `Int`"), (UInt8, "is unsigned"), (BigInt, "wider than `Int`"),
    )
        for f ∈ (
            () -> sYlm(R, IT(4), IT(1)), () -> sYlm(R, 4, IT(1)),
            () -> sYlm(R, IT(4), 1; ℓₘᵢₙ=2), () -> sYlm_matrix(Rs, IT(4), IT(1)),
            () -> sYlm!(Vector{ComplexF64}(undef, Ysize(1, 4)), R, IT(4), IT(1)),
            () -> sYlmCalculator(R, IT(4), IT(1)), () -> sYlmCalculator(R, IT(4), -1:1),
        )
            @test refuses(f, ArgumentError, sentence)
        end
        @test refuses(() -> sYlm(R, 4, 1; ℓₘᵢₙ=IT(2)), ArgumentError, "keyword argument `ℓₘᵢₙ`")
        @test refuses(() -> sYlm(R, 4, 1; ℓₘᵢₙ=IT(2)), ArgumentError, sentence)
    end
    # ... in the calculator forms as well, whether or not a spin weight is given, and for a
    # `HarmonicValues` output too; each spelling of the keyword is named as it was written
    icalc = sYlmCalculator(R, 4, -1:1)
    @test refuses(
        () -> sYlm!(Vector{ComplexF64}(undef, Ysize(1, 4)), icalc, R, Int8(1)), ArgumentError,
        "narrower than `Int`"
    )
    let calc1 = sYlmCalculator(R, 4, 1), Y = Vector{ComplexF64}(undef, Ysize(1, 4))
        Yv, Ym = sYlm(R, 4, 1), sYlm(R, 4, -1:1)
        Ymat = Matrix{ComplexF64}(undef, 3, Ysize(1, 4))
        for (f, name) ∈ (
            (() -> sYlm!(Y, calc1, R; ℓₘᵢₙ=Int32(1)), "ℓₘᵢₙ"),
            (() -> sYlm!(Y, calc1, R; ell_min=Int32(1)), "ell_min"),
            (() -> sYlm!(Ymat, icalc, R; ℓₘᵢₙ=Int32(1)), "ℓₘᵢₙ"),
            (() -> sYlm!(Y, calc1, R, 1; ℓₘᵢₙ=Int32(1)), "ℓₘᵢₙ"),
            (() -> sYlm!(Yv, calc1, R; ell_min=Int32(1)), "ell_min"),
            (() -> sYlm!(Ym, icalc, R; ℓₘᵢₙ=Int16(1)), "ℓₘᵢₙ"),
            (() -> sYlm!(Yv, calc1, R, 1; ell_min=Int32(1)), "ell_min"),
        )
            @test refuses(f, ArgumentError, "got $name = 1::")
            @test refuses(f, ArgumentError, "convert it with `Int`")
        end
        @test refuses(
            () -> sYlm!(Y, calc1, R; ℓₘᵢₙ=Int32(1)), ArgumentError,
            "The indices of this `sYlmCalculator` are integers of type `Int`"
        )
        @test refuses(
            () -> sYlm!(Y, calc1, R; ℓₘᵢₙ=1//2), ArgumentError,
            "The indices of this `sYlmCalculator` are integers of type `Int`"
        )
        # An `Int` is accepted, as it is everywhere
        @test array_view(sYlm!(Y, calc1, R; ℓₘᵢₙ=1)) == array_view(Yv)
        @test sYlm!(Yv, calc1, R; ell_min=1) === Yv
    end
    # A half-odd-integer spelled as a `Rational` of another integer type than `Int` is
    # refused, and says how to write it
    @test refuses(
        () -> sYlm(R, big(7)//2, big(1)//2), ArgumentError,
        "`Rational{BigInt}` is not `Rational{Int}`; write the value with `Int`s, as 7//2"
    )
    @test refuses(
        () -> sYlm_matrix(Rs, Int8(7)//Int8(2), Int8(1)//Int8(2)), ArgumentError,
        "`Rational{Int8}` is not `Rational{Int}`"
    )
end

@testitem "sYlm half-integer orthonormality on the sphere" begin
    import SphericalFunctions: sYlm_matrix, clenshaw_curtis_rings, clenshaw_curtis
    using Quaternionic: Rotor, from_spherical_coordinates
    using LinearAlgebra: Diagonal, I

    # ∫ ₛYₗₘ conj(ₛYₗ′ₘ′) sinθ dθ dϕ = δ δ, by quadrature: Clenshaw–Curtis in θ with N =
    # 2ℓₘₐₓ+1 rings (whose weights include the sinθ), and Nϕ = 2ℓₘₐₓ+1 equally spaced ϕ,
    # both of which are exact at this band limit.  The sample points are the rotors
    # `from_spherical_coordinates(θ, ϕ)`, with ϕ running once around [0, 2π).  Measured
    # worst case over this grid: 1.6e-15 at ℓₘₐₓ = 9/2, so 1e-14 leaves a factor of ≳ 6.
    for ℓₘₐₓ ∈ (1//2, 3//2, 7//2, 9//2)
        N = Int(2ℓₘₐₓ + 1)
        θs, wθ = clenshaw_curtis_rings(N), clenshaw_curtis(N)
        Nϕ = N
        ϕs = [2π * k / Nϕ for k ∈ 0:Nϕ-1]
        Rs = [from_spherical_coordinates(θ, ϕ) for θ ∈ θs for ϕ ∈ ϕs]
        w = [wθ[i] * 2π / Nϕ for i ∈ eachindex(θs) for _ ∈ ϕs]
        for s ∈ -min(ℓₘₐₓ, 3//2):1:min(ℓₘₐₓ, 3//2)
            Y = sYlm_matrix(Rs, ℓₘₐₓ, s)
            G = Y' * Diagonal(w) * Y
            @test maximum(abs, G - I) < 1e-14
        end
    end
end

@testitem "sYlm half-integer antiperiodicity in ϕ" begin
    import SphericalFunctions: sYlm
    using Quaternionic: Rotor, from_spherical_coordinates

    # A circuit in ϕ returns to the antipodal rotor, so for half-integer s the harmonics
    # change sign: ₛY(θ, ϕ+2π) = -ₛY(θ, ϕ).  For integer s they do not.  Both are only
    # reproduced to rounding, because the phases are recomputed from a rotor whose
    # half-angle differs by rounding.  Measured worst case at these angles: 7.8e-16 for the
    # half-integer kind and 8.4e-16 for the integer kind, so 1e-14 leaves a factor of ≳ 12.
    for (θ, ϕ) ∈ ((0.7, 1.2), (2.2, 4.0), (1.0, 0.0), (0.3, 5.9))
        R, R′ = from_spherical_coordinates(θ, ϕ), from_spherical_coordinates(θ, ϕ + 2π)
        for s ∈ (-3//2, -1//2, 1//2, 3//2)
            Y, Y′ = array_view(sYlm(R, 9//2, s)), array_view(sYlm(R′, 9//2, s))
            @test maximum(abs, Y′ + Y) < 1e-14
        end
        for s ∈ (-1, 0, 2)
            Y, Y′ = array_view(sYlm(R, 4, s)), array_view(sYlm(R′, 4, s))
            @test maximum(abs, Y′ - Y) < 1e-14
        end
    end
end

@testitem "sYlm half-integer Float32 rotor" begin
    import SphericalFunctions: sYlm, sYlm_matrix, Ysize
    using Quaternionic: Rotor, from_spherical_coordinates

    R32 = from_spherical_coordinates(0.7f0, 1.2f0)
    R64 = from_spherical_coordinates(0.7, 1.2)
    @test R32 isa Rotor{Float32}
    for s ∈ (-1//2, 1//2, 3//2)
        Y32 = array_view(sYlm(R32, 9//2, s))
        Y64 = array_view(sYlm(R64, 9//2, s))
        @test eltype(Y32) === ComplexF32
        @test length(Y32) == Ysize(abs(s), 9//2)
        @test all(isfinite, Y32)
        # Float32 arithmetic against Float64, at the accuracy Float32 allows
        @test maximum(abs, Y32 .- Y64) < 1e-5
        M32 = sYlm_matrix([R32, -R32], 9//2, s)
        @test M32 isa Matrix{ComplexF32}
        @test M32[1, :] == Y32
    end
end


### Ranges of spin weights.
#
# The calculator is built for the spin weights it will serve, and those may be given either
# singly or as an ascending range.  The items below check that the two spellings are
# accepted in every form the indices take, that a range and the single spin weights
# composing it agree bit for bit — serving several spin weights from one calculator changes
# no arithmetic — and that the flat functions lay a range out the way their docstrings say.

@testitem "sYlmCalculator spin weights, however spelled" setup=[RefusalChecks] begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, spins, spin, HalfOddInteger
    using Quaternionic: Rotor
    using Random

    R = randn(Random.Xoshiro(17), Rotor{Float64})

    # A single spin weight, in each of its three spellings
    for (spelling, value) ∈ ((2, 2), (-2, -2), (3//2, HalfOddInteger(3//2)),
                             (HalfOddInteger(-1//2), HalfOddInteger(-1//2)))
        ℓₘₐₓ = value isa Integer ? 4 : 9//2
        calc = sYlmCalculator(R, ℓₘₐₓ, spelling)
        @test spin(calc) == value
        @test spins(calc) == value:value
        @test length(spins(calc)) == 1
    end

    # A range, in each of its three spellings; the two half-odd ones are the same calculator
    for (spelling, lo, hi) ∈ (
        (-2:2, -2, 2), (1:2, 1, 2), (0:0, 0, 0),
        (-3//2:3//2, HalfOddInteger(-3//2), HalfOddInteger(3//2)),
        (HalfOddInteger(-3//2):HalfOddInteger(3//2), HalfOddInteger(-3//2), HalfOddInteger(3//2)),
        (HalfOddInteger(1//2):HalfOddInteger(5//2), HalfOddInteger(1//2), HalfOddInteger(5//2)),
    )
        ℓₘₐₓ = lo isa Integer ? 4 : 9//2
        calc = sYlmCalculator(R, ℓₘₐₓ, spelling)
        @test spins(calc) == lo:hi
        @test eltype(spins(calc)) === typeof(lo)
        @test length(spins(calc)) == length(spelling)
        # There is no single spin weight to name, so `spin` has no method at all
        @test_throws MethodError spin(calc)
    end

    # A range of an integer type other than `Int` is refused, as a single index of one is,
    # with the spelling to use instead
    for IT ∈ (Int8, Int16, Int32)
        @test refuses(
            () -> sYlmCalculator(R, 4, IT(-1):IT(1)), ArgumentError, "write it with `Int`s, as -1:1"
        )
    end

    # Refusals: a mixture of the two kinds of index, on either side of the colon
    mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
    @test refuses(() -> sYlmCalculator(R, 4, -3//2:3//2), ArgumentError, mixed)
    @test refuses(() -> sYlmCalculator(R, 9//2, -1:1), ArgumentError, mixed)
    @test refuses(() -> sYlmCalculator(R, 4, -5:5), ArgumentError, "exceeds ℓₘₐₓ")
    @test refuses(() -> sYlmCalculator(R, 9//2, -11//2:11//2), ArgumentError, "exceeds ℓₘₐₓ")
end

@testitem "sYlmCalculator ranges agree with single spin weights" begin
    import SphericalFunctions: sYlmCalculator, sYlm, spins, recurrence!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(2024)
    rotors = randn(rng, Rotor{Float64}, 3)
    θs = [0.0, 0.8, 2.6]

    # Serving a range of spin weights changes no arithmetic: a calculator built for a range
    # gives, for each spin weight in it, exactly the bits a calculator built for that one spin
    # weight gives.  Checked for both kinds of index, both shapes of rotor data, and batched
    # as well as single.
    for (ℓₘₐₓ, srange) ∈ ((5, -2:2), (5, 1:2), (9//2, -3//2:3//2), (9//2, 1//2:3//2))
        for data ∈ (rotors[1], rotors, θs[2], θs)
            ranged = sYlmCalculator(data, ℓₘₐₓ, srange)
            batched = data isa AbstractVector
            for s ∈ spins(ranged)
                singly = sYlmCalculator(data, ℓₘₐₓ, s)
                # The slice for one spin weight holds the same numbers a calculator built
                # for that spin weight alone gives
                expected = [collect(b) for (_, b) ∈ singly]
                for (k, (ℓ, b)) ∈ enumerate(ranged)
                    @test collect(batched ? b[:, s, :] : b[s, :]) == expected[k]
                end
            end
        end
    end

    # Below |s| the values are zero, in a range exactly as for a single spin weight
    calc = sYlmCalculator(rotors[1], 4, -2:2)
    blk = recurrence!(calc, 1)
    for s ∈ (-2, 2), m ∈ -1:1
        @test iszero(blk[s, m])
    end
    @test !iszero(blk[0, 0])

    # A range calculator also reproduces the flat `sYlm` for each spin weight in the range
    for s ∈ -2:2
        Y = array_view(sYlm(rotors[1], 4, s; ℓₘᵢₙ=0))
        for (ℓ, b) ∈ sYlmCalculator(rotors[1], 4, -2:2)
            @test all(b[s, m] == Y[SphericalFunctions.Yindex(ℓ, m)] for m ∈ -ℓ:ℓ)
        end
    end
end

@testitem "sYlm flat functions take ranges of spin weights" begin
    import SphericalFunctions: sYlm, sYlm!, sYlm_matrix, Ysize, spins, sYlmCalculator
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(909)
    rotors = randn(rng, Rotor{Float64}, 4)
    R = rotors[1]

    for (ℓₘₐₓ, srange, ℓₘᵢₙ) ∈ (
        (5, -2:2, 0), (5, 1:2, 1), (5, -2:-1, 1),
        (9//2, -3//2:3//2, 1//2), (9//2, 1//2:3//2, 1//2),
    )
        sr = spins(sYlmCalculator(R, ℓₘₐₓ, srange))
        n, nmodes = length(sr), Ysize(ℓₘᵢₙ, ℓₘₐₓ)

        # `sYlm` gives a matrix of spin weights by modes, in the order the range was given,
        # with `ℓₘᵢₙ` defaulting to the smallest |s| in it
        Y = array_view(sYlm(R, ℓₘₐₓ, srange))
        @test Y isa Matrix{ComplexF64}
        @test size(Y) == (n, nmodes)
        for (i, s) ∈ enumerate(sr)
            @test Y[i, :] == array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ))
        end

        # `sYlm!` fills the same thing, and returns it
        Y′ = similar(Y)
        @test parent(array_view(sYlm!(Y′, R, ℓₘₐₓ, srange))) === Y′
        @test Y′ == Y
        # ... as does the calculator form, which needs no spin weight of its own
        calc = sYlmCalculator(R, ℓₘₐₓ, srange)
        fill!(Y′, 0)
        @test parent(array_view(sYlm!(Y′, calc, rotors[2]))) === Y′
        @test Y′ == array_view(sYlm(rotors[2], ℓₘₐₓ, srange))
        # ... and one spin weight of that same calculator still fills a vector
        v = Vector{ComplexF64}(undef, nmodes)
        @test array_view(sYlm!(v, calc, rotors[2], first(sr); ℓₘᵢₙ)) == array_view(sYlm(rotors[2], ℓₘₐₓ, first(sr); ℓₘᵢₙ))

        # `sYlm_matrix` gives a stack of synthesis matrices, indexed [rotor, spin, mode]
        M = sYlm_matrix(rotors, ℓₘₐₓ, srange)
        @test M isa Array{ComplexF64, 3}
        @test size(M) == (length(rotors), n, nmodes)
        for (i, s) ∈ enumerate(sr)
            @test M[:, i, :] == sYlm_matrix(rotors, ℓₘₐₓ, s; ℓₘᵢₙ)
        end
        for (j, Rj) ∈ enumerate(rotors)
            @test M[j, :, :] == array_view(sYlm(Rj, ℓₘₐₓ, srange))
        end
    end

    # An explicit ℓₘᵢₙ is honoured, and the rows below their own |s| are zero
    Y = array_view(sYlm(R, 4, -2:2; ℓₘᵢₙ=0))
    @test size(Y) == (5, Ysize(0, 4))
    @test all(iszero, Y[1, 1:Ysize(0, 1)])   # s = -2 has nothing below ℓ = 2
    @test !all(iszero, Y[3, 1:Ysize(0, 1)])  # ... while s = 0 does

    # Errors: an output of the wrong shape or element type
    @test_throws "Output matrix has size" sYlm!(zeros(ComplexF64, 2, 100), R, 4, -2:2)
    @test_throws "Output matrix has size" sYlm!(zeros(ComplexF64, 5, 3), R, 4, -2:2)
    @test_throws "element type must be Complex{Float64}" sYlm!(
        zeros(ComplexF32, 5, Ysize(0, 4)), R, 4, -2:2
    )
end

@testitem "YlmCalculator is spin weight zero" begin
    using Quaternionic: Rotor
    import SphericalFunctions: spin, spins
    using Random

    rng = Random.Xoshiro(2027)
    ℓₘₐₓ = 5
    R = randn(rng, Rotor{Float64})
    Rs = randn(rng, Rotor{Float64}, 3)

    calc = YlmCalculator(R, ℓₘₐₓ)
    @test calc isa sYlmCalculator
    @test spin(calc) == 0
    @test spins(calc) == 0:0

    # It computes exactly what the spin-weight-zero sYlmCalculator does
    ref = sYlmCalculator(R, ℓₘₐₓ, 0)
    for ℓ ∈ 0:ℓₘₐₓ
        @test recurrence!(calc, ℓ) == recurrence!(ref, ℓ)
    end
    # ... and agrees with the one-shot `Ylm`
    Y = Ylm(R, ℓₘₐₓ)
    for (ℓ, Yˡ) ∈ YlmCalculator(R, ℓₘₐₓ)
        @test Yˡ == Y[ℓ]
    end

    # A collection of rotors gives the batched blocks, as sYlmCalculator does
    batched = YlmCalculator(Rs, ℓₘₐₓ)
    @test axes(recurrence!(batched, 3)) == (1:3, -3:3)

    # Half-integer ℓ has no spin-weight-zero analogue, and says where to look instead
    @test_throws ArgumentError YlmCalculator(R, 7//2)
    @test_throws "`sYlm` and `sYlmCalculator` accept half-integer indices" YlmCalculator(R, 7//2)
end

# The items above check the values.  These cover the entry points and refusals around them:
# the unweighted `Ylm` wrapper in its vector form, reusing a calculator for a different set
# of rotors, and the three ways a call is turned away.

@testitem "sYlm: the `Ylm` wrapper and the calculator refusals" setup=[RefusalChecks] begin
    using Quaternionic: Rotor, RotorF64
    import SphericalFunctions: Nᵣ, spins
    using Random

    rng = Random.Xoshiro(2026)
    ℓₘₐₓ = 5
    R = randn(rng, RotorF64)
    R⃗ = randn(rng, RotorF64, 4)

    # `Ylm` is `sYlm` at spin weight zero, for one rotor and for many
    @test Ylm(R, ℓₘₐₓ) == sYlm(R, ℓₘₐₓ, 0)
    @test Ylm(R⃗, ℓₘₐₓ) == sYlm(R⃗, ℓₘₐₓ, 0)
    @test Ylm(R, ℓₘₐₓ; ℓₘᵢₙ=2) == sYlm(R, ℓₘₐₓ, 0; ℓₘᵢₙ=2)
    @test Ylm(R⃗, ℓₘₐₓ; ℓₘᵢₙ=2) == sYlm(R⃗, ℓₘₐₓ, 0; ℓₘᵢₙ=2)
    @test Nᵣ(Ylm(R⃗, ℓₘₐₓ)) == length(R⃗)
    # ... and one rotor of the batch is the single-rotor answer
    @test array_view(Ylm(R⃗, ℓₘₐₓ))[2, :] == array_view(Ylm(R⃗[2], ℓₘₐₓ))

    # `similar(calc, R)` rebuilds a calculator around new rotor data, but only for the same
    # number of rotors it was built to hold
    c1 = sYlmCalculator(R, ℓₘₐₓ, 0)
    @test Nᵣ(similar(c1, randn(rng, RotorF64))) == 1
    @test refuses(
        () -> similar(c1, randn(rng, RotorF64, 3)), DimensionMismatch, "handles Nᵣ=1"
    )
    c4 = sYlmCalculator(R⃗, ℓₘₐₓ, 0)
    @test Nᵣ(similar(c4, randn(rng, RotorF64, 4))) == 4
    @test refuses(
        () -> similar(c4, randn(rng, RotorF64, 2)), DimensionMismatch, "handles Nᵣ=4"
    )

    # A calculator built for one spin weight refuses another, and names the ones it serves
    cs = sYlmCalculator(R, ℓₘₐₓ, -2)
    @test spins(cs) == -2:-2
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, Ysize(1, ℓₘₐₓ)), cs, R, 1), ArgumentError,
        "serves the spin weights -2, so s=1 is not among them"
    )

    crange = sYlmCalculator(R, ℓₘₐₓ, -2:2)
    @test spins(crange) == -2:2
    # a spin weight inside the range is served, and agrees with computing it afresh
    Y1 = sYlm!(zeros(ComplexF64, Ysize(1, ℓₘₐₓ)), crange, R, 1)
    @test array_view(Y1) ≈ array_view(sYlm(R, ℓₘₐₓ, 1))
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, Ysize(3, ℓₘₐₓ)), crange, R, 3), ArgumentError,
        "not among them"
    )

    # An output vector shorter than the modes it must hold is refused rather than truncated,
    # and one that is longer is labelled over just the modes written
    needed = Ysize(0, ℓₘₐₓ)
    @test length(array_view(sYlm!(zeros(ComplexF64, needed), c1, R, 0))) == needed
    longer = zeros(ComplexF64, needed + 3)
    Y0 = sYlm!(longer, c1, R, 0)
    @test length(array_view(Y0)) == needed && parent(array_view(Y0)) === longer
    @test array_view(Y0) == array_view(sYlm(R, ℓₘₐₓ, 0)) && all(iszero, longer[needed+1:end])
    @test refuses(() -> sYlm!(zeros(ComplexF64, needed - 1), c1, R, 0), DimensionMismatch, "is needed")
    @test refuses(
        () -> sYlm!(zeros(ComplexF64, 3), c1, R, 0), DimensionMismatch, "Output vector has length"
    )

    # A `HarmonicValues` computed for a vector of rotors cannot be refilled from a single
    # rotor, even when the vector held only one, since its storage has a rotor axis
    for Y ∈ (sYlm([R], ℓₘₐₓ, 0), sYlm([R], ℓₘₐₓ, -1:1), sYlm(R⃗, ℓₘₐₓ, 0))
        spin = Y.s
        @test refuses(() -> sYlm!(Y, R, ℓₘₐₓ, spin), ArgumentError, "built for one rotor")
    end
    @test refuses(() -> sYlm!(sYlm([R], ℓₘₐₓ, 0), c1, R), ArgumentError, "built for one rotor")
end


@testitem "sYlmCalculator: each block element is the wedge element `wedge_value` reads" begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, sλlmCalculator, recurrence!, wedge_value,
        zpower, sYlm_coefficient, spins, Nᵣ, ℓₘᵢₙ, ℓₘₐₓ, floattype, number_type
    using Quaternionic: Rotor
    import Random

    # As for the Wigner calculators, the elements H[m, -s] are read in runs along the rows
    # of the wedge, and in one of two orders according to the numbers of rotors and of spin
    # weights.  Each value is compared here, bit for bit, with the one built from
    # `wedge_value` and the same power tables, and each row below ℓ < |s| is zero.
    rng = Random.Xoshiro(20260924)
    R⃗ = randn(rng, Rotor{Float64}, 4)
    θ⃗ = [0.0, 0.7, 2.2, π]
    function expected(calc, H, iᵣ, s, m, ℓ)
        NT, RT = number_type(calc), floattype(calc)
        abs(s) > ℓ && return zero(NT)
        prefactor = √((2ℓ + 1) / (4 * RT(π)))
        value = sYlm_coefficient(NT, RT, 1, m, s, prefactor) * wedge_value(H, iᵣ, m, -s)
        if NT <: Complex && calc.phases[]
            value * (zpower(calc.Z₊, iᵣ, m - s) * zpower(calc.Z₋, iᵣ, m + s))
        else
            value
        end
    end
    count = Ref(0)
    for (ℓmax, spinsets) ∈ ((9, (0, 2, -3, -2:2, 0:3, -9:9)), (17//2, (1//2, -5//2, -3//2:1//2)))
        for s ∈ spinsets, (Ctor, data) ∈ (
            (sYlmCalculator, R⃗), (sYlmCalculator, R⃗[1]), (sYlmCalculator, θ⃗),
            (sλlmCalculator, θ⃗), (sλlmCalculator, θ⃗[2]),
        )
            calc = Ctor(data, ℓmax, s)
            for ℓ ∈ ℓₘᵢₙ(calc):ℓₘₐₓ(calc)
                recurrence!(calc, ℓ)
                H = calc.H.Hˡ
                good = true
                for (i, sᵢ) ∈ enumerate(spins(calc)), m ∈ -ℓ:ℓ, iᵣ ∈ 1:Nᵣ(calc)
                    value = calc.Yˡ[iᵣ, i, Int(m + ℓ) + 1]
                    good &= isequal(value, expected(calc, H, iᵣ, sᵢ, m, ℓ))
                    count[] += 1
                end
                @test good
            end
        end
    end
    @test count[] > 20_000
end


@testitem "sYlm!: one spin weight read from a calculator built for several" begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, sλlmCalculator, sYlm, sλlm, sYlm!, sλlm!,
        recurrence!, spin_row!, Ysize, ℓ, ℓₘᵢₙ, ℓₘₐₓ
    using Quaternionic: Rotor
    import Random

    # Reading one spin weight out of a calculator built for a range assembles only that spin
    # weight's values at each ℓ.  They are the values the calculator would otherwise have
    # given, bit for bit; the calculator then holds no complete block, and says so, and the
    # next full step gives the whole block again.
    rng = Random.Xoshiro(8)
    R, R₂ = randn(rng, Rotor{Float64}, 2)
    for (ℓmax, range, ℓlow) ∈ ((7, -3:3, 0), (15//2, -5//2:3//2, 1//2))
        calc = sYlmCalculator(R, ℓmax, range)
        for s ∈ range
            Y = zeros(ComplexF64, Ysize(ℓlow, ℓmax))
            sYlm!(Y, calc, R₂, s; ℓₘᵢₙ=ℓlow)
            @test isequal(Y, array_view(sYlm(R₂, ℓmax, s; ℓₘᵢₙ=ℓlow)))
            @test ℓ(calc) < ℓₘᵢₙ(calc)
            @test occursin("nothing computed yet", sprint(show, calc))
        end
        fresh = sYlmCalculator(R₂, ℓmax, range)
        @test all(isequal(Array(copy(b)), Array(copy(f))) for ((_, b), (_, f)) ∈ zip(calc, fresh))
        # The matrix form writes every spin weight, so the block of ℓₘₐₓ is held afterwards
        Ym = zeros(ComplexF64, length(range), Ysize(ℓlow, ℓmax))
        sYlm!(Ym, calc, R; ℓₘᵢₙ=ℓlow)
        @test ℓ(calc) == ℓₘₐₓ(calc)
        @test isequal(Ym, array_view(sYlm(R, ℓmax, range; ℓₘᵢₙ=ℓlow)))

        # The same for the real flavor
        cλ = sλlmCalculator(0.4, ℓmax, range)
        for s ∈ range
            Y = zeros(Float64, Ysize(ℓlow, ℓmax))
            sλlm!(Y, cλ, 1.3, s; ℓₘᵢₙ=ℓlow)
            @test isequal(Y, array_view(sλlm(1.3, ℓmax, s; ℓₘᵢₙ=ℓlow)))
        end

        # `spin_row!` steps the calculator and returns the one row, as a block of the shape
        # `recurrence!` would give for that spin weight alone, batched or not
        for data ∈ (R₂, [R, R₂])
            c = sYlmCalculator(data, ℓmax, range)
            f = sYlmCalculator(data, ℓmax, range)
            for ℓ′ ∈ ℓₘᵢₙ(c):ℓₘₐₓ(c), (i, s) ∈ enumerate(range)
                row = spin_row!(c, ℓ′, i)
                ref = recurrence!(f, ℓ′)
                expected = data isa Rotor ? Array(ref)[i, :] : Array(ref)[:, i, :]
                @test isequal(Array(row), expected)
            end
        end
    end
end


@testitem "sYlm_matrix: the rows are sYlm, bit for bit" begin
    import SphericalFunctions: sYlm_matrix, sλlm_matrix, sYlm, sλlm
    using Quaternionic: Rotor
    import Random

    # The values of each ℓ are written straight into the result, for one spin weight or for
    # a range of them, and so are the zeros below |s| of a range that straddles zero.
    rng = Random.Xoshiro(9)
    R⃗ = randn(rng, Rotor{Float64}, 5)
    θ⃗ = [0.1, 0.9, 2.5]
    for (ℓmax, spinsets) ∈ ((8, (0, -2, 3, -2:2, 1:3)), (15//2, (1//2, -3//2, -3//2:5//2)))
        for s ∈ spinsets
            Y = sYlm_matrix(R⃗, ℓmax, s)
            Λ = sλlm_matrix(θ⃗, ℓmax, s)
            if s isa AbstractUnitRange
                @test Y isa Array{ComplexF64, 3} && Λ isa Array{Float64, 3}
                @test all(isequal(Y[i, :, :], array_view(sYlm(R⃗[i], ℓmax, s))) for i ∈ 1:5)
                @test all(isequal(Λ[i, :, :], array_view(sλlm(θ⃗[i], ℓmax, s))) for i ∈ 1:3)
            else
                @test Y isa Matrix{ComplexF64} && Λ isa Matrix{Float64}
                @test all(isequal(Y[i, :], array_view(sYlm(R⃗[i], ℓmax, s))) for i ∈ 1:5)
                @test all(isequal(Λ[i, :], array_view(sλlm(θ⃗[i], ℓmax, s))) for i ∈ 1:3)
            end
            ℓlow = ℓmax isa Integer ? 0 : 1//2
            Y₀ = sYlm_matrix(R⃗, ℓmax, s; ℓₘᵢₙ=ℓlow)
            @test isequal(
                s isa AbstractUnitRange ? Y₀[2, :, :] : Y₀[2, :],
                array_view(sYlm(R⃗[2], ℓmax, s; ℓₘᵢₙ=ℓlow))
            )
        end
    end
end
