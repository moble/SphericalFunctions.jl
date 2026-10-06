# Tests of the spin-weighted spherical harmonics layer: `sYlmCalculator`, `sYlm`, and
# `sYlm_matrix`, against the closed-form expression for ₛYₗₘ and the defining relation to
# the Wigner 𝔇 matrices.

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
    import SphericalFunctions: sYlmCalculator, sYlm, sYlm_matrix, recurrence!, set_R!
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
    )
        @test refuses(f, ArgumentError, "runs downward or is empty")
    end
    calc = sYlmCalculator(R, 4, -2:2)
    # A spin weight the calculator does not serve is out of bounds of the block it returns
    @test_throws BoundsError recurrence!(calc, 2)[3, :]
    @test refuses(() -> recurrence!(calc, 5), ArgumentError, "out of bounds")
    @test refuses(() -> recurrence!(calc, R, -1), ArgumentError, "out of bounds")
    @test refuses(
        () -> recurrence!(calc, 2.0), ArgumentError,
        "The indices of this `sYlmCalculator` are integers of type `Int`, like 3; "
        * "got ℓ = 2.0::Float64"
    )
    @test refuses(() -> recurrence!(calc, Int16(2)), ArgumentError, "narrower than `Int`")
    # A complex "phase" is not a valid rotor for an sYlmCalculator
    @test refuses(() -> recurrence!(calc, cis(0.3), 2), ArgumentError, "rotors")
    @test refuses(() -> recurrence!(calc, [cis(0.3)], 2), ArgumentError, "rotors")
    @test refuses(
        () -> recurrence!(calc, [R, R], 2), DimensionMismatch,
        "This calculator handles Nᵣ=1 rotors, but got 2."
    )
    batched = sYlmCalculator([R, R, R], 4, -2:2)
    @test refuses(
        () -> recurrence!(batched, R, 2), DimensionMismatch,
        "This calculator handles Nᵣ=3 rotors, but a single rotor was given."
    )
    @test refuses(
        () -> recurrence!(batched, [R, R], 2), DimensionMismatch,
        "This calculator handles Nᵣ=3 rotors, but got 2."
    )
    @test refuses(() -> sYlm(R, 2, 3), ArgumentError, "|s|=3 exceeds ℓₘₐₓ=2")
    # A negative ℓₘₐₓ is named as such, rather than as a spin weight too large for it
    @test refuses(() -> sYlm(R, -1, 0), ArgumentError, "ℓₘₐₓ=-1 must be at least 0")
    # The message about ℓₘᵢₙ names the floor of the integer kind, and the bound ℓₘₐₓ
    @test refuses(
        () -> sYlm(R, 2, 0; ℓₘᵢₙ=-1), ArgumentError,
        "ℓₘᵢₙ=-1 must satisfy 0 ≤ ℓₘᵢₙ ≤ ℓₘₐₓ=2."
    )
    @test refuses(() -> sYlm(R, 2, 1; ℓₘᵢₙ=3), ArgumentError, "0 ≤ ℓₘᵢₙ ≤ ℓₘₐₓ=2")
    # A rotor of another float type cannot be pushed through a calculator, which works in
    # the type of the rotor it was built from.  (`check_rotor_type` owns this message.)
    @test refuses(
        () -> set_R!(calc, Rotor{Float32}(R)), ArgumentError,
        "given data would give Float32"
    )
end

@testitem "sYlmCalculator is reused with set_R!" begin
    import SphericalFunctions: sYlmCalculator, sYlm, set_R!, array_view, spins
    using Quaternionic: Rotor
    using Random
    rng = Random.Xoshiro(5)
    ℓₘₐₓ = 6
    Rs = randn(rng, Rotor{Float64}, 4)
    calc = sYlmCalculator(Rs[1], ℓₘₐₓ, -2:2)
    for R ∈ Rs
        # After `set_R!` each block is that of `sYlm` at the new rotor, for every spin
        # weight
        @test set_R!(calc, R) === calc
        Y = sYlm(R, ℓₘₐₓ, -2:2; ℓₘᵢₙ=0)
        Ys = [sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=0) for s ∈ -2:2]
        for (ℓ, block) ∈ calc
            @test block == Y[ℓ]
            @test all(array_view(block[s, :]) == array_view(Ys[s + 3][ℓ]) for s ∈ -2:2)
        end
    end
    # Allocation-free after warm-up.  This is measured inside a function, as a loop reusing
    # the calculator runs.  Before Julia 1.12 a loop that reads the blocks costs a few small
    # allocations in all, rather than any per ℓ; from 1.12 on it costs nothing.
    R = randn(rng, Rotor{Float64})
    function sweep!(calc, R, s)
        set_R!(calc, R)
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[s, ℓ])
        end
        t
    end
    s = first(spins(calc))
    sweep!(calc, R, s)
    @test @allocated(sweep!(calc, R, s)) ≤ (VERSION ≥ v"1.12" ? 0 : 128)
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
                # checked rotors are the plain ones converted, and `Checked` performs each
                # operation in `Float64`, so the values agree bit for bit.
                err = 0.0
                for i ∈ 1:Nᵣ, m ∈ -ℓ:ℓ
                    # A vector of rotors, even of one, gives batched blocks
                    z = blk[i, s, m]
                    zref = refblk[i, s, m]
                    err = max(err, abs(unchecked(real(z)) - real(zref)), abs(unchecked(imag(z)) - imag(zref)))
                end
                @test err == 0
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
    import SphericalFunctions: sYlm, sYlmCalculator, sλlmCalculator, recurrence!
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
        λ = sλlmCalculator(θ, ℓₘₐₓ, s)
        for ℓ ∈ abs(s):ℓₘₐₓ
            block = recurrence!(λ, ℓ)
            for m ∈ -ℓ:ℓ
                @test block[m] ≈ oracle(R, ℓ, m, s) / cispi(Float64(s)) atol=20eps()
            end
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

@testitem "sYlmCalculator half-integer reuse equals sYlm" setup=[RefusalChecks] begin
    import SphericalFunctions: sYlm, sYlmCalculator, set_R!, recurrence!, array_view, spins
    using Quaternionic: Rotor, from_spherical_coordinates

    ℓₘₐₓ = 9//2
    Rs = [from_spherical_coordinates(θ, ϕ) for (θ, ϕ) ∈ ((0.0, 0.0), (0.7, 1.2), (2.2, 4.0), (π, 0.3))]
    calc = sYlmCalculator(Rs[1], ℓₘₐₓ, -3//2:3//2)
    for R ∈ Rs
        # After `set_R!` each spin weight's row of each block is that of `sYlm` at the new
        # rotor, with the indices spelled as `Rational`s
        @test set_R!(calc, R) === calc
        Ys = [sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=1//2) for s ∈ -3//2:3//2]
        for (ℓ, block) ∈ calc, (i, s) ∈ enumerate(-3//2:3//2)
            @test array_view(block[s, :]) == array_view(Ys[i][ℓ])
        end
    end
    # Allocation-free after warm-up, as for the integer kind (and measured inside a function
    # for the same reason)
    R = Rs[2]
    function sweep!(calc, R, s)
        set_R!(calc, R)
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[s, ℓ])
        end
        t
    end
    s = first(spins(calc))
    sweep!(calc, R, s)
    @test @allocated(sweep!(calc, R, s)) ≤ (VERSION ≥ v"1.12" ? 0 : 128)
    # An ℓ of the wrong kind is refused with a message naming the kind of the calculator,
    # rather than with a bare conversion error
    @test refuses(
        () -> recurrence!(calc, 2), ArgumentError,
        "The indices of this `sYlmCalculator` are half-odd-integers, each a "
        * "`HalfOddInteger` or a `Rational{Int}` with denominator 2, like 7//2; got ℓ = 2"
    )
    icalc = sYlmCalculator(R, 4, 1)
    @test refuses(
        () -> recurrence!(icalc, 1//2), ArgumentError,
        "The indices of this `sYlmCalculator` are integers of type `Int`, like 3; "
        * "got ℓ = 1//2"
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
    import SphericalFunctions: sYlm, sYlm_matrix, sYlmCalculator, HalfOddInteger
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
    # ... as is the keyword's ASCII spelling
    @test array_view(sYlm(R, 7//2, 3//2; ell_min=1//2)) ==
        array_view(sYlm(R, ℓₘₐₓ, HalfOddInteger(3//2); ℓₘᵢₙ))

    # A mixture of the two kinds of positional index is refused with a message that lists
    # the indices, names both spellings of a half-odd-integer, and says which kind each
    # index is
    mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
    for f ∈ (
        () -> sYlm(R, 7//2, 1), () -> sYlm(R, 4, 1//2), () -> sYlm(R, 4, HalfOddInteger(1//2)),
        () -> sYlm_matrix(Rs, 7//2, 1), () -> sYlm_matrix(Rs, 4, 1//2),
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
            () -> sYlmCalculator(R, IT(4), IT(1)), () -> sYlmCalculator(R, IT(4), -1:1),
        )
            @test refuses(f, ArgumentError, sentence)
        end
        @test refuses(() -> sYlm(R, 4, 1; ℓₘᵢₙ=IT(2)), ArgumentError, "keyword argument `ℓₘᵢₙ`")
        @test refuses(() -> sYlm(R, 4, 1; ℓₘᵢₙ=IT(2)), ArgumentError, sentence)
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
    import SphericalFunctions: sYlm, sYlm_matrix, Ysize, spins, sYlmCalculator, set_R!,
        array_view
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
        # A calculator built for the range, and given another rotor with `set_R!`, gives the
        # same values one block at a time, and the row of one spin weight in each
        calc = set_R!(sYlmCalculator(R, ℓₘₐₓ, srange), rotors[2])
        Y₂ = sYlm(rotors[2], ℓₘₐₓ, srange)
        Y₁ = sYlm(rotors[2], ℓₘₐₓ, first(sr); ℓₘᵢₙ)
        @test all(b == Y₂[ℓ] for (ℓ, b) ∈ calc if ℓ ∈ keys(Y₂))
        @test all(
            array_view(b[first(sr), :]) == array_view(Y₁[ℓ])
            for (ℓ, b) ∈ calc if ℓ ∈ keys(Y₁)
        )

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
# the unweighted `Ylm` wrapper in its vector form; `similar(calc, R)`, which rebuilds a
# calculator for new rotors only when there are as many as it was built for; and a
# calculator's refusal of mode weights of a spin weight it does not serve.

@testitem "sYlm: the `Ylm` wrapper and the calculator refusals" setup=[RefusalChecks] begin
    using Quaternionic: Rotor, RotorF64
    import SphericalFunctions: Nᵣ, spins, ModeWeights, Ysize
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

    # A calculator built for one spin weight refuses mode weights of another, and names the
    # ones it serves
    cs = sYlmCalculator(R, ℓₘₐₓ, -2)
    @test spins(cs) == -2:-2
    w₁ = ModeWeights(randn(rng, ComplexF64, Ysize(1, ℓₘₐₓ)), 1)
    @test refuses(
        () -> cs * w₁, ArgumentError,
        "spin weight s=1, but this calculator serves only s=-2."
    )

    crange = sYlmCalculator(R, ℓₘₐₓ, -2:2)
    @test spins(crange) == -2:2
    # a spin weight inside the range is served, and agrees with computing it afresh
    @test crange * w₁ ≈ sYlm(R, ℓₘₐₓ, 1) * w₁
    w₃ = ModeWeights(zeros(ComplexF64, Ysize(3, ℓₘₐₓ)), 3)
    @test refuses(() -> crange * w₃, ArgumentError, "this calculator serves only s ∈ -2:2.")
end


@testitem "sYlmCalculator: each block element is the wedge element `wedge_value` reads" begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, sλlmCalculator, recurrence!, wedge_value,
        power_column, sYlm_coefficient, spins, Nᵣ, ℓₘᵢₙ, ℓₘₐₓ, floattype, number_type,
        block_array, nspins
    using Quaternionic: Rotor
    import Random

    # As for the Wigner calculators, the elements H[m, -s] are read in runs along the rows
    # of the wedge, and in one of two orders according to the numbers of rotors and of spin
    # weights, and to whether the phases are computed; the batches of 65 with several spin
    # weights take the order of storage.  Each value is compared here, bit for bit, with the
    # one built from `wedge_value` and the same power tables, and each row below ℓ < |s| is
    # zero.
    rng = Random.Xoshiro(20260924)
    R⃗ = randn(rng, Rotor{Float64}, 4)
    θ⃗ = [0.0, 0.7, 2.2, π]
    R⃗₆₅ = randn(rng, Rotor{Float64}, 65)
    θ⃗₆₅ = [θ⃗; rand(rng, 61) .* π]
    # The power zᵏ of rotor iᵣ in a power table, for k of either sign
    zpower(Z, iᵣ, k) = Z[iᵣ, power_column(Z, k)]
    function expected(calc, H, iᵣ, s, m, ℓ)
        NT, RT = number_type(calc), floattype(calc)
        abs(s) > ℓ && return zero(NT)
        prefactor = √((2ℓ + 1) / (4 * RT(π)))
        value = sYlm_coefficient(NT, RT, 1, m, s, prefactor) * wedge_value(H, iᵣ, m, -s)
        if NT <: Complex && calc.phases[]
            value * (zpower(calc.engine.Z₊, iᵣ, m - s) * zpower(calc.engine.Z₋, iᵣ, m + s))
        else
            value
        end
    end
    count = Ref(0)
    for (ℓmax, spinsets) ∈ ((9, (0, 2, -3, -2:2, 0:3, -9:9)), (17//2, (1//2, -5//2, -3//2:1//2)))
        for s ∈ spinsets, (Ctor, data) ∈ (
            (sYlmCalculator, R⃗), (sYlmCalculator, R⃗[1]), (sYlmCalculator, θ⃗),
            (sλlmCalculator, θ⃗), (sλlmCalculator, θ⃗[2]),
            (sYlmCalculator, R⃗₆₅), (sYlmCalculator, θ⃗₆₅), (sλlmCalculator, θ⃗₆₅),
        )
            calc = Ctor(data, ℓmax, s)
            for ℓ ∈ ℓₘᵢₙ(calc):ℓₘₐₓ(calc)
                recurrence!(calc, ℓ)
                H = calc.engine.H.Hˡ
                Y = block_array(calc, calc.Yˡ, ℓ, Base.OneTo(nspins(calc.s)))  # [iᵣ, s, m]
                good = true
                for (i, sᵢ) ∈ enumerate(spins(calc)), m ∈ -ℓ:ℓ, iᵣ ∈ 1:Nᵣ(calc)
                    value = Y[iᵣ, i, Int(m + ℓ) + 1]
                    good &= isequal(value, expected(calc, H, iᵣ, sᵢ, m, ℓ))
                    count[] += 1
                end
                @test good
            end
        end
    end
    @test count[] > 20_000
end


@testitem "compute_block! writes a harmonic block into any destination" begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, sλlmCalculator, recurrence!, compute_block!,
        block_array, nspins, ℓₘᵢₙ, ℓₘₐₓ
    using Quaternionic: Rotor
    import Random

    # `compute_block!(calc, ℓ, is, A, o)` writes the spin rows `is` of the block of ℓ
    # densely, as [iᵣ, s ∈ is, m], into `A` after its first `o` entries, which is how `sYlm`
    # and `sYlm_matrix` fill their results.  The rows written are bit for bit those of the
    # calculator's own block, whether all of them or one, nothing outside them is touched,
    # and the calculator, whose own buffer was not written, no longer claims to hold a
    # block.
    rng = Random.Xoshiro(20261001)
    R⃗ = randn(rng, Rotor{Float64}, 3)
    θ⃗ = [0.0, 0.7, π]
    for (ℓmax, spinsets) ∈ ((6, (-2, -2:2)), (11//2, (1//2, -3//2:3//2))), s ∈ spinsets,
            (Ctor, data) ∈ (
                (sYlmCalculator, R⃗), (sYlmCalculator, R⃗[1]), (sλlmCalculator, θ⃗),
                (sλlmCalculator, θ⃗[2]),
            )
        calc = Ctor(data, ℓmax, s)
        NT = eltype(calc.Yˡ)
        good = true
        for ℓ ∈ ℓₘᵢₙ(calc):ℓₘₐₓ(calc)
            recurrence!(calc, ℓ)
            full = copy(block_array(calc, calc.Yˡ, ℓ, Base.OneTo(nspins(calc.s))))
            for is ∈ (Base.OneTo(nspins(calc.s)), nspins(calc.s):nspins(calc.s), 1:1)
                rows = vec(full[:, is, :])
                n = length(rows)
                A = fill(NT(NaN), n + 8)
                compute_block!(calc, ℓ, is, A, 3)
                good &= isequal(A[4:(3 + n)], rows)
                good &= all(isnan, A[1:3]) && all(isnan, A[(4 + n):end])
                good &= SphericalFunctions.ℓ(calc) < ℓₘᵢₙ(calc)
            end
        end
        @test good
    end
end


@testitem "sYlmCalculator: one spin weight read from a calculator built for several" begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, sλlmCalculator, sYlm, set_R!, set_θ!,
        recurrence!, spin_row!, ℓ, ℓₘᵢₙ, ℓₘₐₓ
    using Quaternionic: Rotor
    import Random

    # Reading one spin weight out of a calculator built for a range, as `calc * w` does with
    # `spin_row!`, assembles only that spin weight's values at each ℓ.  They are the values
    # the calculator would otherwise have given, bit for bit; the calculator then holds no
    # complete block, and says so, and the next full step gives the whole block again.
    rng = Random.Xoshiro(8)
    R, R₂ = randn(rng, Rotor{Float64}, 2)
    for (ℓmax, range, ℓlow) ∈ ((7, -3:3, 0), (15//2, -5//2:3//2, 1//2))
        calc = set_R!(sYlmCalculator(R, ℓmax, range), R₂)
        for (i, s) ∈ enumerate(range)
            Y = sYlm(R₂, ℓmax, s; ℓₘᵢₙ=ℓlow)
            for ℓ′ ∈ ℓlow:ℓmax
                @test isequal(Array(spin_row!(calc, ℓ′, i)), Array(Y[ℓ′]))
            end
            @test ℓ(calc) < ℓₘᵢₙ(calc)
            @test occursin("nothing computed yet", sprint(show, calc))
        end
        fresh = sYlmCalculator(R₂, ℓmax, range)
        @test all(isequal(Array(copy(b)), Array(copy(f))) for ((_, b), (_, f)) ∈ zip(calc, fresh))
        # Iteration writes every spin weight at each step, so the block of ℓₘₐₓ is held
        # afterwards
        @test ℓ(calc) == ℓₘₐₓ(calc)

        # The same for the real flavor
        cλ = set_θ!(sλlmCalculator(0.4, ℓmax, range), 1.3)
        for (i, s) ∈ enumerate(range)
            single = sλlmCalculator(1.3, ℓmax, s)
            for ℓ′ ∈ ℓₘᵢₙ(cλ):ℓₘₐₓ(cλ)
                @test isequal(Array(spin_row!(cλ, ℓ′, i)), Array(recurrence!(single, ℓ′)))
            end
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
    import SphericalFunctions: sYlm_matrix, sYlm
    using Quaternionic: Rotor
    import Random

    # The values of each ℓ are written straight into the result, for one spin weight or for
    # a range of them, and so are the zeros below |s| of a range that straddles zero.
    rng = Random.Xoshiro(9)
    R⃗ = randn(rng, Rotor{Float64}, 5)
    for (ℓmax, spinsets) ∈ ((8, (0, -2, 3, -2:2, 1:3)), (15//2, (1//2, -3//2, -3//2:5//2)))
        for s ∈ spinsets
            Y = sYlm_matrix(R⃗, ℓmax, s)
            if s isa AbstractUnitRange
                @test Y isa Array{ComplexF64, 3}
                @test all(isequal(Y[i, :, :], array_view(sYlm(R⃗[i], ℓmax, s))) for i ∈ 1:5)
            else
                @test Y isa Matrix{ComplexF64}
                @test all(isequal(Y[i, :], array_view(sYlm(R⃗[i], ℓmax, s))) for i ∈ 1:5)
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
