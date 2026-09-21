# Tests of the spin-weighted spherical harmonics layer: `sYlmCalculator`, `sYlm`, `sYlm!`,
# and `sYlm_matrix`, against the closed-form expression for ₛYₗₘ and the defining relation
# to the Wigner 𝔇 matrices.

@testitem "sYlm vs the closed form, for arbitrary rotors" setup=[Utilities] begin
    # (The Utilities snippet defines a closed-form `sYlm`, so the package function is qualified.)
    import SphericalFunctions
    import SphericalFunctions: Yindex, Ysize
    using Quaternionic: Quaternion, Rotor, components, 𝐢, 𝐣, 𝐤
    using DoubleFloats: Double64
    using Random

    Random.seed!(17)  # `Rrange` draws from the default RNG

    # Reference, category 1 (a closed-form formula, from the conventions pages).  The
    # `Utilities` snippet's `sYlm(s, ℓ, m, θ, ϕ)` is the explicit sum of Ajith et al.,
    # Eqs. (II.7)-(II.8), evaluated on the sphere — that is, at the rotor
    # `from_spherical_coordinates(θ, ϕ)`, whose Euler angles are (ϕ, θ, 0).  A general
    # rotor has a third angle, and ₛYₗₘ depends on it through the spin weight alone:
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
                sYlm(s, ℓ, m, β, α) * cis(-s * γ)
                for ℓ in abs(s):ℓₘₐₓ for m in -ℓ:ℓ
            ]
        end
        Complex{T}.(v)  # in `Yindex(ℓ, m, abs(s))` order
    end

    for T ∈ (Float64, Double64)
        # Measured worst error over everything below: 5.0 eps (Float64), 2.4 eps (Double64).
        ϵ = 20 * eps(T)
        for ℓₘₐₓ ∈ (0, 1, 2, 5, 9)
            for s ∈ -min(3, ℓₘₐₓ):min(3, ℓₘₐₓ)
                for R ∈ Rrange(T, 8)
                    Y = strided(SphericalFunctions.sYlm(R, ℓₘₐₓ, s))
                    @test eltype(Y) === Complex{T}
                    @test length(Y) == Ysize(abs(s), ℓₘₐₓ)
                    # Accumulated and asserted once per rotor: an engine that is wrong
                    # everywhere would otherwise print one failure per (ℓ, m).
                    errY = maximum(abs, Y .- Yref(T, R, s, ℓₘₐₓ))
                    @test errY ≤ ϵ
                    Y₀ = strided(SphericalFunctions.sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=0))
                    @test length(Y₀) == Ysize(0, ℓₘₐₓ)
                    @test all(iszero, Y₀[1:s^2])
                    @test Y₀[s^2+1:end] == Y
                end
            end
        end
    end
end

@testitem "sYlm vs closed form and definition" setup=[Utilities] begin
    # (The Utilities snippet defines a closed-form `sYlm`, so the package function is qualified.)
    import SphericalFunctions
    import SphericalFunctions: Yindex, D
    using Quaternionic
    using Random
    rng = Random.Xoshiro(42)
    ℓₘₐₓ = 5
    for (θ, ϕ) ∈ θϕrange(Float64, 6)
        R = Rotor(from_spherical_coordinates(θ, ϕ))
        for s ∈ -2:2
            Y = strided(SphericalFunctions.sYlm(R, ℓₘₐₓ, s))
            for ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ
                # Closed form on the sphere (from the Utilities snippet)
                @test Y[Yindex(ℓ, m, abs(s))] ≈ sYlm(s, ℓ, m, θ, ϕ) atol=1e-13 rtol=1e-13
            end
        end
    end
    # The definition ₛYₗₘ = (-1)^s √((2ℓ+1)/4π) conj(𝔇ˡₘ,₋ₛ) for random rotors
    for R ∈ randn(rng, Rotor{Float64}, 5)
        𝔇 = D(R, ℓₘₐₓ)
        for s ∈ -2:2
            Y = strided(SphericalFunctions.sYlm(R, ℓₘₐₓ, s))
            for ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ
                @test Y[Yindex(ℓ, m, abs(s))] ≈ (-1)^s * √((2ℓ+1)/(4π)) * conj(𝔇[ℓ][m, -s]) atol=1e-14
            end
        end
    end
end

@testitem "sYlm_matrix" setup=[Utilities] begin
    # (The Utilities snippet defines a closed-form `sYlm`, so the package function is qualified.)
    import SphericalFunctions
    import SphericalFunctions: sYlm_matrix, Ysize
    using Quaternionic: Rotor, Quaternion, components
    using DoubleFloats: Double64
    using Random
    rng = Random.Xoshiro(7)

    # Reference, category 1 (a closed-form formula, from the conventions pages): the same
    # closed form and the same Euler-angle bookkeeping as the "sYlm vs the closed form"
    # item above, evaluated at four times the working precision and rounded.
    function euler_angles(R)
        w, x, y, z = BigFloat.(components(Quaternion(R)))
        ϕₛ, ϕₐ = angle(Complex(w, z)), angle(Complex(y, x))
        (ϕₛ - ϕₐ, 2 * atan(abs(Complex(y, x)), abs(Complex(w, z))), ϕₛ + ϕₐ)
    end
    function Yref(::Type{T}, Rs, s, ℓₘₐₓ) where {T}
        rows = setprecision(BigFloat, 4 * precision(T) + 64) do
            map(Rs) do R
                α, β, γ = euler_angles(R)
                [sYlm(s, ℓ, m, β, α) * cis(-s * γ) for ℓ in abs(s):ℓₘₐₓ for m in -ℓ:ℓ]
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
                @test M[i, :] == strided(SphericalFunctions.sYlm(R, ℓₘₐₓ, s))
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
    @test_throws "exceeds ℓₘₐₓ" sYlm_matrix(randn(rng, Rotor{Float64}, 3), 2, 3)
    @test_throws "exceeds ℓₘₐₓ" SphericalFunctions.sYlm(randn(rng, Rotor{Float64}), 2, 3)
end

@testitem "sYlmCalculator all spins and batches" begin
    import SphericalFunctions
    import SphericalFunctions: sYlmCalculator, sYlm, recurrence!, Yindex
    using Quaternionic: Rotor
    import SphericalFunctions: DegreeBlockBatch, strided
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
    singles = Dict((i, s) => strided(sYlm(Rs[i], ℓₘₐₓ, s; ℓₘᵢₙ=0)) for i ∈ 1:Nᵣ for s ∈ -sₘₐₓ:sₘₐₓ)
    for ℓ ∈ 0:ℓₘₐₓ
        recurrence!(calc, ℓ)
        @test SphericalFunctions.ℓ(calc) == ℓ
        for s ∈ -sₘₐₓ:sₘₐₓ
            blk = calc[ℓ, s]
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
        recurrence!(calc, ℓ)
        for s ∈ (-2, 0, 3)
            blk = calc[ℓ, s]
            for i ∈ 1:Nᵣ, m ∈ -ℓ:ℓ
                @test blk[i, m] == (ℓ < abs(s) ? 0 : singles[(i, s)][Yindex(ℓ, m)])
            end
        end
    end
    # copy keeps the natural axes and survives the next recurrence!; collect is 1-based
    recurrence!(calc, 4)
    v = calc[4, 1]
    c = copy(v)
    a = collect(v)
    @test axes(c) == axes(v)
    @test a isa Matrix{ComplexF64} && size(a) == (Nᵣ, 9)
    recurrence!(calc, 5)
    @test strided(c) == [singles[(i, 1)][Yindex(4, m)] for i ∈ 1:Nᵣ, m ∈ -4:4]
    # Nᵣ == 1: a single rotor and a vector view
    c1 = sYlmCalculator(Rs[2], ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
    recurrence!(c1, 3)
    @test axes(c1[3, -1]) == (-3:3,)
    @test collect(c1[3, -1]) == singles[(2, -1)][Yindex(3, -3):Yindex(3, 3)]
    # Angle input, given at construction, evaluates at (θ, ϕ=0)
    θs = [0.0, 0.7, 1.9, π]
    cθ = sYlmCalculator(θs, 4, -2:2)
    recurrence!(cθ, 2)
    using Quaternionic: from_spherical_coordinates
    for (i, θ) ∈ enumerate(θs), s ∈ -2:2, m ∈ -2:2
        Yref = strided(sYlm(Rotor(from_spherical_coordinates(θ, 0.0)), 4, s; ℓₘᵢₙ=0))[Yindex(2, m)]
        @test cθ[2, s][i, m] ≈ Yref atol=1e-15
        @test imag(cθ[2, s][i, m]) == 0
    end
    # similar and show
    c2 = similar(calc)
    @test SphericalFunctions.Nᵣ(c2) == Nᵣ && SphericalFunctions.ℓₘₐₓ(c2) == ℓₘₐₓ && SphericalFunctions.spins(c2) == -sₘₐₓ:sₘₐₓ
    @test c2.Yˡ !== calc.Yˡ
    str = sprint(show, calc)
    @test occursin("sYlmCalculator", str) && occursin("ℓₘₐₓ=$ℓₘₐₓ", str) && occursin("Nᵣ=$Nᵣ", str)
    @test sprint(show, MIME("text/plain"), c2) isa String
end

@testitem "sYlmCalculator errors" begin
    import SphericalFunctions: sYlmCalculator, sYlm, sYlm!, recurrence!
    using Quaternionic: Rotor
    using Random
    rng = Random.Xoshiro(3)
    R = randn(rng, Rotor{Float64})
    @test_throws "|s| ≤ ℓₘₐₓ" sYlmCalculator(R, 3, 4)
    @test_throws "|s| ≤ ℓₘₐₓ" sYlmCalculator(R, 3, -4:4)
    # A negative spin weight is an ordinary one; what is refused is a range that is not a
    # consecutive run from low to high.
    @test SphericalFunctions.spin(sYlmCalculator(R, 3, -1)) == -1
    @test_throws "step must be 1" sYlmCalculator(R, 3, -2:2:2)
    @test_throws "runs downward or is empty" sYlmCalculator(R, 3, 2:-2)
    @test_throws "runs downward or is empty" sYlmCalculator(R, 3, 2:-1:-2)
    @test_throws "runs downward or is empty" sYlmCalculator(R, 7//2, 3//2:-1:-3//2)
    calc = sYlmCalculator(R, 4, -2:2)
    @test_throws "nothing has been computed" calc[0, 0]
    recurrence!(calc, 2)
    @test_throws "not among them" calc[2, 3]
    @test_throws "currently holds ℓ=2" calc[1, 0]
    @test_throws "out of bounds" calc[5, 0]
    @test_throws "out of bounds" recurrence!(calc, 5)
    @test_throws "out of bounds" recurrence!(calc, R, -1)
    # A complex "phase" is not a valid rotor for an sYlmCalculator
    @test_throws "rotors" recurrence!(calc, cis(0.3), 2)
    @test_throws "rotors" recurrence!(calc, [cis(0.3)], 2)
    @test_throws "Expected 1 rotors" recurrence!(calc, [R, R], 2)
    batched = sYlmCalculator([R, R, R], 4, -2:2)
    @test_throws "expects Nᵣ=3" recurrence!(batched, R, 2)
    @test_throws "Expected 3 rotors" recurrence!(batched, [R, R], 2)
    @test_throws "exceeds ℓₘₐₓ" strided(sYlm(R, 2, 3))
    # The message about ℓₘᵢₙ names the floor of the integer kind
    @test_throws "ℓₘᵢₙ=-1 must satisfy 0 ≤ ℓₘᵢₙ ≤ max(|s|, ℓₘₐₓ)." strided(sYlm(R, 2, 0; ℓₘᵢₙ=-1))
    @test_throws "0 ≤ ℓₘᵢₙ" strided(sYlm(R, 2, 1; ℓₘᵢₙ=3))
    Y = zeros(ComplexF64, 5)
    @test_throws "Output vector has length" sYlm!(Y, R, 3, 0)
    @test_throws "not among them" sYlm!(zeros(ComplexF64, 25), sYlmCalculator(R, 4, 1), R, 2)
    @test_throws "Nᵣ=1" sYlm!(zeros(ComplexF64, 25), batched, R, 1)
    # The output's element type must be the calculator's own; it no longer decides the type
    @test_throws "element type must be Complex{Float64}" sYlm!(zeros(ComplexF32, 25), calc, R, 1; ℓₘᵢₙ=0)
    @test_throws "element type must be Complex{Float64}" sYlm!(zeros(ComplexF32, 25), R, 4, 1; ℓₘᵢₙ=0)
    # ... and a rotor of another float type cannot be pushed through a calculator
    @test_throws "given data would give Float32" sYlm!(zeros(ComplexF64, 25), calc, Rotor{Float32}(R), 1; ℓₘᵢₙ=0)
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
        @test sYlm!(Y, calc, R, s; ℓₘᵢₙ=0) === Y
        @test Y == strided(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=0))
        n = Ysize(abs(s), ℓₘₐₓ)
        sYlm!(Y, calc, R, s)
        @test Y[1:n] == strided(sYlm(R, ℓₘₐₓ, s))
    end
    # Allocation-free after warm-up
    R = randn(rng, Rotor{Float64})
    sYlm!(Y, calc, R, 1; ℓₘᵢₙ=0)
    @test @allocated(sYlm!(Y, calc, R, 1; ℓₘᵢₙ=0)) == 0
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
    # Float32
    Y32 = strided(sYlm(Rotor{Float32}(R64), 20, -1))
    Y64 = strided(sYlm(R64, 20, -1))
    @test eltype(Y32) === ComplexF32
    @test all(isfinite, Y32)
    @test maximum(abs(Y32[i] - Y64[i]) / abs(Y64[i]) for i ∈ eachindex(Y64) if abs(Y64[i]) > 1e-3) < 1e-4
    # BigFloat vs Double64
    YB = strided(sYlm(Rotor{BigFloat}(R64), 4, 2))
    YD = strided(sYlm(Rotor{Double64}(R64), 4, 2))
    @test maximum(abs, YB .- YD) < 1e-30
    # ForwardDiff through the rotor.  ϕ = 0 is included deliberately: there the spinor
    # phase `z₊` is exactly 1, and a `sqrt` of an exact zero inside `complex_powers!` used
    # to make every derivative NaN.  It is also the case the ring-based transforms use.
    θ₀ = 0.8
    for ϕ ∈ (0.0, 0.3)
        f(θ) = real(strided(sYlm(Rotor(from_spherical_coordinates(θ, ϕ)), 3, 1))[Yindex(3, 2, 1)])
        dual = ForwardDiff.derivative(f, θ₀)
        h = 1e-6
        fd = (f(θ₀ + h) - f(θ₀ - h)) / 2h
        @test isfinite(dual)
        @test abs(dual - fd) < 1e-6
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
    for (ℓₘₐₓ, sₘₐₓ, Nᵣ) ∈ ((0, 0, 1), (2, 1, 1), (5, 2, 3), (9, 3, 2))
        Rs = randn(rng, Rotor{Float64}, Nᵣ)
        # A calculator works in the float type of its rotors, so the checked type is applied
        # to them; `ref` keeps the plain-Float64 rotors and is the value to compare against.
        calc = sYlmCalculator(Rotor{NC}.(Rs), ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
        fill!(calc, NaN)  # after construction, which stores the rotor data `fill!` preserves
        ref = sYlmCalculator(Rs, ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
        for ℓ ∈ [0:ℓₘₐₓ; ℓₘₐₓ ÷ 2]
            recurrence!(calc, ℓ)
            recurrence!(ref, ℓ)
            for s ∈ -sₘₐₓ:sₘₐₓ
                blk = calc[ℓ, s]
                refblk = ref[ℓ, s]
                # Every entry must have been written: an untouched one still holds the
                # sentinel NaN, which turns `err` into NaN and fails the comparison.  The
                # values are only approximately equal to the plain-`Float64` run because
                # `@fastmath` has no effect on a wrapper type, so the two round differently.
                err = 0.0
                for i ∈ 1:Nᵣ, m ∈ -ℓ:ℓ
                    z = Nᵣ == 1 ? blk[m] : blk[i, m]
                    zref = Nᵣ == 1 ? refblk[m] : refblk[i, m]
                    err = max(err, abs(unchecked(real(z)) - real(zref)), abs(unchecked(imag(z)) - imag(zref)))
                end
                @test err < 1e-13
            end
        end
        # `fill!` keeps the stored rotor data, as its docstring promises, so the recurrence
        # can be re-run without re-supplying it — and must still write every element.
        fill!(calc, NaN)
        recurrence!(calc, ℓₘₐₓ)
        recurrence!(ref, ℓₘₐₓ)
        for s ∈ -sₘₐₓ:sₘₐₓ
            blk, refblk = calc[ℓₘₐₓ, s], ref[ℓₘₐₓ, s]
            err = 0.0
            for i ∈ 1:Nᵣ, m ∈ -ℓₘₐₓ:ℓₘₐₓ
                z = Nᵣ == 1 ? blk[m] : blk[i, m]
                zref = Nᵣ == 1 ? refblk[m] : refblk[i, m]
                err = max(err, abs(unchecked(real(z)) - real(zref)), abs(unchecked(imag(z)) - imag(zref)))
            end
            @test err < 1e-13
        end
    end
end


### Half-integer indices in the flat functions.
#
# The oracle throughout is `sYlmCalculator`, whose half-integer values are verified against
# the Wigner 𝔇 oracle in `test/wigner/half_integer.jl`; these items check that the flat
# functions lay those values out in the canonical ordering, accept the `Rational` spelling,
# refuse a mixture of the two kinds of index, and keep the two properties the transforms rest
# on — orthonormality on the sphere and antiperiodicity in ϕ.  Angles are fixed wherever a
# value is asserted exactly, so that a failure is reproducible.

@testitem "sYlm half-integer vs sYlmCalculator blocks" begin
    import SphericalFunctions: sYlm, sYlmCalculator, recurrence!, Ysize, Yindex, HalfOddInteger
    using Quaternionic: Rotor, from_euler_angles

    ℓₘₐₓ = 9//2
    Rs = [Rotor(from_euler_angles(α, β, γ)) for (α, β, γ) ∈ ((0.0, 0.0, 0.0), (0.7, 1.1, 2.3), (2.9, 0.4, 5.1), (4.0, 2.2, 0.3))]
    for R ∈ Rs
        calc = sYlmCalculator(R, ℓₘₐₓ, -3//2:3//2)
        for s ∈ (-3//2, -1//2, 1//2, 3//2), ℓₘᵢₙ ∈ (abs(s), 1//2)
            Y = strided(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ))
            @test eltype(Y) === ComplexF64
            @test length(Y) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            # The flat function and the calculator are the same engine, so the values agree
            # exactly, whatever the calculator's own sₘₐₓ.
            for ℓ ∈ ℓₘᵢₙ:1:ℓₘₐₓ
                recurrence!(calc, ℓ)
                blk = calc[ℓ, s]
                for m ∈ -ℓ:ℓ
                    @test Y[Yindex(ℓ, m, ℓₘᵢₙ)] == blk[m]
                end
            end
        end
        # The default ℓₘᵢₙ is |s|, for the half-integer kind as for the integer one.
        @test strided(sYlm(R, ℓₘₐₓ, 3//2)) == strided(sYlm(R, ℓₘₐₓ, 3//2; ℓₘᵢₙ=3//2))
        @test length(strided(sYlm(R, ℓₘₐₓ, 3//2))) == Ysize(3//2, ℓₘₐₓ)
    end
    # For half-integer s the values include the phase i^{2s} = ±i: the ϕ = γ = 0 values, which
    # are real for integer s, are here purely imaginary.
    Rθ = Rotor(from_euler_angles(0.0, 1.1, 0.0))
    for s ∈ (-1//2, 1//2, 3//2)
        Yθ = strided(sYlm(Rθ, ℓₘₐₓ, s))
        @test maximum(abs ∘ real, Yθ) == 0
        @test maximum(abs ∘ imag, Yθ) > 0.1
    end
end

@testitem "sYlm half-integer ℓₘᵢₙ = 1//2 gives zeros below |s|" begin
    import SphericalFunctions: sYlm, sYlm_matrix, Ysize, Yindex
    using Quaternionic: Rotor, from_spherical_coordinates

    ℓₘₐₓ = 9//2
    R = from_spherical_coordinates(0.7, 1.2)
    for s ∈ (-3//2, 3//2, 5//2)
        Y = strided(sYlm(R, ℓₘₐₓ, s))
        Y₀ = strided(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=1//2))
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
    # The floor of ℓₘᵢₙ is 1/2; anything below it, or above max(|s|, ℓₘₐₓ), is refused, with
    # a message that names the floor of this kind of index rather than the integer one.
    @test_throws "must satisfy" strided(sYlm(R, ℓₘₐₓ, 1//2; ℓₘᵢₙ=-1//2))
    @test_throws "must satisfy" strided(sYlm(R, ℓₘₐₓ, 1//2; ℓₘᵢₙ=11//2))
    @test_throws "ℓₘᵢₙ=-1//2 must satisfy 1//2 ≤ ℓₘᵢₙ ≤ max(|s|, ℓₘₐₓ)." strided(sYlm(R, ℓₘₐₓ, 1//2; ℓₘᵢₙ=-1//2))
    @test_throws "1//2 ≤ ℓₘᵢₙ" sYlm_matrix([R, -R], ℓₘₐₓ, 3//2; ℓₘᵢₙ=11//2)
    @test_throws "exceeds ℓₘₐₓ" strided(sYlm(R, 1//2, 3//2))
end

@testitem "sYlm! half-integer, both forms, equals sYlm" begin
    import SphericalFunctions: sYlm, sYlm!, sYlmCalculator, Ysize
    using Quaternionic: Rotor, from_spherical_coordinates

    ℓₘₐₓ = 9//2
    Rs = [from_spherical_coordinates(θ, ϕ) for (θ, ϕ) ∈ ((0.0, 0.0), (0.7, 1.2), (2.2, 4.0), (π, 0.3))]
    calc = sYlmCalculator(Rs[1], ℓₘₐₓ, -3//2:3//2)
    Y = Vector{ComplexF64}(undef, Ysize(1//2, ℓₘₐₓ))
    for R ∈ Rs, s ∈ (-3//2, -1//2, 1//2, 3//2)
        @test sYlm!(Y, calc, R, s; ℓₘᵢₙ=1//2) === Y
        @test Y == strided(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=1//2))
        n = Ysize(abs(s), ℓₘₐₓ)
        sYlm!(Y, calc, R, s)
        @test Y[1:n] == strided(sYlm(R, ℓₘₐₓ, s))
        # The allocating form, with the indices spelled as `Rational`s
        Y′ = Vector{ComplexF64}(undef, n)
        @test sYlm!(Y′, R, ℓₘₐₓ, s) === Y′
        @test Y′ == strided(sYlm(R, ℓₘₐₓ, s))
    end
    # Allocation-free after warm-up, as for the integer kind
    R = Rs[2]
    sYlm!(Y, calc, R, 1//2; ℓₘᵢₙ=1//2)
    @test @allocated(sYlm!(Y, calc, R, 1//2; ℓₘᵢₙ=1//2)) == 0
    # Errors: a spin weight the calculator does not serve, the output length, and a spin
    # weight or ℓₘᵢₙ of the wrong
    # kind for the calculator, which is refused with a message naming the calculator's kind
    # rather than with a bare conversion error.
    @test_throws "not among them" sYlm!(Y, calc, R, 5//2)
    @test_throws "Output vector has length" sYlm!(zeros(ComplexF64, 3), calc, R, 1//2)
    @test_throws "indices are half-odd-integers, like 7//2, so the spin weight s" sYlm!(Y, calc, R, 1)
    @test_throws "indices are half-odd-integers, like 7//2, so ℓₘᵢₙ" sYlm!(Y, calc, R, 1//2; ℓₘᵢₙ=0)
    icalc = sYlmCalculator(R, 4, 1)
    @test_throws "indices are integers, like 3, so the spin weight s" sYlm!(zeros(ComplexF64, 25), icalc, R, 1//2)
    @test_throws "indices are integers, like 3, so ℓₘᵢₙ" sYlm!(zeros(ComplexF64, 25), icalc, R, 1; ℓₘᵢₙ=1//2)
end

@testitem "sYlm_matrix half-integer rows equal sYlm" begin
    import SphericalFunctions: sYlm, sYlm_matrix, Ysize
    using Quaternionic: Rotor, from_spherical_coordinates

    ℓₘₐₓ = 7//2
    Rs = [from_spherical_coordinates(θ, ϕ) for θ ∈ (0.0, 0.7, 2.2, π) for ϕ ∈ (0.0, 1.2, 4.0)]
    for s ∈ (-3//2, -1//2, 1//2, 3//2), ℓₘᵢₙ ∈ (abs(s), 1//2)
        M = sYlm_matrix(Rs, ℓₘₐₓ, s; ℓₘᵢₙ)
        @test M isa Matrix{ComplexF64}
        @test size(M) == (length(Rs), Ysize(ℓₘᵢₙ, ℓₘₐₓ))
        for (i, R) ∈ enumerate(Rs)
            @test M[i, :] == strided(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ))
        end
    end
    @test_throws "exceeds ℓₘₐₓ" sYlm_matrix(Rs, 1//2, 3//2)
end

@testitem "sYlm half-integer spellings and mixed kinds" begin
    import SphericalFunctions: sYlm, sYlm!, sYlm_matrix, sYlmCalculator, Ysize, HalfOddInteger
    using Quaternionic: Rotor, from_spherical_coordinates

    R = from_spherical_coordinates(0.7, 1.2)
    Rs = [R, -R]
    ℓₘₐₓ, s, ℓₘᵢₙ = HalfOddInteger(7//2), HalfOddInteger(1//2), HalfOddInteger(1//2)
    # The `Rational` spelling is normalized at the boundary and gives exactly what the
    # `HalfOddInteger` spelling gives, keyword included.
    @test strided(sYlm(R, 7//2, 1//2)) == strided(sYlm(R, ℓₘₐₓ, s))
    @test strided(sYlm(R, 7//2, 3//2; ℓₘᵢₙ=1//2)) == strided(sYlm(R, ℓₘₐₓ, HalfOddInteger(3//2); ℓₘᵢₙ))
    @test strided(sYlm(R, 7//2, HalfOddInteger(1//2))) == strided(sYlm(R, ℓₘₐₓ, s))
    @test sYlm_matrix(Rs, 7//2, 1//2) == sYlm_matrix(Rs, ℓₘₐₓ, s)
    @test sYlm_matrix(Rs, 7//2, 1//2; ℓₘᵢₙ=1//2) == sYlm_matrix(Rs, ℓₘₐₓ, s; ℓₘᵢₙ)
    Yr = sYlm!(Vector{ComplexF64}(undef, Ysize(1//2, 7//2)), R, 7//2, 1//2)
    Yh = sYlm!(Vector{ComplexF64}(undef, Ysize(1//2, 7//2)), R, ℓₘₐₓ, s)
    @test Yr == Yh
    calc = sYlmCalculator(R, 7//2, 1//2)
    @test sYlm!(similar(Yr), calc, R, 1//2; ℓₘᵢₙ=1//2) == sYlm!(similar(Yr), calc, R, s; ℓₘᵢₙ)
    # A mixture of the two kinds of index is refused with a message naming both spellings.
    msg = "must all be integers, like 3, or all be half-odd-integers, like 7//2"
    @test_throws msg strided(sYlm(R, 7//2, 1))
    @test_throws msg strided(sYlm(R, 4, 1//2))
    @test_throws msg strided(sYlm(R, 7//2, 1//2; ℓₘᵢₙ=0))
    @test_throws msg strided(sYlm(R, 4, 1; ℓₘᵢₙ=1//2))  # the keyword alone of the other kind
    @test_throws msg strided(sYlm(R, 4, HalfOddInteger(1//2)))
    @test_throws msg sYlm_matrix(Rs, 7//2, 1)
    @test_throws msg sYlm_matrix(Rs, 4, 1//2)
    @test_throws msg sYlm!(similar(Yr), R, 7//2, 1)
    @test_throws msg sYlm!(similar(Yr), R, 4, 1//2)
    # A `Rational` that is not a half-odd-integer is refused as such.
    @test_throws "must have denominator 2" strided(sYlm(R, 7//3, 1//3))
    @test_throws "must have denominator 2" strided(sYlm(R, 4//1, 1//1))
    # Integer indices of differing concrete types are unified, and give the `Int` result.
    @test strided(sYlm(R, 4, Int8(1))) == strided(sYlm(R, 4, 1))
    @test strided(sYlm(R, Int8(4), 1; ℓₘᵢₙ=Int16(0))) == strided(sYlm(R, 4, 1; ℓₘᵢₙ=0))
    @test sYlm_matrix(Rs, 4, Int8(1)) == sYlm_matrix(Rs, 4, 1)
    # Indices all of one narrower integer type are kept as that type, all the way into the
    # calculator, and give the `Int` result exactly.
    for IT in (Int8, Int16, Int32)
        @test strided(sYlm(R, IT(4), IT(1))) == strided(sYlm(R, 4, 1))
        @test strided(sYlm(R, IT(4), IT(-1); ℓₘᵢₙ=IT(2))) == strided(sYlm(R, 4, -1; ℓₘᵢₙ=2))
        @test sYlm!(Vector{ComplexF64}(undef, Ysize(1, 4)), R, IT(4), IT(1)) == strided(sYlm(R, 4, 1))
        @test sYlm_matrix(Rs, IT(4), IT(1)) == sYlm_matrix(Rs, 4, 1)
        narrowcalc = sYlmCalculator(R, IT(4), IT(1))
        @test narrowcalc isa sYlmCalculator{IT}
        @test narrowcalc.ℓ isa Base.RefValue{IT}
        @test sYlm!(Vector{ComplexF64}(undef, Ysize(1, 4)), narrowcalc, R, IT(1)) ==
            strided(sYlm(R, 4, 1))
    end
    # A half-odd-integer spelled as a `Rational` of another integer type is the same index.
    @test strided(sYlm(R, big(7)//2, big(1)//2)) == strided(sYlm(R, ℓₘₐₓ, s))
    @test strided(sYlm(R, Int8(7)//Int8(2), Int8(1)//Int8(2))) == strided(sYlm(R, ℓₘₐₓ, s))
    @test sYlm_matrix(Rs, Int8(7)//Int8(2), Int8(1)//Int8(2)) == sYlm_matrix(Rs, ℓₘₐₓ, s)
end

@testitem "sYlm half-integer orthonormality on the sphere" begin
    import SphericalFunctions: sYlm_matrix, clenshaw_curtis_rings, clenshaw_curtis
    using Quaternionic: Rotor, from_spherical_coordinates
    using LinearAlgebra: Diagonal, I

    # ∫ ₛYₗₘ conj(ₛYₗ′ₘ′) sinθ dθ dϕ = δ δ, by quadrature: Clenshaw–Curtis in θ with
    # N = 2ℓₘₐₓ+1 rings (whose weights include the sinθ), and Nϕ = 2ℓₘₐₓ+1 equally spaced ϕ,
    # both of which are exact at this band limit.  The sample points are the rotors
    # `from_spherical_coordinates(θ, ϕ)`, with ϕ running once around [0, 2π).  Measured worst
    # case over this grid: 1.6e-15 at ℓₘₐₓ = 9/2, so 1e-14 leaves a factor of ≳ 6.
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

    # A circuit in ϕ returns to the antipodal rotor, so for half-integer s the harmonics change
    # sign: ₛY(θ, ϕ+2π) = -ₛY(θ, ϕ).  For integer s they do not.  Both are only reproduced to
    # rounding, because the phases are recomputed from a rotor whose half-angle differs by
    # rounding.  Measured worst case at these angles: 7.8e-16 for the half-integer kind and
    # 8.4e-16 for the integer kind, so 1e-14 leaves a factor of ≳ 12.
    for (θ, ϕ) ∈ ((0.7, 1.2), (2.2, 4.0), (1.0, 0.0), (0.3, 5.9))
        R, R′ = from_spherical_coordinates(θ, ϕ), from_spherical_coordinates(θ, ϕ + 2π)
        for s ∈ (-3//2, -1//2, 1//2, 3//2)
            Y, Y′ = strided(sYlm(R, 9//2, s)), strided(sYlm(R′, 9//2, s))
            @test maximum(abs, Y′ + Y) < 1e-14
        end
        for s ∈ (-1, 0, 2)
            Y, Y′ = strided(sYlm(R, 4, s)), strided(sYlm(R′, 4, s))
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
        Y32 = strided(sYlm(R32, 9//2, s))
        Y64 = strided(sYlm(R64, 9//2, s))
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
# singly or as an ascending range.  The items below check that the two spellings are accepted
# in every form the indices take, that a range and the single spin weights composing it agree
# bit for bit — the point being that the change of interface moved no arithmetic — and that
# the flat functions lay a range out the way their docstrings say.

@testitem "sYlmCalculator spin weights, however spelled" begin
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

    # Narrow integer types survive into the calculator, as they do for a single spin weight
    for IT ∈ (Int8, Int16, Int32)
        calc = sYlmCalculator(R, IT(4), IT(-1):IT(1))
        @test calc isa sYlmCalculator{IT}
        @test spins(calc) == -1:1 && eltype(spins(calc)) === IT
    end

    # Refusals: a mixture of the two kinds of index, on either side of the colon
    msg = "must all be integers, like 3, or all be half-odd-integers, like 7//2"
    @test_throws msg sYlmCalculator(R, 4, -3//2:3//2)
    @test_throws msg sYlmCalculator(R, 9//2, -1:1)
    @test_throws "|s| ≤ ℓₘₐₓ" sYlmCalculator(R, 4, -5:5)
    @test_throws "|s| ≤ ℓₘₐₓ" sYlmCalculator(R, 9//2, -11//2:11//2)
end

@testitem "sYlmCalculator ranges agree with single spin weights" begin
    import SphericalFunctions: sYlmCalculator, sYlm, spins, recurrence!, eachℓ
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(2024)
    rotors = randn(rng, Rotor{Float64}, 3)
    θs = [0.0, 0.8, 2.6]

    # The whole point of the change is that it moved no arithmetic: a calculator built for a
    # range gives, for each spin weight in it, exactly the bits a calculator built for that
    # one spin weight gives.  Checked for both kinds of index, both shapes of rotor data, and
    # batched as well as single.
    for (ℓₘₐₓ, srange) ∈ ((5, -2:2), (5, 1:2), (9//2, -3//2:3//2), (9//2, 1//2:3//2))
        for data ∈ (rotors[1], rotors, θs[2], θs)
            ranged = sYlmCalculator(data, ℓₘₐₓ, srange)
            batched = data isa AbstractVector
            for s ∈ spins(ranged)
                singly = sYlmCalculator(data, ℓₘₐₓ, s)
                # `eachℓ(ranged, s)` yields the same kind of block, and the same numbers
                @test [ℓ => copy(b) for (ℓ, b) ∈ eachℓ(ranged, s)] ==
                      [ℓ => copy(b) for (ℓ, b) ∈ singly]
                # ... and so does the corresponding slice of the whole block
                expected = [collect(b) for (_, b) ∈ singly]
                for (k, (ℓ, b)) ∈ enumerate(ranged)
                    @test collect(batched ? b[:, s, :] : b[s, :]) == expected[k]
                end
            end
        end
    end

    # Below |s| the values are zero, in a range exactly as for a single spin weight
    calc = sYlmCalculator(rotors[1], 4, -2:2)
    recurrence!(calc, 1)
    for s ∈ (-2, 2), m ∈ -1:1
        @test iszero(calc[1][s, m])
    end
    @test !iszero(calc[1][0, 0])

    # A range calculator also reproduces the flat `sYlm`, which is the independent oracle
    for s ∈ -2:2
        Y = strided(sYlm(rotors[1], 4, s; ℓₘᵢₙ=0))
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
        Y = strided(sYlm(R, ℓₘₐₓ, srange))
        @test Y isa Matrix{ComplexF64}
        @test size(Y) == (n, nmodes)
        for (i, s) ∈ enumerate(sr)
            @test Y[i, :] == strided(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ))
        end

        # `sYlm!` fills the same thing, and returns it
        Y′ = similar(Y)
        @test sYlm!(Y′, R, ℓₘₐₓ, srange) === Y′
        @test Y′ == Y
        # ... as does the calculator form, which needs no spin weight of its own
        calc = sYlmCalculator(R, ℓₘₐₓ, srange)
        fill!(Y′, 0)
        @test sYlm!(Y′, calc, rotors[2]) === Y′
        @test Y′ == strided(sYlm(rotors[2], ℓₘₐₓ, srange))
        # ... and one spin weight of that same calculator still fills a vector
        v = Vector{ComplexF64}(undef, nmodes)
        @test sYlm!(v, calc, rotors[2], first(sr); ℓₘᵢₙ) == strided(sYlm(rotors[2], ℓₘₐₓ, first(sr); ℓₘᵢₙ))

        # `sYlm_matrix` gives a stack of synthesis matrices, indexed [rotor, spin, mode]
        M = sYlm_matrix(rotors, ℓₘₐₓ, srange)
        @test M isa Array{ComplexF64, 3}
        @test size(M) == (length(rotors), n, nmodes)
        for (i, s) ∈ enumerate(sr)
            @test M[:, i, :] == sYlm_matrix(rotors, ℓₘₐₓ, s; ℓₘᵢₙ)
        end
        for (j, Rj) ∈ enumerate(rotors)
            @test M[j, :, :] == strided(sYlm(Rj, ℓₘₐₓ, srange))
        end
    end

    # An explicit ℓₘᵢₙ is honoured, and the rows below their own |s| are zero
    Y = strided(sYlm(R, 4, -2:2; ℓₘᵢₙ=0))
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
    recurrence!(batched, 3)
    @test axes(batched[3]) == (1:3, -3:3)

    # Half-integer ℓ has no spin-weight-zero analogue
    @test_throws MethodError YlmCalculator(R, 7//2)
end
