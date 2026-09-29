# Tests of 𝔇 and the harmonics at and near the poles β = 0 and β = π, where the
# recurrence's split of the rotor into half angles and phases is singular or
# ill-conditioned.  The values are accurate there, and the derivatives come from the rules
# for automatic differentiation, which form them from the values (see `src/derivatives.jl`),
# so they are accurate there too.  The references are the full polynomial `D_polynomial` of
# the `ExplicitWignerMatrices` module, evaluated in `BigFloat`, which shares no code with
# the package; and, at large ℓ, where the polynomial's coefficients overflow, the product
# 𝔇(R) = 𝔇(R Q⁻¹) 𝔇(Q) with Q = 1 + 𝐣, which moves the pole to β = π/2.
#
# The tolerances all follow one error model: ε times the size √(ℓ+1) (ℓ+1)ᵏ of a k-th
# derivative, at any distance from a pole.  A single constant multiplies that model for
# every number type, so that the tests check how the errors scale with ε rather than fitting
# a tolerance to each type.

@testsnippet PoleTools begin
    import ForwardDiff
    using Quaternionic: Quaternion, Rotor, from_euler_angles
    import SphericalFunctions
    import SphericalFunctions: D, d, sYlm, sλlm, DCalculator, dCalculator, sYlmCalculator,
        recurrence!, array_view, Yindex

    # The direction in which derivatives are taken, which crosses the pole when the path passes
    # through it, and a path through a rotor R₀ in that direction, R₀ exp(t𝐮/2).  The path is
    # written with `cos` and `sin` rather than with Quaternionic's `exp`, which goes through
    # `hypot`, whose derivative at zero is 0/0, so that nested dual numbers give it NaN second
    # derivatives at t = 0 (as of Quaternionic 4.4.1).  The direction is converted to the type
    # of the path before use, so that a path in `BigFloat` follows exactly the same direction as
    # the one it is compared with.
    # A tuple rather than a vector: TestItemRunner on Julia 1.10 evaluates this snippet twice
    # in each item's module, and redefining a constant with an identical value, unlike an
    # equal vector, draws no warning
    const 𝐮 = let u = (0.3, -0.7, 0.2); u ./ sqrt(sum(abs2, u)) end
    qpath(t, u) = Quaternion(cos(t/2), sin(t/2) * u[1], sin(t/2) * u[2], sin(t/2) * u[3])
    through(R₀, ::Type{T}) where {T} = t -> R₀ * qpath(t, T.(𝐮))
    tobig(R) = Quaternion(big(R[1]), big(R[2]), big(R[3]), big(R[4]))

    # The k-th derivative at x, by nested ForwardDiff, of a function returning a real vector
    nthderiv(f, x, k) = k == 0 ? f(x) : ForwardDiff.derivative(t -> nthderiv(f, t, k - 1), x)
    flat(M) = (v = vec(M); vcat(real.(v), imag.(v)))
    maxdiff(a, b) = Float64(maximum(abs, a .- b))

    # Blocks of 𝔇 and rows of the harmonics as ordinary matrices, from the package and from
    # the full polynomial
    block_of(b, n) = [b[m′, m] for m′ ∈ -n:n, m ∈ -n:n]
    Dblock(n, R::Rotor) = block_of(D(R, n)[n], n)
    Dblock(n, R) = Dblock(n, Rotor(R))
    Pblock(n, R) = [ExplicitWignerMatrices.D_polynomial(n, m′, m, R)[1] for m′ ∈ -n:n, m ∈ -n:n]
    Yrow(n, s, R) = (Y = sYlm(Rotor(R), n, s); [Y[n][m] for m ∈ -n:n])
    Prow(n, s, R) = [ExplicitWignerMatrices.sYlm_polynomial(n, m, s, R)[1] for m ∈ -n:n]

    # The distance of a rotor from the nearer pole, and whether that is the north pole
    function pole_distance(R)
        a = Float64(R[1])^2 + Float64(R[4])^2
        b = Float64(R[2])^2 + Float64(R[3])^2
        b ≤ a ? (sqrt(b / (a + b)), true) : (sqrt(a / (a + b)), false)
    end

    # The error model described at the top of this file, times `c`
    tolerance(::Type{T}, n, k; c=10) where {T} = c * √(Float64(n) + 1) * (Float64(n) + 1)^k * Float64(eps(T))
end


@testitem "Poles: the full polynomial agrees with the recurrence" setup=[ExplicitWignerMatrices, PoleTools] begin
    using DoubleFloats: Double64
    using Random
    rng = Random.Xoshiro(17)
    rotors = [randn(rng, 4) for _ ∈ 1:3]

    # Values, against the polynomial in `BigFloat`.  This checks the formula as much as the
    # recurrence.  The largest error measured is 4 ε at ℓ = 20, in each type.
    for T ∈ (Float32, Float64, Double64), q ∈ rotors
        R = Rotor{T}(T.(q)...)
        for n ∈ (0, 1, 2, 4, 8, 12, 20, 1//2, 5//2, 17//2, 33//2)
            𝔇 = Dblock(n, R)
            𝔇ᵖ = Pblock(n, tobig(R))
            @test Float64(maximum(abs, Complex{BigFloat}.(𝔇) .- 𝔇ᵖ)) ≤ 4 * √(n + 1) * eps(T)
        end
        for (n, s) ∈ ((4, -2), (6, 0), (9, 1), (12, 2), (7//2, 3//2), (13//2, -1//2))
            Y = Yrow(n, s, R)
            Yᵖ = Prow(n, s, tobig(R))
            @test Float64(maximum(abs, Complex{BigFloat}.(Y) .- Yᵖ)) ≤ 4 * √(n + 1) * eps(T)
        end
    end

    # The polynomial in the same type, against its own estimate of its rounding error, which is
    # ε times the sum of the magnitudes of its terms.  The largest ratio measured is 2.2.
    for T ∈ (Float32, Float64, Double64), q ∈ rotors, n ∈ (0, 1, 3, 6, 9, 12, 3//2, 13//2)
        R = Rotor{T}(T.(q)...)
        𝔇 = Dblock(n, R)
        ok = true
        for (j′, m′) ∈ enumerate(-n:n), (j, m) ∈ enumerate(-n:n)
            v, scale = ExplicitWignerMatrices.D_polynomial(n, m′, m, R)
            ok &= abs(𝔇[j′, j] - v) ≤ 8 * eps(T) * (scale + √(n + 1))
        end
        @test ok
    end

    # Derivatives of the first three orders, along a path through each rotor, against the
    # polynomial in `BigFloat` dual numbers.
    for T ∈ (Float32, Float64), q ∈ rotors[1:2], n ∈ (1, 2, 4, 6, 3//2, 9//2), k ∈ 1:3
        R = Quaternion(T.(q)...)
        a = nthderiv(t -> flat(Dblock(n, through(R, T)(t))), zero(T), k)
        b = nthderiv(t -> flat(Pblock(n, through(tobig(R), T)(t))), big(zero(T)), k)
        @test maxdiff(a, b) ≤ tolerance(T, n, k; c=4)
    end
end


@testitem "Poles: exact at the poles" setup=[ExplicitWignerMatrices, PoleTools] begin
    north = Quaternion(from_euler_angles(0.3, 0.0, -1.1))
    south = Quaternion(0.0, 0.0, 1.0, 0.0) * Quaternion(cos(0.35), 0.0, 0.0, sin(0.35))
    @test pole_distance(north) == (0.0, true)
    @test pole_distance(south) == (0.0, false)

    # The package's values and derivatives up to fourth order at the poles, for 𝔇 and for the
    # harmonics, against the full polynomial, which is smooth there.  The largest error measured
    # is half the tolerance's unit.
    for R ∈ (north, south), n ∈ (1, 2, 5, 1//2, 7//2), k ∈ 0:4
        a = nthderiv(t -> flat(Dblock(n, through(R, Float64)(t))), 0.0, k)
        b = nthderiv(t -> flat(Pblock(n, through(R, Float64)(t))), 0.0, k)
        @test maxdiff(a, b) ≤ tolerance(Float64, n, k; c=4)
    end
    for R ∈ (north, south), (n, s) ∈ ((3, -1), (4, 2), (5//2, 1//2)), k ∈ 0:2
        a = nthderiv(t -> flat(Yrow(n, s, through(R, Float64)(t))), 0.0, k)
        b = nthderiv(t -> flat(Prow(n, s, through(R, Float64)(t))), 0.0, k)
        @test maxdiff(a, b) ≤ tolerance(Float64, n, k; c=4)
    end

    # The example of issue #67: the derivative with respect to α of an element of 𝔇 at a
    # rotation about the z axis.
    @test ForwardDiff.derivative(α -> imag(D(from_euler_angles(α, 0.0, 0.0), 2)[2][1, 1]), 0.3) ≈
        -cos(0.3) atol=4eps()

    # Gradients and Hessians with respect to all four components, at the identity and at a
    # rotor exactly at the south pole, for integer and half-integer indices.
    for q ∈ ([1.0, 0, 0, 0], [0.0, 0, 1, 0]), (n, m′, m) ∈ ((2, 1, 0), (3, -1, 1), (3//2, 1//2, -1//2))
        F(v) = real(D(Rotor(v...), n)[n][m′, m])
        P(v) = real(ExplicitWignerMatrices.D_polynomial(n, m′, m, Quaternion(v...))[1])
        @test maximum(abs, ForwardDiff.gradient(F, q) - ForwardDiff.gradient(P, q)) ≤ 8eps() * (n + 1)^1.5
        @test maximum(abs, ForwardDiff.hessian(F, q) - ForwardDiff.hessian(P, q)) ≤ 8eps() * (n + 1)^2.5
    end

    # The generators: the derivative of 𝔇(exp(t𝐮/2)) at t = 0 is -i times the matrix of 𝐮·𝐉,
    # with ⟨m|J_z|m⟩ = m and ⟨m±1|J_±|m⟩ = √((ℓ∓m)(ℓ±m+1)), for any unit vector 𝐮.  This shares
    # no code with the package or with the polynomial.
    function J(n, u)
        N = length(-n:n)
        M = zeros(ComplexF64, N, N)
        for (j, m) ∈ enumerate(-n:n)
            M[j, j] += u[3] * m
            m < n && (M[j+1, j] += √((n - m) * (n + m + 1)) * (u[1] - im * u[2]) / 2)
            m > -n && (M[j-1, j] += √((n + m) * (n - m + 1)) * (u[1] + im * u[2]) / 2)
        end
        M
    end
    for n ∈ (1, 2, 5, 1//2, 3//2, 7//2), u ∈ ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0], 𝐮)
        Ḋ = ForwardDiff.derivative(t -> Dblock(n, qpath(t, u)), 0.0)
        @test maximum(abs, Ḋ .+ im .* J(n, u)) ≤ 4eps() * (n + 1)
    end
end


@testitem "Poles: near a pole" setup=[ExplicitWignerMatrices, PoleTools] begin
    using DoubleFloats: Double64
    # Rotors at distances r from 10⁻² down to 10⁻⁸ of either pole, where the recurrence's
    # own k-th derivatives would be wrong by about ε r⁻ᵏ relative to their size.  The
    # package's values and derivatives up to third order agree with the polynomial to within
    # the model, which does not depend on r.
    cases(::Type{Double64}) = ((1, 4, 5//2), (1e-2, 1e-6))
    cases(::Type) = ((1, 2, 4, 8, 3//2, 7//2), (1e-2, 1e-4, 1e-8))
    for T ∈ (Float32, Float64, Double64)
        ns, rs = cases(T)
        for r ∈ rs, isnorth ∈ (true, false), (α, γ) ∈ ((0.3, -1.1), (2.0, 0.4)), n ∈ ns
            δ = 2asin(T(r))
            R = Quaternion(from_euler_angles(T(α), isnorth ? δ : T(π) - δ, T(γ)))
            @test pole_distance(R)[2] == isnorth
            for k ∈ 0:3
                a = nthderiv(t -> flat(Dblock(n, through(R, T)(t))), zero(T), k)
                p = nthderiv(t -> flat(Pblock(n, through(tobig(R), T)(t))), big(zero(T)), k)
                @test maxdiff(a, p) ≤ tolerance(T, n, k)
            end
        end
    end
end


@testitem "Poles: large ℓ" setup=[PoleTools] begin
    # At large ℓ the full polynomial's coefficients overflow, so the reference is the product
    # 𝔇(R) = 𝔇(R Q⁻¹) 𝔇(Q) for Q = 1 + 𝐣, whose 𝔇(Q) is d(π/2) exactly and which moves the pole
    # to β = π/2, where the recurrence is accurate; for the harmonics, ₛY(R) = d(π/2) ₛY(Q⁻¹R).
    # Only 𝔇's norm-independent phases enter, so neither Q nor its inverse 1 - 𝐣 need be
    # normalized.  The rows are restricted to |m′| ≤ 2, which keeps the recurrences cheap.
    Q⁻¹ = Quaternion(1.0, 0.0, -1.0, 0.0)
    for (n, k) ∈ ((1000, 1), (200, 2))
        Δ = let blk = recurrence!(dCalculator(π/2, n), n); [blk[a, b] for a ∈ -n:n, b ∈ -n:n] end
        rows(R) = (blk = recurrence!(DCalculator(Rotor(R), n; m′ₘₐₓ=2), n); [blk[m′, m] for m′ ∈ -2:2, m ∈ -n:n])
        column(R, s) = (Y = sYlm(Rotor(R), n, s); [Y[n][m] for m ∈ -n:n])
        for R₀ ∈ (Quaternion(from_euler_angles(0.3, 0.0, -1.1)), Quaternion(0.0, 0.3, 0.8, 0.0))
            @test pole_distance(R₀)[1] == 0
            a = nthderiv(t -> flat(rows(through(R₀, Float64)(t))), 0.0, k)
            b = nthderiv(t -> flat(rows(through(R₀, Float64)(t) * Q⁻¹) * Δ), 0.0, k)
            @test all(isfinite, a)
            @test maxdiff(a, b) ≤ tolerance(Float64, n, k; c=200)
            a = nthderiv(t -> flat(column(through(R₀, Float64)(t), -2)), 0.0, k)
            b = nthderiv(t -> flat(Δ * column(Q⁻¹ * through(R₀, Float64)(t), -2)), 0.0, k)
            @test maxdiff(a, b) ≤ tolerance(Float64, n, k; c=200)
        end
    end
end


@testitem "Poles: calculators at and near a pole, singly and in batches" setup=[PoleTools] begin
    import SphericalFunctions: set_R!, set_θ!, sYlm_matrix
    using Random
    rng = Random.Xoshiro(41)
    generic = [Rotor(randn(rng, 4)...) for _ ∈ 1:2]
    northR = from_euler_angles(0.3, 0.0, -1.1)
    southR = Rotor(Quaternion(0.0, 0.3, 0.8, 0.0))
    nearR = from_euler_angles(1.2, 1e-7, 0.4)
    batch = [generic[1], northR, southR, nearR, generic[2]]

    # In a batch, the values of every rotor agree to rounding with those of a calculator built
    # for it alone, and the blocks restricted by the keywords are exactly the corresponding
    # parts of the whole blocks.
    same(a, b) = maxdiff(a, b) ≤ 8eps()
    for n ∈ (5, 7//2)
        c = DCalculator(batch, n)
        blk = recurrence!(c, n)
        for (i, R) ∈ enumerate(batch)
            @test same([blk[i, m′, m] for m′ ∈ -n:n, m ∈ -n:n], Dblock(n, R))
        end
        for R ∈ (northR, southR, nearR)
            lo = n isa Integer ? 1 : 1//2
            part = recurrence!(DCalculator(R, n; m′ₘₐₓ=lo+1, m′ₘᵢₙ=-lo, mₘᵢₙ=-lo-1), n)
            whole = recurrence!(DCalculator(R, n), n)
            @test all(part[m′, m] == whole[m′, m] for m′ ∈ -lo:lo+1, m ∈ -lo-1:n)
        end
    end
    for (n, s) ∈ ((5, -2:2), (7//2, -3//2:1//2))
        c = sYlmCalculator(batch, n, s)
        blk = recurrence!(c, n)
        for (i, R) ∈ enumerate(batch), s′ ∈ s
            @test same([blk[i, s′, m] for m ∈ -n:n], Yrow(n, s′, R))
        end
        # The flat interfaces agree with the calculator
        M = sYlm_matrix(batch, n, first(s))
        Y = [sYlm(R, n, first(s)) for R ∈ batch]
        @test all(same(M[i, :], array_view(Y[i])) for i ∈ eachindex(batch))
    end

    # Moving a calculator onto a pole and off it again leaves no trace: the values afterwards
    # are exactly those of a fresh calculator.
    c = DCalculator(generic[1], 6)
    set_R!(c, northR)
    @test block_of(recurrence!(c, 6), 6) == Dblock(6, northR)
    set_R!(c, generic[2])
    @test block_of(recurrence!(c, 6), 6) == Dblock(6, generic[2])

    # An sYlmCalculator given angles after rotors at the poles gives exactly what a calculator
    # built from those angles gives.
    c = sYlmCalculator([northR, southR], 4, -1)
    set_θ!(c, [0.0, 0.7])
    @test array_view(recurrence!(c, 4)) == array_view(recurrence!(sYlmCalculator([0.0, 0.7], 4, -1), 4))

    # `similar` copies the rotor data, and `fill!` leaves it alone.
    c = DCalculator([northR, generic[1], southR], 4)
    c′ = similar(c)
    @test array_view(recurrence!(c′, 4)) == array_view(recurrence!(c, 4))
    fill!(c, NaN)
    @test array_view(recurrence!(c, 4)) == array_view(recurrence!(c′, 4))
    @test !any(isnan, array_view(recurrence!(c, 4)))
end


@testitem "Poles: plain floating-point values" setup=[ExplicitWignerMatrices, PoleTools] begin
    using DoubleFloats: Double64
    import MathChecker: checked, unchecked
    import SphericalFunctions: set_R!
    # At a rotor exactly at a pole, in every floating-point type, the values are within
    # rounding of the exact ones, which include the elements that vanish there.
    for T ∈ (Float32, Float64, Double64, BigFloat), n ∈ (3, 12, 7//2)
        α = T(7) / 10
        for R ∈ (
            Quaternion(cos(α/2), zero(T), zero(T), sin(α/2)),
            Quaternion(zero(T), cos(α/2), sin(α/2), zero(T)),
        )
            𝔇 = Dblock(n, R)
            𝔇ᵖ = Pblock(n, tobig(R))
            @test Float64(maximum(abs, Complex{BigFloat}.(𝔇) .- 𝔇ᵖ)) ≤ 4 * √(n + 1) * eps(T)
        end
    end

    # A rotor whose X² + Y² underflows to zero although X does not is evaluated accurately.
    R = Quaternion(cos(0.35), 1e-170, 0.0, sin(0.35))
    @test R[2]^2 + R[3]^2 == 0
    for n ∈ (4, 9//2)
        @test Float64(maximum(abs, Complex{BigFloat}.(Dblock(n, R)) .- Pblock(n, tobig(R)))) ≤ 4 * √(n + 1) * eps()
    end

    # No NaN arises anywhere, not even in the engine's work for a rotor at a pole:
    # `MathChecker` throws as soon as one takes part in an operation.
    # The checked rotors are the plain ones converted (`Rotor{NC}(R)`, which, unlike
    # `Rotor(NC(R[1]), …)`, does not renormalize them), and `Checked` performs each operation
    # in `Float64`, so the checked values agree with the plain ones bit for bit.
    NC = checked(Float64; precision=false, nan=true, inf=false)
    rotors = [
        from_euler_angles(0.3, 0.0, -1.1), Rotor(Quaternion(0.0, 0.3, 0.8, 0.0)),
        from_euler_angles(1.2, 1e-7, 0.4),
    ]
    for n ∈ (4, 7//2)
        RNC = [Rotor{NC}(R) for R ∈ rotors]
        a = array_view(recurrence!(DCalculator(RNC, n), n))
        b = array_view(recurrence!(DCalculator(rotors, n), n))
        @test complex.(unchecked.(real.(a)), unchecked.(imag.(a))) == b
        s = n isa Integer ? -2 : 1//2
        a = array_view(recurrence!(sYlmCalculator(RNC, n, s), n))
        b = array_view(recurrence!(sYlmCalculator(rotors, n, s), n))
        @test complex.(unchecked.(real.(a)), unchecked.(imag.(a))) == b
    end

    # Once warmed up, moving a calculator onto a pole and computing there allocates nothing.
    c = DCalculator(rotors, 6)
    others = reverse(rotors)
    step!(c, R) = (set_R!(c, R); recurrence!(c, 6); nothing)
    step!(c, rotors); step!(c, others)
    @test @allocated(step!(c, others)) == 0
end


@testitem "Poles: ReverseDiff" setup=[PoleTools] begin
    import ReverseDiff
    # The gradients agree with ForwardDiff's at the poles, near one and away from both, for 𝔇
    # and the harmonics, with integer and half-integer indices.
    functions = (
        v -> real(D(Rotor(v...), 2)[2][1, 0]),
        v -> imag(D(Rotor(v...), 5//2)[5//2][1//2, -3//2]),
        v -> real(array_view(sYlm(Rotor(v...), 3, -1))[Yindex(3, 1, 1)]),
        v -> imag(array_view(sYlm(Rotor(v...), 7//2, 1//2))[Yindex(7//2, -3//2, 1//2)]),
    )
    points = ([1.0, 0, 0, 0], [0.0, 0.0, 1.0, 0.0], [0.6, 0, 0, 0.8], [1.0, 1e-7, 0, 0], [0.3, -0.5, 0.7, 0.2])
    for F ∈ functions, v ∈ points
        g = ReverseDiff.gradient(F, v)
        @test all(isfinite, g)
        @test maximum(abs, g - ForwardDiff.gradient(F, v)) ≤ 8eps()
    end
end


@testitem "Poles: angles differentiate at the poles" setup=[PoleTools] begin
    # The recurrence is smooth in the angle β, so d of an angle or of a phase, and the real
    # harmonics of an angle, need nothing special at the poles.  Their derivatives there are
    # compared with those of 𝔇 at rotors along the corresponding path, which the rules give
    # exactly.
    for β₀ ∈ (0.0, Float64(π)), n ∈ (3, 7//2)
        rotor(β) = from_euler_angles(zero(β), β, zero(β))
        ref = ForwardDiff.derivative(β -> real.(Dblock(n, rotor(β))), β₀)
        fromangle = ForwardDiff.derivative(β -> block_of(d(β, n)[n], n), β₀)
        fromphase = ForwardDiff.derivative(β -> block_of(d(cis(β), n)[n], n), β₀)
        @test maximum(abs, fromangle - ref) ≤ 8eps() * (n + 1)^1.5
        @test maximum(abs, fromphase - ref) ≤ 8eps() * (n + 1)^1.5
        s = n isa Integer ? -1 : 1//2
        # ₛλₗₘ(θ) = ₛYₗₘ(θ, 0) for integer s, and ₛYₗₘ(θ, 0) / i^{2s} for half-odd s
        λ = ForwardDiff.derivative(θ -> [sλlm(θ, n, s)[n][m] for m ∈ -n:n], β₀)
        phase = s isa Integer ? 1 : (1, im, -1, -im)[mod(Int(2s), 4) + 1]
        Y = ForwardDiff.derivative(θ -> Yrow(n, s, rotor(θ)), β₀) ./ phase
        @test maximum(abs, λ - real.(Y)) ≤ 8eps() * (n + 1)^1.5
        @test maximum(abs, imag.(Y)) ≤ 8eps() * (n + 1)^1.5
    end
end
