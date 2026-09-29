# Tests of 𝔇 and the harmonics at and near the poles β = 0 and β = π, where the recurrence's
# split of the rotor into half angles and phases is singular or ill-conditioned, and where the
# calculators evaluate a rotor from the expansion in its Cayley–Klein parameters instead (see
# `src/wigner/poles.jl`).  The references are the full polynomial `D_polynomial` of the
# `ExplicitWignerMatrices` module, evaluated in `BigFloat`, which shares no code with the
# package; the recurrence itself, where it is accurate; and, at large ℓ, where the polynomial's
# coefficients overflow, the product 𝔇(R) = 𝔇(R Q⁻¹) 𝔇(Q) with Q = 1 + 𝐣, which moves the
# pole to β = π/2.
#
# The tolerances all follow one error model, relative to the size √(ℓ+1) (ℓ+1)ᵏ of a k-th
# derivative: at a distance r from a pole (r = sin(β/2) from β = 0, or cos(β/2) from β = π) the
# recurrence's k-th derivatives are wrong by about ε max(1, r⁻ᵏ), and the expansion's by about
# ε + ((ℓ+1) r)^{N+1-k}, where N = `pole_order`.  A single constant multiplies that model for
# every number type, so that the tests check how the errors scale with ε rather than fitting
# a tolerance to each type; the measured errors are below a third of it.

@testsnippet PoleTools begin
    import ForwardDiff
    using Quaternionic: Quaternion, Rotor, from_euler_angles
    import SphericalFunctions
    import SphericalFunctions: D, d, sYlm, sλlm, DCalculator, dCalculator, sYlmCalculator,
        recurrence!, array_view, Yindex, HalfOddInteger
    import SphericalFunctions: pole_element, pole_data, pole_power, pole_order, pole_radius

    # The direction in which derivatives are taken, which crosses the pole when the path passes
    # through it, and a path through a rotor R₀ in that direction, R₀ exp(t𝐮/2).  The path is
    # written with `cos` and `sin` rather than with Quaternionic's `exp`, which goes through
    # `hypot`, whose derivative at zero is 0/0, so that nested dual numbers give it NaN second
    # derivatives at t = 0 (as of Quaternionic 4.4.1).  The direction is converted to the type
    # of the path before use, so that a path in `BigFloat` follows exactly the same direction as
    # the one it is compared with.
    const 𝐮 = let u = [0.3, -0.7, 0.2]; u / sqrt(sum(abs2, u)) end
    qpath(t, u) = Quaternion(cos(t/2), sin(t/2) * u[1], sin(t/2) * u[2], sin(t/2) * u[3])
    through(R₀, ::Type{T}) where {T} = t -> R₀ * qpath(t, T.(𝐮))
    tobig(R) = Quaternion(big(R[1]), big(R[2]), big(R[3]), big(R[4]))

    # The k-th derivative at x, by nested ForwardDiff, of a function returning a real vector
    nthderiv(f, x, k) = k == 0 ? f(x) : ForwardDiff.derivative(t -> nthderiv(f, t, k - 1), x)
    flat(M) = (v = vec(M); vcat(real.(v), imag.(v)))
    maxdiff(a, b) = Float64(maximum(abs, a .- b))

    # Blocks of 𝔇 as ordinary matrices, from the package, from the package's recurrence alone,
    # from the expansion about one pole, and from the full polynomial.  The recurrence's block
    # is taken from a calculator whose record of the rotors near a pole has been cleared, so
    # that `materialize!` leaves the recurrence's values in place; that is meaningful only for
    # a rotor that is not exactly at a pole, whose engine data would be constants.
    index(n) = n isa Integer ? n : HalfOddInteger(n)
    block_of(b, n) = [b[m′, m] for m′ ∈ -n:n, m ∈ -n:n]
    Dblock(n, R::Rotor) = block_of(D(R, n)[n], n)
    Dblock(n, R) = Dblock(n, Rotor(R))
    function recurrence_block(n, R)
        c = DCalculator(Rotor(R), n)
        empty!(c.poles)
        block_of(recurrence!(c, n), n)
    end
    conjpower(z, k) = k ≥ 0 ? conj(pole_power(z, k)) : pole_power(z, -k)  # conj(zᵏ), |z| = 1
    function expansion_block(n, R, north, N=pole_order)
        RT = typeof(R[1] / one(R[1]))
        ζ, κ, z = pole_data(R, RT, north)
        [
            pole_element(
                index(n), index(m′), index(m), north, ζ, κ,
                conjpower(z, Int(north ? m′ + m : m′ - m)), N
            )
            for m′ ∈ -n:n, m ∈ -n:n
        ]
    end
    Pblock(n, R) = [ExplicitWignerMatrices.D_polynomial(n, m′, m, R)[1] for m′ ∈ -n:n, m ∈ -n:n]
    Yrow(n, s, R) = (Y = sYlm(Rotor(R), n, s); [Y[n][m] for m ∈ -n:n])
    Prow(n, s, R) = [ExplicitWignerMatrices.sYlm_polynomial(n, m, s, R)[1] for m ∈ -n:n]

    # The distance of a rotor from the nearer pole, and whether that is the north pole
    function pole_distance(R)
        a = Float64(R[1])^2 + Float64(R[4])^2
        b = Float64(R[2])^2 + Float64(R[3])^2
        b ≤ a ? (sqrt(b / (a + b)), true) : (sqrt(a / (a + b)), false)
    end

    # The error model described at the top of this file, times `c`, at a distance r from a pole.
    # `source` says whether the values being checked come from the recurrence (`:recurrence`),
    # from the expansion (`:expansion`), or from a comparison of the two (`:both`); the default,
    # `:package`, asks for whichever of the two the calculators use there.
    function tolerance(::Type{T}, n, k, r; c=10, source=:package, ℓₘₐₓ=n) where {T}
        L = Float64(n) + 1
        ε = Float64(eps(T))
        if source === :package
            source = r < pole_radius(T, ℓₘₐₓ) ? :expansion : :recurrence
        end
        bound = if source === :recurrence
            ε * max(1, r^(-k))
        elseif source === :expansion
            ε + (L * r)^(pole_order + 1 - k)
        else  # both, as when the recurrence is compared with the expansion
            ε * max(1, r^(-k)) + (L * r)^(pole_order + 1 - k)
        end
        c * √L * L^k * bound
    end
end


@testitem "Poles: the full polynomial agrees with the recurrence" setup=[ExplicitWignerMatrices, PoleTools] begin
    using DoubleFloats: Double64
    using Random
    rng = Random.Xoshiro(17)
    rotors = [randn(rng, 4) for _ ∈ 1:3]

    # Values, against the polynomial in `BigFloat`.  This checks the formula, and with it the
    # expansion that shares its form, as much as the recurrence.  The largest error measured is
    # 4 ε at ℓ = 20, in each type.
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
    # polynomial in `BigFloat` dual numbers.  Random rotors can lie close enough to a pole for
    # the recurrence's derivatives to have lost digits, which the model allows for.
    for T ∈ (Float32, Float64), q ∈ rotors[1:2], n ∈ (1, 2, 4, 6, 3//2, 9//2), k ∈ 1:3
        R = Quaternion(T.(q)...)
        r, _ = pole_distance(R)
        a = nthderiv(t -> flat(Dblock(n, through(R, T)(t))), zero(T), k)
        b = nthderiv(t -> flat(Pblock(n, through(tobig(R), T)(t))), big(zero(T)), k)
        @test maxdiff(a, b) ≤ tolerance(T, n, k, r; c=4)
    end
end


@testitem "Poles: exact at the poles" setup=[ExplicitWignerMatrices, PoleTools] begin
    north = Quaternion(from_euler_angles(0.3, 0.0, -1.1))
    south = Quaternion(0.0, 0.0, 1.0, 0.0) * Quaternion(cos(0.35), 0.0, 0.0, sin(0.35))
    @test pole_distance(north) == (0.0, true)
    @test pole_distance(south) == (0.0, false)

    # The expansion keeps the terms of degree at most `pole_order` in the parameter that
    # vanishes, and the terms it leaves out vanish at the pole together with all of their
    # derivatives of lower order.  So its derivatives there, up to that order, are exactly those
    # of the whole sum — which is to say, bit for bit, since the terms left out contribute
    # exact zeros.
    for (R, isnorth) ∈ ((north, true), (south, false)), n ∈ (1, 3, 6, 1//2, 7//2), k ∈ 1:4
        f(N) = t -> flat(expansion_block(n, through(R, Float64)(t), isnorth, N))
        @test nthderiv(f(pole_order), 0.0, k) == nthderiv(f(10^6), 0.0, k)
    end

    # The package's values and derivatives up to fourth order at the poles, for 𝔇 and for the
    # harmonics, against the full polynomial, which is smooth there.  The largest error measured
    # is half the tolerance's unit.
    for (R, isnorth) ∈ ((north, true), (south, false)), n ∈ (1, 2, 5, 1//2, 7//2), k ∈ 0:4
        a = nthderiv(t -> flat(Dblock(n, through(R, Float64)(t))), 0.0, k)
        b = nthderiv(t -> flat(Pblock(n, through(R, Float64)(t))), 0.0, k)
        @test maxdiff(a, b) ≤ tolerance(Float64, n, k, 0.0; c=4)
    end
    for (R, isnorth) ∈ ((north, true), (south, false)), (n, s) ∈ ((3, -1), (4, 2), (5//2, 1//2)), k ∈ 0:2
        a = nthderiv(t -> flat(Yrow(n, s, through(R, Float64)(t))), 0.0, k)
        b = nthderiv(t -> flat(Prow(n, s, through(R, Float64)(t))), 0.0, k)
        @test maxdiff(a, b) ≤ tolerance(Float64, n, k, 0.0; c=4)
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


@testitem "Poles: the expansion near a pole agrees with the recurrence" setup=[ExplicitWignerMatrices, PoleTools] begin
    using DoubleFloats: Double64
    # Rotors within ⁴√ε of either pole, where the calculators use the expansion.  There its
    # values and derivatives should agree with the recurrence's to within the recurrence's own
    # error, which grows as the pole is approached, like ε r⁻ᵏ for the k-th derivative — at
    # ⁴√ε that is about ε^{1-k/4}, and closer in it is larger still.  Wherever that error is
    # not small compared with the derivative itself, the comparison could not fail, so it is
    # made only where the tolerance is below a hundredth of the largest element being
    # compared; the number of comparisons made for each k is checked at the end.  Against the
    # full polynomial the expansion should be accurate to about ε + ((ℓ+1) r)^{N+1-k}, and
    # that is checked everywhere, up to the fourth derivative.
    #
    # The largest ratios of the measured differences to the model below, for Float32, Float64
    # and Double64 alike, are 2.3 (for the values) and 0.3 (for the derivatives) against the
    # recurrence, and 3.2 against the polynomial.
    cases(::Type{Double64}) = ((1, 4, 5//2), (1, 0.01))
    cases(::Type) = ((1, 2, 4, 8, 3//2, 7//2), (1, 0.1, 0.01))
    compared = zeros(Int, 4)  # comparisons with the recurrence, for k = 0:3
    for T ∈ (Float32, Float64, Double64)
        ns, fractions = cases(T)
        δ₀ = sqrt(sqrt(eps(T)))
        for f ∈ fractions, isnorth ∈ (true, false), (α, γ) ∈ ((0.3, -1.1), (2.0, 0.4)), n ∈ ns
            δ = T(f) * δ₀
            R = Quaternion(from_euler_angles(T(α), isnorth ? δ : T(π) - δ, T(γ)))
            r, n̂ = pole_distance(R)
            @test n̂ == isnorth
            @test r < pole_radius(T, n)  # so that the package does use the expansion here
            for k ∈ 0:4
                e = nthderiv(t -> flat(expansion_block(n, through(R, T)(t), isnorth)), zero(T), k)
                if k ≤ 3 && tolerance(T, n, k, r; source=:both) ≤ maximum(abs, e) / 100
                    h = nthderiv(t -> flat(recurrence_block(n, through(R, T)(t))), zero(T), k)
                    @test maxdiff(e, h) ≤ tolerance(T, n, k, r; source=:both)
                    compared[k+1] += 1
                end
                p = nthderiv(t -> flat(Pblock(n, through(tobig(R), T)(t))), big(zero(T)), k)
                @test maxdiff(e, p) ≤ tolerance(T, n, k, r; source=:expansion)
            end
            # And the package itself, which uses the expansion here, with its phases from the
            # power tables rather than from `pole_power`
            @test maxdiff(Dblock(n, R), expansion_block(n, R, isnorth)) ≤ tolerance(T, n, 0, r; source=:expansion)
        end
    end
    @test all(>(0), compared)
end


@testitem "Poles: the radius of the expansion" setup=[ExplicitWignerMatrices, PoleTools] begin
    # A rotor is evaluated from the expansion exactly when it lies within `pole_radius` of a
    # pole, and the switch is smooth to the accuracy of the side that is worse there: the values
    # agree to a few ulps on either side of it, and the derivatives to the recurrence's
    # accuracy just outside.
    for ℓₘₐₓ ∈ (4, 8, 32), isnorth ∈ (true, false), side ∈ (0.999, 1.001)
        r = side * pole_radius(Float64, ℓₘₐₓ)
        β = isnorth ? 2asin(r) : π - 2asin(r)
        R = Quaternion(from_euler_angles(0.3, β, -1.1))
        c = DCalculator(Rotor(R), ℓₘₐₓ)
        @test length(c.poles) == (side < 1)
        for k ∈ 0:2
            a = nthderiv(t -> flat(Dblock(ℓₘₐₓ, through(R, Float64)(t))), 0.0, k)
            b = nthderiv(t -> flat(Pblock(ℓₘₐₓ, through(tobig(R), Float64)(t))), big(0.0), k)
            @test maxdiff(a, b) ≤ tolerance(Float64, ℓₘₐₓ, k, r)
        end
    end
    # The radius is set by ℓₘₐₓ, so a rotor between the radii of two calculators is evaluated
    # differently by them, but to the same accuracy.
    r = sqrt(pole_radius(Float64, 8) * pole_radius(Float64, 32))
    R = Rotor(from_euler_angles(0.3, 2asin(r), -1.1))
    𝔇₈ = D(R, 8)
    𝔇₃₂ = D(R, 32)
    @test all(maximum(abs, array_view(𝔇₈[n]) - array_view(𝔇₃₂[n])) ≤ 8 * √(n + 1) * eps() for n ∈ 0:8)
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
            @test maxdiff(a, b) ≤ tolerance(Float64, n, k, 0.0; c=200)
            a = nthderiv(t -> flat(column(through(R₀, Float64)(t), -2)), 0.0, k)
            b = nthderiv(t -> flat(Δ * column(Q⁻¹ * through(R₀, Float64)(t), -2)), 0.0, k)
            @test maxdiff(a, b) ≤ tolerance(Float64, n, k, 0.0; c=200)
        end
    end
end


@testitem "Poles: calculators keep track of the rotors near a pole" setup=[PoleTools] begin
    import SphericalFunctions: set_R!, set_θ!, sYlm_matrix
    using Quaternionic: from_spherical_coordinates
    using Random
    rng = Random.Xoshiro(41)
    generic = [Rotor(randn(rng, 4)...) for _ ∈ 1:2]
    northR = from_euler_angles(0.3, 0.0, -1.1)
    southR = Rotor(Quaternion(0.0, 0.3, 0.8, 0.0))
    nearR = from_euler_angles(1.2, 1e-7, 0.4)
    batch = [generic[1], northR, southR, nearR, generic[2]]

    # In a batch, the values of every rotor near a pole are exactly those of a calculator built
    # for it alone, and those of the others agree to rounding; and the blocks restricted by the
    # keywords are exactly the corresponding parts of the whole blocks.
    same(i, a, b) = i ∈ 2:4 ? a == b : maxdiff(a, b) ≤ 8eps()
    for n ∈ (5, 7//2)
        c = DCalculator(batch, n)
        @test sort([p.iᵣ for p ∈ c.poles]) == [2, 3, 4]
        blk = recurrence!(c, n)
        for (i, R) ∈ enumerate(batch)
            @test same(i, [blk[i, m′, m] for m′ ∈ -n:n, m ∈ -n:n], Dblock(n, R))
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
            @test same(i, [blk[i, s′, m] for m ∈ -n:n], Yrow(n, s′, R))
        end
        # The flat interfaces agree with the calculator
        M = sYlm_matrix(batch, n, first(s))
        Y = [sYlm(R, n, first(s)) for R ∈ batch]
        @test all(same(i, M[i, :], array_view(Y[i])) for i ∈ eachindex(batch))
    end

    # Moving a calculator onto a pole and off it again leaves no trace: the record is replaced
    # on every change of the rotors, so the values afterwards are exactly those of a fresh
    # calculator.
    c = DCalculator(generic[1], 6)
    @test isempty(c.poles)
    set_R!(c, northR)
    @test length(c.poles) == 1
    @test block_of(recurrence!(c, 6), 6) == Dblock(6, northR)
    set_R!(c, generic[2])
    @test isempty(c.poles)
    @test block_of(recurrence!(c, 6), 6) == Dblock(6, generic[2])

    # An sYlmCalculator given angles after rotors forgets the rotors near a pole, and gives
    # exactly what a calculator built from those angles gives.
    c = sYlmCalculator([northR, southR], 4, -1)
    @test length(c.poles) == 2
    set_θ!(c, [0.0, 0.7])
    @test isempty(c.poles)
    @test array_view(recurrence!(c, 4)) == array_view(recurrence!(sYlmCalculator([0.0, 0.7], 4, -1), 4))

    # `similar` copies the record, and `fill!` leaves it alone, like the rest of the rotor data.
    c = DCalculator([northR, generic[1], southR], 4)
    c′ = similar(c)
    @test [(p.iᵣ, p.north) for p ∈ c′.poles] == [(p.iᵣ, p.north) for p ∈ c.poles]
    @test array_view(recurrence!(c′, 4)) == array_view(recurrence!(c, 4))
    fill!(c, NaN)
    @test array_view(recurrence!(c, 4)) == array_view(recurrence!(c′, 4))
    @test !any(isnan, array_view(recurrence!(c, 4)))

    # A rotor from spherical coordinates at the north pole is exactly at it.
    @test length(DCalculator(from_spherical_coordinates(0.0, 1.3), 3).poles) == 1
end


@testitem "Poles: plain floating-point values" setup=[ExplicitWignerMatrices, PoleTools] begin
    using DoubleFloats: Double64
    import MathChecker: checked, unchecked
    import SphericalFunctions: set_R!
    # At a rotor exactly at a pole, in every floating-point type, the values are within
    # rounding of the exact ones, and every element that vanishes there is exactly zero.
    for T ∈ (Float32, Float64, Double64, BigFloat), n ∈ (3, 12, 7//2)
        α = T(7) / 10
        for (R, isnorth) ∈ (
            (Quaternion(cos(α/2), zero(T), zero(T), sin(α/2)), true),
            (Quaternion(zero(T), cos(α/2), sin(α/2), zero(T)), false),
        )
            𝔇 = Dblock(n, R)
            𝔇ᵖ = Pblock(n, tobig(R))
            @test Float64(maximum(abs, Complex{BigFloat}.(𝔇) .- 𝔇ᵖ)) ≤ 4 * √(n + 1) * eps(T)
            @test all(
                iszero(𝔇[j′, j]) for (j′, m′) ∈ enumerate(-n:n), (j, m) ∈ enumerate(-n:n)
                if (isnorth ? m′ - m : m′ + m) != 0
            )
        end
    end

    # A rotor whose X² + Y² underflows to zero although X does not is taken to be at the pole,
    # and is evaluated accurately.
    R = Quaternion(cos(0.35), 1e-170, 0.0, sin(0.35))
    @test R[2]^2 + R[3]^2 == 0
    for n ∈ (4, 9//2)
        @test Float64(maximum(abs, Complex{BigFloat}.(Dblock(n, R)) .- Pblock(n, tobig(R)))) ≤ 4 * √(n + 1) * eps()
    end

    # No NaN arises anywhere, not even in the engine's work for a rotor at a pole, whose
    # results are overwritten: `MathChecker` throws as soon as one takes part in an operation.
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
    # Reverse mode runs through every operation it recorded, including the engine's work for a
    # rotor at a pole, whose results are overwritten; that work must therefore take no square
    # root of zero.  The gradients agree with ForwardDiff's at the poles, near one and away
    # from both, for 𝔇 and the harmonics, with integer and half-integer indices.
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
    # compared with those of 𝔇 at rotors along the corresponding path, which the expansion
    # gives exactly.
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
