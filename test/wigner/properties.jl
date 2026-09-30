# Oracle-free (metamorphic) tests of the mathematical properties of the Wigner 𝔇 and d
# matrices produced by the engine.  Nothing here compares the engine with another
# implementation; every check is an identity the recurrence cannot satisfy by accident, so
# these items bear the regression burden alongside the closed-form comparisons in
# `calculators.jl` and `H_calculator.jl`.
#
# Conventions under test (see docs/src/30-conventions):
#   𝔇ˡ_{m′m}(α, β, γ) = e^{-i m′ α} dˡ_{m′m}(β) e^{-i m γ},  R = from_euler_angles(α, β, γ)
#   ₀Yₗₘ(θ, ϕ) = √((2ℓ+1)/4π) conj(𝔇ˡ_{m0}(from_spherical_coordinates(θ, ϕ)))
#
# Tolerances follow the measured Float64 accuracy of the engine: unitarity and the
# representation property at ℓ ≈ 40 are good to ~1.5e-14, so per-block bounds are
# 20·max(1, ℓ)·eps; the character sum accumulates ~ℓ·eps.

@testsnippet WignerPropertyRotors begin
    using Quaternionic: Rotor, 𝐢, 𝐣, 𝐤
    using Random: Xoshiro

    # Rotors at the poles β ∈ {0, π} (where `spinor_phases` takes its special branches and
    # d is exactly diagonal or anti-diagonal) and at β = π/2, in both signs of the double
    # cover.  Ordering: 1, 𝐢, 𝐣, 𝐤, (1+𝐯)/√2, (1-𝐯)/√2 for 𝐯 ∈ (𝐢, 𝐣, 𝐤), then the negatives.
    function axis_rotors(::Type{T}) where {T}
        invsqrt2 = inv(√T(2))
        rotors = [
            Rotor{T}(1);
            [Rotor{T}(𝐯) for 𝐯 ∈ (𝐢, 𝐣, 𝐤)];
            [Rotor{T}(invsqrt2 + invsqrt2 * 𝐯) for 𝐯 ∈ (𝐢, 𝐣, 𝐤)];
            [Rotor{T}(invsqrt2 - invsqrt2 * 𝐯) for 𝐯 ∈ (𝐢, 𝐣, 𝐤)]
        ]
        [rotors; -rotors]
    end

    # The axis rotors followed by `n` seeded random rotors.
    function property_rotors(::Type{T}, n, seed) where {T}
        [axis_rotors(T); randn(Xoshiro(seed), Rotor{T}, n)]
    end
end


@testitem "Wigner D unitarity" setup=[WignerPropertyRotors] begin
    import SphericalFunctions: D
    using Quaternionic: Rotor, 𝐤
    using LinearAlgebra: I, opnorm
    using Random: Xoshiro

    T = Float64

    function test_unitarity(R, ℓₘₐₓ)
        𝔇 = D(R, ℓₘₐₓ)
        @test axes(𝔇, 1) == 0:ℓₘₐₓ
        for ℓ ∈ 0:ℓₘₐₓ
            @test axes(𝔇[ℓ]) == (-ℓ:ℓ, -ℓ:ℓ)
            M = parent(𝔇[ℓ])
            @test M isa Matrix{Complex{T}}
            @test size(M) == (2ℓ + 1, 2ℓ + 1)
            atol = 20 * max(1, ℓ) * eps(T)
            @test opnorm(M * M' - I) ≤ atol
            @test opnorm(M' * M - I) ≤ atol
        end
    end

    # Small and moderate ℓₘₐₓ: every axis rotor plus several random ones
    for (ℓₘₐₓ, nrandom) ∈ ((4, 6), (20, 4))
        for R ∈ property_rotors(T, nrandom, 1000 + ℓₘₐₓ)
            test_unitarity(R, ℓₘₐₓ)
        end
    end

    # Larger ℓₘₐₓ: a pole rotor and a few random ones (the SVDs dominate the runtime)
    let ℓₘₐₓ = 60
        for R ∈ [Rotor{T}(𝐤); randn(Xoshiro(1000 + ℓₘₐₓ), Rotor{T}, 3)]
            test_unitarity(R, ℓₘₐₓ)
        end
    end
end


@testitem "Wigner D representation property" setup=[WignerPropertyRotors] begin
    import SphericalFunctions: D
    using Quaternionic: Rotor, 𝐢, 𝐣, 𝐤
    using LinearAlgebra: opnorm
    using Random: Xoshiro

    T = Float64

    # 𝔇ˡ(R₁ R₂) = 𝔇ˡ(R₁) 𝔇ˡ(R₂),  𝔇ˡ(R⁻¹) = 𝔇ˡ(R)†,  𝔇⁰ = [1]
    function test_representation(rotors, ℓₘₐₓ)
        𝔇s = [D(R, ℓₘₐₓ) for R ∈ rotors]
        for (i₁, R₁) ∈ enumerate(rotors), (i₂, R₂) ∈ enumerate(rotors)
            𝔇₁₂ = D(R₁ * R₂, ℓₘₐₓ)
            for ℓ ∈ 0:ℓₘₐₓ
                atol = 20 * max(1, ℓ) * eps(T)
                @test opnorm(parent(𝔇s[i₁][ℓ]) * parent(𝔇s[i₂][ℓ]) - parent(𝔇₁₂[ℓ])) ≤ atol
            end
        end
        for (i, R) ∈ enumerate(rotors)
            𝔇⁻¹ = D(inv(R), ℓₘₐₓ)
            for ℓ ∈ 0:ℓₘₐₓ
                atol = 20 * max(1, ℓ) * eps(T)
                @test opnorm(parent(𝔇⁻¹[ℓ]) - parent(𝔇s[i][ℓ])') ≤ atol
            end
            @test 𝔇s[i][0][0, 0] == 1
            @test size(parent(𝔇s[i][0])) == (1, 1)
        end
    end

    # Small ℓₘₐₓ: every pair drawn from the axis rotors plus a few random ones
    test_representation(property_rotors(T, 4, 2003), 3)

    # Larger ℓₘₐₓ: a handful of axis rotors (both poles and β = π/2) plus random ones
    let ℓₘₐₓ = 30
        rotors = [
            Rotor{T}(1); Rotor{T}(𝐢); Rotor{T}((1 + 𝐣) / √T(2)); -Rotor{T}(𝐤);
            randn(Xoshiro(2000 + ℓₘₐₓ), Rotor{T}, 4)
        ]
        test_representation(rotors, ℓₘₐₓ)
    end
end


@testitem "Wigner D double cover" setup=[WignerPropertyRotors] begin
    import SphericalFunctions: D
    using Quaternionic: Rotor
    using Random: Xoshiro

    # 𝔇ˡ(-R) = (-1)^{2ℓ} 𝔇ˡ(R).  This holds *exactly*: eⁱᵝ is built from squares of the
    # components, and the phases are powers of z₊ and z₋, which each flip sign under R → -R
    # and appear with a total exponent (m′+m) + (m′-m) = 2m′, whose parity is that of 2ℓ;
    # `complex_powers!` rotates z and -z to the same phase in the half-open sector
    # -π/4 < arg z ≤ π/4, so (-z)^k is computed as exactly (-1)^k z^k.  That includes the
    # phases on the diagonals, which are the edges of the sector:
    # z₊ = (R[1] + i R[4]) / √(R[1]² + R[4]²) is on one when |R[1]| = |R[4]| ≠ 0, and
    # z₋ = (R[3] - i R[2]) / √(R[2]² + R[3]²) when |R[2]| = |R[3]| ≠ 0.  Four of the axis
    # rotors are of that kind, and so are the rotors of `diagonal`, which include turns by
    # 2π/3 about the diagonals of a cube and by π about the diagonals of its faces.
    T = Float64
    ℓₘₐₓ = 8
    diagonal = [
        Rotor(T[1, 1, 1, 1]...), Rotor(T[1, 1, 1, -1]...), Rotor(T[1, -1, -1, 1]...),
        Rotor(T[0, 1, 1, 0]...), Rotor(T[0, 1, 1, 1]...),
        Rotor(T[0.5, 0.3, 0.2, -0.5]...), Rotor(T[0.3, 0.4, 0.4, 0.8]...),
    ]
    @test all(R -> abs(R[1]) == abs(R[4]) ≠ 0 || abs(R[2]) == abs(R[3]) ≠ 0, diagonal)
    for R ∈ [randn(Xoshiro(3000), Rotor{T}, 8); axis_rotors(T); diagonal]
        𝔇₊ = D(R, ℓₘₐₓ)
        𝔇₋ = D(-R, ℓₘₐₓ)
        for ℓ ∈ 0:ℓₘₐₓ
            @test parent(𝔇₋[ℓ]) == parent(𝔇₊[ℓ])
        end
        𝔇₊ = D(R, ℓₘₐₓ - 1//2)
        𝔇₋ = D(-R, ℓₘₐₓ - 1//2)
        for ℓ ∈ 1//2:ℓₘₐₓ - 1//2
            @test parent(𝔇₋[ℓ]) == -parent(𝔇₊[ℓ])
        end
    end
end


@testitem "Wigner d symmetries" begin
    import SphericalFunctions: d
    using LinearAlgebra: I
    using Random: Xoshiro

    T = Float64
    ℓₘₐₓ = 8
    rng = Xoshiro(4000)
    βs = T[T(π) .* rand(rng, 8); 0; nextfloat(zero(T)); T(π) / 2; prevfloat(T(π)); T(π)]

    for β ∈ βs
        d₊ = d(β, ℓₘₐₓ)
        d₋ = d(-β, ℓₘₐₓ)
        for ℓ ∈ 0:ℓₘₐₓ
            @test axes(d₊[ℓ]) == (-ℓ:ℓ, -ℓ:ℓ)
            @test parent(d₊[ℓ]) isa Matrix{T}
            for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                # d_{m′m}(β) = (-1)^{m′-m} d_{mm′}(β) = d_{-m,-m′}(β)
                @test d₊[ℓ][m′, m] ≈ (-1)^(m′ - m) * d₊[ℓ][m, m′] atol=4eps(T)
                @test d₊[ℓ][m′, m] ≈ d₊[ℓ][-m, -m′] atol=4eps(T)
                # d(-β) = d(β)ᵀ, exactly: cis(-β) = conj(cis(β)) and every element has a
                # definite parity in sin β
                @test d₋[ℓ][m, m′] == d₊[ℓ][m′, m]
            end
        end
    end

    # β may be any Real (or the phase eⁱᵝ); the result is computed in float(typeof(β))
    @test d(1, ℓₘₐₓ) == d(one(T), ℓₘₐₓ)
    @test d(cis(one(T)), ℓₘₐₓ) == d(one(T), ℓₘₐₓ)
    @test eltype(d(one(Float32), ℓₘₐₓ)[ℓₘₐₓ]) === Float32
    @test eltype(d(one(BigFloat), ℓₘₐₓ)[ℓₘₐₓ]) === BigFloat

    # β = 0: the identity.  (Not exactly: the recurrence leaves a sub-eps residue, e.g. an
    # anti-diagonal of ±1.6e-16 at ℓ = 5, so the whole matrix is checked to a few eps.)
    d₀ = d(zero(T), ℓₘₐₓ)
    for ℓ ∈ 0:ℓₘₐₓ
        @test maximum(abs, parent(d₀[ℓ]) - I) ≤ 4eps(T)
    end

    # β = π: anti-diagonal, dˡ_{m′m}(π) = (-1)^{ℓ-m} δ_{m′,-m}.  The sign follows from the
    # Varshalovich closed form Eq. 4.3.1(2): at β = π only cos(β/2)⁰ survives, which forces
    # m′ = -m and leaves the single term (-1)^{ℓ-m}.  (That closed form is also applied
    # directly at β = π by the "dCalculator vs closed form" item.)
    # With β = π as a Float64, sin β ≈ 1.2e-16 rather than 0, so the off-anti-diagonal
    # elements are O(ℓ eps) rather than exactly zero; with eⁱᵝ = -1 given exactly they
    # vanish exactly.
    dπ = d(T(π), ℓₘₐₓ)
    dπ_exact = d(complex(-one(T)), ℓₘₐₓ)
    for ℓ ∈ 0:ℓₘₐₓ
        for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
            expected = m′ == -m ? T((-1)^(ℓ - m)) : zero(T)
            @test dπ[ℓ][m′, m] ≈ expected atol=2(ℓ + 1) * eps(T)
            @test dπ_exact[ℓ][m′, m] ≈ expected atol=4eps(T)
            if m′ != -m
                @test dπ_exact[ℓ][m′, m] == 0
            end
        end
    end
end


@testitem "Wigner D Euler-angle factorization" setup=[Utilities] begin
    import SphericalFunctions: D, d
    import .Utilities: sYlm_closed_form
    using Quaternionic: from_euler_angles, from_spherical_coordinates
    using Random: Xoshiro

    T = Float64
    ℓₘₐₓ = 6
    rng = Xoshiro(5000)

    # 𝔇ˡ_{m′m}(R_{αβγ}) = e^{-i m′ α} dˡ_{m′m}(β) e^{-i m γ}, including at the poles β ∈ {0, π}
    # where α and γ are individually undefined but the combination is not.
    angles = [(2T(π) * rand(rng), T(π) * rand(rng), 2T(π) * rand(rng)) for _ ∈ 1:10]
    append!(
        angles,
        [
            (T(0.9), zero(T), T(2.6)), (T(0.9), T(π), T(2.6)),
            (zero(T), zero(T), zero(T)), (T(π), T(π) / 2, T(π)),
            (2T(π) - T(1e-9), T(1e-9), T(0.3))
        ]
    )
    for (α, β, γ) ∈ angles
        𝔇 = D(from_euler_angles(α, β, γ), ℓₘₐₓ)
        dβ = d(β, ℓₘₐₓ)
        for ℓ ∈ 0:ℓₘₐₓ, m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
            @test 𝔇[ℓ][m′, m] ≈ cis(-m′ * α) * dβ[ℓ][m′, m] * cis(-m * γ) atol=40eps(T)
        end
    end

    # Spherical-harmonic special case, which pins the overall convention (issue #42):
    #   ₀Yₗₘ(θ, ϕ) = √((2ℓ+1)/4π) conj(𝔇ˡ_{m0}(R)),  R = from_spherical_coordinates(θ, ϕ)
    # with `sYlm_closed_form` the closed-form expression from the `Utilities` module.
    θϕs = [(T(π) * rand(rng), 2T(π) * rand(rng)) for _ ∈ 1:6]
    append!(θϕs, [(zero(T), one(T)), (T(π), T(2)), (T(π) / 2, zero(T)), (T(1e-9), T(0.5))])
    for (θ, ϕ) ∈ θϕs
        𝔇 = D(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ)
        for ℓ ∈ 0:ℓₘₐₓ, m ∈ -ℓ:ℓ
            @test 𝔇[ℓ][m, 0] ≈ √(4T(π) / (2ℓ + 1)) * conj(sYlm_closed_form(0, ℓ, m, θ, ϕ)) atol=1e-13
        end
    end
end


@testitem "Wigner d character identity" begin
    import SphericalFunctions: dCalculator, recurrence!, d, array_view

    # χˡ(β) = Σₘ dˡₘₘ(β) = sin((2ℓ+1)β/2) / sin(β/2)
    T = Float64
    βs = T[0.3, 1.1, 2.0, 2.9]
    ℓs = (5, 30, 100, 300)
    χ(ℓ, β) = sin((2ℓ + 1) * β / 2) / sin(β / 2)
    # Float64 accumulates ~ℓ eps in the trace; measured: at most 2.4ℓ eps
    atol(ℓ) = 20ℓ * eps(T)

    # Only the block of the largest ℓ is needed, so a calculator is stepped up to it, rather
    # than `d(β, ℓ)` materializing every block (which, at ℓ = 300, allocates 300 MB)
    single(β, ℓ) = recurrence!(dCalculator(β, ℓ), ℓ)
    for ℓ ∈ ℓs, β ∈ βs
        dβ = single(β, ℓ)
        @test sum(dβ[m, m] for m ∈ -ℓ:ℓ) ≈ χ(ℓ, β) atol=atol(ℓ)
    end

    # The same identity through the batched calculator (one recurrence for all β at once),
    # whose blocks must agree exactly with those of a calculator for each β alone, and with
    # the convenience function `d`, which is checked at the smallest ℓ
    calc = dCalculator(βs, maximum(ℓs))
    for ℓ ∈ 1:maximum(ℓs)
        block = recurrence!(calc, ℓ)
        if ℓ ∈ ℓs
            @test axes(block) == (1:length(βs), -ℓ:ℓ, -ℓ:ℓ)
            for (iᵣ, β) ∈ enumerate(βs)
                @test sum(block[iᵣ, m, m] for m ∈ -ℓ:ℓ) ≈ χ(ℓ, β) atol=atol(ℓ)
                @test array_view(block)[iᵣ, :, :] == array_view(single(β, ℓ))
                ℓ == minimum(ℓs) && @test array_view(block)[iᵣ, :, :] == array_view(d(β, ℓ)[ℓ])
            end
        end
    end
end


@testitem "Wigner d Float64 vs BigFloat" begin
    import SphericalFunctions: dCalculator, recurrence!, array_view
    using DoubleFloats: Double64

    # The block of the largest ℓ alone, from a calculator stepped up to it, rather than every
    # block from `d(β, ℓ)`
    block(β, ℓ) = array_view(recurrence!(dCalculator(β, ℓ), ℓ))

    # Float64 against BigFloat (the same β in both, so only the arithmetic differs)
    for ℓ ∈ (25, 100), β ∈ (0.4, 1.3, 2.7)
        d₆₄ = block(β, ℓ)
        dᵦᵢ = block(big(β), ℓ)
        @test eltype(d₆₄) === Float64
        @test eltype(dᵦᵢ) === BigFloat
        Δ = abs.(d₆₄ .- dᵦᵢ)
        @test maximum(Δ) ≤ 5e-15
        significant = abs.(dᵦᵢ) .> 1e-3
        @test maximum(Δ[significant] ./ abs.(dᵦᵢ[significant])) ≤ 1e-12
    end

    # Double64 against BigFloat
    let ℓ = 40
        for β ∈ (0.4, 1.3, 2.7)
            d₆₄₆₄ = block(Double64(β), ℓ)
            dᵦᵢ = block(big(β), ℓ)
            @test eltype(d₆₄₆₄) === Double64
            @test maximum(abs.(BigFloat.(d₆₄₆₄) .- dᵦᵢ)) ≤ 1e-30
        end
    end
end


@testitem "Wigner D Float64 vs BigFloat, with phases just below an axis" begin
    import SphericalFunctions: DCalculator, recurrence!, array_view
    using Quaternionic: Rotor, from_euler_angles

    # 𝔇 multiplies d by powers of the phases z₊ = e^{i(α+γ)/2} and z₋ = e^{i(α-γ)/2}, so its
    # error is that of d plus that of the powers.  These rotors put both phases at a small
    # distance δ below an axis, and `complex_powers!` rotates each of them to e^{-iδ}, near 1.
    # At ℓ = 32 the error in the block, against the same rotor in BigFloat, is at most 3.6 eps,
    # and the bound is 6 eps.  Were the phases rotated to e^{i(π/2 - δ)} instead, where Re z is
    # near 0, the errors would be 6.2 to 14.7 eps for these rotors.
    block(R, ℓ) = array_view(recurrence!(DCalculator(R, ℓ), ℓ))
    ℓ = 32
    for (α, β, γ) ∈ (
        (π - 2e-3, 0.9, 0.0), (-2e-12, 0.9, 0.0),
        (2π - 0.035, 1.3, 0.0), (π/2 - 0.035, 0.9, -π/2),
    )
        R = from_euler_angles(α, β, γ)
        @test maximum(abs, block(R, ℓ) .- block(Rotor{BigFloat}(R), ℓ)) ≤ 6eps(Float64)
    end
end

@testitem "Wigner D large ℓ sanity" begin
    import SphericalFunctions: SphericalFunctions, DCalculator, recurrence!
    using Quaternionic: Rotor
    using LinearAlgebra: I, opnorm
    using Random: Xoshiro

    # One rotor, ℓₘₐₓ = 512, stepping the calculator through every ℓ (D(R, 512) would
    # allocate 513 matrices).  No NaN/Inf anywhere, and the last block is still unitary.
    T = Float64
    L = 512
    R = randn(Xoshiro(8000), Rotor{T})
    calc = DCalculator(R, L)
    @test SphericalFunctions.ℓₘₐₓ(calc) == L
    @test SphericalFunctions.Nᵣ(calc) == 1

    @test recurrence!(calc, 0)[0, 0] == 1
    @test SphericalFunctions.ℓ(calc) == 0
    nonfinite_ℓs = Int[]  # mutated, not reassigned, so no soft-scope ambiguity in the loop
    for ℓ ∈ 1:L
        all(isfinite, recurrence!(calc, ℓ)) || push!(nonfinite_ℓs, ℓ)
    end
    @test isempty(nonfinite_ℓs)
    @test SphericalFunctions.ℓ(calc) == L
    @test_throws MethodError calc[L - 1]

    M = parent(recurrence!(calc, L))
    @test size(M) == (2L + 1, 2L + 1)
    @test opnorm(M * M' - I) ≤ 1e-11
    @test opnorm(M' * M - I) ≤ 1e-11
end
