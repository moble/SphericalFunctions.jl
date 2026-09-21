# Robustness tests for the v3 Wigner engine: uninitialized memory, generic number types,
# differentiability, and interface details (show, allocation, independence of calculators).

@testitem "Wigner calculators never read uninitialized memory" begin
    import SphericalFunctions: HCalculator, DCalculator, dCalculator,
        recurrence!, wedge_value
    import SphericalFunctions
    import MathChecker: checked, unchecked, NaNError
    import Quaternionic: Rotor, from_euler_angles

    # A `Float64` that throws as soon as a NaN takes part in an operation.  Only the NaN
    # check is enabled: the `Precision` check would object to the plain `Float64` values
    # this test mixes in deliberately, and the others flag conditions it is not about.
    NC = checked(Float64; precision=false, nan=true, inf=false)

    # Sanity checks on the detector itself
    @test_throws NaNError NC(NaN) + NC(1.0)
    @test_throws NaNError NC(NaN) * NC(2.0)
    @test NC(1.0) + NC(2.0) == NC(3.0)

    # Rotor data: for Nᵣ == 1 a generic rotor; for Nᵣ == 3 a generic rotor plus both poles.
    function rotor_data(Nᵣ)
        αβγ = if Nᵣ == 1
            [(0.3, 0.7, 1.1)]
        else
            [(0.3, 0.7, 1.1), (0.0, 0.0, 0.0), (1.2, Float64(π), -0.4)]
        end
        single(v) = length(v) == 1 ? v[1] : v
        βs = [t[2] for t in αβγ]
        (
            β=single(βs),
            βNC=single(NC.(βs)),
            eⁱᵝ=single(cis.(βs)),
            eⁱᵝNC=single(cis.(NC.(βs))),
            R=single([from_euler_angles(t...) for t in αβγ]),
            RNC=single([from_euler_angles(NC.(t)...) for t in αβγ]),
        )
    end

    # Every ℓ in order, then a backwards jump (which restarts the recurrence), then a forward
    # jump that skips the intermediate values.
    schedule(ℓₘₐₓ) = [collect(0:ℓₘₐₓ); ℓₘₐₓ ÷ 2; ℓₘₐₓ]

    # Symmetric limits m′ₘₐₓ ∈ (0, 1, ℓₘₐₓ), plus a couple of asymmetric blocks.
    function block_limits(ℓₘₐₓ)
        limits = [
            (m′ₘₐₓ=k, m′ₘᵢₙ=-k, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-ℓₘₐₓ)
            for k in unique((0, 1, ℓₘₐₓ)) if k ≤ ℓₘₐₓ
        ]
        if ℓₘₐₓ ≥ 2
            push!(limits, (m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-min(3, ℓₘₐₓ)))
            push!(limits, (m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=0, mₘₐₓ=1, mₘᵢₙ=-ℓₘₐₓ))
        end
        limits
    end

    # Read every element of the current block of the checked calculator (any read of a NaN
    # throws) and compare to the plain Float64 calculator.  Exact equality is not possible:
    # `@fastmath` in `recurrence_step4!` and `complex_powers!` lets the Float64 path contract
    # multiply-adds into FMAs, while `Checked` arithmetic goes through the generic methods.
    function check_block(calc, calcF, ℓ, atol)
        # `strided` is what turns a labelled block into a plain 1-based array; the containers
        # deliberately have no linear indexing of their own, so `eachindex` goes through it.
        A = strided(calc[ℓ])
        B = strided(calcF[ℓ])
        axes(A) == axes(B) || return false
        for i in eachindex(A)
            abs(unchecked(A[i]) - B[i]) ≤ atol || return false
        end
        true
    end

    # Same for the H wedge, reading every (m′, m) with |m′| ≤ m′ₘₐₓ through `wedge_value`.
    function check_wedge(calc, calcF, ℓ, atol)
        H = calc.Hˡ
        HF = calcF.Hˡ
        SphericalFunctions.ℓ(H) == ℓ || return false
        m′ₘ = SphericalFunctions.m′ₘₐₓ(H)
        m′ₘ == min(ℓ, SphericalFunctions.m′ₘₐₓ(calc)) || return false
        for iᵣ in 1:SphericalFunctions.Nᵣ(H), m′ in -m′ₘ:m′ₘ, m in -ℓ:ℓ
            a = unchecked(wedge_value(H, iᵣ, m′, m))
            b = wedge_value(HF, iᵣ, m′, m)
            abs(a - b) ≤ atol || return false
        end
        true
    end

    for ℓₘₐₓ in (0, 1, 2, 5, 9), Nᵣ in (1, 3)
        data = rotor_data(Nᵣ)
        atol = 4 * max(1, ℓₘₐₓ) * eps(Float64)

        # The raw H engine, driven by β
        for m′ₘₐₓ in unique((0, 1, ℓₘₐₓ))
            m′ₘₐₓ ≤ ℓₘₐₓ || continue
            calc = HCalculator(data.βNC, ℓₘₐₓ; m′ₘₐₓ)
            calcF = HCalculator(data.β, ℓₘₐₓ; m′ₘₐₓ)
            fill!(calc, NaN)
            @test all(isnan, parent(calc.Hˡ))
            @test all(isnan, parent(calc.h⃗ᵃ))
            @test all(isnan, parent(calc.h⃗ᵇ))
            recurrence!(calc, 0)
            recurrence!(calcF, 0)
            if ℓₘₐₓ > 0
                # The detector works: storage for the larger ℓ values is still NaN, and
                # touching it throws.
                @test_throws NaNError parent(calc.Hˡ)[end] + NC(1.0)
            end
            for ℓ in schedule(ℓₘₐₓ)
                recurrence!(calc, ℓ)
                recurrence!(calcF, ℓ)
                @test check_wedge(calc, calcF, ℓ, atol)
            end
            # Start over from NaN and jump straight to ℓₘₐₓ
            fill!(calc, NaN)
            recurrence!(calc, data.βNC, ℓₘₐₓ)
            recurrence!(calcF, ℓₘₐₓ)
            @test check_wedge(calc, calcF, ℓₘₐₓ, atol)
        end

        # The d and 𝔇 calculators.  d is driven by eⁱᵝ (Nᵣ == 1) or by Rotors (Nᵣ == 3), 𝔇
        # by Rotors.
        dNC, dF = Nᵣ == 1 ? (data.eⁱᵝNC, data.eⁱᵝ) : (data.RNC, data.R)
        for lim in block_limits(ℓₘₐₓ)
            for (Calc, RNC, RF) in (
                (dCalculator, dNC, dF), (DCalculator, data.RNC, data.R)
            )
                calc = Calc(RNC, ℓₘₐₓ; lim...)
                calcF = Calc(RF, ℓₘₐₓ; lim...)
                # `fill!` poisons everything the recurrence is responsible for writing.
                # It deliberately leaves the rotor data (`eⁱᵝ`, the half angles, and the
                # phase powers `Z₊`, `Z₋`) alone — those are `set_rotors!`'s job, and its
                # docstring promises they survive — so they are not asserted NaN here.
                fill!(calc, NaN)
                @test all(isnan, parent(calc.H.Hˡ))
                @test all(x -> isnan(real(x)), calc.Wˡ)
                for ℓ in schedule(ℓₘₐₓ)
                    recurrence!(calc, ℓ)
                    recurrence!(calcF, ℓ)
                    @test check_block(calc, calcF, ℓ, atol)
                end
                # Start over from NaN and jump straight to ℓₘₐₓ
                fill!(calc, NaN)
                recurrence!(calc, RNC, ℓₘₐₓ)
                recurrence!(calcF, ℓₘₐₓ)
                @test check_block(calc, calcF, ℓₘₐₓ, atol)
                # And again *without* re-supplying the rotors, which `fill!` promises to
                # keep.  This is the strongest form of the check: every element the block
                # needs must be rewritten by the recurrence from the surviving rotor data
                # alone, or a signaling NaN is read.
                fill!(calc, NaN)
                recurrence!(calc, ℓₘₐₓ)
                @test check_block(calc, calcF, ℓₘₐₓ, atol)
            end
        end
    end
end


@testitem "Wigner calculators with ForwardDiff duals" begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D, d
    import ForwardDiff
    import Quaternionic: Rotor, from_euler_angles

    ℓₘₐₓ = 4
    β = 0.7
    h = 1e-6

    # d with a Dual β: one Dual computation gives every d^ℓ_{m′m} and its β-derivative.
    βd = ForwardDiff.Dual(β, one(β))
    dd = d(βd, ℓₘₐₓ)
    @test eltype(dd[ℓₘₐₓ]) <: ForwardDiff.Dual
    @test axes(dd) == (0:ℓₘₐₓ,)
    value(x) = ForwardDiff.value(x)
    deriv(x) = ForwardDiff.partials(x, 1)

    # Analytic derivatives of the ℓ = 1 elements
    @test value(dd[1][0, 0]) ≈ cos(β) atol=4eps()
    @test deriv(dd[1][0, 0]) ≈ -sin(β) atol=4eps()
    @test value(dd[1][1, 1]) ≈ (1 + cos(β)) / 2 atol=4eps()
    @test deriv(dd[1][1, 1]) ≈ -sin(β) / 2 atol=4eps()
    @test value(dd[1][1, 0]) ≈ -sin(β) / √2 atol=4eps()
    @test deriv(dd[1][1, 0]) ≈ -cos(β) / √2 atol=4eps()

    # Values agree with the Float64 path and derivatives with central finite differences
    d₀ = d(β, ℓₘₐₓ)
    d₊ = d(β + h, ℓₘₐₓ)
    d₋ = d(β - h, ℓₘₐₓ)
    for ℓ in 0:ℓₘₐₓ
        @test axes(dd[ℓ]) == axes(d₀[ℓ])
        for m′ in -ℓ:ℓ, m in -ℓ:ℓ
            @test value(dd[ℓ][m′, m]) ≈ d₀[ℓ][m′, m] atol=4*max(1, ℓ)*eps()
            fd = (d₊[ℓ][m′, m] - d₋[ℓ][m′, m]) / 2h
            @test deriv(dd[ℓ][m′, m]) ≈ fd atol=1e-6
        end
    end

    # The same through an explicit calculator, with Nᵣ > 1.  The vector of `Dual` angles is
    # itself what makes the calculator a `Dual` one.
    calc = dCalculator([βd, ForwardDiff.Dual(2β, one(β))], ℓₘₐₓ)
    recurrence!(calc, ℓₘₐₓ)
    dd2 = d(ForwardDiff.Dual(2β, one(β)), ℓₘₐₓ)
    @test calc[ℓₘₐₓ][1] == dd[ℓₘₐₓ]
    @test calc[ℓₘₐₓ][2] == dd2[ℓₘₐₓ]

    # 𝔇 with a Rotor{Dual}: derivative with respect to β vs finite differences
    𝔇(α, θ, γ) = D(from_euler_angles(α, θ, γ), 3)[3]
    g = ForwardDiff.derivative(θ -> real(𝔇(0.3, θ, 1.1)[2, -1]), β)
    fd = (real(𝔇(0.3, β + h, 1.1)[2, -1]) - real(𝔇(0.3, β - h, 1.1)[2, -1])) / 2h
    @test g ≈ fd atol=1e-6
    for (m′, m) in ((2, -1), (-3, 3), (0, 1), (1, 0), (3, 3))
        gr = ForwardDiff.derivative(θ -> real(𝔇(0.3, θ, 1.1)[m′, m]), β)
        gi = ForwardDiff.derivative(θ -> imag(𝔇(0.3, θ, 1.1)[m′, m]), β)
        fdr = (real(𝔇(0.3, β + h, 1.1)[m′, m]) - real(𝔇(0.3, β - h, 1.1)[m′, m])) / 2h
        fdi = (imag(𝔇(0.3, β + h, 1.1)[m′, m]) - imag(𝔇(0.3, β - h, 1.1)[m′, m])) / 2h
        @test gr ≈ fdr atol=1e-6
        @test gi ≈ fdi atol=1e-6
        # Derivatives with respect to α and γ are analytic in the settled convention
        # 𝔇 = e^{-im′α} d e^{-imγ}: ∂α𝔇 = -im′ 𝔇 and ∂γ𝔇 = -im 𝔇.
        𝔇₀ = 𝔇(0.3, β, 1.1)[m′, m]
        gα = ForwardDiff.derivative(α -> real(𝔇(α, β, 1.1)[m′, m]), 0.3)
        @test gα ≈ m′ * imag(𝔇₀) atol=1e-13
        gγ = ForwardDiff.derivative(γ -> imag(𝔇(0.3, β, γ)[m′, m]), 1.1)
        @test gγ ≈ -m * real(𝔇₀) atol=1e-13
    end

    # And through an explicit 𝔇 calculator holding Rotor{Dual} data.  The dual β promotes
    # into the rotor, and the calculator's type follows the rotor's — there is nothing else
    # left to tell it.
    Rd = from_euler_angles(0.3, βd, 1.1)
    @test Rd isa Rotor{<:ForwardDiff.Dual}
    calcD = DCalculator(Rd, 3)
    recurrence!(calcD, 3)
    @test eltype(calcD[3]) <: Complex{<:ForwardDiff.Dual}
    @test value(real(calcD[3][2, -1])) ≈ real(𝔇(0.3, β, 1.1)[2, -1]) atol=40eps()
    @test deriv(real(calcD[3][2, -1])) ≈ fd atol=1e-6
end


@testitem "Wigner calculators with Float32 and Float16" begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D, d
    import Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(3216)

    # Float32: the recurrence's absolute error grows like ℓ·eps, so on the small elements
    # a relative criterion alone is not meaningful; use the atol rule with an rtol on top.
    let ℓₘₐₓ = 30, T = Float32
        atol = 20 * ℓₘₐₓ * eps(T)
        rtol = 1e-4
        for _ in 1:3
            R = randn(rng, Rotor{Float64})
            D32 = D(Rotor{T}(R), ℓₘₐₓ)
            D64 = D(R, ℓₘₐₓ)
            @test eltype(D32[ℓₘₐₓ]) === Complex{T}
            @test axes(D32) == (0:ℓₘₐₓ,)
            for ℓ in 0:ℓₘₐₓ
                @test axes(D32[ℓ]) == axes(D64[ℓ])
                @test !any(isnan, D32[ℓ])
                @test all(isapprox.(D32[ℓ], D64[ℓ]; atol, rtol))
                # On the well-conditioned elements the relative error is small.  The mask
                # and the selection it drives are 1-based, so both go through `strided`:
                # the containers have no logical indexing of their own.
                A32, A64 = strided(D32[ℓ]), strided(D64[ℓ])
                big = abs.(A64) .> 0.1
                @test all(abs.(A32[big] .- A64[big]) .≤ rtol .* abs.(A64[big]))
            end

            β = rand(rng) * π
            d32 = d(T(β), ℓₘₐₓ)
            d64 = d(β, ℓₘₐₓ)
            @test eltype(d32[ℓₘₐₓ]) === T
            for ℓ in 0:ℓₘₐₓ
                @test !any(isnan, d32[ℓ])
                @test all(isapprox.(d32[ℓ], d64[ℓ]; atol, rtol))
                a32, a64 = strided(d32[ℓ]), strided(d64[ℓ])
                big = abs.(a64) .> 0.1
                @test all(abs.(a32[big] .- a64[big]) .≤ rtol .* abs.(a64[big]))
            end
            # The phase form of the input gives the same type and values
            d32ϕ = d(cis(T(β)), ℓₘₐₓ)
            @test eltype(d32ϕ[ℓₘₐₓ]) === T
            @test all(isapprox.(d32ϕ[ℓₘₐₓ], d64[ℓₘₐₓ]; atol, rtol))
        end

        # Batched calculators in Float32
        Rs = randn(rng, Rotor{Float64}, 4)
        calc = DCalculator(Rotor{T}.(Rs), ℓₘₐₓ)
        recurrence!(calc, ℓₘₐₓ)
        @test eltype(calc[ℓₘₐₓ]) === Complex{T}
        for (i, R) in enumerate(Rs)
            @test all(isapprox.(strided(calc[ℓₘₐₓ][i]), strided(D(R, ℓₘₐₓ)[ℓₘₐₓ]); atol, rtol))
        end
    end

    # Float16: finite and roughly right
    let ℓₘₐₓ = 8, T = Float16
        atol = 1e-2
        for _ in 1:3
            R = randn(rng, Rotor{Float64})
            D16 = D(Rotor{T}(R), ℓₘₐₓ)
            D64 = D(R, ℓₘₐₓ)
            @test eltype(D16[ℓₘₐₓ]) === Complex{T}
            for ℓ in 0:ℓₘₐₓ
                @test all(isfinite, D16[ℓ])
                @test all(isapprox.(D16[ℓ], D64[ℓ]; atol))
            end

            β = rand(rng) * π
            d16 = d(T(β), ℓₘₐₓ)
            d64 = d(β, ℓₘₐₓ)
            @test eltype(d16[ℓₘₐₓ]) === T
            for ℓ in 0:ℓₘₐₓ
                @test all(isfinite, d16[ℓ])
                @test all(isapprox.(d16[ℓ], d64[ℓ]; atol))
            end
        end
        βs = [T(0.4), T(1.9), T(3.0)]
        calc = dCalculator(βs, ℓₘₐₓ)
        recurrence!(calc, ℓₘₐₓ)
        @test eltype(calc[ℓₘₐₓ]) === T
        @test all(isfinite, calc[ℓₘₐₓ])
        for (i, β) in enumerate(βs)
            @test all(isapprox.(strided(calc[ℓₘₐₓ][i]), strided(d(Float64(β), ℓₘₐₓ)[ℓₘₐₓ]); atol))
        end
    end
end


@testitem "Wigner calculators show" begin
    import SphericalFunctions: HCalculator, DCalculator, dCalculator,
        recurrence!
    import Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(11)
    ℓₘₐₓ = 7
    for Nᵣ in (1, 3)
        Rs = randn(rng, Rotor{Float64}, Nᵣ)
        βs = rand(rng, Nᵣ) .* π
        βs32 = Float32.(βs)  # the element type shown is the data's own
        for (calc, name, R) in (
            (DCalculator(Rs, ℓₘₐₓ), "DCalculator", Rs),
            (dCalculator(βs, ℓₘₐₓ), "dCalculator", βs),
            (HCalculator(βs, ℓₘₐₓ), "HCalculator", βs),
        )
            s = sprint(show, MIME("text/plain"), calc)
            @test occursin(name, s)
            @test occursin("ℓₘₐₓ=$ℓₘₐₓ", s)
            @test occursin("Nᵣ=$Nᵣ", s)
            @test occursin("Float64", s)
            @test !occursin("ℓ=4", s)
            recurrence!(calc, Nᵣ == 1 ? R[1] : R, 4)
            s = sprint(show, MIME("text/plain"), calc)
            @test occursin(name, s)
            @test occursin("ℓₘₐₓ=$ℓₘₐₓ", s)
            @test occursin("Nᵣ=$Nᵣ", s)
            @test occursin("ℓ=4", s)
            # A different number type and non-default limits also show up
            if calc isa HCalculator
                @test occursin(
                    "m′ₘₐₓ=2",
                    sprint(show, MIME("text/plain"), HCalculator(βs32, ℓₘₐₓ; m′ₘₐₓ=2))
                )
                @test occursin(
                    "Float32", sprint(show, MIME("text/plain"), HCalculator(βs32, ℓₘₐₓ))
                )
                # The wedge itself can be displayed
                @test sprint(show, MIME("text/plain"), calc.Hˡ) isa String
            else
                s = sprint(
                    show, MIME("text/plain"),
                    dCalculator(βs32, ℓₘₐₓ; m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=3)
                )
                @test occursin("Float32", s)
                @test occursin("m′=-1:2", s)
                @test occursin("m=-3:3", s)
                # The returned block can be displayed
                @test sprint(show, MIME("text/plain"), calc[4]) isa String
            end
        end
    end
end


@testitem "Wigner calculators two-argument show" begin
    # `show(io, x)` without a MIME is what `repr`, `print`, `@show`, string interpolation,
    # and the display of a container holding a calculator fall back to.  It must not throw.
    import SphericalFunctions: HCalculator, DCalculator, dCalculator,
        recurrence!
    import Quaternionic: from_euler_angles

    R = from_euler_angles(0.3, 0.7, 1.1)
    β = 0.7
    for calc in (DCalculator(R, 3), dCalculator(β, 3), HCalculator(β, 3))
        @test sprint(show, calc) isa String
        @test repr(calc) isa String
        @test sprint(show, [calc]) isa String
    end
    calc = HCalculator(β, 3)
    @test sprint(show, calc.Hˡ) isa String
    @test sprint(show, calc.h⃗ᵃ) isa String
end


@testitem "Wigner calculators allocation" begin
    import SphericalFunctions: DCalculator, recurrence!
    import Quaternionic: Rotor
    import Random

    # Measure inside functions so that the calculator's type is concrete at the call site.
    alloc_H(calc, ℓ) = @allocated recurrence!(calc.H, ℓ)
    alloc_D(calc, ℓ) = @allocated recurrence!(calc, ℓ)
    alloc_set(calc, R, ℓ) = @allocated recurrence!(calc, R, ℓ)
    alloc_index(calc, ℓ) = @allocated calc[ℓ]

    rng = Random.Xoshiro(64)
    ℓₘₐₓ = 64
    Nᵣ = 8
    ℓ = 40
    Rs = randn(rng, Rotor{Float64}, Nᵣ)
    calc = DCalculator(Rs, ℓₘₐₓ)
    recurrence!(calc, Rs, ℓ)  # warm-up (compiles everything)

    # The H engine: recomputing the current ℓ, and stepping ℓ-1 → ℓ (which also runs step 2)
    alloc_H(calc, ℓ)
    aH = alloc_H(calc, ℓ)
    @test aH == 0
    recurrence!(calc.H, ℓ - 1)  # backwards: restarts from 0 and runs to ℓ-1
    aH_step = alloc_H(calc, ℓ)
    @test aH_step == 0
    # A full restart from ℓ = 0 up to ℓ
    recurrence!(calc.H, ℓ)
    alloc_restart(calc, ℓ) = @allocated (recurrence!(calc.H, 0); recurrence!(calc.H, ℓ))
    alloc_restart(calc, ℓ)
    aH_restart = alloc_restart(calc, ℓ)
    @test aH_restart == 0

    # The 𝔇 calculator (recurrence + materialization)
    alloc_D(calc, ℓ)
    aD = alloc_D(calc, ℓ)
    @test aD ≤ 512
    recurrence!(calc, ℓ - 1)
    aD_step = alloc_D(calc, ℓ)
    @test aD_step ≤ 512

    # Setting the rotors, and indexing (a view, so small)
    alloc_set(calc, Rs, ℓ)
    aS = alloc_set(calc, Rs, ℓ)
    @test aS ≤ 512
    alloc_index(calc, ℓ)
    aI = alloc_index(calc, ℓ)
    @test aI ≤ 512
    @info "Allocation (bytes)" aH aH_step aH_restart aD aD_step aS aI
end


@testitem "Wigner calculators thread safety via similar" begin
    import SphericalFunctions: DCalculator, dCalculator, HCalculator,
        recurrence!, wedge_value
    import SphericalFunctions
    import Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(8)
    ℓₘₐₓ = 20
    Rs = randn(rng, Rotor{Float64}, 8)

    # `similar` gives an independent calculator with the same sizes and types, holding a copy
    # of the same rotor data with nothing computed.  The d and H calculators keep only β, so
    # the rotors handed to them here serve only to fix Nᵣ = 3 and 4 — and, for the Float32 d
    # calculator, the element type it works in.
    for calc in (
        DCalculator(Rs[1:2], 5; m′ₘₐₓ=3, m′ₘᵢₙ=-2, mₘₐₓ=5, mₘᵢₙ=-4),
        dCalculator(Rotor{Float32}.(Rs[1:3]), 5; m′ₘₐₓ=1),
        HCalculator(Rs[1:4], 5; m′ₘₐₓ=2),
    )
        c = similar(calc)
        @test typeof(c) === typeof(calc)
        @test c !== calc
        for f in (
            SphericalFunctions.ℓₘₐₓ, SphericalFunctions.m′ₘₐₓ, SphericalFunctions.m′ₘᵢₙ,
            SphericalFunctions.Nᵣ,
        )
            @test f(c) == f(calc)
        end
        if calc isa HCalculator
            @test parent(c.Hˡ) !== parent(calc.Hˡ)
            @test parent(c.h⃗ᵃ) !== parent(calc.h⃗ᵃ)
            @test parent(c.h⃗ᵇ) !== parent(calc.h⃗ᵇ)
            @test c.eⁱᵝ !== calc.eⁱᵝ
            @test c.eⁱᵝ == calc.eⁱᵝ
        else
            @test SphericalFunctions.mₘₐₓ(c) == SphericalFunctions.mₘₐₓ(calc)
            @test SphericalFunctions.mₘᵢₙ(c) == SphericalFunctions.mₘᵢₙ(calc)
            @test parent(c.H.Hˡ) !== parent(calc.H.Hˡ)
            @test c.Wˡ !== calc.Wˡ
            @test c.Z₊ !== calc.Z₊
            @test c.Z₋ !== calc.Z₋
            @test c.Z₊ == calc.Z₊
            @test c.Z₋ == calc.Z₋
            @test_throws "currently holds" c[0]
        end
    end

    # 𝔇 for 8 rotors on 8 tasks, each with its own calculator, vs serial results
    calc = DCalculator(Rs[1], ℓₘₐₓ)
    serial = [[copy(recurrence!(calc, R, ℓ)[ℓ]) for ℓ in 0:ℓₘₐₓ] for R in Rs]
    recurrence!(calc, Rs[1], ℓₘₐₓ)  # leave the template holding data while the tasks run
    tasks = map(Rs) do R
        Threads.@spawn begin
            c = similar(calc)
            [copy(recurrence!(c, R, ℓ)[ℓ]) for ℓ in 0:ℓₘₐₓ]
        end
    end
    parallel = fetch.(tasks)
    @test parallel == serial
    @test calc[ℓₘₐₓ] == serial[1][end]  # the template was not disturbed

    # The same with a batched d calculator and interleaved ℓ orders
    βs = [rand(rng, 2) .* π for _ in 1:8]
    calcd = dCalculator(βs[1], ℓₘₐₓ)
    seriald = [[copy(recurrence!(calcd, β, ℓ)[ℓ]) for ℓ in 0:ℓₘₐₓ] for β in βs]
    tasksd = map(enumerate(βs)) do (i, β)
        Threads.@spawn begin
            c = similar(calcd)
            # Odd tasks go up, even tasks start at the top (forcing restarts) and go down
            order = isodd(i) ? (0:ℓₘₐₓ) : (ℓₘₐₓ:-1:0)
            out = Vector{Any}(undef, ℓₘₐₓ + 1)
            for ℓ in order
                out[ℓ + 1] = copy(recurrence!(c, β, ℓ)[ℓ])
            end
            out
        end
    end
    paralleld = fetch.(tasksd)
    @test paralleld == seriald

    # Raw H engines in parallel, compared through wedge_value
    calcH = HCalculator(βs[1][1], ℓₘₐₓ; m′ₘₐₓ=6)
    wedge(c, ℓ) = [wedge_value(c.Hˡ, 1, m′, m) for m′ in -min(ℓ, 6):min(ℓ, 6), m in -ℓ:ℓ]
    serialH = [[wedge(recurrence!(calcH, β[1], ℓ), ℓ) for ℓ in 0:ℓₘₐₓ] for β in βs]
    tasksH = map(βs) do β
        Threads.@spawn begin
            c = similar(calcH)
            [wedge(recurrence!(c, β[1], ℓ), ℓ) for ℓ in 0:ℓₘₐₓ]
        end
    end
    @test fetch.(tasksH) == serialH
end
