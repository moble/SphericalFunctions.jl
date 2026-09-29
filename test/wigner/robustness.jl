# Robustness tests for the Wigner engine: uninitialized memory, generic number types,
# differentiability, and interface details (show, allocation, independence of calculators).

# A refusal is recognized by the type of the exception and by a fragment of its message, so
# that an unrelated error raised later, for another reason, cannot stand in for it.  (A
# function given to `@test_throws` is applied to the message rather than to the exception, so
# it cannot check the type.)  This is used by the test items of the Wigner and harmonic
# calculators and of their containers, here and in `test/wigner/*.jl`, `test/sYlm/*.jl`,
# `test/hwedge.jl`, `test/haxis.jl` and `test/mode_weights/mode_weights.jl`.
@testsnippet RefusalChecks begin
    function refuses(f, ::Type{T}, fragment::AbstractString) where {T<:Exception}
        try
            f()
        catch e
            return e isa T && occursin(fragment, sprint(showerror, e))
        end
        false
    end
end

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
            eⁱᵝNC=single(Complex{NC}.(cis.(βs))),
            R=single([from_euler_angles(t...) for t in αβγ]),
            RNC=single([Rotor{NC}(from_euler_angles(t...)) for t in αβγ]),
        )
    end

    # The smallest index of the kind of ℓₘₐₓ, which is 0 for integers and 1/2 for
    # half-integers, given here as `Rational`s
    lowest(ℓₘₐₓ) = ℓₘₐₓ isa Integer ? 0 : 1//2

    # Every ℓ in order, then a backwards jump (which restarts the recurrence), then a forward
    # jump that skips the intermediate values.
    function schedule(ℓₘₐₓ)
        lo = lowest(ℓₘₐₓ)
        [collect(lo:1:ℓₘₐₓ); lo + (ℓₘₐₓ - lo) ÷ 2; ℓₘₐₓ]
    end

    # Symmetric limits m′ₘₐₓ ∈ (lo, lo+1, ℓₘₐₓ), plus a couple of asymmetric blocks.  For
    # half-integers the rows m′ = ±1/2 must both be present, so the asymmetric windows stop
    # at -1/2 rather than at 0.
    function block_limits(ℓₘₐₓ)
        lo = lowest(ℓₘₐₓ)
        limits = [
            (m′ₘₐₓ=k, m′ₘᵢₙ=-k, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-ℓₘₐₓ)
            for k in unique((lo, lo + 1, ℓₘₐₓ)) if k ≤ ℓₘₐₓ
        ]
        if ℓₘₐₓ ≥ lo + 2
            push!(limits, (m′ₘₐₓ=lo+2, m′ₘᵢₙ=-lo-1, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-min(lo+3, ℓₘₐₓ)))
            push!(limits, (m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-lo, mₘₐₓ=lo+1, mₘᵢₙ=-ℓₘₐₓ))
        end
        limits
    end

    # Read every element of the current block of the checked calculator (any read of a NaN
    # throws) and compare to the plain Float64 calculator.  The checked inputs are the plain
    # ones converted, not recomputed, and `Checked` performs each operation in `Float64`, so
    # the two calculators agree bit for bit, and `atol` is zero.
    function check_block(block, blockF, atol)
        # `array_view` is what turns a labelled block into a plain 1-based array; the containers
        # deliberately have no linear indexing of their own, so `eachindex` goes through it.
        A = array_view(block)
        B = array_view(blockF)
        axes(A) == axes(B) || return false
        for i in eachindex(A)
            abs(unchecked(A[i]) - B[i]) ≤ atol || return false
        end
        true
    end

    # Same for the H wedge, reading every (m′, m) with |m′| ≤ m′ₘₐₓ through `wedge_value`.
    function check_wedge(calc, calcF, ℓ, atol)
        ℓ = ℓ isa Rational ? SphericalFunctions.HalfOddInteger(ℓ) : ℓ  # the wedge's index type
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

    # The half-integer ℓₘₐₓ exercise the paths of their own: the half-angle buffers, the seed
    # of the rows m′ = ±1/2, whose last column is peeled because the axis slot it would read
    # holds another order's data (or, after `fill!`, a NaN that a zero coefficient would not
    # annihilate), and the half-integer phases.
    for ℓₘₐₓ in (0, 1, 2, 5, 9, 1//2, 5//2, 9//2), Nᵣ in (1, 3)
        data = rotor_data(Nᵣ)
        atol = 0.0
        lo = lowest(ℓₘₐₓ)

        # The raw H engine, driven by β
        for m′ₘₐₓ in unique((lo, lo + 1, ℓₘₐₓ))
            m′ₘₐₓ ≤ ℓₘₐₓ || continue
            calc = HCalculator(data.βNC, ℓₘₐₓ; m′ₘₐₓ)
            calcF = HCalculator(data.β, ℓₘₐₓ; m′ₘₐₓ)
            fill!(calc, NaN)
            @test all(isnan, parent(calc.Hˡ))
            @test all(isnan, parent(calc.h⃗ᵃ))
            @test all(isnan, parent(calc.h⃗ᵇ))
            recurrence!(calc, lo)
            recurrence!(calcF, lo)
            if ℓₘₐₓ > lo
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
                    @test check_block(recurrence!(calc, ℓ), recurrence!(calcF, ℓ), atol)
                end
                # Start over from NaN and jump straight to ℓₘₐₓ
                fill!(calc, NaN)
                @test check_block(
                    recurrence!(calc, RNC, ℓₘₐₓ), recurrence!(calcF, ℓₘₐₓ), atol
                )
                # And again *without* re-supplying the rotors, which `fill!` promises to
                # keep.  This is the strongest form of the check: every element the block
                # needs must be rewritten by the recurrence from the surviving rotor data
                # alone, or a signaling NaN is read.
                fill!(calc, NaN)
                @test check_block(recurrence!(calc, ℓₘₐₓ), recurrence!(calcF, ℓₘₐₓ), atol)
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

    # Values agree with the Float64 path and derivatives with central finite differences,
    # which with h = 1e-6 are accurate to about eps/h ≈ 2e-10 (measured at most 1.5e-10 here)
    d₀ = d(β, ℓₘₐₓ)
    d₊ = d(β + h, ℓₘₐₓ)
    d₋ = d(β - h, ℓₘₐₓ)
    for ℓ in 0:ℓₘₐₓ
        @test axes(dd[ℓ]) == axes(d₀[ℓ])
        for m′ in -ℓ:ℓ, m in -ℓ:ℓ
            @test value(dd[ℓ][m′, m]) ≈ d₀[ℓ][m′, m] atol=4*max(1, ℓ)*eps()
            fd = (d₊[ℓ][m′, m] - d₋[ℓ][m′, m]) / 2h
            @test deriv(dd[ℓ][m′, m]) ≈ fd atol=1e-8
        end
    end

    # The same through an explicit calculator, with Nᵣ > 1.  The vector of `Dual` angles is
    # itself what makes the calculator a `Dual` one.
    calc = dCalculator([βd, ForwardDiff.Dual(2β, one(β))], ℓₘₐₓ)
    blk = recurrence!(calc, ℓₘₐₓ)
    dd2 = d(ForwardDiff.Dual(2β, one(β)), ℓₘₐₓ)
    @test blk[1] == dd[ℓₘₐₓ]
    @test blk[2] == dd2[ℓₘₐₓ]

    # 𝔇 with a Rotor{Dual}: derivative with respect to β vs finite differences
    𝔇(α, θ, γ) = D(from_euler_angles(α, θ, γ), 3)[3]
    g = ForwardDiff.derivative(θ -> real(𝔇(0.3, θ, 1.1)[2, -1]), β)
    fd = (real(𝔇(0.3, β + h, 1.1)[2, -1]) - real(𝔇(0.3, β - h, 1.1)[2, -1])) / 2h
    @test g ≈ fd atol=1e-8
    for (m′, m) in ((2, -1), (-3, 3), (0, 1), (1, 0), (3, 3))
        gr = ForwardDiff.derivative(θ -> real(𝔇(0.3, θ, 1.1)[m′, m]), β)
        gi = ForwardDiff.derivative(θ -> imag(𝔇(0.3, θ, 1.1)[m′, m]), β)
        fdr = (real(𝔇(0.3, β + h, 1.1)[m′, m]) - real(𝔇(0.3, β - h, 1.1)[m′, m])) / 2h
        fdi = (imag(𝔇(0.3, β + h, 1.1)[m′, m]) - imag(𝔇(0.3, β - h, 1.1)[m′, m])) / 2h
        @test gr ≈ fdr atol=1e-8
        @test gi ≈ fdi atol=1e-8
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
    blkD = recurrence!(calcD, 3)
    @test eltype(blkD) <: Complex{<:ForwardDiff.Dual}
    @test value(real(blkD[2, -1])) ≈ real(𝔇(0.3, β, 1.1)[2, -1]) atol=40eps()
    @test deriv(real(blkD[2, -1])) ≈ fd atol=1e-8
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
                # and the selection it drives are 1-based, so both go through `array_view`:
                # the containers have no logical indexing of their own.
                A32, A64 = array_view(D32[ℓ]), array_view(D64[ℓ])
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
                a32, a64 = array_view(d32[ℓ]), array_view(d64[ℓ])
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
        blk = recurrence!(calc, ℓₘₐₓ)
        @test eltype(blk) === Complex{T}
        for (i, R) in enumerate(Rs)
            @test all(isapprox.(array_view(blk[i]), array_view(D(R, ℓₘₐₓ)[ℓₘₐₓ]); atol, rtol))
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
        blk = recurrence!(calc, ℓₘₐₓ)
        @test eltype(blk) === T
        @test all(isfinite, blk)
        for (i, β) in enumerate(βs)
            @test all(isapprox.(array_view(blk[i]), array_view(d(Float64(β), ℓₘₐₓ)[ℓₘₐₓ]); atol))
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
                @test occursin("HWedge", sprint(show, MIME("text/plain"), calc.Hˡ))
            else
                s = sprint(
                    show, MIME("text/plain"),
                    dCalculator(βs32, ℓₘₐₓ; m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=3)
                )
                @test occursin("Float32", s)
                @test occursin("m′=-1:2", s)
                @test occursin("m=-3:3", s)
                # The returned block can be displayed
                @test occursin("for ℓ=4", sprint(show, MIME("text/plain"), recurrence!(calc, 4)))
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
    for (calc, name) in (
        (DCalculator(R, 3), "DCalculator"), (dCalculator(β, 3), "dCalculator"),
        (HCalculator(β, 3), "HCalculator"),
    )
        @test occursin(name, sprint(show, calc))
        @test repr(calc) == sprint(show, calc)
        @test occursin(name, sprint(show, [calc]))
    end
    calc = HCalculator(β, 3)
    @test occursin("HWedge", sprint(show, calc.Hˡ))
    @test occursin("HAxis", sprint(show, calc.h⃗ᵃ))
end


@testitem "Wigner calculators allocation" begin
    import SphericalFunctions: DCalculator, recurrence!
    import Quaternionic: Rotor
    import Random

    # Measure inside functions so that the calculator's type is concrete at the call site.
    alloc_H(calc, ℓ) = @allocated recurrence!(calc.H, ℓ)
    alloc_D(calc, ℓ) = @allocated recurrence!(calc, ℓ)
    alloc_set(calc, R, ℓ) = @allocated recurrence!(calc, R, ℓ)

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
    @test aD == 0
    recurrence!(calc, ℓ - 1)
    aD_step = alloc_D(calc, ℓ)
    @test aD_step == 0

    # Setting the rotors.  (The block `recurrence!` returns is a view, and is counted in
    # `aD` above; there is no separate indexing step to measure.)
    alloc_set(calc, Rs, ℓ)
    aS = alloc_set(calc, Rs, ℓ)
    @test aS == 0
    @info "Allocation (bytes)" aH aH_step aH_restart aD aD_step aS
end


@testitem "Wigner calculators thread safety via similar" begin
    import SphericalFunctions: DCalculator, dCalculator, HCalculator, sYlmCalculator,
        sλlmCalculator, recurrence!, wedge_value
    import SphericalFunctions
    import Quaternionic: Rotor
    import Random

    rng = Random.Xoshiro(8)
    ℓₘₐₓ = 20
    Rs = randn(rng, Rotor{Float64}, 8)

    # Whether two objects share no mutable part, at any depth.  Arrays are the leaves: each
    # must be a different object in the two, except that an empty one is exempt, because Julia
    # makes every empty `Memory` of a type one shared object, and there is nothing in them to
    # share (the half-angle buffers of an integer-index calculator are empty).  Every other
    # mutable object, such as an `HWedge`, an `HAxis` or a `Ref`, must be a different object
    # too, and is searched in turn, as is every immutable one.
    function unshared(a, b)
        all(fieldnames(typeof(a))) do f
            x, y = getfield(a, f), getfield(b, f)
            if isbitstype(typeof(x)) || x isa Union{Type, Symbol, Module}
                true
            elseif x isa AbstractArray
                x !== y || isempty(x)
            elseif ismutable(x)
                x !== y && unshared(x, y)
            else
                unshared(x, y)
            end
        end
    end
    # The helper does find storage shared at depth, inside distinct mutable wrappers
    let v = [1.0]
        @test !unshared((w=Ref(v),), (w=Ref(v),))
        @test unshared((w=Ref(v),), (w=Ref(copy(v)),))
    end

    # `similar` gives an independent calculator with the same sizes and types, holding a copy
    # of the same rotor data with nothing computed.  The d and H calculators keep only β, so
    # the rotors handed to them here serve only to fix Nᵣ = 3 and 4 — and, for the Float32 d
    # calculator, the element type it works in.
    # The half-integer calculators are included because only they use the half-angle
    # buffers `cβ½` and `sβ½`, and the harmonic calculators built from angles because only
    # they hold the `phases` flag false.
    for calc in (
        DCalculator(Rs[1:2], 5; m′ₘₐₓ=3, m′ₘᵢₙ=-2, mₘₐₓ=5, mₘᵢₙ=-4),
        dCalculator(Rotor{Float32}.(Rs[1:3]), 5; m′ₘₐₓ=1),
        HCalculator(Rs[1:4], 5; m′ₘₐₓ=2),
        DCalculator(Rs[1:2], 7//2; m′ₘₐₓ=3//2),
        HCalculator(Rs[1:4], 7//2),
        sYlmCalculator(Rs[1:3], 5, -2:2),
        sYlmCalculator(Rs[1], 9//2, -3//2:3//2),
        sYlmCalculator([0.3, 1.1], 5, 1),
        sλlmCalculator([0.3, 1.1], 9//2, 1//2),
    )
        c = similar(calc)
        @test typeof(c) === typeof(calc)
        @test c !== calc
        # Nothing that can change is shared — not the buffers, nor the `Ref`s such as `ℓ` and
        # `axes_valid` that record what has been computed — at any depth
        @test unshared(c, calc)
        for f in (SphericalFunctions.ℓₘₐₓ, SphericalFunctions.Nᵣ, SphericalFunctions.ℓₘᵢₙ)
            @test f(c) == f(calc)
        end
        H, Hc = calc isa HCalculator ? (calc, c) : (calc.H, c.H)
        @test SphericalFunctions.m′ₘₐₓ(Hc) == SphericalFunctions.m′ₘₐₓ(H)
        @test parent(Hc.Hˡ) !== parent(H.Hˡ)
        @test Hc.Hˡ.row_index !== H.Hˡ.row_index
        @test parent(Hc.h⃗ᵃ) !== parent(H.h⃗ᵃ)
        @test parent(Hc.h⃗ᵇ) !== parent(H.h⃗ᵇ)
        @test Hc.eⁱᵝ !== H.eⁱᵝ
        @test Hc.eⁱᵝ == H.eⁱᵝ
        @test Hc.cβ½ == H.cβ½ && Hc.sβ½ == H.sβ½
        if !(calc isa HCalculator)
            @test c.Z₊ !== calc.Z₊
            @test c.Z₋ !== calc.Z₋
            # A calculator built from angles never fills its phase buffers, which then hold
            # whatever the allocation held, NaN included, so the copies are compared with
            # `isequal`, under which a NaN equals itself
            @test isequal(c.Z₊, calc.Z₊)
            @test isequal(c.Z₋, calc.Z₋)
            # `similar` copies no results: `ℓ` reports that nothing has been computed
            @test SphericalFunctions.ℓ(c) == SphericalFunctions.ℓₘᵢₙ(c) - 1
            # ... and stepping it gives what the original gives, bit for bit
            @test all(
                copy(recurrence!(c, ℓ)) == copy(recurrence!(calc, ℓ))
                for ℓ ∈ SphericalFunctions.ℓₘᵢₙ(c):SphericalFunctions.ℓₘₐₓ(c)
            )
        end
        if calc isa SphericalFunctions.WignerCalculator
            @test c.Wˡ !== calc.Wˡ
            for f in (SphericalFunctions.m′ₘₐₓ, SphericalFunctions.m′ₘᵢₙ,
                    SphericalFunctions.mₘₐₓ, SphericalFunctions.mₘᵢₙ)
                @test f(c) == f(calc)
            end
        elseif calc isa SphericalFunctions.HarmonicCalculator
            @test c.Yˡ !== calc.Yˡ
            @test c.phases[] == calc.phases[] && c.phases !== calc.phases
            @test SphericalFunctions.spins(c) == SphericalFunctions.spins(calc)
        end
    end

    # 𝔇 for 8 rotors on 8 tasks, each with its own calculator, vs serial results
    calc = DCalculator(Rs[1], ℓₘₐₓ)
    serial = [[copy(recurrence!(calc, R, ℓ)) for ℓ in 0:ℓₘₐₓ] for R in Rs]
    recurrence!(calc, Rs[1], ℓₘₐₓ)  # leave the template holding data while the tasks run

    # Two calculators stepped alternately on one task interleave as thoroughly as threads
    # could, and deterministically.  The spawned tasks below prove nothing when
    # `Threads.nthreads() == 1`, since they then run one after another; this does.
    c₁, c₂ = similar(calc), similar(calc)
    @test all(
        recurrence!(c₁, Rs[2], ℓ) == serial[2][ℓ+1] && recurrence!(c₂, Rs[3], ℓ) == serial[3][ℓ+1]
        for ℓ in 0:ℓₘₐₓ
    )
    @test recurrence!(calc, ℓₘₐₓ) == serial[1][end]
    tasks = map(Rs) do R
        Threads.@spawn begin
            c = similar(calc)
            [copy(recurrence!(c, R, ℓ)) for ℓ in 0:ℓₘₐₓ]
        end
    end
    parallel = fetch.(tasks)
    @test parallel == serial
    @test recurrence!(calc, ℓₘₐₓ) == serial[1][end]  # the template was not disturbed

    # The same with a batched d calculator and interleaved ℓ orders
    βs = [rand(rng, 2) .* π for _ in 1:8]
    calcd = dCalculator(βs[1], ℓₘₐₓ)
    seriald = [[copy(recurrence!(calcd, β, ℓ)) for ℓ in 0:ℓₘₐₓ] for β in βs]
    tasksd = map(enumerate(βs)) do (i, β)
        Threads.@spawn begin
            c = similar(calcd)
            # Odd tasks go up, even tasks start at the top (forcing restarts) and go down
            order = isodd(i) ? (0:ℓₘₐₓ) : (ℓₘₐₓ:-1:0)
            out = Vector{Any}(undef, ℓₘₐₓ + 1)
            for ℓ in order
                out[ℓ + 1] = copy(recurrence!(c, β, ℓ))
            end
            out
        end
    end
    paralleld = fetch.(tasksd)
    @test paralleld == seriald

    # Raw H engines in parallel, compared through wedge_value
    calcH = HCalculator(βs[1][1], ℓₘₐₓ; m′ₘₐₓ=6)
    wedge(H, ℓ) = [wedge_value(H, 1, m′, m) for m′ in -min(ℓ, 6):min(ℓ, 6), m in -ℓ:ℓ]
    serialH = [[wedge(recurrence!(calcH, β[1], ℓ), ℓ) for ℓ in 0:ℓₘₐₓ] for β in βs]
    tasksH = map(βs) do β
        Threads.@spawn begin
            c = similar(calcH)
            [wedge(recurrence!(c, β[1], ℓ), ℓ) for ℓ in 0:ℓₘₐₓ]
        end
    end
    @test fetch.(tasksH) == serialH
end

@testitem "Offset rotor and output arrays are refused before any write" begin
    import SphericalFunctions: DCalculator, dCalculator, HCalculator, sYlmCalculator,
        sλlmCalculator, sYlm, sYlm!, sλlm!, sYlm_matrix, set_R!, set_β!, set_θ!, array_view,
        WignerMatrix, WignerMatrixBatch, DegreeBlock, DegreeBlockBatch, SpinMatrix,
        SpinMatrixBatch, WignerSeries, relabel
    import Quaternionic: from_euler_angles
    import OffsetArrays: OffsetVector, OffsetArray

    # Each of these writes the calculator's 1-based buffers (or the caller's output) at the
    # input's own indices, under `@inbounds`, so that an offset array would write outside
    # them.  Each must therefore be refused, with the ArgumentError from
    # `Base.require_one_based_indexing`, before anything is written.
    offset = "offset arrays are not supported"
    R = [from_euler_angles(0.1i, 0.2i, 0.3i) for i ∈ 1:2]
    β = [0.3, 0.4]
    Ro, βo = OffsetVector(R, 0:1), OffsetVector(β, 0:1)

    # Rotor data, at construction and when reset
    @test_throws offset DCalculator(Ro, 2)
    @test_throws offset dCalculator(βo, 2)
    @test_throws offset dCalculator(cis.(βo), 2)
    @test_throws offset dCalculator(βo, 7//2)
    @test_throws offset HCalculator(βo, 2)
    @test_throws offset sYlmCalculator(Ro, 2, 1)
    @test_throws offset sλlmCalculator(βo, 2, 1)
    @test_throws offset sYlm(Ro, 2, 1)
    @test_throws offset sYlm_matrix(Ro, 2, 1)
    @test_throws offset set_R!(DCalculator(R, 2), Ro)
    @test_throws offset set_R!(sYlmCalculator(R, 2, 1), Ro)
    @test_throws offset set_β!(dCalculator(β, 2), βo)
    @test_throws offset set_θ!(sλlmCalculator(β, 2, 1), βo)

    # A refused reset leaves the calculator as it was
    c = DCalculator(R, 2)
    @test_throws offset set_R!(c, Ro)
    @test c.H.eⁱᵝ == DCalculator(R, 2).H.eⁱᵝ

    # Output arrays, for one spin weight and for a range of them
    Y = OffsetVector(zeros(ComplexF64, 8), 0:7)
    @test_throws offset sYlm!(Y, R[1], 2, 1)
    @test_throws offset sYlm!(Y, sYlmCalculator(R[1], 2, 1), R[1])
    @test_throws offset sYlm!(OffsetArray(zeros(ComplexF64, 3, 9), 0:2, 0:8), R[1], 2, -1:1)
    @test_throws offset sλlm!(OffsetVector(zeros(8), 0:7), sλlmCalculator(0.3, 2, 1), 0.3)

    # The block containers and `WignerSeries`, which index their storage as 1-based, refuse
    # an offset parent however they are built: directly, by `relabel`, or as a series
    @test_throws offset WignerMatrix(OffsetArray(zeros(ComplexF64, 5, 5), -2:2, -2:2), 2)
    @test_throws offset WignerMatrixBatch(OffsetArray(zeros(ComplexF64, 2, 5, 5), 0:1, -2:2, -2:2), 2)
    @test_throws offset DegreeBlock(OffsetVector(zeros(5), -2:2), 2)
    @test_throws offset DegreeBlockBatch(OffsetArray(zeros(2, 5), 0:1, -2:2), 2)
    @test_throws offset SpinMatrix(OffsetArray(zeros(3, 5), -1:1, -2:2), 2; sₘₐₓ=1, sₘᵢₙ=-1)
    @test_throws offset SpinMatrixBatch(OffsetArray(zeros(2, 3, 5), 0:1, -1:1, -2:2), 2; sₘₐₓ=1, sₘᵢₙ=-1)
    w = WignerMatrix(zeros(ComplexF64, 5, 5), 2)
    @test_throws offset relabel(w, OffsetArray(zeros(ComplexF64, 5, 5), -2:2, -2:2))
    blocks = [WignerMatrix(zeros(ComplexF64, 2ℓ+1, 2ℓ+1), ℓ) for ℓ ∈ 0:2]
    @test_throws offset WignerSeries(OffsetVector(blocks, 0:2), 0, 2)
    @test WignerSeries(blocks, 0, 2)[2] === blocks[3]

    # The ordinary 1-based forms are unaffected
    @test array_view(sYlm!(zeros(ComplexF64, 8), R[1], 2, 1)) == array_view(sYlm(R[1], 2, 1))
end
