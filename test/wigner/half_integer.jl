# Stage-2 tests of the half-integer Wigner engine (design memo §6): `calc[ℓ][m′, m]`,
# `D` and `d` against the stage-1 oracle of `test/wigner/half_integer_oracle.jl`, the
# oracle-free metamorphic identities that stay valid at large `J`, the validation of
# `Rational` index ranges, and the half-integer containers.
#
# Every comparison accumulates a single worst-case error over a whole block (or a whole
# sweep) and asserts once.  A half-integer block has O(ℓ²) elements, and asserting element
# by element would turn one wrong sign into minutes of printed failures.  The accumulators
# live inside `let` blocks because a test item is evaluated at top level, where a bare
# `for` loop cannot assign to an enclosing variable.
#
# All tolerances below are quoted from measurements made against this tree; the measured
# worst case and the resulting margin are stated in a comment at the point of use.


@testitem "Half-integer oracle: the two references agree" setup=[HalfIntegerOracle] begin
    import Quaternionic: Rotor, from_euler_angles
    import .HalfIntegerOracle: WignerDElement, D_oracle, d_oracle, d_table, βvalues

    # Stage 1 of memo §6.  Varshalovich's closed form (Eq. 4.3.1(2), exact factorials) and
    # the quaternionic algorithm of Boyle (2016) are independent implementations; here they
    # are checked against each other and against Varshalovich's transcribed tables, so that
    # the oracle every item below leans on is anchored before it is used.  The two Literate
    # comparison pages make the same checks against the originals; repeating them here
    # catches a drift between those originals and the copies in `half_integer_oracle.jl`.

    αβγs = [(0.0, 0.0, 0.0), (0.7, 1.1, 2.3), (2.9, 0.4, 5.1), (4.0, 2.2, 0.3), (1.3, 2.9, 4.7)]

    # Boyle (2016) vs Varshalovich, both read in this package's convention, for J ≤ 15/2.
    # Measured worst case over this grid: 3.2e-15 (14.6 eps) at J = 15/2, dominated by the
    # Float64 arithmetic of `WignerDElement`; 1e-13 leaves a factor of ≳ 31.
    err = let e = 0.0
        for (α, β, γ) ∈ αβγs
            R = Rotor(from_euler_angles(α, β, γ))
            for J ∈ 1//2:1:15//2, m′ ∈ -J:J, m ∈ -J:J
                varsh = cis(-m′ * α) * d_oracle(J, m′, m, β) * cis(-m * γ)
                e = max(e, abs(varsh - D_oracle(R, J, m′, m)))
            end
        end
        e
    end
    @test err < 1e-13

    # The Boyle (2016) reference respects the double cover: 𝔇(-R) = -𝔇(R).  Unlike the
    # package's own 𝔇 (asserted bit-exactly in the metamorphic item below), the reference
    # forms its phases with `angle`/`cis` and so only reproduces this to rounding: measured
    # worst case 5.4e-15 over this grid, so 1e-13 leaves a factor of ≳ 18.
    errdc = let e = 0.0
        for (α, β, γ) ∈ αβγs
            R = Rotor(from_euler_angles(α, β, γ))
            for J ∈ 1//2:1:15//2, m′ ∈ -J:J, m ∈ -J:J
                e = max(e, abs(D_oracle(-R, J, m′, m) + D_oracle(R, J, m′, m)))
            end
        end
        e
    end
    @test errdc < 1e-13

    # Every entry Varshalovich prints in Tables 4.3-4.12 agrees with his closed form.  The
    # transcription returns `nothing` for the entries the book omits; the counts are
    # asserted so that a transcription which silently became all-`nothing` would fail here
    # rather than quietly gut the table comparison in the item below.
    errtab, nchecked, nomitted = let e = 0.0, nc = 0, no = 0
        for β ∈ βvalues, J ∈ 1//2:1:9//2, M ∈ -J:J, M′ ∈ -J:J
            t = d_table(J, M, M′, big(β))
            if t === nothing
                no += 1
                continue
            end
            nc += 1
            e = max(e, abs(Float64(t) - Float64(d_oracle(J, M, M′, big(β)))))
        end
        (e, nc, no)
    end
    @test nchecked == 144 * length(βvalues)   # 144 distinct transcribed entries
    @test nomitted == 76 * length(βvalues)
    @test errtab < 1e-70   # measured 0.0: both sides are evaluated in BigFloat
    @info "Stage-1 oracle cross-checks" err errdc errtab nchecked nomitted

    # The tables are transcribed only for 1/2 ≤ J ≤ 9/2, and only for half-integers
    @test_throws "Only J = 1/2" d_table(11//2, 1//2, 1//2, 0.3)
    @test_throws "Only half-integer" d_table(2//1, 1//1, 1//1, 0.3)
    @test_throws "must be ≤ J" d_table(3//2, 5//2, 1//2, 0.3)

    # `WignerDElement` refuses the range where its closed form loses accuracy
    @test_throws "maximum supported ℓ" WignerDElement(Rotor{Float64}(1), 17//2, 1//2, 1//2)
end


@testitem "Half-integer d vs the Varshalovich closed form" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: d, dCalculator, recurrence!
    import .HalfIntegerOracle: d_oracle, βvalues

    # `d(β, ℓₘₐₓ)` and `dCalculator` against Varshalovich Eq. 4.3.1(2), evaluated in
    # `BigFloat` (exact factorials, so the oracle itself contributes nothing at Float64
    # resolution).
    #
    # Measured worst case over J ≤ 15/2 and all nine β: 4.4e-16 (2.0 eps) for `d`, and
    # 4.4e-16 (2.0 eps) through the calculator with restricted m′ blocks and Nᵣ = 3.  A
    # wider scan over 39 β reaches 6.7e-16 (3.0 eps).  `40eps()` = 8.9e-15 leaves a factor
    # of ≳ 13 on that wider scan.
    ϵ = 40eps()

    errbyJ = Dict{Rational{Int}, Float64}()
    for J ∈ 1//2:1:15//2
        e = 0.0
        for β ∈ βvalues
            𝔡 = d(β, J)
            for ℓ ∈ 1//2:1:J, m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                e = max(e, abs(𝔡[ℓ][m′, m] - Float64(d_oracle(ℓ, m′, m, big(β)))))
            end
        end
        errbyJ[J] = e
        @test e < ϵ
    end
    @info "Half-integer d vs the Varshalovich closed form" errbyJ

    # The calculator, stepping ℓ by ℓ, with a restricted m′ block and several β at a time
    errcalc = let e = 0.0
        for Jₘₐₓ ∈ (1//2, 7//2, 15//2), m′ₘₐₓ ∈ (1//2, 3//2, Jₘₐₓ)
            m′ₘₐₓ > Jₘₐₓ && continue
            βs = βvalues[1:3]
            calc = dCalculator(βs, Jₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ)
            for ℓ ∈ 1//2:1:Jₘₐₓ
                recurrence!(calc, ℓ)
                block = calc[ℓ]
                for iᵣ ∈ eachindex(βs), m′ ∈ max(-ℓ, -m′ₘₐₓ):min(ℓ, m′ₘₐₓ), m ∈ -ℓ:ℓ
                    e = max(
                        e, abs(block[iᵣ, m′, m] - Float64(d_oracle(ℓ, m′, m, big(βs[iᵣ]))))
                    )
                end
            end
        end
        e
    end
    @test errcalc < ϵ
    @info "Half-integer d through dCalculator" errcalc

    # A batch of β must agree bitwise with the same β values computed one at a time
    βs = βvalues[3:5]
    batch = dCalculator(βs, 9//2)
    recurrence!(batch, 9//2)
    for (iᵣ, β) ∈ enumerate(βs)
        single = dCalculator(β, 9//2)
        recurrence!(single, 9//2)
        @test collect(batch[9//2][iᵣ]) == collect(single[9//2])
    end
end


@testitem "Half-integer 𝔇 vs Boyle (2016)" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: D, DCalculator, recurrence!
    import .HalfIntegerOracle: D_oracle, rotors

    # `D(R, ℓₘₐₓ)` and `DCalculator` against the independent quaternionic reference of
    # Boyle (2016), whose convention is the complex conjugate of this package's.
    #
    # Measured worst case over J ≤ 15/2 and the ten rotors below: 2.7e-15 (12.3 eps) for
    # `D`, and 1.9e-15 (8.3 eps) through the calculator with restricted blocks and Nᵣ = 2.
    # A wider scan over 45 rotors reaches 4.1e-15 (18.4 eps) at J = 13/2, so the tolerance
    # is set from that rather than from the fixed list: `200eps()` = 4.4e-14 leaves a factor
    # of ≳ 11 on the wide scan and ≳ 16 on the list actually used here.  Part of the error
    # is the Float64 reference's own.
    ϵ = 200eps()

    Rs = rotors()

    errbyJ = Dict{Rational{Int}, Float64}()
    for J ∈ 1//2:1:15//2
        e = 0.0
        for R ∈ Rs
            𝔇 = D(R, J)
            for ℓ ∈ 1//2:1:J, m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                e = max(e, abs(𝔇[ℓ][m′, m] - D_oracle(R, ℓ, m′, m)))
            end
        end
        errbyJ[J] = e
        @test e < ϵ
    end
    @info "Half-integer 𝔇 vs Boyle (2016)" errbyJ

    # The calculator, with restricted m′ blocks and batched rotors
    errcalc = let e = 0.0
        for Jₘₐₓ ∈ (1//2, 7//2, 15//2), m′ₘₐₓ ∈ (1//2, 3//2, Jₘₐₓ)
            m′ₘₐₓ > Jₘₐₓ && continue
            Rbatch = Rs[5:6]
            calc = DCalculator(Rbatch, Jₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ)
            for ℓ ∈ 1//2:1:Jₘₐₓ
                recurrence!(calc, ℓ)
                block = calc[ℓ]
                for iᵣ ∈ 1:2, m′ ∈ max(-ℓ, -m′ₘₐₓ):min(ℓ, m′ₘₐₓ), m ∈ -ℓ:ℓ
                    e = max(e, abs(block[iᵣ, m′, m] - D_oracle(Rbatch[iᵣ], ℓ, m′, m)))
                end
            end
        end
        e
    end
    @test errcalc < ϵ
    @info "Half-integer 𝔇 through DCalculator" errcalc

    # Jumping around in ℓ must give exactly what a fresh, sequential calculator gives
    calc = DCalculator(Rs[6], 9//2)
    for ℓ ∈ (9//2, 1//2, 5//2, 7//2, 3//2)
        recurrence!(calc, ℓ)
        fresh = DCalculator(Rs[6], 9//2)
        recurrence!(fresh, ℓ)
        @test collect(calc[ℓ]) == collect(fresh[ℓ])
    end
end


@testitem "Half-integer d vs the Varshalovich tables" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: d
    import .HalfIntegerOracle: d_table, βvalues

    # Varshalovich's Tables 4.3-4.12 print d^J_{MM′} explicitly for J ≤ 9/2, but only for
    # the rows M ≥ 1/2 and, within those, only the entries not obtainable from ones already
    # given; `d_table` returns `nothing` for the rest, so every caller must guard on that.
    # The counts are asserted so that the guard cannot silently swallow the whole table.
    #
    # Measured worst case over all 144 transcribed entries and all nine β: 3.3e-16
    # (1.5 eps).  `40eps()` = 8.9e-15 leaves a factor of ≳ 26.
    ϵ = 40eps()

    err, nchecked, nomitted = let e = 0.0, nc = 0, no = 0
        for β ∈ βvalues
            for J ∈ 1//2:1:9//2
                block = d(β, J)[J]
                for M ∈ -J:J, M′ ∈ -J:J
                    t = d_table(J, M, M′, big(β))
                    if t === nothing
                        no += 1
                        continue
                    end
                    nc += 1
                    e = max(e, abs(block[M, M′] - Float64(t)))
                end
            end
        end
        (e, nc, no)
    end
    @test nchecked == 144 * length(βvalues)
    @test nomitted == 76 * length(βvalues)
    @test err < ϵ
    @info "Half-integer d vs the Varshalovich tables" err nchecked nomitted
end


@testitem "Half-integer metamorphic identities" setup=[HalfIntegerOracle] begin
    import LinearAlgebra: I
    import SphericalFunctions: D, d
    import .HalfIntegerOracle: rotors

    # Oracle-free identities, which stay meaningful far past the J where any closed form is
    # usable.  Each is accumulated to a single worst case per J.
    Js = (1//2, 7//2, 15//2, 21//2, 41//2, 61//2, 101//2)
    Rs = rotors()
    R₁, R₂ = Rs[6], Rs[7]

    # 1. The character identity Σₘ d^J_{mm}(β) = sin((2J+1)β/2) / sin(β/2).  β is kept away
    #    from 0, where the right-hand side approaches its maximum 2J+1 and the absolute
    #    error grows with it.  Measured worst case over β ∈ {0.3, 1.1, 2.0, 2.9} and all J
    #    above: 1.1e-14 at J = 61/2, growing roughly like J·eps.  `2e-13` gives ≳ 18.
    errchar = Dict{Rational{Int}, Float64}()
    for J ∈ Js
        e = 0.0
        for β ∈ (0.3, 1.1, 2.0, 2.9)
            block = d(β, J)[J]
            χ = sum(block[m, m] for m ∈ -J:J)
            e = max(e, abs(χ - sin((2J + 1) * β / 2) / sin(β / 2)))
        end
        errchar[J] = e
        @test e < 2e-13
    end

    # 2. 𝔇(-R) = -𝔇(R), exactly.  The double-cover sign lives entirely in the α and γ
    #    phases, which are formed from the rotor components without a branch, so this is
    #    bit-exact and is asserted as such rather than with a tolerance.
    for J ∈ Js
        @test Matrix(D(R₁, J)[J]) == -Matrix(D(-R₁, J)[J])
    end

    # 3. Unitarity ‖𝔇†𝔇 - I‖∞, and 4. the representation property 𝔇(R₁R₂) = 𝔇(R₁)𝔇(R₂).
    #    Measured worst cases with these two rotors: unitarity 3.1e-15, representation
    #    1.0e-15, both at J = 61/2.  A wider scan (six random rotor pairs per J) reaches
    #    1.3e-14 and 7.7e-15, so `2e-13` leaves a factor of ≳ 15 even there.
    errunit = Dict{Rational{Int}, Float64}()
    errrep = Dict{Rational{Int}, Float64}()
    for J ∈ Js
        A = Matrix(D(R₁, J)[J])
        B = Matrix(D(R₂, J)[J])
        AB = Matrix(D(R₁ * R₂, J)[J])
        errunit[J] = maximum(abs, A' * A - I)
        errrep[J] = maximum(abs, AB - A * B)
        @test errunit[J] < 2e-13
        @test errrep[J] < 2e-13
    end

    @info "Half-integer metamorphic identities" errchar errunit errrep
end


@testitem "Half-integer Float64 vs BigFloat" setup=[HalfIntegerOracle] begin
    import Quaternionic: Rotor
    import SphericalFunctions: d, D
    import .HalfIntegerOracle: rotors

    # The same β, once in Float64 and once in BigFloat, at the three J of memo §5.5.
    # BigFloat's default 256-bit precision makes its result exact at Float64 resolution, so
    # this measures the engine's own error growth with J where no closed form is available.
    #
    # Measured worst cases (absolute; and relative, restricted to |d| > 1e-3, where a
    # relative error is meaningful):
    #     J = 21/2 : 4.4e-16 (2.0 eps), 1.3e-13
    #     J = 41/2 : 8.9e-16 (4.0 eps), 8.0e-14
    #     J = 61/2 : 1.3e-15 (6.0 eps), 9.8e-14
    # `40eps()` = 8.9e-15 and `1e-11` leave factors of ≳ 6.7 and ≳ 78.
    ϵabs = 40eps()
    ϵrel = 1e-11

    errabs = Dict{Rational{Int}, Float64}()
    errrel = Dict{Rational{Int}, Float64}()
    for J ∈ (21//2, 41//2, 61//2)
        ea = 0.0
        er = 0.0
        for β ∈ (0.3, 1.1, 2.0, 2.9)
            f = Matrix(d(β, J)[J])
            b = Matrix(d(big(β), J)[J])
            for i ∈ eachindex(f)
                exact = Float64(b[i])
                δ = abs(f[i] - exact)
                ea = max(ea, δ)
                abs(exact) > 1e-3 && (er = max(er, δ / abs(exact)))
            end
        end
        errabs[J] = ea
        errrel[J] = er
        @test ea < ϵabs
        @test er < ϵrel
    end
    @info "Half-integer Float64 vs BigFloat" errabs errrel

    # The BigFloat path must also be internally consistent: unitarity of 𝔇 at BigFloat
    # precision is far tighter than anything Float64 can reach.  Measured 6.9e-77 at the
    # default 256-bit precision; `1e-70` leaves a factor of 1.4e6.
    Rbig = Rotor{BigFloat}(rotors()[6])
    A = Matrix(D(Rbig, 21//2)[21//2])
    @test maximum(abs, A' * A - one(A)) < 1e-70
end


@testitem "Half-integer β conventions" setup=[HalfIntegerOracle] begin
    import Quaternionic: Rotor, from_euler_angles
    import SphericalFunctions: d, D

    # Half-integer d has period 4π in β, not 2π, so the angle must be tracked all the way
    # into the half-angle pair (cos(β/2), sin(β/2)) rather than through cos β and sin β.
    #
    # Measured worst cases over J ∈ {1/2, 3/2, 7/2} and β ∈ {0.3, 1.1, 2.0, 2.9}:
    #     d(β+2π) + d(β) : 1.1e-15     d(β+4π) - d(β) : 1.6e-15
    # `40eps()` = 8.9e-15 leaves a factor of ≳ 5.7.
    ϵ = 40eps()

    e2π, e4π = let a2 = 0.0, a4 = 0.0
        for J ∈ (1//2, 3//2, 7//2), β ∈ (0.3, 1.1, 2.0, 2.9)
            a = Matrix(d(β, J)[J])
            a2 = max(a2, maximum(abs, Matrix(d(β + 2π, J)[J]) + a))
            a4 = max(a4, maximum(abs, Matrix(d(β + 4π, J)[J]) - a))
        end
        (a2, a4)
    end
    @test e2π < ϵ
    @test e4π < ϵ
    @info "Half-integer β periodicity" e2π e4π

    # An integer ℓ has period 2π, for contrast
    for ℓ ∈ (1, 3), β ∈ (0.3, 1.1)
        @test collect(d(β + 2π, ℓ)[ℓ]) ≈ collect(d(β, ℓ)[ℓ]) atol=ϵ
    end

    # A bare phase e^{iβ} fixes β only modulo 2π, so the branch β ∈ (-π, π] is used.  For β
    # in that branch the phase and the angle give bitwise identical results, and the same
    # phase built from β + 2π gives the same answer again -- i.e. the double-cover sign is
    # lost, exactly as the `dCalculator` docstring says.
    for J ∈ (1//2, 3//2, 7//2), β ∈ (0.3, 1.1, 2.0, 2.9)
        # `cis(β)` for β in the branch: bitwise identical to passing the angle itself
        @test Matrix(d(cis(β), J)[J]) == Matrix(d(β, J)[J])
        # `cis(β + 2π)` is the same point of the circle to within rounding, and gives the
        # same d -- with the opposite sign to `d(β + 2π, J)`, which is the information the
        # bare phase cannot preserve.  Measured difference 4.4e-16; `40eps()` gives ≳ 20.
        @test Matrix(d(cis(β + 2π), J)[J]) ≈ Matrix(d(β, J)[J]) atol=ϵ
        @test Matrix(d(cis(β + 2π), J)[J]) ≈ -Matrix(d(β + 2π, J)[J]) atol=ϵ
    end

    # A `Rotor` restricts β to [0, π] only, so d sees the same β for R and -R; the
    # double-cover sign of a half-integer 𝔇 lives entirely in the α and γ phases.
    for J ∈ (1//2, 5//2), (α, β, γ) ∈ ((0.7, 1.1, 2.3), (2.9, 0.4, 5.1))
        R = Rotor(from_euler_angles(α, β, γ))
        @test Matrix(d(R, J)[J]) == Matrix(d(-R, J)[J])
        @test Matrix(d(R, J)[J]) ≈ Matrix(d(β, J)[J]) atol=ϵ
        @test Matrix(D(-R, J)[J]) == -Matrix(D(R, J)[J])
    end
end


@testitem "Half-integer index validation and error messages" begin
    import Quaternionic: Rotor
    import SphericalFunctions: D, d, DCalculator, dCalculator,
        HCalculator, WignerMatrix, WignerDMatrix, recurrence!,
        HalfOddInteger, half_integer

    # Half-integer indices may be spelled as `Rational`s with denominator exactly 2, which
    # the public entry points convert to `HalfOddInteger`.  Everything else must fail loudly,
    # and in particular must never be silently floored to a neighbouring index.
    𝟙 = Rotor{Float64}(1)

    # ℓₘₐₓ must be a half-integer, not an integer-valued Rational and not a Float
    @test_throws "must have denominator 2" D(𝟙, 3//1)
    @test_throws "must have denominator 2" d(1.1, 4//2)
    @test_throws "must have denominator 2" DCalculator(𝟙, 7//3)
    @test_throws MethodError D(𝟙, 3.5)
    @test_throws "must be non-negative" DCalculator(𝟙, -1//2)

    # The four block limits must be half-integers of the same type as ℓₘₐₓ
    @test_throws TypeError DCalculator(𝟙, 7//2; m′ₘₐₓ=1)
    @test_throws "must have denominator 2" DCalculator(𝟙, 7//2; m′ₘₐₓ=2//1, m′ₘᵢₙ=-2//1)

    # Both rows m′ = ±1/2 are needed to seed the half-integer ladder, so the m′ and m
    # windows must bracket ±ℓₘᵢₙ
    @test_throws "too large for this index type" DCalculator(𝟙, 7//2; m′ₘₐₓ=3//2, m′ₘᵢₙ=1//2)
    @test_throws "too small for this index type" DCalculator(𝟙, 7//2; m′ₘₐₓ=-1//2, m′ₘᵢₙ=-3//2)
    @test_throws "too large for this index type" DCalculator(𝟙, 7//2; mₘₐₓ=7//2, mₘᵢₙ=1//2)
    @test_throws "is too large for ℓₘₐₓ" DCalculator(𝟙, 7//2; m′ₘₐₓ=9//2)

    # ...but a legal narrow window is fine, including the narrowest one
    @test DCalculator(𝟙, 7//2; m′ₘₐₓ=1//2, m′ₘᵢₙ=-1//2) isa DCalculator
    @test dCalculator(1.1, 1//2) isa dCalculator
    @test HCalculator(1.1, 1//2; m′ₘₐₓ=1//2) isa HCalculator

    # `recurrence!` and `calc[ℓ]` reject the wrong parity of ℓ, and ℓ out of range
    calc = DCalculator(𝟙, 5//2)
    @test_throws InexactError recurrence!(calc, 𝟙, 2)
    @test_throws "out of bounds" recurrence!(calc, 𝟙, 7//2)
    recurrence!(calc, 𝟙, 5//2)
    @test_throws "out of bounds" calc[9//2]
    @test_throws "not ℓ=" calc[3//2]

    # An integer ℓ on a half-integer `WignerSeries` must say so, rather than throwing a bare
    # `InexactError` out of the index arithmetic (or, under `@inbounds`, quietly returning a
    # neighbouring block)
    𝔇 = D(𝟙, 7//2)
    @test_throws "is not one of the ℓ values" 𝔇[2]
    @test_throws "is not one of the ℓ values" 𝔇[3//1]
    @test_throws BoundsError 𝔇[9//2]
    @test 𝔇[5//2] === 𝔇[5//2]

    # Index-type mismatches on a block are a `MethodError`, not a silent conversion
    block = 𝔇[7//2]
    @test_throws MethodError block[1, 1]
    @test_throws BoundsError block[9//2, 1//2]
    @test_throws BoundsError block[1//2, 9//2]

    # The containers apply the same rules as the calculators
    @test_throws "must have denominator 2" WignerMatrix(zeros(ComplexF64, 3, 3), 1//1)
    @test_throws TypeError WignerMatrix(zeros(ComplexF64, 3, 3), half_integer(1//2); m′ₘₐₓ=1)
    @test_throws "must have denominator 2" WignerDMatrix(ComplexF64, 5//3)
    @test WignerDMatrix(ComplexF64, 5//2) isa WignerMatrix{HalfOddInteger}
end


@testitem "Half-integer containers" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: D, d, DCalculator, recurrence!, sYlmCalculator,
        WignerMatrix, WignerMatrixBatch, DegreeBlock, DegreeBlockBatch, WignerSeries,
        WignerDMatrix, WignerdMatrix, WignerRange, HalfOddInteger,
        SpinMatrix, SpinMatrixBatch,
        ℓ, ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, sₘₐₓ, sₘᵢₙ
    import .HalfIntegerOracle: rotors

    R = rotors()[6]
    J = 5//2

    @testset "WignerSeries" begin
        𝔇 = D(R, J)
        @test 𝔇 isa WignerSeries
        @test ℓₘᵢₙ(𝔇) == 1//2 && ℓₘₐₓ(𝔇) == J
        @test length(𝔇) == 3
        @test size(𝔇) == (3,) && size(𝔇, 1) == 3 && size(𝔇, 2) == 1
        @test axes(𝔇) == (1//2:J,)
        @test keys(𝔇) == 1//2:1:J
        @test firstindex(𝔇) == 1//2 && lastindex(𝔇) == J
        @test collect(𝔇) == [𝔇[ℓ] for ℓ ∈ 1//2:1:J]   # iteration
        @test eltype(𝔇) <: WignerDMatrix

        # copy is deep: mutating the copy must not touch the original
        𝔇c = copy(𝔇)
        @test 𝔇c == 𝔇
        @test 𝔇c[J] !== 𝔇[J]
        𝔇c[J][1//2, 1//2] += 1
        @test 𝔇c != 𝔇
        @test D(R, J) == 𝔇

        # similar keeps the shape and the ℓ range, with fresh storage
        𝔇s = similar(𝔇)
        @test 𝔇s isa WignerSeries
        @test length(𝔇s) == length(𝔇) && axes(𝔇s) == axes(𝔇)
        @test axes(𝔇s[J]) == axes(𝔇[J])

        # d gives the real sibling
        𝔡 = d(1.1, J)
        @test 𝔡 isa WignerSeries
        @test eltype(𝔡) <: WignerdMatrix
        @test occursin("WignerSeries", sprint(show, 𝔇))
        @test occursin("ℓ ∈ 1//2:5//2", sprint(show, 𝔇))
        @test occursin("ℓ = 5//2", sprint(show, MIME("text/plain"), 𝔇))
    end

    @testset "WignerMatrix" begin
        w = D(R, J)[J]
        @test w isa WignerMatrix && w isa WignerDMatrix
        @test !(w isa AbstractMatrix)   # half-integer axes cannot satisfy that interface
        @test ndims(w) == 2
        @test ℓ(w) == J && ℓₘᵢₙ(w) == 1//2
        @test m′ₘₐₓ(w) == J && m′ₘᵢₙ(w) == -J && mₘₐₓ(w) == J && mₘᵢₙ(w) == -J
        @test axes(w) == (-J:J, -J:J)
        @test axes(w, 1) == -J:J && axes(w, 3) == Base.OneTo(1)
        @test size(w) == (6, 6) && size(w, 2) == 6 && size(w, 3) == 1
        @test length(w) == 36
        @test eltype(w) === ComplexF64

        # getindex/setindex! use the natural indices; the 1-based parent is the same data
        @test w[-J, -J] == parent(w)[1, 1]
        @test w[J, J] == parent(w)[6, 6]
        @test w[1//2, -3//2] == parent(w)[4, 2]

        # collect / Matrix / Array all give a plain 1-based matrix in increasing (m′, m)
        M = Matrix(w)
        @test M isa Matrix{ComplexF64} && size(M) == (6, 6)
        @test collect(w) == M && Array(w) == M
        @test M == [w[m′, m] for m′ ∈ -J:J, m ∈ -J:J]
        @test collect(Iterators.take(w, length(w))) == vec(M)   # iteration order

        # copy is independent of the original and keeps the indices
        wc = copy(w)
        @test wc == w && axes(wc) == axes(w) && ℓ(wc) == ℓ(w)
        wc[1//2, 1//2] += 1
        @test wc != w

        # similar has the same indices and fresh, block-sized, plain storage
        ws = similar(w)
        @test ws isa WignerMatrix
        @test axes(ws) == axes(w) && ℓ(ws) == ℓ(w) && eltype(ws) === eltype(w)
        @test parent(ws) isa Matrix{ComplexF64} && size(parent(ws)) == size(w)
        @test eltype(similar(w, Float64)) === Float64
        fill!(parent(ws), 0)
        for m′ ∈ -J:J, m ∈ -J:J     # the storage must be writable through the natural indices
            ws[m′, m] = w[m′, m]
        end
        @test ws == w

        @test occursin("(-5//2:5//2)×(-5//2:5//2)", sprint(show, w))
        @test occursin("for ℓ=5//2", sprint(show, w))
        @test occursin("for ℓ=5//2", sprint(show, MIME("text/plain"), w))
    end

    @testset "WignerMatrixBatch" begin
        Rs = rotors()[6:8]
        calc = DCalculator(Rs, J)
        recurrence!(calc, J)
        b = calc[J]
        @test b isa WignerMatrixBatch
        @test ndims(b) == 3
        @test size(b) == (3, 6, 6) && size(b, 1) == 3 && size(b, 4) == 1
        @test length(b) == 108
        @test axes(b) == (1:3, -J:J, -J:J)
        @test axes(b, 2) == -J:J && axes(b, 4) == Base.OneTo(1)

        # `b[iᵣ]` is the WignerMatrix view of one rotor, and agrees with a single-rotor run
        for iᵣ ∈ 1:3
            @test b[iᵣ] isa WignerMatrix
            @test collect(b[iᵣ]) == collect(D(Rs[iᵣ], J)[J])
        end
        @test all(b[iᵣ, m′, m] == b[iᵣ][m′, m] for iᵣ ∈ 1:3, m′ ∈ -J:J, m ∈ -J:J)

        A = Array(b)
        @test A isa Array{ComplexF64, 3} && size(A) == (3, 6, 6)
        @test collect(b) == A
        @test_throws "3-dimensional" Matrix(b)

        # iteration follows `Array(b)`, so the reducers work as they do on the integer path
        @test collect(Iterators.take(b, length(b))) == vec(A)
        # unitary blocks, so Σ|𝔇|² = Nᵣ(2J+1); measured residual 3.6e-15
        @test sum(abs2, b) ≈ 3 * (2J + 1) atol=1e-13

        bc = copy(b)
        @test bc == b && axes(bc) == axes(b)
        bc[1, 1//2, 1//2] += 1
        @test bc != b

        bs = similar(b)
        @test bs isa WignerMatrixBatch
        @test axes(bs) == axes(b) && eltype(bs) === eltype(b)
        @test parent(bs) isa Array{ComplexF64, 3} && size(parent(bs)) == size(b)

        @test occursin("(1:3)×(-5//2:5//2)×(-5//2:5//2)", sprint(show, b))
        @test occursin("WignerMatrixBatch", sprint(show, MIME("text/plain"), b))
    end

    @testset "DegreeBlock and DegreeBlockBatch" begin
        # `sYlmCalculator` is the half-integer producer of the 1-dimensional containers
        calc = sYlmCalculator(R, J, -3//2:3//2)
        recurrence!(calc, J)
        v = calc[J, 1//2]
        @test v isa DegreeBlock
        @test ndims(v) == 1
        @test size(v) == (6,) && size(v, 1) == 6 && size(v, 2) == 1
        @test length(v) == 6
        @test axes(v) == (-J:J,) && axes(v, 1) == -J:J && axes(v, 2) == Base.OneTo(1)
        @test firstindex(v) == -J && lastindex(v) == J && keys(v) == -J:1:J
        @test v[-J] == parent(v)[1] && v[J] == parent(v)[6]
        @test collect(v) == Vector(v) == Array(v) == [v[m] for m ∈ -J:J]
        @test collect(Iterators.take(v, length(v))) == collect(v)

        vc = copy(v)
        @test vc == v
        vc[1//2] += 1
        @test vc != v

        vs = similar(v)
        @test vs isa DegreeBlock
        @test axes(vs) == axes(v) && eltype(vs) === eltype(v)
        @test parent(vs) isa Vector{ComplexF64} && length(parent(vs)) == length(v)
        @test occursin("(-5//2:5//2) DegreeBlock", sprint(show, v))

        # The axis itself is a `WignerRange`, and it has the unit step and the length that
        # `Base` would otherwise try to build from `oneunit` and `zero`, both of which are
        # deliberately absent for `HalfOddInteger`; so it can be shown, measured and collected
        r = axes(v, 1)
        @test r isa WignerRange{HalfOddInteger}
        @test step(r) === 1
        @test length(r) == 6 && lastindex(r) == 6
        @test collect(r) == [h for h ∈ -J:J]
        @test sprint(show, r) == "-5//2:1:5//2"
        @test sprint(show, MIME("text/plain"), r) == "-5//2:1:5//2"
        @test sprint(show, axes(v)) == "(-5//2:1:5//2,)"
        # ... while the integer axis is exactly what `Base` gives it
        ri = WignerRange(-2:2)
        @test step(ri) === 1 && length(ri) == 5
        @test sprint(show, ri) == "-2:1:2"
        @test sprint(show, (ri, ri)) == "(-2:1:2, -2:1:2)"
        # Both are indexed by position, so their own axes are 1-based, and broadcasting a
        # function over one maps its values rather than mis-reading them as positions
        @test axes(ri) == (Base.OneTo(5),) && axes(r) == (Base.OneTo(6),)
        @test Rational.(ri) == [-2//1, -1//1, 0//1, 1//1, 2//1]
        @test Rational.(r) == [k//2 for k ∈ -5:2:5]

        calcb = sYlmCalculator([R, R], J, -3//2:3//2)
        recurrence!(calcb, J)
        vb = calcb[J, 1//2]
        @test vb isa DegreeBlockBatch
        @test ndims(vb) == 2
        @test size(vb) == (2, 6) && size(vb, 1) == 2 && size(vb, 3) == 1
        @test axes(vb) == (1:2, -J:J) && axes(vb, 3) == Base.OneTo(1)
        @test vb[1] isa DegreeBlock
        @test collect(vb[1]) == collect(v)
        @test collect(vb) == Matrix(vb) == [vb[iᵣ, m] for iᵣ ∈ 1:2, m ∈ -J:J]
        @test collect(Iterators.take(vb, length(vb))) == vec(Matrix(vb))

        vbc = copy(vb)
        @test vbc == vb
        vbc[1, 1//2] += 1
        @test vbc != vb

        vbs = similar(vb)
        @test vbs isa DegreeBlockBatch
        @test axes(vbs) == axes(vb) && eltype(vbs) === eltype(vb)
        @test parent(vbs) isa Matrix{ComplexF64} && size(parent(vbs)) == size(vb)
        @test occursin("(1:2)×(-5//2:5//2) DegreeBlockBatch", sprint(show, vb))
    end

    @testset "SpinMatrix and SpinMatrixBatch" begin
        # A calculator built for a range of spin weights is the half-integer producer of the
        # 2- and 3-dimensional (s, m) containers, as one built for a single spin weight is of
        # the 1- and 2-dimensional ones above
        sr = -3//2:3//2
        calc = sYlmCalculator(R, J, sr)
        recurrence!(calc, J)
        b = calc[J]
        @test b isa SpinMatrix
        @test ndims(b) == 2
        @test size(b) == (4, 6) && size(b, 1) == 4 && size(b, 3) == 1
        @test length(b) == 24
        @test axes(b) == (sr, -J:J) && axes(b, 1) == sr && axes(b, 3) == Base.OneTo(1)
        @test sₘₐₓ(b) == 3//2 && sₘᵢₙ(b) == -3//2 && mₘₐₓ(b) == J && mₘᵢₙ(b) == -J
        @test ℓ(b) == J && ℓₘᵢₙ(b) == 1//2 && keys(b) == sr
        @test b[1//2, -J] == parent(b)[3, 1]
        # Each row is the very block a single-spin calculator would give
        for s ∈ sr
            @test b[s, :] isa DegreeBlock
            @test collect(b[s, :]) == collect(calc[J, s])
            @test collect(b[s, :]) == collect(recurrence!(sYlmCalculator(R, J, s), J)[J])
        end
        @test collect(b) == Matrix(b) == Array(b) == [b[s, m] for s ∈ sr, m ∈ -J:J]
        @test collect(Iterators.take(b, length(b))) == vec(Matrix(b))

        bc = copy(b)
        @test bc == b
        bc[1//2, 1//2] += 1
        @test bc != b

        bs = similar(b)
        @test bs isa SpinMatrix
        @test axes(bs) == axes(b) && eltype(bs) === eltype(b)
        @test parent(bs) isa Matrix{ComplexF64} && size(parent(bs)) == size(b)
        @test occursin("(-3//2:3//2)×(-5//2:5//2) SpinMatrix", sprint(show, b))
        @test occursin("SpinMatrix", sprint(show, MIME("text/plain"), b))

        calcb = sYlmCalculator([R, R], J, sr)
        recurrence!(calcb, J)
        bb = calcb[J]
        @test bb isa SpinMatrixBatch
        @test ndims(bb) == 3
        @test size(bb) == (2, 4, 6) && size(bb, 1) == 2 && size(bb, 4) == 1
        @test axes(bb) == (1:2, sr, -J:J) && axes(bb, 4) == Base.OneTo(1)
        @test sₘₐₓ(bb) == 3//2 && sₘᵢₙ(bb) == -3//2 && ℓ(bb) == J
        @test bb[1] isa SpinMatrix
        @test collect(bb[1]) == collect(b)
        @test bb[:, 1//2, :] isa DegreeBlockBatch
        @test collect(bb[:, 1//2, :]) == collect(calcb[J, 1//2])
        @test collect(bb) == Array(bb) == [bb[iᵣ, s, m] for iᵣ ∈ 1:2, s ∈ sr, m ∈ -J:J]
        @test collect(Iterators.take(bb, length(bb))) == vec(Array(bb))

        bbc = copy(bb)
        @test bbc == bb
        bbc[1, 1//2, 1//2] += 1
        @test bbc != bb

        bbs = similar(bb)
        @test bbs isa SpinMatrixBatch
        @test axes(bbs) == axes(bb) && eltype(bbs) === eltype(bb)
        @test parent(bbs) isa Array{ComplexF64, 3} && size(parent(bbs)) == size(bb)
        @test occursin("(1:2)×(-3//2:3//2)×(-5//2:5//2) SpinMatrixBatch", sprint(show, bb))
        @test occursin("SpinMatrixBatch", sprint(show, MIME("text/plain"), bb))
    end
end
