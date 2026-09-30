# Tests of the parts of `src/mode_weights/operators.jl` that compute no numbers.
#
# `test/operators.jl` checks what the operators compute, with the help of the explicit
# derivatives in `test/utilities/explicit_operators.jl`.  The items here check what each
# operator *is*: a singleton value with its own name and display, its shift of the spin weight,
# the band of the matrix it occupies, its ASCII alias, and its documentation.  Those are the
# small dispatch tables at the top of the file, and they are what a user sees when an operator
# turns up in an error message or at the REPL.

@testitem "Differential operators: names, display and spin shift" begin
    import SphericalFunctions: DifferentialOperator, Δspin, bandstructure, coefftype
    import SphericalFunctions: DiagonalBand, SubdiagonalBand, SuperdiagonalBand, TridiagonalBand

    # Every exported operator is a value, not a function, and knows its own name
    operators = (L², R², Lz, L₊, L₋, Lx, Ly, Rz, R₊, R₋, ð, ð̄)
    names     = (:L², :R², :Lz, :L₊, :L₋, :Lx, :Ly, :Rz, :R₊, :R₋, :ð, :ð̄)

    for (op, nm) ∈ zip(operators, names)
        @test op isa DifferentialOperator
        # Zero-size singletons, so dispatching on one costs nothing and every trait folds away
        @test Base.issingletontype(typeof(op))
        @test sizeof(op) == 0
        @test nameof(op) == nm
        # ... bound to the name it reports
        @test getfield(SphericalFunctions, nm) === op
        # `show` prints the name rather than the internal struct, so that an operator in an
        # error message reads as `ð` and not as `SphericalFunctions.SpinRaising()`
        @test sprint(show, op) == string(nm)
        @test repr(op) == string(nm)
        @test !occursin("SphericalFunctions", sprint(show, op))
    end

    # All twelve names are distinct, so none of the `nameof` methods shadows another, and no
    # two operators share a type, so no two can share a trait by accident
    @test length(unique(nameof.(operators))) == length(operators)
    @test length(unique(typeof.(operators))) == length(operators)

    # The spin weight is changed only by the right-handed ladder and the eth operators
    for op ∈ (L², R², Lz, L₊, L₋, Lx, Ly, Rz)
        @test Δspin(op) == 0
    end
    @test Δspin(R₊) == 1
    @test Δspin(ð) == 1
    @test Δspin(R₋) == -1
    @test Δspin(ð̄) == -1

    # Which band each operator occupies decides which builder and kernel apply
    for op ∈ (L², R², Lz, Rz, R₊, R₋, ð, ð̄)
        @test bandstructure(op) isa DiagonalBand
    end
    @test bandstructure(L₊) isa SubdiagonalBand
    @test bandstructure(L₋) isa SuperdiagonalBand
    @test bandstructure(Lx) isa TridiagonalBand
    @test bandstructure(Ly) isa TridiagonalBand

    # `Ly` is the only operator whose matrix entries are complex
    for op ∈ (L², R², Lz, L₊, L₋, Lx, Rz, R₊, R₋, ð, ð̄)
        @test coefftype(op, Float64) == Float64
        @test coefftype(op, Float32) == Float32
    end
    @test coefftype(Ly, Float64) == ComplexF64
    @test coefftype(Ly, Float32) == ComplexF32
end

@testitem "Differential operators: the coefficients vanish below their support" begin
    import SphericalFunctions: support_ℓ, diagonal_coefficient
    import SphericalFunctions: subdiagonal_coefficient, superdiagonal_coefficient
    import SphericalFunctions: Casimir, LeftZ, RightZ, RightRaising, RightLowering
    import SphericalFunctions: SpinRaising, SpinLowering, LeftRaising, LeftLowering, LeftX

    # A spin-weighted function has no modes with ℓ < |s|, and the spin-raising operators
    # additionally have none below |s+1|, because their support is set by the *output* spin.
    @test support_ℓ(Casimir(), 2) == 2
    @test support_ℓ(Casimir(), -2) == 2
    @test support_ℓ(RightRaising(), 2) == 3
    @test support_ℓ(SpinRaising(), 2) == 3
    @test support_ℓ(RightRaising(), -2) == 2
    @test support_ℓ(SpinRaising(), -3) == 3

    # Below the support every coefficient is exactly zero, and of the requested type
    for T ∈ (Float64, Float32)
        z = zero(T)
        @test diagonal_coefficient(Casimir(), T, 3, 2, 0) === z
        @test diagonal_coefficient(LeftZ(), T, 3, 2, 0) === z
        @test diagonal_coefficient(RightZ(), T, 3, 2, 0) === z
        @test diagonal_coefficient(RightRaising(), T, 3, 2, 0) === z
        @test diagonal_coefficient(RightLowering(), T, 3, 2, 0) === z
        @test diagonal_coefficient(SpinRaising(), T, 3, 2, 0) === z
        @test subdiagonal_coefficient(LeftRaising(), T, 3, 2, 0) === z
        @test superdiagonal_coefficient(LeftLowering(), T, 3, 2, 0) === z
        @test subdiagonal_coefficient(LeftX(), T, 3, 2, 0) === z

        # The spin-lowering coefficient negates the masked value, so its vanishing entries
        # are `-0.0` and even `isequal` agrees with what the matrix holds
        @test diagonal_coefficient(SpinLowering(), T, 3, 2, 0) === -z
        @test isequal(diagonal_coefficient(SpinLowering(), T, 3, 2, 0), -z)

        # Above the support they are the familiar square roots, and the ladder coefficients
        # vanish exactly at the edge of each ℓ block, which is what stops ℓ blocks coupling
        @test diagonal_coefficient(Casimir(), T, 0, 3, 0) == T(12)
        @test diagonal_coefficient(LeftZ(), T, 0, 3, 2) == T(2)
        @test diagonal_coefficient(RightZ(), T, 1, 3, 2) == T(1)
        @test diagonal_coefficient(RightRaising(), T, 1, 3, 0) == √T((3-1)*(3+1+1))
        @test diagonal_coefficient(RightLowering(), T, 1, 3, 0) == √T((3+1)*(3-1+1))
        @test subdiagonal_coefficient(LeftRaising(), T, 0, 3, -3) == z
        @test superdiagonal_coefficient(LeftLowering(), T, 0, 3, 3) == z
        @test subdiagonal_coefficient(LeftX(), T, 0, 3, 2) ==
            subdiagonal_coefficient(LeftRaising(), T, 0, 3, 2) / 2

        # `ð` and `R₊` share a coefficient, as do `ð̄` and `-R₋`
        @test diagonal_coefficient(SpinRaising(), T, 1, 3, 0) ==
            diagonal_coefficient(RightRaising(), T, 1, 3, 0)
        @test diagonal_coefficient(SpinLowering(), T, 1, 3, 0) ==
            -diagonal_coefficient(RightLowering(), T, 1, 3, 0)
    end
end

@testitem "Differential operators: ASCII aliases and documentation" begin
    import SphericalFunctions
    import SphericalFunctions: DifferentialOperator, Δspin, Deltaspin
    import SphericalFunctions: L2, Lplus, Lminus, R2, Rplus, Rminus, eth, ethbar

    # Each alias is the operator itself, so it shares its identity, its name and its display
    for (alias, op) ∈ (
        (L2, L²), (Lplus, L₊), (Lminus, L₋), (R2, R²), (Rplus, R₊), (Rminus, R₋),
        (eth, ð), (ethbar, ð̄),
    )
        @test alias === op
        @test repr(alias) == repr(op)
        @test alias(1, 3) == op(1, 3)
    end
    @test Deltaspin === Δspin
    @test Deltaspin(ethbar) == -1

    # The aliases are not exported, since names such as `L2` are generic enough to clash with
    # a user's own
    for name ∈ (:L2, :Lplus, :Lminus, :R2, :Rplus, :Rminus, :eth, :ethbar, :Deltaspin)
        @test !Base.isexported(SphericalFunctions, name)
    end

    # The abstract type, `Δspin`, and each alias have docstrings of their own.  An alias of a
    # value, unlike one of a function, does not lead to the value's docstring, so each alias
    # needs its own for `?eth` to say anything.  (The registry is read directly, because
    # `Base.Docs.hasdoc` exists only from Julia 1.11 on.)
    documented = Base.Docs.meta(SphericalFunctions)
    for name ∈ (:DifferentialOperator, :Δspin, :L2, :Lplus, :Lminus, :R2, :Rplus, :Rminus, :eth, :ethbar)
        @test haskey(documented, Base.Docs.Binding(SphericalFunctions, name))
    end
    @test occursin("ASCII alias of the operator", string(@doc eth))
    @test occursin("`ethbar` may be used in place of `ð̄`", string(@doc ð̄))
    @test occursin("Deltaspin", string(@doc Δspin))
    @test occursin("not an extension point", string(@doc DifferentialOperator))
    # ... and the note shared by the operators' docstrings says what `T` is and how an operator
    # applies to mode weights
    @test occursin("real floating-point type", string(@doc L²))
    @test occursin("without building the matrix", string(@doc ð))
end
