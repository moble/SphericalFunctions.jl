# Tests of the parts of `src/utilities/operators.jl` that carry no numbers.
#
# `test/operators.jl` and `test/utilities/explicit_operators.jl` check what the operators
# compute.  What neither checks is that each operator can say what it *is*: its name, how it
# prints, how it shifts the spin weight, and which band of the matrix it occupies.  Those are
# the small dispatch tables at the top of the file, and they are what a user sees when an
# operator turns up in an error message or at the REPL.

@testitem "Differential operators: names, display and spin shift" begin
    import SphericalFunctions: DifferentialOperator, Δspin, bandstructure, coefftype
    import SphericalFunctions: DiagonalBand, SubdiagonalBand, SuperdiagonalBand, TridiagonalBand

    # Every exported operator is a value, not a function, and knows its own name
    operators = (L², R², Lz, L₊, L₋, Lx, Ly, Rz, R₊, R₋, ð, ð̄)
    names     = (:L², :R², :Lz, :L₊, :L₋, :Lx, :Ly, :Rz, :R₊, :R₋, :ð, :ð̄)

    for (op, nm) ∈ zip(operators, names)
        @test op isa DifferentialOperator
        @test nameof(op) == nm
        # `show` prints the name rather than the internal struct, so that an operator in an
        # error message reads as `ð` and not as `SphericalFunctions.SpinRaising()`
        @test sprint(show, op) == string(nm)
        @test !occursin("SphericalFunctions", sprint(show, op))
    end

    # All twelve names are distinct, so none of the `nameof` methods shadows another
    @test length(unique(nameof.(operators))) == length(operators)

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
