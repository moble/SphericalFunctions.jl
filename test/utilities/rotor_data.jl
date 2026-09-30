# Tests of the functions in `src/calculators/rotors.jl` that read the element type and the
# number of rotors from rotor data.

@testitem "Utilities: the element type and the number of rotors of rotor data" begin
    import SphericalFunctions: floattype, nrotors, check_rotor_type, DCalculator,
        dCalculator, HCalculator, sYlmCalculator, sλlmCalculator
    using Quaternionic: Rotor, Quaternion, QuatVec
    using DoubleFloats: Double64

    # The element type is the `float` of the component type of the data, for a single rotor,
    # angle or phase, and for a vector of any one of those, and it is the same for the data
    # and for its type
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        R = Rotor{T}(1, 0, 0, 0)
        for data ∈ (R, [R, R], T(0.3), T[0.3, 0.4], cis(T(0.3)), [cis(T(0.3))])
            @test floattype(data) === T
            @test floattype(typeof(data)) === T
        end
    end
    @test floattype(3) === Float64
    @test floattype([1, 2]) === Float64
    @test floattype(1 + 0im) === Float64
    @test floattype(Int) === Float64
    @test floattype(Rotor{Int}) === Float64
    @test floattype(Rotor{Int}(1, 0, 0, 0)) === Float64

    # A calculator's is the type it works in, for the calculator and for its type
    for (calc, T) ∈ (
        (DCalculator(Rotor{Float32}(1, 0, 0, 0), 2), Float32), (dCalculator(big(0.3), 2), BigFloat),
        (HCalculator([0.3, 0.4], 2), Float64), (sλlmCalculator(0.3f0, 2, 0), Float32),
        (sYlmCalculator(Rotor{Double64}(1, 0, 0, 0), 2, 0), Double64),
    )
        @test floattype(calc) === floattype(typeof(calc)) === T
    end

    # Data whose component type is abstract does not say what type to work in, and is
    # refused rather than answered with a guess, which would be `float(Real) === Float64`
    for data ∈ (
        Complex{Real}(1, 2.0), Complex{Real}[1 + 0im, 2.0 + 0im], Complex{AbstractFloat}[1.0 + 0im],
        Rotor{Real}(1.0, 0, 0, 0), Rotor{Real}[Rotor(1.0)], Real[1.0, 2.0],
        Union{Float64, Float32}[1.0, 2.0f0],
    )
        @test_throws ArgumentError floattype(data)
        @test_throws "which does not say what floating-point type to work in" floattype(data)
    end
    @test_throws "has components of type Real" floattype(Rotor{Real}[Rotor(1.0)])
    @test_throws ArgumentError DCalculator(Rotor{Real}[Rotor(1.0)], 2)

    # Data that is not rotor data at all is refused with the forms that are accepted
    for data ∈ (Any[1.0, 2.0], Rotor[Rotor(1.0)], "β", [1.0 2.0; 3.0 4.0], nothing)
        @test_throws ArgumentError floattype(data)
        @test_throws "Cannot build a calculator from rotor data of type" floattype(data)
        @test_throws "the accepted forms are a Rotor" floattype(data)
    end
    # ... a `Quaternion` is the rotation of its normalization ...
    @test floattype(Quaternion(1.0, 0, 0, 0)) === Float64
    @test floattype([Quaternion(1.0f0, 0, 0, 0)]) === Float32
    # ... and a `QuatVec` is told how to make one
    for data ∈ (QuatVec(0.0, 1, 0, 0), [QuatVec(0.0, 1, 0, 0)])
        @test_throws ArgumentError floattype(data)
        @test_throws "Rotations are taken as `Rotor`s" floattype(data)
    end

    # The number of rotors is one for a single one, and the length of a non-empty vector
    @test nrotors(Rotor(1.0)) == nrotors(0.3) == nrotors(cis(0.3)) == 1
    @test nrotors([0.3]) == 1
    @test nrotors(fill(Rotor(1.0), 5)) == 5
    @test_throws ArgumentError nrotors(Float64[])
    @test_throws "needs at least one rotor, but got an empty Vector{Float64}" nrotors(Float64[])
    @test_throws ArgumentError nrotors("β")
    @test_throws "Cannot build a calculator from rotor data of type String" nrotors("β")

    # A calculator's element type is fixed by the data it was built from
    calc = DCalculator(Rotor(1.0), 2)
    @test check_rotor_type(calc, Rotor(0.0, 1.0, 0.0, 0.0)) === nothing
    @test check_rotor_type(calc, [0.3, 0.4]) === nothing
    @test_throws ArgumentError check_rotor_type(calc, Rotor{Float32}(1, 0, 0, 0))
    @test_throws "This calculator works in Float64, but the given data would give Float32" check_rotor_type(calc, Float32(0.3))

    # The element type is known from the type of the data alone
    for data ∈ (Rotor(1.0), [Rotor(1.0f0)], 0.3, Double64[0.3], cis(0.3), [cis(big(0.3))])
        @test Base.return_types(floattype, (typeof(data),)) == [Type{floattype(data)}]
        @test Base.return_types(floattype, (Type{typeof(data)},)) == [Type{floattype(data)}]
    end
end
