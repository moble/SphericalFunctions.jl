# Tests of `@index_methods`, the macro that writes the boundary between the indices a caller
# writes and the `Int` or `HalfOddInteger` indices the package computes with.
#
# The functions are defined in a test module rather than in the items, because the items run
# inside a `@testset`, where neither a `struct` nor a docstring can be defined, and because
# `Test.detect_ambiguities` wants a module to search.  Each definition exercises one shape of
# signature that the package's own boundary methods use.

@testmodule IndexMethodsExamples begin
    using SphericalFunctions: @index_methods, IndexType, IndexRange, IndexOrRange,
        HalfOddInteger

    # Two index arguments, one of them with a positional default that refers to the other,
    # and the ASCII-alias pattern for a keyword index, of which only the Unicode name is used.
    """
        pair(ℓ, m=ℓ; mp_max=ℓ, m′ₘₐₓ=mp_max)

    A docstring written above the macro call.
    """
    @index_methods function pair(
        ℓ::IndexType, m::IndexType=ℓ; mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max
    )
        (ℓ, m, m′ₘₐₓ)
    end

    # A type variable bounded by a marker, used in the body, and a positional element type
    # with a default, as the pixelizations and operators take it.
    @index_methods typed(s::IT, ℓ::IT, ::Type{T}=Float64) where {IT<:IndexType, T} =
        (IT, s, ℓ, T)

    # A spin weight that may be one index or a range of them, beside a non-index argument, and
    # a keyword default computed from the converted positional arguments.
    @index_methods function spins(
        R, ℓₘₐₓ::IndexType, s::IndexOrRange; ℓₘᵢₙ::IndexType=lowest(s)
    )
        (R, ℓₘₐₓ, s, ℓₘᵢₙ)
    end
    lowest(s::Union{Int, HalfOddInteger}) = abs(s)
    lowest(s::AbstractUnitRange) =
        first(s) ≤ 0 ≤ last(s) ? zero_of(first(s)) : min(abs(first(s)), abs(last(s)))
    zero_of(::Int) = 0
    zero_of(::HalfOddInteger) = HalfOddInteger(1//2)

    # An argument that must be a range.
    @index_methods only_range(r::IndexRange) = r

    # Integer indices only, with a positional default, as the equiangular grids.
    @index_methods integer_only function grid(
        s::IndexType, ℓₘₐₓ::IndexType, ::Type{T}=Float64
    ) where {T}
        (s, ℓₘₐₓ, T)
    end

    # Integer indices only, with a hint saying what to use instead: one written over several
    # lines as a concatenation, and one as a single literal.
    @index_methods integer_only (
        "Functions of half-integer spin may be sampled with `spiral`, "
        * "which accepts half-integer indices."
    ) function hinted(s::IndexType, ℓₘₐₓ::IndexType, ::Type{T}=Float64) where {T}
        (s, ℓₘₐₓ, T)
    end
    @index_methods integer_only "Use `pair` for half-integer indices." single_hint(
        ℓ::IndexType
    ) = ℓ

    # An integer-only definition beside one of the same function whose other argument is less
    # specific.
    @index_methods layered(x::AbstractVector, s::IndexOrRange) = (:any, s)
    @index_methods integer_only layered(x::Vector{Float64}, s::IndexType) = (:integer, s)

    # Two arities of one function, as `Ysize(ℓₘₐₓ)` and `Ysize(ℓₘᵢₙ, ℓₘₐₓ)`, with the body
    # branching on the index type, which is resolved at compile time.
    @index_methods function count_modes(ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {IT<:IndexType}
        IT === Int ? (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2 : ((2ℓₘₐₓ + 2)^2 - (2ℓₘᵢₙ)^2) ÷ 4
    end
    @index_methods count_modes(ℓₘₐₓ::IT) where {IT<:IndexType} =
        count_modes(IT === Int ? 0 : HalfOddInteger(1//2), ℓₘₐₓ)

    # A callable object, with the two call shapes of the operators.
    struct Operator end
    Base.show(io::IO, ::Operator) = print(io, "Ω")
    @index_methods function (op::Operator)(
        s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, ::Type{T}=Float64
    ) where {IT<:IndexType, T}
        (op, s, ℓₘᵢₙ, ℓₘₐₓ, T)
    end
    @index_methods (op::Operator)(s::IndexType, ℓₘₐₓ::IndexType, ::Type{T}=Float64) where {T} =
        op(s, abs(s), ℓₘₐₓ, T)

    # An inner constructor, as `ModeWeights` validates its indices, and a parameterized outer
    # constructor, as `ModeWeights{T}(undef, s, ℓₘₐₓ)`.
    struct Labelled{T, IT}
        data::Vector{T}
        ℓ::IT
        @index_methods function Labelled(data::Vector{T}, ℓ::IT) where {T, IT<:IndexType}
            ℓ < 0 && throw(ArgumentError("ℓ=$ℓ must be non-negative."))
            new{T, IT}(data, ℓ)
        end
    end
    @index_methods Labelled{T}(::UndefInitializer, ℓ::IndexType) where {T} =
        Labelled(Vector{T}(undef, 2), ℓ)

    # A keyword that may be `nothing`, a keyword typed by the type variable of the positional
    # indices, keywords that are not indices, and a keyword splat.
    @index_methods optional(ℓ::IndexType; ℓₘᵢₙ::Union{Nothing, IndexType}=nothing) = (ℓ, ℓₘᵢₙ)
    @index_methods same_kind(ℓ::IT; m::IT=ℓ) where {IT<:IndexType} = (IT, ℓ, m)
    @index_methods passing(ℓ::IndexType; m::IndexType=ℓ, Nᵣ::Int=1, rest...) =
        (ℓ, m, Nᵣ, values(rest))

    # A `where` parameter that only a keyword uses, such as the element type of a keyword
    # `T::Type{TT}=Float64`, which the methods without that keyword must not declare.
    @index_methods element_keyword(ℓ::IndexType; T::Type{TT}=Float64) where {TT} = (ℓ, TT)

    # A keyword index captured by a closure in the body, which must not be boxed.
    @index_methods captured(ℓ::IndexType; m::IndexType=ℓ) = map(x -> x + m, (1, 2))

    # A keyword that may be `nothing` or of the positional type variable, and a keyword that
    # may be a range.
    @index_methods function optional_same_kind(
        ℓ::IT; m::Union{Nothing, IT}=nothing
    ) where {IT<:IndexType}
        (ℓ, m)
    end
    @index_methods ranged(ℓ::IndexType; s::IndexOrRange=ℓ) = (ℓ, s)

    # A keyword default that combines an earlier keyword index with a positional one, and so
    # needs the earlier keyword normalized.
    @index_methods function clipped(
        ℓ::IndexType; m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max, mₘᵢₙ::IndexType=max(-mₘₐₓ, -ℓ)
    )
        (ℓ, mₘₐₓ, mₘᵢₙ)
    end

    # Positional defaults that follow an index argument, which are evaluated where the indices
    # have been converted: one computed from an index, and one that uses the index type
    # variable.
    @index_methods sampled(ℓₘₐₓ::IndexType, Nϕ::Int=2ℓₘₐₓ+1) = (ℓₘₐₓ, Nϕ)
    @index_methods from_smallest(ℓ::IT, m::IT=smallest(IT)) where {IT<:IndexType} = (ℓ, m)
    smallest(::Type{Int}) = 0
    smallest(::Type{HalfOddInteger}) = HalfOddInteger(1//2)

    # An optional first index argument, as in `ModeWeights(data, s=0)`, followed by an
    # optional element type, and a keyword default computed from the index.
    @index_methods function completed(
        data::AbstractVector, s::IndexType=0, ::Type{T}=Float64;
        ell_min::IndexType=abs(s), ℓₘᵢₙ::IndexType=ell_min
    ) where {T}
        (s, T, ℓₘᵢₙ)
    end

    # An optional index argument before a vararg, and a `where` parameter needed only through
    # the bound of another.
    @index_methods optional_then_rest(ℓ::IndexType, m::IndexType=ℓ, rest...) = (ℓ, m, rest)
    @index_methods function chained(
        x::V, ℓ::IndexType, ::Type{T}=Float64
    ) where {S, V<:AbstractVector{S}, T}
        (S, ℓ, T)
    end

    # The defaults of a call that is refused are not evaluated.
    const evaluated_defaults = Any[]
    @index_methods witnessed(ℓ::IndexType, x=push!(evaluated_defaults, ℓ)) = ℓ

    # A non-index vararg, and an applied macro.
    @index_methods trailing(ℓ::IndexType, rest...) = (ℓ, rest)
    @index_methods @inline inlined(ℓ::IndexType) = ℓ + 1

    # A docstring given with `@doc`, as the pixelizations give theirs.
    @doc raw"""
        raw_documented(ℓ)

    A raw docstring, with a backslash: \alpha.
    """
    @index_methods raw_documented(ℓ::IndexType) = ℓ

    # The docstrings of a name in this module, as pairs of their text and the signature they
    # are attached to.  These are read from the module's own table of docstrings rather than
    # through `@doc`, whose result depends on whether Markdown is loaded.
    docstrings(name::Symbol) = [
        join(string.(d.text)) => d.data[:typesig]
        for d ∈ values(Base.Docs.meta(@__MODULE__)[Base.Docs.Binding(@__MODULE__, name)].docs)
    ]

    # The exception a call throws, or `nothing`.
    caught(thunk) = try thunk(); nothing catch e; e end
    message(thunk) = let e = caught(thunk)
        e isa ArgumentError ? e.msg : error("Expected an ArgumentError; got $e.")
    end

    # Allocation counts of calls through the conversion method, measured inside a function so
    # that nothing about the call site allocates.
    allocations_pair(ℓ, m) = @allocated pair(ℓ, m)
    allocations_pair_keyword(ℓ, x) = @allocated pair(ℓ; m′ₘₐₓ=x)
    allocations_typed(s, ℓ) = @allocated typed(s, ℓ)
    allocations_spins(ℓ, s, x) = @allocated spins(nothing, ℓ, s; ℓₘᵢₙ=x)
    allocations_operator(op, s, ℓ) = @allocated op(s, ℓ)
    allocations_clipped(ℓ, x) = @allocated clipped(ℓ; m_max=x)
    allocations_completed(data, s) = @allocated completed(data, s)
    allocations_sampled(ℓ) = @allocated sampled(ℓ)
end


@testitem "@index_methods: dispatch" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: HalfOddInteger
    import .IndexMethodsExamples: pair, typed, spins, only_range, count_modes, caught
    h(x) = HalfOddInteger(x)

    # Integers of type `Int` reach the `Int` work method untouched.
    @test pair(3) === (3, 3, 3)
    @test pair(3, 2) === (3, 2, 3)
    @test typed(1, 3) === (Int, 1, 3, Float64)

    # A `Rational{Int}` with denominator 2 and a `HalfOddInteger` are the same index, and may
    # be mixed in one call; either way the body sees `HalfOddInteger`s.
    @test pair(7//2) === (h(7//2), h(7//2), h(7//2))
    @test pair(7//2, 3//2) === (h(7//2), h(3//2), h(7//2))
    @test pair(h(7//2), 3//2) === (h(7//2), h(3//2), h(7//2))
    @test pair(7//2, h(3//2)) === (h(7//2), h(3//2), h(7//2))
    @test pair(h(7//2), h(3//2)) === (h(7//2), h(3//2), h(7//2))
    @test typed(1//2, 3//2) === (HalfOddInteger, h(1//2), h(3//2), Float64)
    @test typed(1//2, h(3//2), Float32) === (HalfOddInteger, h(1//2), h(3//2), Float32)
    @test -7//2 isa Rational{Int} && pair(-7//2) === (h(-7//2), h(-7//2), h(-7//2))

    # Every other value of an index type, and every mixture of the two kinds, reaches the
    # fallback, which throws an `ArgumentError`.
    for bad ∈ (
        (0, 7//2), (0, h(7//2)), (7//2, 0), (Int8(3), 2), (Int32(3), 2), (UInt(3), 2),
        (big(3), 2), (Int128(3), 2), (true, 2), (Int8(7)//Int8(2), 3//2), (3//1, 2),
        (3//1, 3//1), (1//3, 2), (Int8(3), Int8(2)),
    )
        @test caught(() -> pair(bad...)) isa ArgumentError
        @test caught(() -> typed(bad...)) isa ArgumentError
    end
    @test caught(() -> pair(Int8(3))) isa ArgumentError
    @test caught(() -> pair(7//4)) isa ArgumentError

    # A non-index value is not an index argument at all, and remains a `MethodError`, as does
    # a call with the arguments in the wrong order.
    @test caught(() -> pair(3.0)) isa MethodError
    @test caught(() -> typed(Float64, 1, 3)) isa MethodError

    # Ranges: a `UnitRange{Int}` and a `UnitRange{HalfOddInteger}` are used as they are, and a
    # `UnitRange{Rational{Int}}` is converted endpoint by endpoint.
    @test spins(:R, 3, -2:2) === (:R, 3, -2:2, 0)
    @test spins(:R, 3, 2) === (:R, 3, 2, 2)
    @test spins(:R, 7//2, -3//2:3//2) === (:R, h(7//2), h(-3//2):h(3//2), h(1//2))
    @test spins(:R, 7//2, h(-3//2):h(3//2)) === (:R, h(7//2), h(-3//2):h(3//2), h(1//2))
    @test spins(:R, h(7//2), 1//2) === (:R, h(7//2), h(1//2), h(1//2))
    @test only_range(-3//2:3//2) === h(-3//2):h(3//2)
    @test only_range(1//2:-1//2) === h(1//2):h(-1//2)
    @test isempty(only_range(1//2:-1//2))
    @test only_range(1:3) === 1:3
    for bad ∈ (
        3:-1:-3, -3:1:3, h(3//2):-1:h(-3//2), 3//2:-1:-3//2, Base.OneTo(3), Int8(1):Int8(3),
        1//1:3//1, 1//3:4//3, (Int8(1)//Int8(2)):(Int8(5)//Int8(2)),
    )
        @test caught(() -> only_range(bad)) isa ArgumentError
        @test caught(() -> spins(:R, 3, bad)) isa ArgumentError
    end
    @test caught(() -> spins(:R, 3, -3//2:3//2)) isa ArgumentError
    @test caught(() -> spins(:R, 7//2, -2:2)) isa ArgumentError
    @test caught(() -> only_range(1.0:3.0)) isa MethodError

    # Two arities of one function, whose bodies branch on the index type.
    @test count_modes(3) === 16
    @test count_modes(1, 3) === 15
    @test count_modes(7//2) === count_modes(1//2, 7//2) === 20
    @test count_modes(h(3//2), 7//2) === 18
    @test caught(() -> count_modes(0, 7//2)) isa ArgumentError
end


@testitem "@index_methods: ambiguities" setup=[IndexMethodsExamples] begin
    using Test: detect_ambiguities, detect_unbound_args
    import SphericalFunctions
    @test isempty(detect_ambiguities(IndexMethodsExamples))
    @test isempty(detect_unbound_args(IndexMethodsExamples))

    # The package's own methods, which include those that the macro generates for its boundary
    # methods beside those written by hand.  Aqua runs the same checks, but in a process of its
    # own, whose report is harder to read.  Here the other items of the suite may already have
    # loaded packages and test modules into the process, such as a test module whose `*`
    # takes a `Rational` first.  What is found would then depend on which items happened to
    # run before this one, so only the pairs among the methods of the package, of Base, of
    # the package's direct dependencies, and of DoubleFloats and the package's extension for
    # it are looked at.  DoubleFloats is loaded here, so that its constructors, which a
    # conversion of a `HalfOddInteger` to every float type at once would be ambiguous with,
    # are always among them.
    using DoubleFloats: DoubleFloats
    function relevant(m::Method)
        root = Base.moduleroot(m.module)
        root ∈ (SphericalFunctions, Base, Core, DoubleFloats) ||
            nameof(root) === :SphericalFunctionsDoubleFloatsExt ||
            Base.identify_package(
                Base.PkgId(SphericalFunctions), String(nameof(root))
            ) == Base.PkgId(root)
    end
    extension = Base.get_extension(SphericalFunctions, :SphericalFunctionsDoubleFloatsExt)
    @test extension isa Module
    modules = extension isa Module ? (SphericalFunctions, extension) : (SphericalFunctions,)
    ambiguities = detect_ambiguities(modules...; recursive=true)
    @test isempty(filter(pair -> all(relevant, pair), ambiguities))
    @test isempty(detect_unbound_args(modules...; recursive=true))
end


@testitem "@index_methods: error messages" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: HalfOddInteger
    import .IndexMethodsExamples: pair, spins, only_range, message
    h(x) = HalfOddInteger(x)

    # The message names the function and states the rule, and lists every index argument with
    # its value and type, including those that are acceptable on their own.
    let m = message(() -> pair(0, 7//2))
        @test occursin("The indices of one call to `pair` must all be integers of type `Int`", m)
        @test occursin("or all be half-odd-integers", m)
        @test occursin("`Rational{Int}` with denominator 2, like 7//2", m)
        @test occursin("ℓ = 0::Int64", m)
        @test occursin("m = 7//2::Rational{Int64}", m)
        @test occursin("mixes integers (ℓ) with half-odd-integers (m)", m)
    end
    @test occursin("mixes integers (ℓ) with half-odd-integers (m)", message(() -> pair(0, h(7//2))))

    # Each kind of value that could not be an index on its own has its own sentence.
    @test occursin("`Int8` is narrower than `Int`", message(() -> pair(Int8(3), 2)))
    @test occursin("`Int32` is narrower than `Int`", message(() -> pair(Int32(3), 2)))
    @test occursin("`UInt64` is unsigned", message(() -> pair(UInt(3), 2)))
    @test occursin("`BigInt` is wider than `Int`", message(() -> pair(big(3), 2)))
    @test occursin("`Int128` is wider than `Int`", message(() -> pair(Int128(3), 2)))
    @test occursin("A `Bool` is not an index", message(() -> pair(true, 2)))
    @test occursin(
        "`Rational{Int8}` is not `Rational{Int}`; write the value with `Int`s, as 7//2",
        message(() -> pair(Int8(7)//Int8(2), 3//2))
    )
    @test occursin("3//1 is a whole number; write it as the integer 3", message(() -> pair(3//1, 2)))
    @test occursin("1//3 is neither an integer nor a half-odd-integer", message(() -> pair(1//3, 2)))
    let m = message(() -> spins(:R, 3, 3:-1:-3))
        @test occursin("s = (3:-1:-3)::StepRange{Int64, Int64}", m)
        @test occursin("A range of indices must be a unit range, running upward in steps of 1", m)
    end
    @test occursin("and may be written -3:3", message(() -> spins(:R, 3, -3:1:3)))
    @test occursin("write this `Base.OneTo{Int64}` as 1:3", message(() -> only_range(Base.OneTo(3))))
    @test occursin("are whole numbers; write it as 1:3", message(() -> only_range(1//1:3//1)))
    @test occursin("must have denominator 2", message(() -> only_range(1//3:4//3)))
    @test occursin("of type `Int8`", message(() -> only_range(Int8(1):Int8(3))))
    @test occursin(
        "of type `Bool`; write it with `Int`s, as 0:1", message(() -> only_range(false:true))
    )

    # A value that is acceptable on its own has no sentence of its own.
    let m = message(() -> pair(Int8(3), 7//2))
        @test occursin("m = 7//2::Rational{Int64}\n", m * "\n")
        @test !occursin("mixes", m)
    end
end


@testitem "@index_methods: keyword indices" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: HalfOddInteger
    import .IndexMethodsExamples: pair, spins, optional, same_kind, passing, captured,
        optional_same_kind, ranged, clipped, element_keyword, message, caught
    h(x) = HalfOddInteger(x)

    # A keyword index of the kind of the positional ones is used as it is, and a `Rational`
    # is converted where that kind is half-odd, whether the positional indices were given as
    # `Rational`s or as `HalfOddInteger`s.
    @test pair(3; m′ₘₐₓ=1) === (3, 3, 1)
    @test pair(7//2; m′ₘₐₓ=1//2) === (h(7//2), h(7//2), h(1//2))
    @test pair(h(7//2); m′ₘₐₓ=1//2) === (h(7//2), h(7//2), h(1//2))
    @test pair(7//2; m′ₘₐₓ=h(1//2)) === (h(7//2), h(7//2), h(1//2))
    @test spins(:R, h(7//2), -3//2:3//2; ℓₘᵢₙ=5//2) === (:R, h(7//2), h(-3//2):h(3//2), h(5//2))

    # The ASCII alias sets the Unicode keyword, of either kind, and the Unicode name wins
    # where both are given.
    @test pair(3; mp_max=2) === (3, 3, 2)
    @test pair(7//2; mp_max=1//2) === (h(7//2), h(7//2), h(1//2))
    @test pair(7//2; mp_max=h(1//2)) === (h(7//2), h(7//2), h(1//2))
    @test pair(7//2; mp_max=5//2, m′ₘₐₓ=1//2) === (h(7//2), h(7//2), h(1//2))

    # A keyword of the other kind, or not an index at all, is an `ArgumentError` naming it.
    let m = message(() -> pair(3; m′ₘₐₓ=1//2))
        @test occursin("The keyword argument `m′ₘₐₓ` of `pair`", m)
        @test occursin("which are integers of type `Int`, like 3", m)
        @test occursin("got m′ₘₐₓ = 1//2::Rational{Int64}", m)
    end
    let m = message(() -> pair(7//2; m′ₘₐₓ=1))
        @test occursin("which are half-odd-integers", m)
        @test occursin("got m′ₘₐₓ = 1::Int64", m)
    end
    @test occursin("The keyword argument `mp_max`", message(() -> pair(7//2; mp_max=1)))
    @test occursin("`Int8` is narrower", message(() -> pair(3; m′ₘₐₓ=Int8(1))))
    @test occursin("`UInt64` is unsigned", message(() -> pair(3; mp_max=UInt(1))))
    @test occursin("3//1 is a whole number", message(() -> pair(7//2; m′ₘₐₓ=3//1)))
    @test occursin(
        "`Rational{Int8}` is not `Rational{Int}`",
        message(() -> pair(7//2; m′ₘₐₓ=Int8(1)//Int8(2)))
    )
    @test caught(() -> pair(3; m′ₘₐₓ=1.0)) isa TypeError
    # The positional indices are checked first, so a bad positional index is what is named.
    @test occursin("The indices of one call to `pair`", message(() -> pair(0, 7//2; m′ₘₐₓ=1)))

    # A keyword that may be `nothing`, a keyword typed by the positional type variable, and
    # keywords that are not indices, which pass through untouched.
    @test optional(3) === (3, nothing)
    @test optional(7//2) === (h(7//2), nothing)
    @test optional(7//2; ℓₘᵢₙ=1//2) === (h(7//2), h(1//2))
    @test caught(() -> optional(7//2; ℓₘᵢₙ=1)) isa ArgumentError
    @test same_kind(7//2) === (HalfOddInteger, h(7//2), h(7//2))
    @test same_kind(7//2; m=1//2) === (HalfOddInteger, h(7//2), h(1//2))
    @test same_kind(3; m=1) === (Int, 3, 1)
    @test caught(() -> same_kind(3; m=1//2)) isa ArgumentError
    @test passing(7//2; m=1//2, Nᵣ=2, x=1//2) === (h(7//2), h(1//2), 2, (x=1//2,))
    @test passing(3; x=:y) === (3, 3, 1, (x=:y,))
    @test element_keyword(7//2) === (h(7//2), Float64)
    @test element_keyword(3; T=Float32) === (3, Float32)

    # A keyword that may be `nothing` or of the positional type variable, and a keyword range,
    # whose step must be 1 as that of a positional range must.
    @test optional_same_kind(3) === (3, nothing)
    @test optional_same_kind(7//2; m=1//2) === (h(7//2), h(1//2))
    @test occursin(
        "got m = 1//2::Rational{Int64}", message(() -> optional_same_kind(3; m=1//2))
    )
    @test ranged(7//2; s=-3//2:3//2) === (h(7//2), h(-3//2):h(3//2))
    @test ranged(3; s=-1:1) === (3, -1:1)
    @test occursin("must be a unit range", message(() -> ranged(3; s=1:-1:-1)))

    # A default that refers to an earlier keyword index sees it normalized.
    @test clipped(7//2) === (h(7//2), h(7//2), h(-7//2))
    @test clipped(7//2; m_max=1//2) === (h(7//2), h(1//2), h(-1//2))
    @test clipped(h(7//2); mₘₐₓ=1//2) === (h(7//2), h(1//2), h(-1//2))
    @test clipped(3; m_max=1) === (3, 1, -1)
    @test occursin(
        "The keyword argument `m_max` of `clipped`", message(() -> clipped(3; m_max=1//2))
    )

    # A keyword captured by a closure in the body is bound once, and so is not boxed.
    @test (@inferred captured(3; m=1)) === (2, 3)
    @test (@inferred captured(7//2; m=1//2)) === (h(3//2), h(5//2))

    # Each keyword is listed once by `methods` and in the hints of a `MethodError`, even where
    # a later default sees it normalized; so too for the package's own aliased keywords.
    package = (SphericalFunctions.D, SphericalFunctions.sYlm, SphericalFunctions.ModeWeights)
    for f ∈ (pair, clipped, package...), m ∈ methods(f)
        @test allunique(Base.kwarg_decl(m))
    end
    @test Base.kwarg_decl(which(clipped, (Int,))) == [:m_max, :mₘₐₓ, :mₘᵢₙ]
end


@testitem "@index_methods: positional defaults" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: HalfOddInteger
    import .IndexMethodsExamples: sampled, from_smallest, completed, optional_then_rest,
        chained, witnessed, evaluated_defaults, message
    h(x) = HalfOddInteger(x)

    # A default that follows an index argument sees the index converted, and may use the index
    # type variable.
    @test sampled(7//2) === (h(7//2), 8)
    @test sampled(h(7//2)) === (h(7//2), 8)
    @test sampled(3) === (3, 7)
    @test sampled(7//2, 5) === (h(7//2), 5)
    @test from_smallest(7//2) === (h(7//2), h(1//2))
    @test from_smallest(3) === (3, 0)
    @test from_smallest(7//2, 3//2) === (h(7//2), h(3//2))

    # An optional first index argument is filled in before the indices are dispatched on, and
    # the optional arguments after it once the index has been converted.
    @test completed([1.0]) === (0, Float64, 0)
    @test completed([1.0]; ℓₘᵢₙ=2) === (0, Float64, 2)
    @test completed([1.0], 1//2) === (h(1//2), Float64, h(1//2))
    @test completed([1.0], -1//2, Float32; ell_min=3//2) === (h(-1//2), Float32, h(3//2))
    @test occursin("got ℓₘᵢₙ = 1//2", message(() -> completed([1.0]; ℓₘᵢₙ=1//2)))
    @test occursin("s = 3::Int8", message(() -> completed([1.0], Int8(3))))

    # An optional index argument before a vararg, and a `where` parameter needed only by the
    # bound of another.
    @test optional_then_rest(7//2) === (h(7//2), h(7//2), ())
    @test optional_then_rest(7//2, 1//2, :x) === (h(7//2), h(1//2), (:x,))
    @test chained([1], 7//2) === (Int, h(7//2), Float64)
    @test chained([1], 3, Float32) === (Int, 3, Float32)

    # A call that is refused is refused before its defaults are evaluated, and the message
    # lists the indices it was given.
    empty!(evaluated_defaults)
    @test occursin("ℓ = 3::Int8", message(() -> witnessed(Int8(3))))
    @test isempty(evaluated_defaults)
    @test witnessed(7//2) === h(7//2)
    @test evaluated_defaults == [h(7//2)]
    let m = message(() -> sampled(Int8(3)))
        @test occursin("ℓₘₐₓ = 3::Int8", m)
        @test !occursin("Nϕ", m)
    end
end


@testitem "@index_methods: integer_only" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: HalfOddInteger
    import .IndexMethodsExamples: grid, hinted, single_hint, layered, message, caught

    @test grid(0, 3) === (0, 3, Float64)
    @test grid(0, 3, Float32) === (0, 3, Float32)
    let m = message(() -> grid(1//2, 7//2))
        @test occursin("The indices of `grid` must be integers of type `Int`, like 3", m)
        @test occursin("`grid` does not accept half-odd-integers", m)
        @test occursin("s = 1//2::Rational{Int64}\n        1//2 is a half-odd-integer.", m)
    end
    @test occursin("is a half-odd-integer", message(() -> grid(HalfOddInteger(1//2), 7//2)))
    @test occursin("3//1 is a whole number", message(() -> grid(3//1, 3)))
    @test occursin("`Int8` is narrower", message(() -> grid(Int8(0), 3)))
    @test occursin("does not accept half-odd-integers", message(() -> grid(0, 7//2)))
    # The half-odd method refuses as the conversion method does.
    let m = which(grid, Tuple{HalfOddInteger, HalfOddInteger, Type{Float64}})
        @test Base.unwrap_unionall(m.sig).parameters[2] === HalfOddInteger
    end
    @test Base.return_types(grid, (HalfOddInteger, HalfOddInteger)) == [Union{}]
    # Without a hint, the message ends with the list of arguments.
    @test endswith(message(() -> grid(1//2, 7//2)), "7//2 is a half-odd-integer.")

    # A hint is appended, on a line of its own, to the message of a call with a half-odd-integer
    # index, however it reaches the refusal, and is left out of that of a call refused for
    # another reason.
    hint = (
        "\nFunctions of half-integer spin may be sampled with `spiral`, which accepts "
        * "half-integer indices."
    )
    @test hinted(0, 3) === (0, 3, Float64)
    @test hinted(0, 3, Float32) === (0, 3, Float32)
    let m = message(() -> hinted(1//2, 7//2))
        @test occursin("`hinted` does not accept half-odd-integers.  This call has", m)
        @test occursin("ℓₘₐₓ = 7//2::Rational{Int64}\n        7//2 is a half-odd-integer.", m)
        @test endswith(m, "7//2 is a half-odd-integer." * hint)
    end
    @test endswith(message(() -> hinted(HalfOddInteger(1//2), HalfOddInteger(7//2))), hint)
    @test endswith(message(() -> hinted(0, 7//2, Float32)), hint)
    @test endswith(message(() -> hinted(HalfOddInteger(1//2), 3)), hint)
    let m = message(() -> hinted(Int8(0), 3))
        @test occursin("`Int8` is narrower than `Int`", m)
        @test !occursin("spiral", m)
    end
    @test !occursin("spiral", message(() -> hinted(3//1, 3)))
    @test endswith(message(() -> single_hint(1//2)), "\nUse `pair` for half-integer indices.")
    @test endswith(
        message(() -> single_hint(HalfOddInteger(1//2))), "for half-integer indices."
    )
    @test single_hint(2) === 2

    # An integer-only definition beside one of the same function whose other argument is less
    # specific: each call reaches the definition it matches, and none is ambiguous.
    @test layered(Float32[1], 1//2) === (:any, HalfOddInteger(1//2))
    @test layered(Float32[1], -1//2:1//2) === (:any, HalfOddInteger(-1//2):HalfOddInteger(1//2))
    @test layered([1.0], 1) === (:integer, 1)
    @test occursin("does not accept half-odd-integers", message(() -> layered([1.0], 1//2)))
    @test occursin(
        "does not accept half-odd-integers", message(() -> layered([1.0], HalfOddInteger(1//2)))
    )
end


@testitem "@index_methods: other definition shapes" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: HalfOddInteger
    import .IndexMethodsExamples: Operator, Labelled, trailing, inlined, message, caught
    h(x) = HalfOddInteger(x)
    Ω = Operator()

    # A callable object, with both call shapes and the positional element type.
    @test Ω(1, 1, 3) === (Ω, 1, 1, 3, Float64)
    @test Ω(1, 3) === (Ω, 1, 1, 3, Float64)
    @test Ω(-1//2, 7//2, Float32) === (Ω, h(-1//2), h(1//2), h(7//2), Float32)
    @test Ω(1//2, 1//2, 7//2) === (Ω, h(1//2), h(1//2), h(7//2), Float64)
    @test occursin("The indices of one call to `Ω`", message(() -> Ω(1//2, 3)))
    @test occursin("The indices of one call to `Ω`", message(() -> Ω(1//2, 1, 3, Float32)))

    # An inner constructor, and a parameterized outer one.
    @test Labelled([1.0, 2.0], 3) isa Labelled{Float64, Int}
    @test Labelled([1.0, 2.0], 7//2).ℓ === h(7//2)
    @test Labelled{Float32}(undef, 7//2) isa Labelled{Float32, HalfOddInteger}
    @test Labelled{Float32}(undef, 3) isa Labelled{Float32, Int}
    @test caught(() -> Labelled([1.0], -1)) isa ArgumentError
    @test occursin("ℓ=-1 must be non-negative", message(() -> Labelled([1.0], -1)))
    @test caught(() -> Labelled([1.0], Int8(3))) isa ArgumentError
    @test occursin("Labelled{Float32}", message(() -> Labelled{Float32}(undef, UInt(3))))

    # A non-index vararg is passed on, and an applied macro is applied to each method.
    @test trailing(7//2, :a, 2) === (h(7//2), (:a, 2))
    @test trailing(3) === (3, ())
    @test inlined(3) === 4
    @test inlined(7//2) === h(9//2)
end


@testitem "@index_methods: docstrings" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: IndexType
    import .IndexMethodsExamples: docstrings

    # The docstring written above the macro call is attached once, to the fallback that takes
    # every positional argument, whose index arguments are typed with their markers.  This
    # holds whether the docstring is a plain string or given with `@doc`.
    let d = docstrings(:pair)
        @test length(d) == 1
        @test occursin("A docstring written above the macro call.", d[1].first)
        @test d[1].second == Tuple{IndexType, IndexType}
    end
    let d = docstrings(:raw_documented)
        @test length(d) == 1
        @test occursin("A raw docstring, with a backslash: \\alpha.", d[1].first)
        @test d[1].second == Tuple{IndexType}
    end
end


@testitem "@index_methods: inference and allocation" setup=[IndexMethodsExamples] begin
    using SphericalFunctions: HalfOddInteger
    import .IndexMethodsExamples: pair, typed, spins, count_modes, Operator, clipped, sampled,
        completed, allocations_pair, allocations_pair_keyword, allocations_typed,
        allocations_spins, allocations_operator, allocations_clipped, allocations_completed,
        allocations_sampled
    Ω = Operator()

    # The conversion method returns exactly what the work method does.
    @inferred pair(7//2, 3//2)
    @inferred pair(7//2)
    @inferred pair(7//2; m′ₘₐₓ=1//2)
    @inferred pair(7//2; mp_max=1//2)
    @inferred typed(1//2, 3//2)
    @inferred typed(1//2, 3//2, Float32)
    @inferred spins(:R, 7//2, -3//2:3//2; ℓₘᵢₙ=1//2)
    @inferred count_modes(7//2)
    @inferred Ω(1//2, 7//2)
    @inferred Ω(1//2, 1//2, 7//2, Float32)
    @inferred clipped(7//2; m_max=1//2)
    @inferred sampled(7//2)
    @inferred completed([1.0])
    @inferred completed([1.0], 1//2)

    # A call that can only be refused is inferred to throw.
    @test Base.return_types(pair, (Int8,)) == [Union{}]
    @test Base.return_types(pair, (Int, Rational{Int})) == [Union{}]
    @test Base.return_types(sampled, (Int8,)) == [Union{}]

    # And it allocates nothing, once compiled.
    for _ ∈ 1:2
        allocations_pair(7//2, 3//2)
        allocations_pair_keyword(7//2, 1//2)
        allocations_typed(1//2, 3//2)
        allocations_spins(7//2, -3//2:3//2, 1//2)
        allocations_operator(Ω, 1//2, 7//2)
        allocations_clipped(7//2, 1//2)
        allocations_completed([1.0], 1//2)
        allocations_sampled(7//2)
    end
    @test allocations_pair(7//2, 3//2) == 0
    @test allocations_pair(3, 2) == 0
    @test allocations_pair_keyword(7//2, 1//2) == 0
    @test allocations_typed(1//2, 3//2) == 0
    @test allocations_spins(7//2, -3//2:3//2, 1//2) == 0
    @test allocations_operator(Ω, 1//2, 7//2) == 0
    @test allocations_clipped(7//2, 1//2) == 0
    @test allocations_completed([1.0], 1//2) == 0
    @test allocations_sampled(7//2) == 0
end


@testitem "@index_methods: misuse" setup=[IndexMethodsExamples] begin
    import .IndexMethodsExamples

    # The error raised while expanding a definition, which some versions of Julia wrap in a
    # `LoadError`.
    function expansion_error(ex)
        try
            macroexpand(IndexMethodsExamples, ex)
            nothing
        catch e
            e isa LoadError ? e.error : e
        end
    end
    function refused(ex, text)
        e = expansion_error(ex)
        e isa ArgumentError && occursin(text, e.msg)
    end

    @test refused(:(@index_methods f(x) = x), "has no positional index argument")
    @test refused(:(@index_methods f(x; m::IndexType=1) = x), "has no positional index argument")
    @test refused(:(@index_methods f(::IndexType) = 1), "must be named")
    @test refused(:(@index_methods f(ℓ::IndexType...) = ℓ), "may not be a vararg")
    @test refused(
        :(@index_methods f(ℓ::IT, x::Vector{IT}) where {IT<:IndexType} = ℓ),
        "may annotate index arguments only"
    )
    # A default may use the index type variable once an index argument precedes it.
    @test expansion_error(
        :(@index_methods f(ℓ::IT, m::IT=smallest(IT)) where {IT<:IndexType} = ℓ)
    ) === nothing
    @test refused(
        :(@index_methods f(ℓ::IT=zero(IT)) where {IT<:IndexType} = ℓ),
        "no index argument precedes that argument"
    )
    @test refused(
        :(@index_methods f(x::Int=1, ℓ::IndexType) = ℓ),
        "must come last, or just before a vararg"
    )
    @test refused(
        :(@index_methods f(s::IT, t::IT) where {IT<:IndexOrRange} = s),
        "would all have to be indices or all be ranges"
    )
    @test refused(
        :(@index_methods Box{IT}(ℓ::IT) where {IT<:IndexType} = ℓ),
        "`Box{IT}`, the function defined, uses a type variable bounded by an index marker"
    )
    @test refused(
        :(@index_methods f(ℓ::IT)::IT where {IT<:IndexType} = ℓ),
        "applies to its return type"
    )
    @test refused(:(@index_methods @generated f(ℓ::IndexType) = :ℓ), "`@generated` function")
    @test refused(:(@index_methods const x = 1), "applies to a function definition")
    @test refused(:(@index_methods half_integer f(ℓ::IndexType) = ℓ), "takes one function definition")
    @test refused(
        :(@index_methods integer_only "a" "b" f(ℓ::IndexType) = ℓ),
        "takes one function definition"
    )
    # The hint is a string literal, or a concatenation of them, and follows `integer_only`.
    @test expansion_error(
        :(@index_methods integer_only ("a" * "b") f(ℓ::IndexType) = ℓ)
    ) === nothing
    @test refused(
        :(@index_methods integer_only 3 f(ℓ::IndexType) = ℓ), "must be a string literal"
    )
    @test refused(
        :(@index_methods integer_only ("a" * b) f(ℓ::IndexType) = ℓ), "must be a string literal"
    )
    @test refused(
        :(@index_methods integer_only "a $b" f(ℓ::IndexType) = ℓ), "must be a string literal"
    )
    @test refused(:(@index_methods "hint" f(ℓ::IndexType) = ℓ), "only after `integer_only`")
end
