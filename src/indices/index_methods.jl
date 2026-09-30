# The boundary between the indices a caller writes and the indices the package computes
# with.
#
# A caller may write an integer index as `3` and a half-odd-integer index as `7//2` or as a
# `HalfOddInteger`, but the code behind every public function sees only `Int`s or only
# `HalfOddInteger`s, never a `Rational`, a narrow or unsigned integer, or a mixture of the
# two kinds.  The macro `@index_methods` writes that boundary for one function definition:
# the definition is written once, with its index arguments annotated by one of the markers
# below, and the macro generates the methods that convert, refuse and re-dispatch.
# Containers' `getindex` and `setindex!`, `recurrence!`, and `wedge_value` are not written
# this way; they check their indices, with `checked_index` below, against the kind of the
# object they are given, which dispatch on the arguments alone cannot determine.


### Markers

"""
    IndexType

`Union{Integer, Rational, HalfOddInteger}`, the marker for an index argument of a function
defined with [`@index_methods`](@ref).  These are the types a caller may pass, but not all
of them are accepted: the indices of one call must all be `Int`s or all be
half-odd-integers, given as [`HalfOddInteger`](@ref)s or as `Rational{Int}`s with
denominator 2.  Any other value of one of these types reaches a method that throws an
`ArgumentError` explaining why.

See also [`IndexRange`](@ref) and [`IndexOrRange`](@ref).
"""
const IndexType = Union{Integer, Rational, HalfOddInteger}

"""
    IndexRange

`AbstractUnitRange{<:IndexType}`, the marker for an argument of a function defined with
[`@index_methods`](@ref) that is a range of indices, such as the spin weights `-3//2:3//2`.
The range must be a `UnitRange`, written `a:b`, and the body of the function sees a
`UnitRange{Int}` or a `UnitRange{HalfOddInteger}`; a `UnitRange{Rational{Int}}` is converted
endpoint by endpoint.  Any other range is refused with an `ArgumentError` that says how to
write it: one whose step is not 1, such as `0:2:4`, and also one with the right elements but
of another type, such as `0:1:4` or `Base.OneTo(4)`, which is to be written `0:4` or `1:4`.

See also [`IndexType`](@ref) and [`IndexOrRange`](@ref).
"""
const IndexRange = AbstractUnitRange{<:IndexType}

"""
    IndexOrRange

`Union{IndexType, IndexRange}`, the marker for an argument of a function defined with
[`@index_methods`](@ref) that may be either one index or a range of them.  The spin weight
of [`sYlm`](@ref) is such an argument.

See also [`IndexType`](@ref) and [`IndexRange`](@ref).
"""
const IndexOrRange = Union{IndexType, IndexRange}

# The argument types of the generated methods, for each marker: those of the fallback, of
# the conversion method, and of the `Int` and `HalfOddInteger` work methods.  The fallback
# admits any range of indices rather than only a unit range, so that a range whose step is
# not 1, or any other range that is not a `UnitRange`, reaches the explanation rather than a
# `MethodError`.  A keyword index is annotated with the fallback's type for the same reason.
const index_method_types = (
    IndexType = (
        IndexType,
        Union{Rational, HalfOddInteger},
        Int,
        HalfOddInteger,
    ),
    IndexRange = (
        AbstractRange{<:IndexType},
        AbstractUnitRange{<:Union{Rational, HalfOddInteger}},
        UnitRange{Int},
        UnitRange{HalfOddInteger},
    ),
    IndexOrRange = (
        Union{IndexType, AbstractRange{<:IndexType}},
        Union{Rational, HalfOddInteger, AbstractUnitRange{<:Union{Rational, HalfOddInteger}}},
        Union{Int, UnitRange{Int}},
        Union{HalfOddInteger, UnitRange{HalfOddInteger}},
    ),
)
const index_marker_names = (:IndexType, :IndexRange, :IndexOrRange)


### Conversion
#
# The conversion method checks every index argument before converting any, so that a call
# which is refused is refused with all of its arguments in the message.  The check and the
# conversion are separate functions for the same reason, and because the conversion then
# cannot fail: a `Rational{Int}` that has passed the check has an odd numerator.  A range is
# checked and converted by its endpoints; `last` of a `UnitRange{Rational{Int}}` is its
# first element plus a whole number, so the two endpoints always have the same denominator.

@inline is_half_odd_index(::HalfOddInteger) = true
@inline is_half_odd_index(x::Rational{Int}) = denominator(x) == 2
@inline is_half_odd_index(::AbstractUnitRange{HalfOddInteger}) = true
@inline is_half_odd_index(r::AbstractUnitRange{Rational{Int}}) =
    is_half_odd_index(first(r)) && is_half_odd_index(last(r))
@inline is_half_odd_index(x) = false

@inline half_odd_index(x::HalfOddInteger) = x
@inline half_odd_index(x::Rational{Int}) = unsafe_half_odd_integer(numerator(x))
@inline half_odd_index(r::UnitRange{HalfOddInteger}) = r
@inline half_odd_index(r::AbstractUnitRange{HalfOddInteger}) = first(r):last(r)
@inline half_odd_index(r::AbstractUnitRange{Rational{Int}}) =
    half_odd_index(first(r)):half_odd_index(last(r))

# Keyword arguments take no part in dispatch, so the work methods normalize each keyword
# annotated with a marker against the kind of the positional indices, `Int` or
# `HalfOddInteger`, which dispatch has already settled.  A value of that kind is returned
# unchanged, a `Rational{Int}` with denominator 2 is converted where the kind is
# `HalfOddInteger`, and `nothing` passes through for a keyword annotated `Union{Nothing,
# IndexType}`.  Anything else is refused, naming the keyword.
@inline index_keyword(f, name, x::Nothing, ::Type) = x
@inline index_keyword(f, name, x::Union{Int, UnitRange{Int}}, ::Type{Int}) = x
@inline index_keyword(
    f, name, x::Union{HalfOddInteger, UnitRange{HalfOddInteger}}, ::Type{HalfOddInteger}
) = x
@inline function index_keyword(
    f, name, x::Union{Rational, AbstractUnitRange{<:Rational}}, ::Type{HalfOddInteger}
)
    is_half_odd_index(x) || throw(index_keyword_error(f, name, x, HalfOddInteger))
    half_odd_index(x)
end
index_keyword(f, name, x, ::Type{K}) where {K} = throw(index_keyword_error(f, name, x, K))

# The index `x` of an object whose kind of index is fixed — a container, a calculator, or a
# wedge, whose index type is `IT` — converted to that type.  Wherever a call names such an
# object and its indices, as indexing a container, `recurrence!(calc, ℓ)`, and `wedge_value`
# do, the kind of index is that of the object rather than one that dispatch on the arguments
# could choose, so these calls are written by hand rather than with `@index_methods`; but
# they accept what the index methods accept: an `Int` for an object with integer indices,
# and a `HalfOddInteger` or a `Rational{Int}` with denominator 2 for one with half-integer
# indices.  Anything else is refused with the reason the index methods give; a narrow
# integer, in particular, is told to be converted, because the arithmetic that finds an
# element's position is not closed under it.  The message names the index by `name` and the
# object `owner` by `container_name`.  An index already of the object's own type is returned
# as it is, which is the path that every loop over the elements takes.
@inline checked_index(::Type{Int}, x::Int, owner, name) = x
@inline checked_index(::Type{HalfOddInteger}, x::HalfOddInteger, owner, name) = x
@inline function checked_index(::Type{HalfOddInteger}, x::Rational{Int}, owner, name)
    is_half_odd_index(x) || throw(checked_index_error(HalfOddInteger, x, owner, name))
    half_odd_index(x)
end
checked_index(::Type{IT}, x, owner, name) where {IT} =
    throw(checked_index_error(IT, x, owner, name))


### Errors

# The kind of a value that is acceptable as an index on its own — `:integer` or `:half` — or
# `nothing` for a value that is not.
index_kind(::Int) = :integer
index_kind(::UnitRange{Int}) = :integer
index_kind(::HalfOddInteger) = :half
index_kind(::UnitRange{HalfOddInteger}) = :half
index_kind(x::Union{Rational{Int}, AbstractUnitRange{Rational{Int}}}) =
    is_half_odd_index(x) ? :half : nothing
index_kind(x) = nothing

# What is wrong with a value that is not acceptable as an index, as a sentence addressed to
# the caller.  The integer types are refused rather than converted because the index
# arithmetic is not closed under them: `ℓ^2` overflows a narrow type, `-m` and `ℓ - 1` wrap
# around in an unsigned one, and a wider type would reach code that stores `Int`s.
function index_problem(x::Integer)
    T = typeof(x)
    if x isa Bool
        "A `Bool` is not an index; write 0 or 1."
    elseif x isa Unsigned
        (
            "`$T` is unsigned, and index arithmetic such as `-m` and `ℓ - 1` wraps around in "
            * "it; convert it with `Int`."
        )
    elseif x isa Union{Int8, Int16, Int32, Int64, Int128} && sizeof(T) < sizeof(Int)
        (
            "`$T` is narrower than `Int`, and index arithmetic such as `ℓ^2` overflows in it; "
            * "convert it with `Int`."
        )
    elseif x isa Union{Int8, Int16, Int32, Int64, Int128, BigInt}
        (
            "`$T` is wider than `Int`, which is the type in which indices are stored and "
            * "computed; convert it with `Int`."
        )
    else
        "`$T` is not `Int`; convert it with `Int`."
    end
end
function index_problem(x::Rational)
    T = typeof(x)
    if denominator(x) == 1
        "$(repr(x)) is a whole number; write it as the integer $(numerator(x))."
    elseif denominator(x) != 2
        "$(repr(x)) is neither an integer nor a half-odd-integer."
    else
        "`$T` is not `Rational{Int}`; write the value with `Int`s, as $(numerator(x))//2."
    end
end
function index_problem(r::AbstractRange)
    lo, hi = repr(first(r)), repr(last(r))
    if !(r isa AbstractUnitRange)
        (
            "A range of indices must be a unit range, running upward in steps of 1 and "
            * "written a:b; $(repr(r)) is a `$(typeof(r))`"
            * (step(r) == 1 ? ", and may be written $lo:$hi." : ".")
        )
    elseif eltype(r) <: Union{Int, HalfOddInteger}
        (
            "A range of indices must be a `UnitRange`, written a:b; write this "
            * "`$(typeof(r))` as $lo:$hi."
        )
    elseif eltype(r) === Rational{Int}
        if denominator(first(r)) == 1
            (
                "The endpoints of $(repr(r)) are whole numbers; write it as "
                * "$(numerator(first(r))):$(numerator(last(r)))."
            )
        else
            (
                "The endpoints of a range of half-odd-integers must have denominator 2, as "
                * "in -3//2:3//2; those of $(repr(r)) do not."
            )
        end
    elseif eltype(r) <: Integer
        (
            "The elements of $(repr(r)) are of type `$(eltype(r))`; write it with `Int`s, as "
            * "$(Int(first(r))):$(Int(last(r)))."
        )
    else
        (
            "The elements of $(repr(r)) are of type `$(eltype(r))`; write it with `Int`s, as "
            * "$lo:$hi."
        )
    end
end
index_problem(x) = "`$(typeof(x))` is not an index type."

# A value as it is printed in a message, with its type, written as Julia would parse it.
typed_repr(x) = "$(repr(x))::$(typeof(x))"
typed_repr(r::AbstractRange) = "($(repr(r)))::$(typeof(r))"

# The sentence, if any, to print below one argument in the message.  An argument of either
# kind is acceptable on its own, and needs no sentence, except that a half-odd-integer needs
# one where only integers are accepted.
function index_note(x, integer_only::Bool)
    kind = index_kind(x)
    if kind === :integer
        nothing
    elseif kind === :half
        !integer_only ? nothing :
            x isa AbstractRange ? "$(repr(x)) is a range of half-odd-integers." :
            "$(repr(x)) is a half-odd-integer."
    else
        index_problem(x)
    end
end

"""
    index_argument_error(f, names, values, integer_only=false, hint="")

The `ArgumentError` thrown when the index arguments of a call to `f`, a function defined
with [`@index_methods`](@ref), are not all `Int`s or all half-odd-integers.  The arguments
`names` and `values` are tuples of the names and the values of every index argument of the
call, and each is listed in the message with its type.  Below each argument that could not
be an index on its own — a narrow, unsigned or wide integer, a `Bool`, a `Rational` whose
integer type is not `Int` or whose denominator is not 2, or a range that is not a
`UnitRange` written `a:b`, such as `0:2:4`, `0:1:4` or `Base.OneTo(4)` — a sentence says
what is wrong with it and how to write it instead.  A call whose arguments are each
acceptable but of different kinds is told which are integers and which are
half-odd-integers.

With `integer_only` the message says that `f` accepts integer indices only, and a
half-odd-integer argument is marked as the problem.  A non-empty `hint` is then appended, on
a line of its own, to the message of a call that has a half-odd-integer argument; it is
meant to say what to use instead, such as the functions that do accept half-integer indices.

`f` is printed with `print`, so it may be the name of the function as a `String`, or the
function or callable object itself.
"""
@noinline function index_argument_error(
    f, names::Tuple, values::Tuple, integer_only::Bool=false, hint::AbstractString=""
)
    io = IOBuffer()
    if integer_only
        print(io,
            "The indices of `", f, "` must be integers of type `Int`, like 3; `", f,
            "` does not accept half-odd-integers.  This call has"
        )
    else
        print(io,
            "The indices of one call to `", f, "` must all be integers of type `Int`, like ",
            "3, or all be half-odd-integers, each a `HalfOddInteger` or a `Rational{Int}` ",
            "with denominator 2, like 7//2.  This call has"
        )
    end
    for (name, value) ∈ zip(names, values)
        print(io, "\n    ", name, " = ", typed_repr(value))
        let note = index_note(value, integer_only)
            note === nothing || print(io, "\n        ", note)
        end
    end
    if integer_only
        if !isempty(hint) && any(value -> index_kind(value) === :half, values)
            print(io, "\n", hint)
        end
    else
        kinds = map(index_kind, values)
        integers = [name for (name, kind) ∈ zip(names, kinds) if kind === :integer]
        halves = [name for (name, kind) ∈ zip(names, kinds) if kind === :half]
        if !isempty(integers) && !isempty(halves)
            print(io,
                "\nand so mixes integers (", join(integers, ", "), ") with half-odd-integers (",
                join(halves, ", "), ")."
            )
        end
    end
    ArgumentError(String(take!(io)))
end

# The error for a keyword index whose kind is not that of the positional indices, or which is
# not an index at all.
@noinline function index_keyword_error(f, name::Symbol, x, ::Type{K}) where {K}
    kind = K === Int ? "integers of type `Int`, like 3" : (
        "half-odd-integers, each a `HalfOddInteger` or a `Rational{Int}` with denominator 2, "
        * "like 7//2"
    )
    message = (
        "The keyword argument `$name` of `$f` must be an index of the same kind as the "
        * "positional indices of the call, which are $kind; got $name = $(typed_repr(x))."
    )
    if index_kind(x) === nothing
        message *= "  " * index_problem(x)
    end
    ArgumentError(message)
end

# The error for an index that is not of the kind of the object it is given with, or which is
# not an index at all; see `checked_index`.
@noinline function checked_index_error(::Type{IT}, x, owner, name) where {IT}
    kind = IT === HalfOddInteger ? (
        "half-odd-integers, each a `HalfOddInteger` or a `Rational{Int}` with denominator 2, "
        * "like 7//2"
    ) : "integers of type `Int`, like 3"
    message = (
        "The indices of this `$(container_name(owner))` are $kind; got $name = "
        * "$(typed_repr(x))."
    )
    index_kind(x) === nothing && (message *= "  " * index_problem(x))
    ArgumentError(message)
end

# What a message calls the owner of an index: the name of its type, or, for a calculator
# whose struct serves two constructors, the name of the constructor that the caller wrote,
# which the methods beside those calculators give.
container_name(owner) = nameof(typeof(owner))

# The kinds of index of two objects, which must agree when one is applied to the other.
# `HalfOddInteger` and `Integer` deliberately do not promote — but `≤` between them is well
# defined and returns an ordinary `Bool`, so a mismatch of kinds sails straight through a
# test of a range of ℓ, and surfaces much later as an `InexactError` from `convert`, or as a
# complaint from inside a loop.  The products of the labelled containers therefore compare
# the kind of the mode weights `w`, which is the type of their degrees, with the kind `IT`
# of the series or the calculator applied to them, which `what` names, before anything else.
index_kind_name(::Type{<:Integer}) = "integers"
index_kind_name(::Type{HalfOddInteger}) = "half-odd-integers"

function check_same_kind(::Type{IT}, w, what) where {IT<:IntegerHalf}
    JT = typeof(ℓₘᵢₙ(w))
    (IT <: Integer) === (JT <: Integer) && return nothing
    throw(ArgumentError(
        "These mode weights are indexed by $(index_kind_name(JT)) — "
        * "ℓ ∈ $(ℓₘᵢₙ(w)):$(ℓₘₐₓ(w)) — but $what is indexed by $(index_kind_name(IT)); "
        * "the two must be of one kind."
    ))
end


### The macro

"""
    @index_methods function f(args...; kwargs...) ... end
    @index_methods integer_only function f(args...; kwargs...) ... end
    @index_methods integer_only "hint" function f(args...; kwargs...) ... end

Define `f` for index arguments of any type a caller may reasonably write, while the body
sees only indices of type `Int` or only [`HalfOddInteger`](@ref)s.

The index arguments of the definition are annotated with one of the markers
[`IndexType`](@ref), [`IndexRange`](@ref) and [`IndexOrRange`](@ref), either directly, as in
`ℓₘₐₓ::IndexType`, or through a type variable bounded by a marker, as in `ℓₘₐₓ::IT, s::IT`
with `where {IT<:IndexType}`, in which case the body may use `IT`.  Other arguments,
positional defaults and `where` clauses are written as usual.  The macro generates four
methods of `f`:

1. A fallback, whose index arguments are typed with their markers, and which throws the
   `ArgumentError` of [`index_argument_error`](@ref), naming every index argument with its
   value and type and saying what is wrong.  This is the method that a call reaches if its
   indices are not all of one kind, or if any of them is an integer of a type other than
   `Int`, a `Bool`, a `Rational` whose integer type is not `Int` or whose denominator is not
   2, or a range that is not a `UnitRange` written `a:b`: one whose step is not 1, or one
   such as `0:1:4` or `Base.OneTo(4)`, whose message says to write it `0:4` or `1:4`.  The
   docstring written above the macro call is attached to this method.
2. A conversion method, whose index arguments are typed `Union{Rational, HalfOddInteger}`,
   or the unit ranges of those, which converts each `Rational{Int}` with denominator 2 to a
   `HalfOddInteger`, and a range endpoint by endpoint, and calls `f` again with the keyword
   arguments forwarded untouched.  Any other `Rational` is refused with the fallback's
   error.
3. A work method whose index arguments are `Int`s or `UnitRange{Int}`s (with `IT<:Int` for a
   type variable), whose keyword arguments are exactly those written, and whose body is the
   body written.
4. The same with `HalfOddInteger` and `UnitRange{HalfOddInteger}`.

Defaults, positional and keyword, are evaluated in the work methods, where the positional
indices have already been converted, so that a default computed from an index, such as
`Nϕ::Int=2ℓₘₐₓ+1`, sees an `Int` or a `HalfOddInteger`, and a default may use `IT`.  The
exception is the default of an argument that no index argument precedes, including that of
the first index argument itself: no index is known where it is evaluated, so it may not use
`IT`, and a call that omits it is completed by one more generated method, which evaluates it
and calls `f` again.  A definition with optional positional arguments has a fallback and a
conversion method for each number of arguments with which it may be called, and the
docstring is attached to the fallback that takes them all.  A call that is refused is
refused before any default is evaluated.

Keyword arguments take no part in dispatch, so each keyword annotated with a marker, or with
an index type variable, which then stands for its marker, or with the `Union` of either with
`Nothing`, is normalized at the start of each work method, after the defaults have been
evaluated.  A keyword of the kind of the positional indices is used as it is, a
`Rational{Int}` with denominator 2 is converted where that kind is `HalfOddInteger`,
`nothing` is passed through where the annotation allows it, and anything else is an
`ArgumentError` naming the keyword.  A default that refers to an earlier keyword index sees
that index normalized in the same way.  This is what makes the ASCII-alias pattern

```julia
@index_methods function f(ℓ::IndexType; mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max)
    # only m′ₘₐₓ is used here
end
```

accept either spelling, of either kind.  A keyword that is not annotated with a marker is
passed through untouched.

With `integer_only`, the conversion method and the fourth method both throw an
`ArgumentError` saying that `f` accepts integer indices only.  The fourth method is kept,
rather than omitted, so that the methods of `f` remain free of ambiguities when another
definition of `f` has less specific non-index arguments.  A string written after
`integer_only` is a hint, which is appended to that message when the call has a
half-odd-integer index, to say what to use instead:

```julia
@index_methods integer_only (
    "Functions of half-integer spin may be sampled with `golden_ratio_spiral_pixels`, "
    * "`leja_pixels` or `sorted_ring_pixels`."
) function driscoll_healy_pixels(s::IndexType, ℓₘₐₓ::IndexType, ::Type{T}=Float64) where {T}
    ...
end
```

The hint is a string literal, or a concatenation of string literals with `*` as here, so
that a long one may be written over several lines; it is fixed when the definition is
expanded.

The definition must have at least one positional index argument, since that is what the work
methods dispatch on, and every index argument must be named.  A type variable bounded by a
marker may annotate index arguments only, and one bounded by `IndexOrRange` may annotate
only one of them, because the arguments it annotated would then all have to be of one type,
all indices or all ranges.  A short-form definition with both a return type and a `where` is
written with the signature in parentheses, as `(f(ℓ::IT)::IT) where {IT<:IndexType} = ...`,
since otherwise the `where` applies to the return type.  Macros applied to the definition,
such as `@inline` but not `@generated`, are written between `@index_methods` and `function`,
and are applied to each generated method.

# Example

```julia
\"\"\"
    Ysize(ℓₘᵢₙ, ℓₘₐₓ)

...
\"\"\"
@index_methods function Ysize(ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {IT<:IndexType}
    ...
end
```
"""
macro index_methods(arguments...)
    index_methods_expansion(__source__, arguments...)
end

# The marker a type annotation names, as a `Symbol`, or `nothing` if it names none.  A
# marker may be written bare or qualified by its module.
function index_marker(T)
    name = if T isa Symbol
        T
    elseif T isa GlobalRef
        T.name
    elseif Meta.isexpr(T, :., 2) && T.args[2] isa QuoteNode
        T.args[2].value
    else
        nothing
    end
    name ∈ index_marker_names ? name : nothing
end

# The marker that a keyword's annotation names directly or through a type variable bounded
# by a marker, or `nothing`.
function index_marker(T, index_variables)
    marker = index_marker(T)
    marker !== nothing && return marker
    T isa Symbol && haskey(index_variables, T) ? index_variables[T] : nothing
end

# The marker of a keyword's annotation, and whether the keyword may also be `nothing`.
function keyword_marker(T, index_variables)
    marker = index_marker(T, index_variables)
    marker !== nothing && return (marker, false)
    if Meta.isexpr(T, :curly, 3) &&
            (T.args[1] === :Union || T.args[1] == GlobalRef(Core, :Union))
        a, b = T.args[2], T.args[3]
        if a === :Nothing && index_marker(b, index_variables) !== nothing
            return (index_marker(b, index_variables), true)
        elseif b === :Nothing && index_marker(a, index_variables) !== nothing
            return (index_marker(a, index_variables), true)
        end
    end
    (nothing, false)
end

# The value of a string written as a literal or as a concatenation of literals with `*`, or
# `nothing` for any other expression.
literal_string(ex::String) = ex
function literal_string(ex::Expr)
    (Meta.isexpr(ex, :call) && ex.args[1] === :*) || return nothing
    pieces = map(literal_string, ex.args[2:end])
    any(isnothing, pieces) ? nothing : join(pieces)
end
literal_string(ex) = nothing

# Whether an expression mentions any of the given symbols.
expression_mentions(ex::Symbol, names) = ex ∈ names
expression_mentions(ex::Expr, names) = any(arg -> expression_mentions(arg, names), ex.args)
expression_mentions(ex, names) = false

# A copy of an expression in which each symbol that is a key of `names` is replaced by its
# value.  The name of a keyword argument in a call is left as it is, a keyword written by
# name alone, as in `f(; m)`, is given its value explicitly, and nothing quoted is changed.
replace_symbols(ex::Symbol, names) = get(names, ex, ex)
function replace_symbols(ex::Expr, names)
    Meta.isexpr(ex, (:quote, :inert)) && return ex
    Meta.isexpr(ex, :kw, 2) && return Expr(:kw, ex.args[1], replace_symbols(ex.args[2], names))
    if Meta.isexpr(ex, :parameters)
        return Expr(:parameters, map(ex.args) do arg
            arg isa Symbol && haskey(names, arg) ? Expr(:kw, arg, names[arg]) :
                replace_symbols(arg, names)
        end...)
    end
    Expr(ex.head, map(arg -> replace_symbols(arg, names), ex.args)...)
end
replace_symbols(ex, names) = ex

# The name of the type variable that a `where` parameter introduces, or `nothing` for a form
# not recognized here.
parameter_name(p::Symbol) = p
parameter_name(p::Expr) =
    Meta.isexpr(p, (:<:, :>:), 2) ? p.args[1] :
    Meta.isexpr(p, :comparison, 5) ? p.args[3] :
    nothing
parameter_name(p) = nothing

# The `where` parameters that a method with the given argument types needs, in their order:
# those that the types mention, and those that the bounds of the parameters needed mention
# in turn.  A bound may refer only to a parameter listed before its own, so one pass from
# the last parameter to the first is enough.  A parameter of a form not recognized is kept.
function needed_parameters(parameters, types)
    needed = falses(length(parameters))
    for i ∈ reverse(eachindex(parameters))
        name = parameter_name(parameters[i])
        inner = i+1:length(parameters)
        needed[i] = name === nothing ||
            any(T -> expression_mentions(T, (name,)), types) ||
            any(j -> needed[j] && expression_mentions(parameters[j], (name,)), inner)
    end
    parameters[needed]
end

function index_methods_expansion(source::LineNumberNode, arguments...)
    usage = (
        "Use `@index_methods function f(...) ... end`, optionally with `integer_only`, or "
        * "`integer_only` followed by a hint string, before `function`."
    )
    integer_only, hint = false, ""
    if length(arguments) == 1
        definition = arguments[1]
    elseif 2 ≤ length(arguments) ≤ 3 && arguments[1] === :integer_only
        integer_only, definition = true, arguments[end]
        if length(arguments) == 3
            hint = literal_string(arguments[2])
            hint === nothing && throw(ArgumentError(
                "The hint given after `integer_only` must be a string literal, or a "
                * "concatenation of string literals with `*`; `$(arguments[2])` is neither."
            ))
        end
    elseif length(arguments) == 2 && literal_string(arguments[1]) !== nothing
        throw(ArgumentError(
            "A hint may be given to `@index_methods` only after `integer_only`, since it is "
            * "appended to the message refusing half-odd-integer indices.  $usage"
        ))
    else
        throw(ArgumentError("`@index_methods` takes one function definition.  $usage"))
    end

    # Macros applied to the definition, outermost first, to be applied to every method.  The
    # body of a `@generated` function returns the code to run, which the generated methods'
    # bodies are not.
    wrappers = Any[]
    while Meta.isexpr(definition, :macrocall)
        applied = definition.args[1]
        applied = applied isa GlobalRef ? applied.name :
            Meta.isexpr(applied, :., 2) && applied.args[2] isa QuoteNode ?
            applied.args[2].value : applied
        applied === Symbol("@generated") && throw(ArgumentError(
            "`@index_methods` does not apply to a `@generated` function, because the "
            * "methods it generates run their bodies rather than return code."
        ))
        push!(wrappers, definition.args[1:end-1])
        definition = definition.args[end]
    end
    rewrap(method) = foldr((w, m) -> Expr(:macrocall, w..., m), wrappers; init=method)

    # The long form `function f(...) ... end` or the short form `f(...) = ...`.
    short_form = Meta.isexpr(definition, :(=), 2) &&
        Meta.isexpr(definition.args[1], (:call, :where, :(::)))
    if Meta.isexpr(definition, :function, 2) || short_form
        signature, body = definition.args
    else
        throw(ArgumentError("`@index_methods` applies to a function definition.  $usage"))
    end
    body = Meta.isexpr(body, :block) ? body : Expr(:block, source, body)
    first_line = something(findfirst(x -> x isa LineNumberNode, body.args), 0)
    line = first_line == 0 ? source : body.args[first_line]

    # The `where` parameters, outermost first, which is the order in which a single `where`
    # must list them; the return type; and the call.
    where_parameters = Any[]
    while Meta.isexpr(signature, :where)
        append!(where_parameters, signature.args[2:end])
        signature = signature.args[1]
    end
    return_type = nothing
    if Meta.isexpr(signature, :(::), 2)
        signature, return_type = signature.args
    end
    # In the short form, `f(ℓ::IT)::IT where {IT<:IndexType} = ...` parses with the `where`
    # applied to the return type rather than to the method.
    if Meta.isexpr(return_type, :where) &&
            any(p -> Meta.isexpr(p, :<:, 2) && index_marker(p.args[2]) !== nothing,
                return_type.args[2:end])
        throw(ArgumentError(
            "The `where` of this definition applies to its return type, "
            * "`$(return_type.args[1])`, rather than to the method.  Write the signature in "
            * "parentheses, as in `(f(ℓ::IT)::IT) where {IT<:IndexType} = ...`."
        ))
    end
    Meta.isexpr(signature, :call) ||
        throw(ArgumentError("`@index_methods` could not read the signature.  $usage"))

    # The type variables bounded by a marker.
    index_variables = Dict{Symbol, Symbol}()
    for p ∈ where_parameters
        if Meta.isexpr(p, :<:, 2) && p.args[1] isa Symbol && index_marker(p.args[2]) !== nothing
            index_variables[p.args[1]] = index_marker(p.args[2])
        end
    end
    variable_names = collect(keys(index_variables))
    other_parameters = filter(where_parameters) do p
        !(Meta.isexpr(p, :<:, 2) && haskey(index_variables, p.args[1]))
    end

    # The function being defined, the expression that calls it, and what the error messages
    # call it: its name, or for a callable object or a parameterized type, the object
    # itself.  A parameterized type is shown in the messages of a call as it was written,
    # with the values of its parameters, such as `ModeWeights{Float64}`, rather than as the
    # `UnionAll` it evaluates to, which may print with every parameter and bound.
    callee = signature.args[1]
    expression_mentions(callee, variable_names) && throw(ArgumentError(
        "`$callee`, the function defined, uses a type variable bounded by an index marker, "
        * "which may annotate index arguments only."
    ))
    if callee isa Symbol || Meta.isexpr(callee, :.)
        caller, name = callee, string(callee)
        shown = name
    elseif Meta.isexpr(callee, :curly)
        caller, name = callee, callee
        type_parameters = callee.args[2:end]
        shown = Expr(:string, string(callee.args[1]), "{",
            [i == 1 ? p : Expr(:string, ", ", p) for (i, p) ∈ enumerate(type_parameters)]...,
            "}")
    elseif Meta.isexpr(callee, :(::), 2)
        caller = name = shown = callee.args[1]
    elseif Meta.isexpr(callee, :(::), 1)
        caller = name = shown = gensym(:callable)
        callee = Expr(:(::), caller, callee.args[1])
    else
        throw(ArgumentError("`@index_methods` could not read the name of the function."))
    end
    for p ∈ other_parameters
        expression_mentions(p, variable_names) && throw(ArgumentError(
            "The `where` parameter `$p` of `$name` uses a type variable bounded by an index "
            * "marker, which may annotate index arguments only."
        ))
    end

    # The keyword arguments.  A keyword index is annotated with the type of the fallback's
    # argument for its marker, so that a range that is not a `UnitRange` reaches the
    # explanation.
    arguments = signature.args[2:end]
    keyword_items = Any[]
    if !isempty(arguments) && Meta.isexpr(arguments[1], :parameters)
        keyword_items = arguments[1].args
        arguments = arguments[2:end]
    end
    keywords = map(keyword_items) do item
        if Meta.isexpr(item, :...)
            return (; item, name=nothing, argument=item, default=nothing, has_default=false,
                marker=nothing, annotation=nothing)
        end
        has_default = Meta.isexpr(item, :kw, 2)
        argument, default = has_default ? (item.args[1], item.args[2]) : (item, nothing)
        keyword_name, T = argument isa Symbol ? (argument, nothing) :
            Meta.isexpr(argument, :(::), 2) ? (argument.args[1], argument.args[2]) :
            throw(ArgumentError(
                "`@index_methods` could not read the keyword argument `$item`."
            ))
        marker, optional = T === nothing ? (nothing, false) : keyword_marker(T, index_variables)
        annotation = marker === nothing ? nothing :
            optional ? Union{Nothing, index_method_types[marker][1]} :
            index_method_types[marker][1]
        (; item, name=keyword_name, argument, default, has_default, marker, annotation)
    end
    index_keywords = [k for k ∈ keywords if k.marker !== nothing]

    # The positional arguments.
    positional = map(enumerate(arguments)) do (position, argument)
        argument, default, has_default = Meta.isexpr(argument, :kw, 2) ?
            (argument.args[1], argument.args[2], true) : (argument, nothing, false)
        vararg = Meta.isexpr(argument, :..., 1)
        vararg && (argument = argument.args[1])
        argument_name, T = argument isa Symbol ? (argument, nothing) :
            Meta.isexpr(argument, :(::), 2) ? (argument.args[1], argument.args[2]) :
            Meta.isexpr(argument, :(::), 1) ? (nothing, argument.args[1]) :
            throw(ArgumentError("`@index_methods` could not read the argument `$argument`."))
        marker = T === nothing ? nothing : index_marker(T)
        variable = nothing
        if marker === nothing && T isa Symbol && haskey(index_variables, T)
            variable, marker = T, index_variables[T]
        end
        if marker !== nothing
            argument_name === nothing && throw(ArgumentError(
                "Every index argument of `$name` must be named; `::$T` is not."
            ))
            vararg && throw(ArgumentError(
                "The index argument `$argument_name` of `$name` may not be a vararg."
            ))
        elseif expression_mentions(T, variable_names)
            throw(ArgumentError(
                "The type `$T` of an argument of `$name` uses a type variable bounded by an "
                * "index marker, which may annotate index arguments only."
            ))
        end
        (; name=argument_name, type=T, default, has_default, vararg, marker, variable, position)
    end
    index_arguments = [a for a ∈ positional if a.marker !== nothing]
    isempty(index_arguments) && throw(ArgumentError(
        "`$name` has no positional index argument, which is what the methods generated by "
        * "`@index_methods` dispatch on; annotate one with `IndexType`, `IndexRange` or "
        * "`IndexOrRange`."
    ))
    first_index = first(index_arguments).position
    for (variable, marker) ∈ index_variables
        annotated = [a.name for a ∈ positional if a.variable === variable]
        marker === :IndexOrRange && length(annotated) > 1 && throw(ArgumentError(
            "The type variable `$variable` of `$name` is bounded by `IndexOrRange`, and so "
            * "the arguments it annotates ($(join(annotated, ", "))) would all have to be "
            * "indices or all be ranges; annotate each of them with `IndexOrRange` instead."
        ))
    end

    # The optional positional arguments, which Julia requires to come last or just before a
    # vararg, and the numbers of positional arguments of the calls that omit some of them.
    n = length(positional)
    defaulted = [a.position for a ∈ positional if a.has_default]
    if !isempty(defaulted)
        k, q = first(defaulted), last(defaulted)
        (defaulted == k:q && (q == n || (q == n - 1 && positional[n].vararg))) ||
            throw(ArgumentError(
                "The optional positional arguments of `$name` must come last, or just before "
                * "a vararg."
            ))
    end
    shorter = isempty(defaulted) ? (1:0) : (first(defaulted) - 1):(last(defaulted) - 1)
    for a ∈ positional
        if a.has_default && a.position ≤ first_index &&
                expression_mentions(a.default, variable_names)
            throw(ArgumentError(
                "The default `$(a.default)` of an argument of `$name` uses a type variable "
                * "bounded by an index marker, but no index argument precedes that argument, "
                * "so the default is evaluated where the variable does not exist."
            ))
        end
    end

    # The pieces that the methods assemble in different ways.  An anonymous argument needs a
    # name to be passed on.
    as_written(a, argument_name=a.name) = let ex =
            argument_name === nothing ? Expr(:(::), a.type) :
            a.type === nothing ? argument_name : Expr(:(::), argument_name, a.type)
        a.vararg ? Expr(:..., ex) : ex
    end
    typed(a, i) = Expr(:(::), a.name, index_method_types[a.marker][i])
    forwarded_names = [a.name === nothing ? gensym(:argument) : a.name for a ∈ positional]
    forwarded(a) = let n = forwarded_names[a.position]
        a.vararg ? Expr(:..., n) : n
    end
    # What the `where` parameters of a method for the given arguments may be needed by: the
    # types of its arguments, and the function itself, which may be a parameterized type.
    types_of(given) = Any[
        callee; [a.type for a ∈ given if a.marker === nothing && a.type !== nothing]
    ]
    function method(parameters, method_arguments, where_list, returns, method_body)
        call = Expr(:call, callee, parameters..., method_arguments...)
        call = returns === nothing ? call : Expr(:(::), call, returns)
        call = isempty(where_list) ? call : Expr(:where, call, where_list...)
        Expr(:function, call, method_body)
    end
    helper(f) = GlobalRef(@__MODULE__, f)
    # The collected keywords are called `kwargs`, as they would be written by hand, unless
    # that is the name of a positional argument.  The same expression collects them in a
    # signature and passes them on in a call.
    kwargs = any(a -> a.name === :kwargs, positional) ? gensym(:kwargs) : :kwargs
    collected = isempty(keyword_items) ? Any[] : Any[Expr(:parameters, Expr(:..., kwargs))]
    function refusal(given)
        indices = [a for a ∈ given if a.marker !== nothing]
        index_names = Expr(:tuple, [QuoteNode(a.name) for a ∈ indices]...)
        index_values = Expr(:tuple, [a.name for a ∈ indices]...)
        :($(GlobalRef(Core, :throw))($(helper(:index_argument_error))(
            $shown, $index_names, $index_values, $integer_only, $hint
        )))
    end
    normalized_keyword(keyword, K) =
        :($(helper(:index_keyword))($shown, $(QuoteNode(keyword)), $keyword, $K))
    normalization(keyword, K) = Expr(:(=), keyword, normalized_keyword(keyword, K))

    # 1. The fallback, for the first `m` positional arguments, with the index arguments
    # typed as those of the fallback (`i = 1`), or, where only integers are accepted, as
    # those of the `HalfOddInteger` work method that it replaces (`i = 4`).
    function refusing(m, i)
        given = positional[1:m]
        method(
            collected,
            [a.marker === nothing ? as_written(a) : typed(a, i) for a ∈ given],
            needed_parameters(other_parameters, types_of(given)), nothing,
            Expr(:block, source, refusal(given))
        )
    end
    fallback(m) = refusing(m, 1)

    # The completion of a call that gives no index argument, by the next default.
    function completion(m)
        given = positional[1:m]
        method(
            collected,
            [as_written(a, forwarded_names[a.position]) for a ∈ given],
            needed_parameters(other_parameters, types_of(given)), nothing,
            Expr(:block, source,
                Expr(:call, caller, collected..., forwarded.(given)..., positional[m+1].default)
            )
        )
    end

    # 2. The conversion, for the first `m` positional arguments.
    function conversion(m)
        given = positional[1:m]
        conversion_body = if integer_only
            Expr(:block, source, refusal(given))
        else
            checks = [
                :($(helper(:is_half_odd_index))($(a.name)))
                for a ∈ given if a.marker !== nothing
            ]
            Expr(:block, source,
                Expr(:||,
                    length(checks) == 1 ? checks[1] : Expr(:&&, checks...), refusal(given)
                ),
                Expr(:call, caller, collected...,
                    [
                        a.marker === nothing ? forwarded(a) :
                        :($(helper(:half_odd_index))($(a.name)))
                        for a ∈ given
                    ]...
                )
            )
        end
        method(
            collected,
            [
                a.marker === nothing ? as_written(a, forwarded_names[a.position]) : typed(a, 2)
                for a ∈ given
            ],
            needed_parameters(other_parameters, types_of(given)), nothing, conversion_body
        )
    end

    # 3 and 4. The work methods, with the positional defaults that follow an index argument,
    # the keywords as written, and the body behind a prologue that normalizes the keyword
    # indices.  The prologue is a `let`, so that each keyword is bound once in the body, and
    # a closure there that uses one captures it without boxing.  A keyword default that
    # refers to an earlier keyword index is given that index normalized, by a `let` of its
    # own, which binds it to a generated name that replaces the keyword in a copy of the
    # default.  A `let` that bound the keyword's own name would add a second local variable
    # of that name to the method, and `methods(f)` and the hints of a `MethodError` would
    # list the keyword twice; a generated name contains `#`, and is not listed.
    function work_method(i, K)
        work_parameters = map(where_parameters) do p
            if Meta.isexpr(p, :<:, 2) && haskey(index_variables, p.args[1])
                Expr(:<:, p.args[1], index_method_types[index_variables[p.args[1]]][i])
            else
                p
            end
        end
        work_arguments = map(positional) do a
            ex = a.marker === nothing ? as_written(a) :
                a.variable === nothing ? typed(a, i) : Expr(:(::), a.name, a.variable)
            a.has_default && a.position > first_index ? Expr(:kw, ex, a.default) : ex
        end
        normalized = Symbol[]
        work_keywords = map(keywords) do k
            k.name === nothing && return k.item
            argument = k.marker === nothing ? k.argument : Expr(:(::), k.name, k.annotation)
            item = if k.has_default
                earlier = filter(e -> expression_mentions(k.default, (e,)), normalized)
                renamed = Dict(e => gensym(e) for e ∈ earlier)
                Expr(:kw, argument,
                    isempty(earlier) ? k.default : Expr(:let,
                        Expr(:block,
                            [Expr(:(=), renamed[e], normalized_keyword(e, K)) for e ∈ earlier]...
                        ),
                        replace_symbols(k.default, renamed)
                    )
                )
            else
                argument
            end
            k.marker === nothing || push!(normalized, k.name)
            item
        end
        work_body = if isempty(index_keywords)
            body
        else
            Expr(:block, line,
                Expr(:let,
                    Expr(:block, [normalization(k.name, K) for k ∈ index_keywords]...),
                    body
                )
            )
        end
        method(
            isempty(keyword_items) ? Any[] : Any[Expr(:parameters, work_keywords...)],
            work_arguments, work_parameters, return_type, work_body
        )
    end

    # Each method is a separate copy, because macro expansion of the bodies may happen in
    # place, and the body written is shared by the two work methods.  Where only integers
    # are accepted, the `HalfOddInteger` work method is replaced by one that refuses, rather
    # than omitted: another definition of `f` whose other arguments are less specific would
    # otherwise leave its own `HalfOddInteger` method ambiguous with this conversion method.
    methods = Any[
        Expr(:macrocall,
            GlobalRef(Base, Symbol("@__doc__")), source, deepcopy(rewrap(fallback(n)))
        ),
    ]
    for m ∈ shorter
        push!(methods, deepcopy(rewrap(m < first_index ? completion(m) : fallback(m))))
    end
    for m ∈ [n; [m for m ∈ shorter if m ≥ first_index]]
        push!(methods, deepcopy(rewrap(conversion(m))))
        integer_only && push!(methods, deepcopy(rewrap(refusing(m, 4))))
    end
    push!(methods, deepcopy(rewrap(work_method(3, Int))))
    integer_only || push!(methods, deepcopy(rewrap(work_method(4, HalfOddInteger))))
    esc(Expr(:block, methods...))
end
