# Angular-momentum operators, as matrices acting on mode weights and as maps from one set of
# mode weights to another.
#
# Each operator is a value: one of twelve zero-size singletons, each of its own subtype of
# `DifferentialOperator`.  Called with indices, in either of the shapes
#
#     op(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])   or   op(s, ℓₘₐₓ, [T=Float64])  (with ℓₘᵢₙ = abs(s))
#
# it returns a (sparse) matrix that acts on a vector of mode weights ordered as
#
#     [ f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ ]
#
# with spin weight `s`.  Entries with ℓ < |s| are zero.  The conventions are those of the
# "Conventions" section of the documentation: L and R are the left and right Lie derivatives
#
#     L_𝐮 f(𝐑) =  i d/dϵ f(e^{-ϵ𝐮/2} 𝐑),        R_𝐮 f(𝐑) = -i d/dϵ f(𝐑 e^{-ϵ𝐮/2}),
#
# with L_± = L_x ± i L_y and R_± = R_x ± i R_y, so that [L_z, L_±] = ±L_± and [R_z, R_±] =
# ±R_±.  Spin-weighted functions satisfy R_z η = s η, and ð = R_+, ð̄ = -R_-.
#
# The two call shapes are written once, for all twelve operators, with `@index_methods`, so
# that their bodies see three indices of one kind — `Int`s or `HalfOddInteger`s — and a call
# that mixes the kinds, or passes an index of another type, is refused with a message naming
# the operator as it prints (`L²`, `ð`).  The bodies are generic over the two kinds, because
# every quantity they form — ℓ ± m, ℓ ± s, 2ℓ — is an `Int` for either; the two exceptions,
# ℓ(ℓ+1) and the conversion of an index to the matrix element type, go through
# `casimir_eigenvalue` below and `index_value` in `half_odd_integer.jl`.
#
# What distinguishes one operator from another is a handful of traits — the change it makes
# to the spin weight (`Δspin`), the band of the matrix it occupies (`bandstructure`), the
# element type of that matrix (`coefftype`) and the ℓ below which it vanishes (`support_ℓ`)
# — together with its matrix elements, the coefficient functions below.  The coefficients
# are evaluated both by `operator_matrix`, which builds the matrix, and by
# `apply_operator!`, which applies the operator to the storage of a `ModeWeights` without
# building one; `op * w`, `op(w)` and `mul!(w′, op, w)` are written with the latter.

const operator_signature_note = """
The argument `ℓₘᵢₙ` may be omitted, in which case it defaults to `abs(s)`.  The result acts
on a vector of mode weights ordered as `[f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]`; any
entries with ``ℓ < |s|`` are mapped to zero.  The indices `s`, `ℓₘᵢₙ` and `ℓₘₐₓ` may be
integers of type `Int` or half-odd-integers, the latter given as [`HalfOddInteger`](@ref)s
or as `Rational{Int}`s with denominator 2 — as in `L²(1//2, 7//2)` — in which case every
``ℓ`` and ``m`` of the ordering is a half-odd-integer.  The indices in one call must all be
of one kind; a call that mixes them, such as `L²(1//2, 0, 7//2)`, is refused, as is an index
of any other type (see [`IndexType`](@ref)).  The matrices of operators that change the spin
weight can be multiplied together only when they share the range of ``ℓ``, so `ℓₘᵢₙ` must
then be given explicitly, as in `ð̄(1, 0, ℓₘₐₓ) * ð(0, 0, ℓₘₐₓ)`.

The argument `T` is the real floating-point type in which the entries are computed; it
defaults to `Float64`.  The matrix is returned as a `Diagonal`, `Bidiagonal` or
`Tridiagonal` matrix of element type `T`, except that the entries of [`Ly`](@ref) are of
type `Complex{T}`.

Applied to a [`ModeWeights`](@ref) `w` instead, as `op(w)` or `op * w`, the operator returns
new mode weights over the same range of ``ℓ``, with the spin weight changed by
[`Δspin`](@ref), and it does so without building the matrix; `mul!(w′, op, w)` writes the
result into `w′`.  A product such as `ð̄ * ð * w` applies the operators from right to left.
See [`DifferentialOperator`](@ref) for the details.
"""

# The docstrings below are `raw` strings, because they are full of LaTeX backslashes, and a
# `raw` string does not interpolate, so `$(operator_signature_note)` written inside one is
# literal text.  This replaces that text with the note; `@doc` accepts any expression that
# evaluates to a string.
splice_signature_note(s) =
    replace(s, "\$(operator_signature_note)" => operator_signature_note)

# The eigenvalue ℓ(ℓ+1) of L² and R², as the value the typed comprehension in those
# functions converts to `T`.  For an integer ℓ it is the integer ℓ(ℓ+1).  For a half-odd ℓ
# the product is a quarter-integer, which the index type cannot hold and which must not be
# formed as a `Rational`; instead the `Int` (2ℓ)(2ℓ+2) is converted to `T` and divided by 4.
# Division by 4 is exact in every binary floating-point type, so this path is as exact as
# the integer one.
@inline casimir_eigenvalue(::Type{T}, ℓ::Integer) where {T} = ℓ*(ℓ+1)
@inline casimir_eigenvalue(::Type{T}, ℓ::HalfOddInteger) where {T} = T((2ℓ)*(2ℓ+2)) / 4


### The operators as objects, and their matrix elements.
#
# Each operator is a zero-size singleton, so that dispatching on it costs nothing and every
# trait below folds away at compile time.  The twelve exported names are instances of these
# types.
#
# The matrix elements are written *once*, here, and evaluated both by the comprehensions
# that build the operator matrices below and by the loops that apply an operator to a
# `ModeWeights` without building one.  That is what makes the two agree: not a coincidence
# to be tested for, but the same expression evaluated twice.  Note in particular that the `ℓ
# < …` mask lives *inside* the coefficient, so that neither caller can forget it.

"""
    DifferentialOperator

The abstract type of the twelve angular-momentum operators, each of which is a singleton
instance of its own subtype.

The operators are [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`Lx`](@ref) and
[`Ly`](@ref), built from the left Lie derivative, and [`R²`](@ref), [`Rz`](@ref),
[`R₊`](@ref), [`R₋`](@ref), [`ð`](@ref) and [`ð̄`](@ref), built from the right one.  Each of
them `op` can be used in these ways:

  * `op(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])` and `op(s, ℓₘₐₓ, [T=Float64])` return the matrix of the
    operator acting on mode weights of spin weight `s`, as described in the docstring of
    each operator.
  * `op(w)` and `op * w`, for a [`ModeWeights`](@ref) `w`, return the mode weights of the
    result, without building the matrix.  They are labelled with the spin weight `spin(w) +
    Δspin(op)` (see [`Δspin`](@ref)) and with the range of ``ℓ`` of `w`; modes below the new
    ``|s|`` are zero rather than dropped.  The element type is that of the product of the
    matrix entries and the weights, so that [`Ly`](@ref) applied to real weights gives
    complex ones.
  * `mul!(w′, op, w)` writes the same result into `w′`, which must be labelled with that
    spin weight and range of ``ℓ`` and must not share storage with `w`, since the banded
    operators read neighboring modes.  `w′` may also be a bare vector at least as long as
    the result, in which case the result is returned as a `ModeWeights` over its first
    entries.
  * `a * b * w` for operators `a` and `b` is `a * (b * w)`, so a product of operators is
    applied from right to left, as the notation means.  No product of two operators on their
    own is defined.
  * In a broadcast, an operator is a scalar, so that `ð .* ws` applies `ð` to each element
    of a collection `ws` of mode weights, as `ð.(ws)` does.

For finite data, `op * w` agrees with the product of the matrix and the weights under `==`.
Where the weights are not finite the two can differ, because the stored zero diagonal of a
banded matrix turns an infinite weight into a `NaN` in the row of that mode, which `op * w`
does not.

This type is provided for dispatch.  It is not an extension point: the behavior of each
operator is set by internal traits and coefficient functions, which are not part of the
public interface, so a new subtype defined elsewhere would not work.
"""
abstract type DifferentialOperator end

struct Casimir       <: DifferentialOperator end
struct RightCasimir  <: DifferentialOperator end
struct LeftZ         <: DifferentialOperator end
struct LeftRaising   <: DifferentialOperator end
struct LeftLowering  <: DifferentialOperator end
struct LeftX         <: DifferentialOperator end
struct LeftY         <: DifferentialOperator end
struct RightZ        <: DifferentialOperator end
struct RightRaising  <: DifferentialOperator end
struct RightLowering <: DifferentialOperator end
struct SpinRaising   <: DifferentialOperator end
struct SpinLowering  <: DifferentialOperator end

"""
    Δspin(op)
    Deltaspin(op)

The change in spin weight produced by the operator `op`: `1` for [`ð`](@ref) and
[`R₊`](@ref), `-1` for [`ð̄`](@ref) and [`R₋`](@ref), and `0` for the other operators.

This is the amount by which the spin weight of `op * w` differs from that of the mode
weights `w`, and so it is what a destination for `mul!` must be labelled with:

```julia
w′ = ModeWeights(
    similar(parent(w)), spin(w) + Δspin(ð),
    SphericalFunctions.ℓₘᵢₙ(w), SphericalFunctions.ℓₘₐₓ(w)
)
mul!(w′, ð, w)
```

The ASCII alias `Deltaspin` may be used in place of `Δspin`.  Neither name is exported.
"""
function Δspin end
@inline Δspin(::DifferentialOperator) = 0
@inline Δspin(::Union{RightRaising,  SpinRaising})  =  1
@inline Δspin(::Union{RightLowering, SpinLowering}) = -1

# Which band of the matrix the operator occupies, and hence which builder and which kernel
# apply.  A sub-diagonal entry takes `f[ℓ, m-1]` into `out[ℓ, m]`; a super-diagonal one
# takes `f[ℓ, m+1]`.
abstract type BandStructure end
struct DiagonalBand      <: BandStructure end
struct SubdiagonalBand   <: BandStructure end
struct SuperdiagonalBand <: BandStructure end
struct TridiagonalBand   <: BandStructure end

@inline bandstructure(::DifferentialOperator) = DiagonalBand()
@inline bandstructure(::LeftRaising)  = SubdiagonalBand()
@inline bandstructure(::LeftLowering) = SuperdiagonalBand()
@inline bandstructure(::Union{LeftX, LeftY}) = TridiagonalBand()

# The element type of the matrix, given the real type it was asked for.  `Ly` is the only
# one whose entries are complex — they are ∓i/2 times those of `L₊` and `L₋`.
@inline coefftype(::DifferentialOperator, ::Type{T}) where {T} = T
@inline coefftype(::LeftY, ::Type{T}) where {T} = Complex{T}

# Each operator says its own name: `nameof` because code (and tests) reach for it, and
# `show` so that an operator prints as `ð` rather than as
# `SphericalFunctions.SpinRaising()`, in particular in the messages of the errors that name
# it.
Base.nameof(::Casimir)       = :L²
Base.nameof(::RightCasimir)  = :R²
Base.nameof(::LeftZ)         = :Lz
Base.nameof(::LeftRaising)   = :L₊
Base.nameof(::LeftLowering)  = :L₋
Base.nameof(::LeftX)         = :Lx
Base.nameof(::LeftY)         = :Ly
Base.nameof(::RightZ)        = :Rz
Base.nameof(::RightRaising)  = :R₊
Base.nameof(::RightLowering) = :R₋
Base.nameof(::SpinRaising)   = :ð
Base.nameof(::SpinLowering)  = :ð̄
Base.show(io::IO, op::DifferentialOperator) = print(io, nameof(op))

# An operator is a single value, so a broadcast such as `ð .* ws` treats it as a scalar
# rather than trying to iterate over it.
Base.broadcastable(op::DifferentialOperator) = Ref(op)

# A product of operators applied to mode weights, `ð̄ * ð * w`, is parsed as `*(ð̄, ð, w)`,
# which `Base` would fold from the left into `(ð̄ * ð) * w`.  No product of two operators is
# defined, so the fold is taken from the right instead, as `ð̄ * (ð * w)`, which applies the
# operators in the order the notation means.  Longer products recurse through this method,
# as `L₊ * L₋ * Lz * w` is `L₊ * (L₋ * (Lz * w))`.
Base.:*(a::DifferentialOperator, b::DifferentialOperator, c, xs...) = a * *(b, c, xs...)

# The ℓ below which the result vanishes.  For the spin-changing operators the cutoff is set
# by the *output* spin weight, which is the `s′` of the matrix builders.
@inline support_ℓ(::DifferentialOperator, s) = abs(s)
@inline support_ℓ(::Union{RightRaising,  SpinRaising},  s) = max(abs(s), abs(s + 1))
@inline support_ℓ(::Union{RightLowering, SpinLowering}, s) = max(abs(s), abs(s - 1))

# `casimir_eigenvalue` and `index_value` return an `Int` for an integer index, so the
# coefficients convert to `T` explicitly; without that, the loops that share them would be
# type-unstable on the integer path.
@inline function diagonal_coefficient(op::Casimir, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : convert(T, casimir_eigenvalue(T, ℓ))
end
@inline diagonal_coefficient(::RightCasimir, ::Type{T}, s, ℓ, m) where {T} =
    diagonal_coefficient(Casimir(), T, s, ℓ, m)
@inline function diagonal_coefficient(op::LeftZ, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : convert(T, index_value(T, m))
end
@inline function diagonal_coefficient(op::RightZ, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : convert(T, index_value(T, s))
end
@inline function diagonal_coefficient(op::RightRaising, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : √T((ℓ-s)*(ℓ+s+1))
end
@inline function diagonal_coefficient(op::RightLowering, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : √T((ℓ+s)*(ℓ-s+1))
end
@inline diagonal_coefficient(::SpinRaising, ::Type{T}, s, ℓ, m) where {T} =
    diagonal_coefficient(RightRaising(), T, s, ℓ, m)
# The negation is applied to the *masked* value, as `-R₋` does, so that the vanishing
# entries are `-0.0` here too and even `isequal` agrees with the matrix.
@inline diagonal_coefficient(::SpinLowering, ::Type{T}, s, ℓ, m) where {T} =
    -diagonal_coefficient(RightLowering(), T, s, ℓ, m)

# The coefficients a_m = √((ℓ-m+1)(ℓ+m)) and b_m = √((ℓ+m+1)(ℓ-m)) of the ladder operators,
# in the real type `T`.  Both vanish exactly where the neighbor they multiply lies outside
# -ℓ:ℓ, and are then given as zero rather than as the square root of zero, whose derivative
# is infinite: when `T` is a dual number, the zero partials of the constant would be
# multiplied by that infinity, and give `NaN`.  They are computed in `recurrence_type(T)`,
# the floating-point type underneath any dual numbers (see `src/derivatives/lifting.jl`),
# since they are constants.  The derivatives in `src/derivatives/kernels.jl` use them too.
@inline function ladder_coefficient(n::Int, ::Type{T}) where {T}
    let F = recurrence_type(T)
        n == 0 ? zero(F) : √F(n)
    end
end
@inline ladder_down(ℓ, m, ::Type{T}) where {T} = ladder_coefficient(Int(ℓ - m + 1) * Int(ℓ + m), T)
@inline ladder_up(ℓ, m, ::Type{T}) where {T} = ladder_coefficient(Int(ℓ + m + 1) * Int(ℓ - m), T)

# The matrix elements of L₊ and L₋, indexed by the *output* mode `(ℓ, m)`, are a_m and b_m.
# Both vanish exactly at the edge of their ℓ block — at `m = -ℓ` for the raising one and at
# `m = +ℓ` for the lowering one — which is what keeps an ℓ block from coupling to its
# neighbors.
@inline function subdiagonal_coefficient(op::LeftRaising, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : convert(T, ladder_down(ℓ, m, T))
end
@inline function superdiagonal_coefficient(op::LeftLowering, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : convert(T, ladder_up(ℓ, m, T))
end
@inline subdiagonal_coefficient(::LeftX, ::Type{T}, s, ℓ, m) where {T} =
    subdiagonal_coefficient(LeftRaising(), T, s, ℓ, m) / 2
@inline superdiagonal_coefficient(::LeftX, ::Type{T}, s, ℓ, m) where {T} =
    superdiagonal_coefficient(LeftLowering(), T, s, ℓ, m) / 2
# The entries of `Ly` are built as complex numbers with a zero real part, rather than by
# dividing those of `L₊` and `L₋` by `2im`: that is a complex division, which for some float
# types does not give exactly `∓i x/2`, whereas this form does, so that `2im .* Ly` equals
# `L₊ - L₋` exactly.
@inline subdiagonal_coefficient(::LeftY, ::Type{T}, s, ℓ, m) where {T} =
    Complex{T}(zero(T), -subdiagonal_coefficient(LeftRaising(), T, s, ℓ, m) / 2)
@inline superdiagonal_coefficient(::LeftY, ::Type{T}, s, ℓ, m) where {T} =
    Complex{T}(zero(T), superdiagonal_coefficient(LeftLowering(), T, s, ℓ, m) / 2)


### The two call shapes, written once for every operator.
#
# The first builds the matrix from the band structure and the coefficients above; the second
# supplies the default ℓₘᵢₙ = |s|, computed from the spin weight after it has been
# converted.

@index_methods function (op::DifferentialOperator)(
    s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, ::Type{T}=Float64
) where {IT<:IndexType, T}
    # The range of ℓ is validated here, for every band structure alike: the diagonal builder
    # never calls `Ysize`, and would otherwise read a negative ℓₘᵢₙ as 0 and an inverted
    # range as an empty one.
    Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    operator_matrix(op, bandstructure(op), s, ℓₘᵢₙ, ℓₘₐₓ, T)
end
@index_methods (op::DifferentialOperator)(
    s::IndexType, ℓₘₐₓ::IndexType, ::Type{T}=Float64
) where {T} = op(s, abs(s), ℓₘₐₓ, T)

# One builder per band structure.  The `ifelse` in the ladder ranges drops the one mode that
# has no band entry: the very first for a sub-diagonal, the very last for a super-diagonal.
function operator_matrix(op, ::DiagonalBand, s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, ::Type{T}) where {IT<:IntegerHalf, T}
    Diagonal(
        coefftype(op, T)[
            diagonal_coefficient(op, T, s, ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ
        ]
    )
end
function operator_matrix(op, ::SubdiagonalBand, s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, ::Type{T}) where {IT<:IntegerHalf, T}
    Bidiagonal(
        zeros(coefftype(op, T), Ysize(ℓₘᵢₙ, ℓₘₐₓ)),
        coefftype(op, T)[
            subdiagonal_coefficient(op, T, s, ℓ, m)
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ ifelse(ℓ==ℓₘᵢₙ,-ℓ+1,-ℓ):ℓ
        ],
        :L
    )
end
function operator_matrix(op, ::SuperdiagonalBand, s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, ::Type{T}) where {IT<:IntegerHalf, T}
    Bidiagonal(
        zeros(coefftype(op, T), Ysize(ℓₘᵢₙ, ℓₘₐₓ)),
        coefftype(op, T)[
            superdiagonal_coefficient(op, T, s, ℓ, m)
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ifelse(ℓ==ℓₘₐₓ,ℓ-1,ℓ)
        ],
        :U
    )
end
function operator_matrix(op, ::TridiagonalBand, s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, ::Type{T}) where {IT<:IntegerHalf, T}
    CT = coefftype(op, T)
    Tridiagonal(
        CT[
            subdiagonal_coefficient(op, T, s, ℓ, m)
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ ifelse(ℓ==ℓₘᵢₙ,-ℓ+1,-ℓ):ℓ
        ],
        zeros(CT, Ysize(ℓₘᵢₙ, ℓₘₐₓ)),
        CT[
            superdiagonal_coefficient(op, T, s, ℓ, m)
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ifelse(ℓ==ℓₘₐₓ,ℓ-1,ℓ)
        ],
    )
end


### Applying an operator without building its matrix.
#
# The matrix builders above and the loops below evaluate the *same* coefficient functions,
# so for finite data the two agree under `==` rather than merely to within rounding, which
# is what the tests of `op * w` assert.  Reproducing that exactly is why these loops are
# written plainly: no `@simd`, no `@fastmath`, no `muladd`, and the two terms of a
# tridiagonal row summed left to right, as `LinearAlgebra`'s own `l[i-1]*b₋ + d[i]*b₀ +
# u[i]*b₊` does with `d` identically zero.  The two differ in what that zero diagonal does:
# the matrix multiplies it by the weight, which can give `-0.0` where the loop gives `0.0`,
# and `NaN` where the weight is infinite or `NaN`, whereas the loop never reads the diagonal
# at all.
#
# The ladder coefficients vanish *exactly* at the edge of each ℓ block — at `m = -ℓ` for
# the raising one and at `m = +ℓ` for the lowering one — so no block ever couples to its
# neighbor and the loops need no per-block special case.  Only the very first and very last
# position in the whole vector need a branch, because there the matrix has no band entry at
# all.

function apply_operator!(
    out, op, ::DiagonalBand, in, s::IT, ℓ₀::IT, ℓ₁::IT, ::Type{T}
) where {IT<:IntegerHalf, T}
    @inbounds for ℓ ∈ ℓ₀:ℓ₁
        i = Yindex(ℓ, -ℓ, ℓ₀)
        for m ∈ -ℓ:ℓ
            out[i] = diagonal_coefficient(op, T, s, ℓ, m) * in[i]
            i += 1
        end
    end
    out
end

function apply_operator!(
    out, op, ::SubdiagonalBand, in, s::IT, ℓ₀::IT, ℓ₁::IT, ::Type{T}
) where {IT<:IntegerHalf, T}
    Z = zero(eltype(out))
    @inbounds for ℓ ∈ ℓ₀:ℓ₁
        i = Yindex(ℓ, -ℓ, ℓ₀)
        for m ∈ -ℓ:ℓ
            out[i] = i == 1 ? Z : subdiagonal_coefficient(op, T, s, ℓ, m) * in[i-1]
            i += 1
        end
    end
    out
end

function apply_operator!(
    out, op, ::SuperdiagonalBand, in, s::IT, ℓ₀::IT, ℓ₁::IT, ::Type{T}
) where {IT<:IntegerHalf, T}
    N = length(out)
    Z = zero(eltype(out))
    @inbounds for ℓ ∈ ℓ₀:ℓ₁
        i = Yindex(ℓ, -ℓ, ℓ₀)
        for m ∈ -ℓ:ℓ
            out[i] = i == N ? Z : superdiagonal_coefficient(op, T, s, ℓ, m) * in[i+1]
            i += 1
        end
    end
    out
end

function apply_operator!(
    out, op, ::TridiagonalBand, in, s::IT, ℓ₀::IT, ℓ₁::IT, ::Type{T}
) where {IT<:IntegerHalf, T}
    N = length(out)
    Z = zero(eltype(out))
    @inbounds for ℓ ∈ ℓ₀:ℓ₁
        i = Yindex(ℓ, -ℓ, ℓ₀)
        for m ∈ -ℓ:ℓ
            lo = i == 1 ? Z : subdiagonal_coefficient(op, T, s, ℓ, m)   * in[i-1]
            hi = i == N ? Z : superdiagonal_coefficient(op, T, s, ℓ, m) * in[i+1]
            out[i] = lo + hi
            i += 1
        end
    end
    out
end


### Operators on mode weights
#
# One method covers all twelve: the operator is a value, so it says its own effect on the
# spin weight through `Δspin`, and the container already holds three indices of one kind,
# `Int` or `HalfOddInteger`.  The range of ℓ is unchanged even where the spin weight moves.
# The entries of the result below the new |s| belong to no harmonic; each is the product of
# an entry of the input with a coefficient that vanishes there, so it is zero where the
# input is finite, and is not dropped.  `ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` is what changes the
# range, and drops those entries.
function Base.:*(op::DifferentialOperator, w::ModeWeights{T}) where {T}
    # The result is allocated at the length of the input's storage, so this one check covers
    # both of the vectors that the kernel indexes.
    check_storage_length(w)
    Treal = real(float(T))
    out = similar(w.data, Base.promote_op(*, coefftype(op, Treal), T))
    apply_operator!(out, op, bandstructure(op), w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ, Treal)
    ModeWeights(out, w.s + Δspin(op), w.ℓₘᵢₙ, w.ℓₘₐₓ)
end
(op::DifferentialOperator)(w::ModeWeights) = op * w

# The in-place form, for a loop over many sets of weights.  Aliasing is refused for the
# banded operators, whose kernels read a neighbor that an in-place write may already have
# clobbered; it would be safe for the diagonal ones, but allowing it there only would be a
# trap.
function LinearAlgebra.mul!(
    w′::ModeWeights, op::DifferentialOperator, w::ModeWeights{T}
) where {T}
    if spin(w′) != w.s + Δspin(op) || ℓₘᵢₙ(w′) != w.ℓₘᵢₙ || ℓₘₐₓ(w′) != w.ℓₘₐₓ
        throw(operator_output_error(w′, op, w))
    end
    check_storage_length(w)
    check_storage_length(w′)
    if Base.mightalias(w′.data, w.data)
        throw(ArgumentError(
            "The output aliases the input.  $(nameof(op)) reads neighboring modes, so it "
            * "cannot be applied in place; pass a separate destination, such as `similar(w)`."
        ))
    end
    apply_operator!(
        w′.data, op, bandstructure(op), w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ, real(float(T))
    )
    w′
end
# The operators keep the range of ℓ of their input, so a destination allocated with the
# default range of its own spin weight, `abs(s′):ℓₘₐₓ`, is refused whenever |s′| ≠ |s|; the
# message says how to allocate one, and how to change the range of the result afterwards.
@noinline function operator_output_error(w′::ModeWeights, op, w::ModeWeights{T}) where {T}
    s′ = w.s + Δspin(op)
    ArgumentError(
        "The output has s=$(spin(w′)) and ℓ ∈ $(ℓₘᵢₙ(w′)):$(ℓₘₐₓ(w′)), but $(nameof(op)) "
        * "applied to these weights gives s=$s′ and ℓ ∈ $(w.ℓₘᵢₙ):$(w.ℓₘₐₓ), since the "
        * "operators keep the range of ℓ of their input.  Allocate the output with "
        * "`ModeWeights{$T}(undef, $s′, $(w.ℓₘᵢₙ), $(w.ℓₘₐₓ))`, and use "
        * "`ModeWeights(w′; ℓₘᵢₙ, ℓₘₐₓ)` to copy the result into another range of ℓ."
    )
end

# Bare storage, at least as long as the result, is accepted as the output too, and the
# result comes back labelled, as a `ModeWeights` over it (see `mode_weights_view`).
function LinearAlgebra.mul!(w′::AbstractVector, op::DifferentialOperator, w::ModeWeights)
    mul!(mode_weights_view(w′, w.s + Δspin(op), w.ℓₘᵢₙ, w.ℓₘₐₓ), op, w)
end

@doc splice_signature_note(raw"""
    L²(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    L²(s, ℓₘₐₓ, [T=Float64])
    L² * w

Compute the total angular-momentum operator (the Casimir operator) for spin weight `s`.

This is the standard ``L^2`` operator, familiar from basic physics, extended to work with
SWSHs.  It is equal to
```math
L^2 = L_x^2 + L_y^2 + L_z^2 = \frac{L_+L_- + L_-L_+ + 2L_zL_z}{2}.
```
Note that these are the left Lie derivatives, but ``L^2 = R^2``, where ``R`` is the right
Lie derivative.  See the [conventions summary](@ref summary_L_R_definitions) or
[Boyle](@cite Boyle_2016) for more details.

In terms of the SWSHs, we can write the action of ``L^2`` as
```math
L^2 {}_{s}Y_{ℓ,m} = ℓ\,(ℓ+1) {}_{s}Y_{ℓ,m}.
```

$(operator_signature_note)

The ASCII alias `L2` may be used in place of `L²`; it is public but not exported.

See also [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const L² = Casimir()


@doc splice_signature_note(raw"""
    Lz(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    Lz(s, ℓₘₐₓ, [T=Float64])
    Lz * w

Compute the angular-momentum operator associated with the ``z`` direction.  This is the
standard ``L_z`` operator, familiar from basic physics, extended to work with SWSHs.  Note
that this is the left Lie derivative; see [`Rz`](@ref) for the equivalent right Lie
derivative.  See the [conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite
Boyle_2016) for more details.

In terms of the SWSHs, we can write the action of ``L_z`` as
```math
L_z {}_{s}Y_{ℓ,m} = m\, {}_{s}Y_{ℓ,m}.
```

$(operator_signature_note)

See also [`L²`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const Lz = LeftZ()


@doc splice_signature_note(raw"""
    L₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    L₊(s, ℓₘₐₓ, [T=Float64])
    L₊ * w

Compute the angular-momentum raising operator.  This is the standard ``L_+`` operator,
familiar from basic physics, extended to work with SWSHs.  Note that this is the left Lie
derivative; see [`R₊`](@ref) for the equivalent right Lie derivative.  See the [conventions
summary](@ref summary_L_R_definitions) or [Boyle](@cite Boyle_2016) for more details.

We define ``L_+`` to be the raising operator for the left Lie derivative with respect to
rotation about ``z``: ``L_z``.  By definition, this implies the commutator relation ``[L_z,
L_+] = L_+``, which allows us to derive ``L_+ = L_x + i\, L_y.``

In terms of the SWSHs, we can write the action of ``L_+`` as
```math
L_+ {}_{s}Y_{ℓ,m} = \sqrt{(ℓ-m)(ℓ+m+1)}\, {}_{s}Y_{ℓ,m+1}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{L_+(f)\right\}_{s,ℓ,m} = \sqrt{(ℓ+m)(ℓ-m+1)}\,\left\{f\right\}_{s,ℓ,m-1}.
```

$(operator_signature_note)

The ASCII alias `Lplus` may be used in place of `L₊`; it is public but not exported.

See also [`L²`](@ref), [`Lz`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const L₊ = LeftRaising()


@doc splice_signature_note(raw"""
    L₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    L₋(s, ℓₘₐₓ, [T=Float64])
    L₋ * w

Compute the angular-momentum lowering operator.  This is the standard ``L_-`` operator,
familiar from basic physics, extended to work with SWSHs.  Note that this is the left Lie
derivative; see [`R₋`](@ref) for the equivalent right Lie derivative.  See the [conventions
summary](@ref summary_L_R_definitions) or [Boyle](@cite Boyle_2016) for more details.

We define ``L_-`` to be the lowering operator for the left Lie derivative with respect to
rotation about ``z``: ``L_z``.  By definition, this implies the commutator relation ``[L_z,
L_-] = -L_-``, which allows us to derive ``L_- = L_x - i\, L_y.``

In terms of the SWSHs, we can write the action of ``L_-`` as
```math
L_- {}_{s}Y_{ℓ,m} = \sqrt{(ℓ+m)(ℓ-m+1)}\, {}_{s}Y_{ℓ,m-1}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{L_-(f)\right\}_{s,ℓ,m} = \sqrt{(ℓ-m)(ℓ+m+1)}\,\left\{f\right\}_{s,ℓ,m+1}.
```

$(operator_signature_note)

The ASCII alias `Lminus` may be used in place of `L₋`; it is public but not exported.

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const L₋ = LeftLowering()


@doc splice_signature_note(raw"""
    Lx(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    Lx(s, ℓₘₐₓ, [T=Float64])
    Lx * w

Compute the ``x`` component of the left angular-momentum operator, ``L_x = (L_+ + L_-)/2``.

This is the standard ``L_x`` operator, familiar from basic physics, extended to work with
SWSHs.  See the [conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite
Boyle_2016) for more details.  The matrix is real, symmetric and tridiagonal.

$(operator_signature_note)

Note that there are no corresponding `Rx` and `Ry` functions.  The right raising and
lowering operators change the spin weight, so ``R_x = (R_+ + R_-)/2`` would map a function
of spin weight ``s`` to a sum of a function of spin weight ``s+1`` and one of spin weight
``s-1``.  That is not a linear map on the mode weights of a single spin weight, and so it
cannot be represented by a matrix of the kind these functions return.

See also [`Ly`](@ref), [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref),
[`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const Lx = LeftX()


@doc splice_signature_note(raw"""
    Ly(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    Ly(s, ℓₘₐₓ, [T=Float64])
    Ly * w

Compute the ``y`` component of the left angular-momentum operator, ``L_y = (L_+ -
L_-)/(2i)``.

This is the standard ``L_y`` operator, familiar from basic physics, extended to work with
SWSHs.  See the [conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite
Boyle_2016) for more details.  The matrix is tridiagonal and purely imaginary, so its
element type is `Complex{T}` rather than `T`.

$(operator_signature_note)

There are no corresponding `Rx` and `Ry` functions; see [`Lx`](@ref) for why.

See also [`Lx`](@ref), [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref),
[`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const Ly = LeftY()


@doc splice_signature_note(raw"""
    R²(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    R²(s, ℓₘₐₓ, [T=Float64])
    R² * w

Compute the total angular-momentum operator (the Casimir operator) for spin weight `s`, in
terms of the right Lie derivative.

This is the ``R^2`` operator, much like the ``L^2`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  It is equal to
```math
R^2 = R_x^2 + R_y^2 + R_z^2 = \frac{R_+R_- + R_-R_+ + 2R_zR_z}{2}.
```
Note that these are the right Lie derivatives, but ``L^2 = R^2``, where ``L`` is the left
Lie derivative, so this function returns the same matrix as [`L²`](@ref).  See the
[conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite Boyle_2016) for more
details.

In terms of the SWSHs, we can write the action of ``R^2`` as
```math
R^2 {}_{s}Y_{ℓ,m} = ℓ\,(ℓ+1) {}_{s}Y_{ℓ,m}.
```

$(operator_signature_note)

The ASCII alias `R2` may be used in place of `R²`; it is public but not exported.

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const R² = RightCasimir()


@doc splice_signature_note(raw"""
    Rz(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    Rz(s, ℓₘₐₓ, [T=Float64])
    Rz * w

Compute the *right* angular-momentum operator associated with the ``z`` direction.

This is the ``R_z`` operator, much like the ``L_z`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  See [`Lz`](@ref)
for the equivalent left Lie derivative.  See the [conventions summary](@ref
summary_spin_weight) or [Boyle](@cite Boyle_2016) for more details.

Spin-weighted functions are precisely the eigenfunctions of ``R_z``, with eigenvalue equal
to the spin weight.  In particular,
```math
R_z {}_{s}Y_{ℓ,m} = s\, {}_{s}Y_{ℓ,m}.
```

$(operator_signature_note)

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const Rz = RightZ()


@doc splice_signature_note(raw"""
    R₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    R₊(s, ℓₘₐₓ, [T=Float64])
    R₊ * w

Compute the *right* angular-momentum raising operator.

This is the ``R_+`` operator, much like the ``L_+`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  See [`L₊`](@ref)
for the equivalent left Lie derivative.  See the [conventions summary](@ref
summary_L_R_definitions) or [Boyle](@cite Boyle_2016) for more details.

We define ``R_+`` to be the raising operator for the right Lie derivative with respect to
rotation about ``z``: ``R_z``.  By definition, this implies the commutator relation ``[R_z,
R_+] = R_+``, which allows us to derive ``R_+ = R_x + i\, R_y.``  Because the eigenvalue of
``R_z`` on a spin-weighted function is the spin weight ``s``, this operator raises the spin
weight by one; it is identical to the spin-raising operator [`ð`](@ref).

In terms of the SWSHs, we can write the action of ``R_+`` as
```math
R_+ {}_{s}Y_{ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓ,m}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{R_+(f)\right\}_{s+1,ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\,\left\{f\right\}_{s,ℓ,m},
```
where the argument `s` of this function is the spin weight of the *input*.

$(operator_signature_note)

The ASCII alias `Rplus` may be used in place of `R₊`; it is public but not exported.

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const R₊ = RightRaising()


@doc splice_signature_note(raw"""
    R₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    R₋(s, ℓₘₐₓ, [T=Float64])
    R₋ * w

Compute the *right* angular-momentum lowering operator.

This is the ``R_-`` operator, much like the ``L_-`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  See [`L₋`](@ref)
for the equivalent left Lie derivative.  See the [conventions summary](@ref
summary_L_R_definitions) or [Boyle](@cite Boyle_2016) for more details.

We define ``R_-`` to be the lowering operator for the right Lie derivative with respect to
rotation about ``z``: ``R_z``.  By definition, this implies the commutator relation ``[R_z,
R_-] = -R_-``, which allows us to derive ``R_- = R_x - i\, R_y.``  Because the eigenvalue of
``R_z`` on a spin-weighted function is the spin weight ``s``, this operator lowers the spin
weight by one; it is the negative of the spin-lowering operator [`ð̄`](@ref).

In terms of the SWSHs, we can write the action of ``R_-`` as
```math
R_- {}_{s}Y_{ℓ,m} = \sqrt{(ℓ+s)(ℓ-s+1)}\, {}_{s-1}Y_{ℓ,m}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{R_-(f)\right\}_{s-1,ℓ,m} = \sqrt{(ℓ+s)(ℓ-s+1)}\,\left\{f\right\}_{s,ℓ,m},
```
where the argument `s` of this function is the spin weight of the *input*.

$(operator_signature_note)

The ASCII alias `Rminus` may be used in place of `R₋`; it is public but not exported.

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`Lx`](@ref), [`Ly`](@ref),
[`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const R₋ = RightLowering()


@doc splice_signature_note(raw"""
    ð(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    ð(s, ℓₘₐₓ, [T=Float64])
    ð * w

Compute the spin-raising operator ``\eth``.

This operator was originally defined by [Newman and Penrose](@cite Newman_1966), but is more
completely defined by [Boyle](@cite Boyle_2016).  It is identical to [`R₊`](@ref); see the
[conventions summary](@ref summary_spin_weight).

By definition, the spin-raising operator satisfies the commutator relation ``[R_z, \eth] =
\eth`` (recall that ``R_z`` multiplies a spin-weighted function by its spin weight).  In
terms of the SWSHs, we can write the action of ``\eth`` as
```math
\eth {}_{s}Y_{ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓ,m}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{\eth f\right\}_{s+1,ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\,\left\{f\right\}_{s,ℓ,m},
```
where the argument `s` of this function is the spin weight of the *input*.

$(operator_signature_note)

The ASCII alias `eth` may be used in place of `ð`; it is public but not exported.

See also [`ð̄`](@ref), [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref),
[`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref).
""")
const ð = SpinRaising()


@doc splice_signature_note(raw"""
    ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])
    ð̄(s, ℓₘₐₓ, [T=Float64])
    ð̄ * w

Compute the spin-lowering operator ``\bar{\eth}``.

This operator was originally defined by [Newman and Penrose](@cite Newman_1966), but is more
completely defined by [Boyle](@cite Boyle_2016).  It is the negative of [`R₋`](@ref):
``\bar{\eth} = -R_-``; the sign is Newman and Penrose's.  See the [conventions summary](@ref
summary_spin_weight).

By definition, the spin-lowering operator satisfies the commutator relation ``[R_z,
\bar{\eth}] = -\bar{\eth}`` (recall that ``R_z`` multiplies a spin-weighted function by its
spin weight).  In terms of the SWSHs, we can write the action of ``\bar{\eth}`` as
```math
\bar{\eth} {}_{s}Y_{ℓ,m} = -\sqrt{(ℓ+s)(ℓ-s+1)}\, {}_{s-1}Y_{ℓ,m}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{\bar{\eth} f\right\}_{s-1,ℓ,m} = -\sqrt{(ℓ+s)(ℓ-s+1)}\,\left\{f\right\}_{s,ℓ,m},
```
where the argument `s` of this function is the spin weight of the *input*.

$(operator_signature_note)

The ASCII alias `ethbar` may be used in place of `ð̄`; it is public but not exported.

See also [`ð`](@ref), [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref),
[`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref).
""")
const ð̄ = SpinLowering()

