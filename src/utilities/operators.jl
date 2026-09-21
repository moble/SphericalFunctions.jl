# Angular-momentum operators as matrices acting on mode weights.
#
# Every function here has the signature
#
#     f(s, ℓₘᵢₙ, ℓₘₐₓ, [T=Float64])   or   f(s, ℓₘₐₓ, [T=Float64])  (with ℓₘᵢₙ = abs(s))
#
# and returns a (sparse) matrix that acts on a vector of mode weights ordered as
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
# Each operator is one function with two kinds of method: a *boundary* method typed
# `IndexSpelling`, which normalizes the indices and re-dispatches, and a *worker* method typed
# `where {IT<:IntegerHalf}`, which does the arithmetic.  The worker is reached by dispatch
# rather than by a separate underscore-prefixed name: `IT<:IntegerHalf` with one `IT` for all
# three indices is strictly more specific than three independent `IndexSpelling`s, so the
# worker always wins once the indices agree, and `unify_indices` guarantees that they do.
# Being methods of the exported name, the workers are simply undocumented rather than hidden.
#
# The indices `s`, `ℓₘᵢₙ` and `ℓₘₐₓ` may be integers or half-odd-integers, the latter spelled
# as `Rational`s with denominator 2 or as `HalfOddInteger`s.  Each public function is a
# boundary method, typed `IndexSpelling` on its indices, which does nothing but normalize the
# three with `unify_indices` and re-dispatch to a private worker — `L²` for `L²`, and so on —
# whose signature is `where {IT<:IntegerHalf, T}`.  The worker therefore sees three indices of
# one concrete type and never a `Rational`, and a call that mixes the two kinds of index is
# refused at the boundary.  The three-argument forms normalize first as well, so that the
# default `ℓₘᵢₙ = abs(s)` is computed from the normalized spin weight.  The worker bodies are
# generic over the two kinds, because every quantity they form — ℓ ± m, ℓ ± s, 2ℓ — is an
# `Int` for either; the two exceptions, ℓ(ℓ+1) and the conversion of an index to the matrix
# element type, go through `casimir_eigenvalue` below and `index_value` in
# `half_odd_integer.jl`.  The operator methods on a `ModeWeights`, in `mode_weights.jl`, call
# the workers directly, since the container already holds normalized indices.

const _operator_signature_note = """
The argument `ℓₘᵢₙ` may be omitted, in which case it defaults to `abs(s)`.  The result acts
on a vector of mode weights ordered as `[f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]`; any
entries with ``ℓ < |s|`` are mapped to zero.  The indices `s`, `ℓₘᵢₙ` and `ℓₘₐₓ` may be
integers or half-odd-integers, the latter spelled as `Rational`s with denominator 2 — as in
`L²(1//2, 7//2)` — in which case every ``ℓ`` and ``m`` of the ordering is a half-odd-integer.
The indices in one call must all be of one kind; a call that mixes them, such as
`L²(1//2, 0, 7//2)`, is an error.
"""

# The docstrings below are `raw` strings, because they are full of LaTeX backslashes, and a
# `raw` string does not interpolate — so `$(_operator_signature_note)` written inside one
# stays there as literal text rather than being replaced by the note.  (It used to, and the
# note was never rendered.)  This splices it in explicitly; `@doc` accepts any expression
# that evaluates to a string.
splice_signature_note(s) =
    replace(s, "\$(_operator_signature_note)" => _operator_signature_note)

# The eigenvalue ℓ(ℓ+1) of L² and R², as the value the typed comprehension in those functions
# converts to `T`.  For an integer ℓ it is the integer ℓ(ℓ+1), exactly as it always was.  For
# a half-odd ℓ the product is a quarter-integer, which the index type cannot hold and which
# must not be formed as a `Rational`; instead the `Int` (2ℓ)(2ℓ+2) is converted to `T` and
# divided by 4.  Division by 4 is exact in every binary floating-point type, so this path is
# as exact as the integer one.
@inline casimir_eigenvalue(::Type{T}, ℓ::Integer) where {T} = ℓ*(ℓ+1)
@inline casimir_eigenvalue(::Type{T}, ℓ::HalfOddInteger) where {T} = T((2ℓ)*(2ℓ+2)) / 4


### The operators as objects, and their matrix elements.
#
# Each operator is a zero-size singleton, so that dispatching on it costs nothing and every
# trait below folds away at compile time.  The twelve exported names are instances of these
# types.
#
# The matrix elements are written *once*, here, and evaluated both by the comprehensions that
# build the operator matrices below and by the loops that apply an operator to a `ModeWeights`
# without building one.  That is what makes the two agree bit for bit — not a coincidence to
# be tested for, but the same expression evaluated twice.  Note in particular that the
# `ℓ < …` mask lives *inside* the coefficient, so that neither caller can forget it.

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

# How each operator changes the spin weight of what it acts on.
@inline Δspin(::DifferentialOperator) = 0
@inline Δspin(::Union{RightRaising,  SpinRaising})  =  1
@inline Δspin(::Union{RightLowering, SpinLowering}) = -1

# Which band of the matrix the operator occupies, and hence which builder and which kernel
# apply.  A sub-diagonal entry takes `f[ℓ, m-1]` into `out[ℓ, m]`; a super-diagonal one takes
# `f[ℓ, m+1]`.
abstract type BandStructure end
struct DiagonalBand      <: BandStructure end
struct SubdiagonalBand   <: BandStructure end
struct SuperdiagonalBand <: BandStructure end
struct TridiagonalBand   <: BandStructure end

@inline bandstructure(::DifferentialOperator) = DiagonalBand()
@inline bandstructure(::LeftRaising)  = SubdiagonalBand()
@inline bandstructure(::LeftLowering) = SuperdiagonalBand()
@inline bandstructure(::Union{LeftX, LeftY}) = TridiagonalBand()

# The element type of the matrix, given the real type it was asked for.  `Ly` is the only one
# whose entries are complex — they are ∓i/2 times those of `L₊` and `L₋`.
@inline coefftype(::DifferentialOperator, ::Type{T}) where {T} = T
@inline coefftype(::LeftY, ::Type{T}) where {T} = Complex{T}

# The operators are values now rather than functions, so they need to say their own names:
# `nameof` because code (and tests) reach for it, and `show` so that one prints as `ð` rather
# than as `SphericalFunctions.SpinRaising()`.
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

# The ℓ below which the result vanishes.  For the spin-changing operators the cutoff is set by
# the *output* spin weight, which is the `s′` of the matrix builders.
@inline support_ℓ(::DifferentialOperator, s) = abs(s)
@inline support_ℓ(::Union{RightRaising,  SpinRaising},  s) = max(abs(s), abs(s + 1))
@inline support_ℓ(::Union{RightLowering, SpinLowering}, s) = max(abs(s), abs(s - 1))

# `casimir_eigenvalue` and `index_value` return an `Int` for an integer index; the typed
# comprehensions used to do the conversion, so the coefficients must do it explicitly or the
# loops that share them go type-unstable on the integer path.
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
# The negation is applied to the *masked* value, as `-R₋` does, so that the vanishing entries
# are `-0.0` here too and even `isequal` agrees with the matrix.
@inline diagonal_coefficient(::SpinLowering, ::Type{T}, s, ℓ, m) where {T} =
    -diagonal_coefficient(RightLowering(), T, s, ℓ, m)

# The ladder coefficients, indexed by the *output* mode `(ℓ, m)`.  Both vanish exactly at the
# edge of their ℓ block — `√0` at `m = -ℓ` for the raising one and at `m = +ℓ` for the lowering
# one — which is what keeps an ℓ block from coupling to its neighbours.
@inline function subdiagonal_coefficient(op::LeftRaising, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : √T((ℓ+m)*(ℓ-m+1))
end
@inline function superdiagonal_coefficient(op::LeftLowering, ::Type{T}, s, ℓ, m) where {T}
    ℓ < support_ℓ(op, s) ? zero(T) : √T((ℓ-m)*(ℓ+m+1))
end
@inline subdiagonal_coefficient(::LeftX, ::Type{T}, s, ℓ, m) where {T} =
    subdiagonal_coefficient(LeftRaising(), T, s, ℓ, m) / 2
@inline superdiagonal_coefficient(::LeftX, ::Type{T}, s, ℓ, m) where {T} =
    superdiagonal_coefficient(LeftLowering(), T, s, ℓ, m) / 2
# Built as a complex number with a zero real part rather than divided by `2im`, for the reason
# given at `Ly` below.
@inline subdiagonal_coefficient(::LeftY, ::Type{T}, s, ℓ, m) where {T} =
    Complex{T}(zero(T), -subdiagonal_coefficient(LeftRaising(), T, s, ℓ, m) / 2)
@inline superdiagonal_coefficient(::LeftY, ::Type{T}, s, ℓ, m) where {T} =
    Complex{T}(zero(T), superdiagonal_coefficient(LeftLowering(), T, s, ℓ, m) / 2)


### The three call shapes, written once for every operator.
#
# These replace the thirty-six methods — three per operator — that the twelve names used to
# carry between them.  The first two are the `IndexSpelling` boundaries, which normalize and
# re-dispatch; the third is the worker, reached once the three indices agree in kind, and it
# builds the matrix from the band structure and the coefficients above.

function (op::DifferentialOperator)(
    s::IndexSpelling, ℓₘᵢₙ::IndexSpelling, ℓₘₐₓ::IndexSpelling, ::Type{T}=Float64
) where T
    op(unify_indices(s, ℓₘᵢₙ, ℓₘₐₓ)..., T)
end
function (op::DifferentialOperator)(s::IndexSpelling, ℓₘₐₓ::IndexSpelling, ::Type{T}=Float64) where T
    s, ℓₘₐₓ = unify_indices(s, ℓₘₐₓ)
    op(s, abs(s), ℓₘₐₓ, T)
end
function (op::DifferentialOperator)(s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, ::Type{T}) where {IT<:IntegerHalf, T}
    operator_matrix(op, bandstructure(op), s, ℓₘᵢₙ, ℓₘₐₓ, T)
end

# One builder per band structure.  The `ifelse` in the ladder ranges drops the one mode that
# has no band entry — the very first for a sub-diagonal, the very last for a super-diagonal —
# exactly as the hand-written builders did.
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
# The matrix builders above and the loops below evaluate the *same* coefficient functions, so
# the two agree bit for bit rather than merely to within rounding — which is what the existing
# tests in `test/mode_weights/mode_weights.jl` assert, with `==` rather than `≈`.  Reproducing
# that exactly is why these loops are written plainly: no `@simd`, no `@fastmath`, no `muladd`,
# and the two terms of a tridiagonal row summed left to right, as `LinearAlgebra`'s own
# `l[i-1]*b₋ + d[i]*b₀ + u[i]*b₊` does with `d` identically zero.
#
# The ladder coefficients vanish *exactly* at the edge of each ℓ block — √0 at `m = -ℓ` for the
# raising one and at `m = +ℓ` for the lowering one — so no block ever couples to its neighbour
# and the loops need no per-block special case.  Only the very first and very last position in
# the whole vector need a branch, because there the matrix has no band entry at all.

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


@doc splice_signature_note(raw"""
    L²(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    L²(s, ℓₘₐₓ, [T])

Compute the total angular-momentum operator (the Casimir operator) for spin weight `s`.

This is the standard ``L^2`` operator, familiar from basic physics, extended to work with
SWSHs.  It is equal to
```math
L^2 = L_x^2 + L_y^2 + L_z^2 = \frac{L_+L_- + L_-L_+ + 2L_zL_z}{2}.
```
Note that these are the left Lie derivatives, but ``L^2 = R^2``, where ``R`` is the right Lie
derivative.  See the [conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite
Boyle_2016) for more details.

In terms of the SWSHs, we can write the action of ``L^2`` as
```math
L^2 {}_{s}Y_{ℓ,m} = ℓ\,(ℓ+1) {}_{s}Y_{ℓ,m}.
```

$(_operator_signature_note)

See also [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const L² = Casimir()


@doc splice_signature_note(raw"""
    Lz(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    Lz(s, ℓₘₐₓ, [T])

Compute the angular-momentum operator associated with the ``z`` direction.  This is the
standard ``L_z`` operator, familiar from basic physics, extended to work with SWSHs.  Note
that this is the left Lie derivative; see [`Rz`](@ref) for the equivalent right Lie
derivative.  See the [conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite
Boyle_2016) for more details.

In terms of the SWSHs, we can write the action of ``L_z`` as
```math
L_z {}_{s}Y_{ℓ,m} = m\, {}_{s}Y_{ℓ,m}.
```

$(_operator_signature_note)

See also [`L²`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const Lz = LeftZ()


@doc splice_signature_note(raw"""
    L₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    L₊(s, ℓₘₐₓ, [T])

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

$(_operator_signature_note)

See also [`L²`](@ref), [`Lz`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const L₊ = LeftRaising()


@doc splice_signature_note(raw"""
    L₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    L₋(s, ℓₘₐₓ, [T])

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

$(_operator_signature_note)

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const L₋ = LeftLowering()


@doc splice_signature_note(raw"""
    Lx(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    Lx(s, ℓₘₐₓ, [T])

Compute the ``x`` component of the left angular-momentum operator, ``L_x = (L_+ + L_-)/2``.

This is the standard ``L_x`` operator, familiar from basic physics, extended to work with
SWSHs.  See the [conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite
Boyle_2016) for more details.  The matrix is real, symmetric and tridiagonal.

$(_operator_signature_note)

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
    Ly(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    Ly(s, ℓₘₐₓ, [T])

Compute the ``y`` component of the left angular-momentum operator,
``L_y = (L_+ - L_-)/(2i)``.

This is the standard ``L_y`` operator, familiar from basic physics, extended to work with
SWSHs.  See the [conventions summary](@ref summary_L_R_definitions) or [Boyle](@cite
Boyle_2016) for more details.  The matrix is tridiagonal and purely imaginary, so its
element type is `Complex{T}` rather than `T`.

$(_operator_signature_note)

There are no corresponding `Rx` and `Ry` functions; see [`Lx`](@ref) for why.

See also [`Lx`](@ref), [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref),
[`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const Ly = LeftY()


@doc splice_signature_note(raw"""
    R²(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    R²(s, ℓₘₐₓ, [T])

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

$(_operator_signature_note)

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const R² = RightCasimir()


@doc splice_signature_note(raw"""
    Rz(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    Rz(s, ℓₘₐₓ, [T])

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

$(_operator_signature_note)

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const Rz = RightZ()


@doc splice_signature_note(raw"""
    R₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    R₊(s, ℓₘₐₓ, [T])

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

$(_operator_signature_note)

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const R₊ = RightRaising()


@doc splice_signature_note(raw"""
    R₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    R₋(s, ℓₘₐₓ, [T])

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

$(_operator_signature_note)

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`Lx`](@ref), [`Ly`](@ref),
[`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref), [`ð`](@ref), [`ð̄`](@ref).
""")
const R₋ = RightLowering()


@doc splice_signature_note(raw"""
    ð(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    ð(s, ℓₘₐₓ, [T])

Compute coefficients for the spin-raising operator ``\eth``.

This operator was originally defined by [Newman and Penrose](@cite Newman_1966), but is
more completely defined by [Boyle](@cite Boyle_2016).  It is identical to [`R₊`](@ref); see
the [conventions summary](@ref summary_spin_weight).

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

$(_operator_signature_note)

See also [`ð̄`](@ref), [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref),
[`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref).
""")
const ð = SpinRaising()


@doc splice_signature_note(raw"""
    ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, [T])
    ð̄(s, ℓₘₐₓ, [T])

Compute coefficients for the spin-lowering operator ``\bar{\eth}``.

This operator was originally defined by [Newman and Penrose](@cite Newman_1966), but is
more completely defined by [Boyle](@cite Boyle_2016).  It is the negative of [`R₋`](@ref):
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

$(_operator_signature_note)

See also [`ð`](@ref), [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref),
[`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref).
""")
const ð̄ = SpinLowering()
