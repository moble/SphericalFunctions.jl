# [Half-integer indices](@id interface_half_integers)

```@meta
CurrentModule = SphericalFunctions
```

Wigner's matrices and the spin-weighted spherical harmonics are
defined for half-integer ``ℓ, m', m`` and ``s`` as readily as for
integer ones.  These are the irreducible representations of
``𝐒𝐩𝐢𝐧(3) ≅ 𝐒𝐔(2)`` that do not descend to ``𝐒𝐎(3)``, and they
are supported throughout this package — by ``𝔇``, ``d``,
``{}_{s}Y_{ℓ,m}``, and ``Y_{ℓ,m}``, by the mode weights and the
operators acting on them, and by the transforms — at the same accuracy
and essentially the same speed as integer indices.

Because half-integers appear in every part of the interface, the way
they are input, the types they produce, and the few places where they
behave differently are all collected on this page.

## Creating a half-integer

A half-integer index can be input as a [`Rational`](@extref
Base.Rational) whose denominator is exactly 2:
```julia
𝔇 = D(R, 7//2)
𝔡 = d(β, 7//2)
calculator = WignerDCalculator(rotors, 7//2)
Y = sYlm(R, 7//2, -3//2)
𝒯 = SSHT(1//2, 7//2)
```
Note that `Rational` numbers must be written with the double-slash
`//` operator, not a single slash; `7//2` is a `Rational`, but Julia
immediately converts `7/2` to `3.5` — a `Float64` — which will produce
and error.  Optional keyword limits such as `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ`
and `mₘᵢₙ` may also occur in some of these functions.

The indices within a single call must all be of one kind.  An integer
``ℓ`` goes with an integer ``m`` and an integer spin weight; a
half-odd ``ℓ`` with half-odd ones.  There is no mixed case, because no
single representation contains both.  A call that mixes them is
refused with a specific, explanatory error message — rather than with
the bare `MethodError` that dispatch alone would produce:
```julia
julia> sYlm(R, 7//2, 1)
ERROR: ArgumentError: The indices in one call must all be integers or
       all be half-odd-integers; got 7//2, 1.
```

## The `HalfOddInteger` type

These indices are not kept as `Rational`s internally.  At the boundary
of every public function they are converted to [`HalfOddInteger`](@ref
SphericalFunctions.HalfOddInteger), a type that stores the odd
numerator of ``n//2`` and nothing else.  It exists for two reasons.

  1. It makes the *kind* of index a property of the type.  Every
     container and calculator is parameterized by its index type, so
     the integer and half-integer cases dispatch to separate,
     separately compiled code, rather than branching on a denominator
     at run time.
  2. It keeps the recurrences in integer arithmetic.  The sum or
     difference of two `HalfOddInteger`s is an `Int`, and so is
     ``2ℓ``, so a recurrence coefficient such as ``(ℓ-m)(ℓ+m+1)`` is
     computed entirely in integers while being written exactly as the
     references write it.  The same expression in `Rational`
     arithmetic is roughly thirty times slower.

The type is deliberately spare.  It defines only the operations the
recurrences and the public interface actually need, and in particular
defines no `promote_rule`, so that an unanticipated operation is a
`MethodError` rather than a silent — and invisibly slow — fall-back to
`Rational`.  Because a half-odd-integer is never zero and never one,
`zero`, `one` and `oneunit` throw rather than returning a value of
some other type.

Together with the `Integer`s, these make up [`HalfInteger`](@ref
SphericalFunctions.HalfInteger), the union of index types over which
the recurrences are defined.  Neither name is exported, and neither
*needs* to be constructed by hand: because the conversion happens at
the interface boundary, `7//2` remains a natural thing to write
everywhere — though it is also acceptable to construct a
`HalfOddInteger` directly and pass that.  What does deserve noticing
is that the values coming *out* — from [`ℓ`](@ref
SphericalFunctions.ℓ), [`m′ₘₐₓ`](@ref SphericalFunctions.m′ₘₐₓ),
`axes`, and the other accessors — are `HalfOddInteger`s, though they
display simply as `7//2` and compare equal to it.  Call `Rational(i)`
to convert one back.

## What the calls return

Nothing, which is the point.  Through version 2 the integer path
returned [`OffsetArray`](@ref OffsetArrays.OffsetArray)s, so that a
block could be indexed by its natural ``m`` and ``m'``; half-odd
indices could not be handled that way, because `OffsetArray`s cannot
use non-integer axes, and this package grew its own family of
containers for them.  Since version 3 those containers are what *every*
call returns, whichever kind of index is in play, so there is no table
of correspondences to learn:

| call | what comes back |
|---|---|
| `D(R, ℓₘₐₓ)`, `d(β, ℓₘₐₓ)` | a [`WignerSeries`](@ref) of blocks |
| a block of `D` | [`WignerDMatrix`](@ref) |
| a block of `d` | [`WignerdMatrix`](@ref) |
| either, over a batch of rotors | [`WignerMatrixBatch`](@ref) |
| `sYlm(R, ℓₘₐₓ, s)` | [`HarmonicValues`](@ref) |
| a block of an `sYlmCalculator`, one spin weight | [`DegreeBlock`](@ref) |
| the same, over a batch of rotors | [`DegreeBlockBatch`](@ref) |
| a block of a range of spin weights | [`SpinMatrix`](@ref) |
| the same, over a batch of rotors | [`SpinMatrixBatch`](@ref) |

The reason for retiring the `OffsetArray`s from the integer path is not
uniformity but safety, and is set out under [Containers](@ref
interface_containers): an `OffsetArray` with non-trivial offsets accepts
`*` and `mul!` and returns silently wrong answers.  [`strided`](@ref) is
the explicit route to a plain 1-based array, for either kind of index,
and [`relabel`](@ref) is the way back.

``ℓ`` runs over `1//2, 3//2, …, ℓₘₐₓ`, so `ℓₘᵢₙ` is `1//2` rather than
`0`, and that is where a loop over a calculator begins.  Indexing a
series with something that is not one of those values — `𝔇[3]`,
`𝔇[1//1]` — is an error.

These containers support `[m′, m]` (or `[iᵣ, m′, m]` where relevant),
`axes`, `size`, `size(w, d)`, `length`, `ndims`, `eltype`, `parent`,
`copy`, `similar`, `collect`, `Array`, `Matrix`, `==`, iteration and
`show`, and a `WignerMatrixBatch` additionally gives `w[iᵣ]` — a
[`WignerMatrix`](@ref) view of one rotor's block.  They are
deliberately **not** `AbstractMatrix`es, so linear algebra does not
apply to them directly: operations like `B'`, `B * C`, `lu(B)` and the
like are `MethodError`s, and broadcasting (`B .+ 1`) returns an
ordinary 1-based `Matrix`, dropping the natural indices.  Call
[`strided`](@ref) for a 1-based view of the same storage, on which
BLAS works at full speed, or `Matrix(B)` for an independent copy.

## ``β``, the double cover, and the ``H`` wedge

Half-integer ``d`` has period ``4π`` in ``β`` rather than ``2π``:
``d^{(ℓ)}(β + 2π) = -d^{(ℓ)}(β)``.  An angle or a `Rotor` therefore
pins the result down unambiguously, but a bare phase ``e^{iβ}`` fixes
``β`` only modulo ``2π``; in that case the branch ``β ∈ (-π, π]`` is
used, so `d(cis(β), ℓₘₐₓ)` agrees with `d(β, ℓₘₐₓ)` bitwise inside
that branch and differs by ``(-1)^{2ℓ}`` outside it.

The double cover shows up in ``𝔇`` as ``𝔇^{(ℓ)}(-R) =
-𝔇^{(ℓ)}(R)``, which holds *exactly* — the two results are bitwise
negatives of each other, with no rounding.  A `Rotor` passed to `d`
contributes only ``β ∈ [0, π]``, so `d(-R, ℓₘₐₓ) == d(R, ℓₘₐₓ)`: the
sign lives entirely in the ``α`` and ``γ`` phases.

One caveat applies to the ``H`` wedge underneath.  For half-integer
indices ``H`` is symmetric only up to the sign ``σ =
\mathrm{sgn}(m)\,\mathrm{sgn}(m')``, so code that reads `calc.Hˡ` and
transposes it by hand will get the wrong sign for half of the
elements; use [`wedge_value`](@ref SphericalFunctions.wedge_value),
which applies ``σ`` for you.  See the notes on the [``H``
recursion](@ref "Algorithm for computing ``H`` (redesigned)").

## Half-integer spin weight

Two things differ mathematically when the spin weight itself is
half-odd.  The prefactor ``(-1)^s`` in the definition of
``{}_{s}Y_{ℓ,m}`` is then ``\pm i``; this package uses the principal
branch ``(-1)^s ≡ e^{iπs} = i^{2s}``, settled in the [conventions
summary](@ref summary_swsh).  Consequently ``{}_{s}λ_{ℓ,m}(θ)`` — the
value at ``(θ, ϕ=0)``, which is real for integer ``s`` — is not real
here.

The other difference is in what a function of half-integer spin weight
is a function *of*.  The two rotors ``±𝐑`` above a point of the
sphere give values of opposite sign, so "the value at ``(θ, ϕ)``" is
defined only once the rotor is named.  That has consequences for the
sampling that the transforms do, and is discussed on the
[transformations page](@ref transformations_half_integer).


## What is not supported

Half-integer indices are available throughout most of the package,
except where otherwise noted.  Three things do not accept them.
[`Ylm`](@ref) is integer-only, because the ordinary spherical
harmonics are the spin-weight-zero members of the family, and no
half-integer series contains ``s = 0``.  The two equiangular grids,
`driscoll_healy_pixels` and `mcewen_wiaux_pixels`, take an integer
band limit, since they ignore the spin weight and no transform uses
them by default.  And the `"Minimal"` transform method is
integer-only, because its bookkeeping of rings and aliased modes is
written for integer indices and has not *yet* been extended;
`SSHT(1//2, 7//2; method="Minimal")` is refused with a message naming
the `"RS"` and `"Matrix"` methods, which do accept half-integer
indices.


## Docstrings

```@docs
SphericalFunctions.HalfOddInteger
SphericalFunctions.HalfInteger
```
