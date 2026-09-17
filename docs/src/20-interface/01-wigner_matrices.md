# [Wigner's ``𝔇`` and ``d`` matrices](@id interface_wigner_matrices)

```@meta
CurrentModule = SphericalFunctions
```

Wigner's ``𝔇`` matrices — and to a lesser extent, the related ``d`` matrices —
are extremely important in the theory of rotations.  Each element is, itself, a
special function of the rotation group: in particular, an eigenfunction of [the
left- and right-Lie derivatives](@ref "Differential operators"), and
thus a spin-weighted spherical function.  Collectively, they describe how
spin-weighted spherical functions transform under rotation.  But their accurate
and efficient computation is surprisingly subtle.  This package implements the
current state-of-the-art techniques for their fast and accurate computation,
based on the [``H`` recursion](@ref "Algorithm for computing ``H`` (redesigned)").

The convention used here is that
```math
𝔇^{(ℓ)}_{m',m}\left( e^{α𝐤/2} 𝐐 e^{γ𝐤/2} \right)
=
e^{-i m' α}\, 𝔇^{(ℓ)}_{m',m}\left( 𝐐 \right)\, e^{-i m γ}.
```
Note that this equation *does not use Euler angles* per se; ``𝔇`` is
not defined in terms of Euler angles in this package, and the
parameters ``α`` and ``γ`` are simply unrelated parameters in this
expression.  Nonetheless, we can relate our expression to the more
common Euler-angle form: our convention *agrees with*
```math
𝔇^{(ℓ)}_{m',m}(α, β, γ) = e^{-i m' α}\, d^{(ℓ)}_{m',m}(β)\, e^{-i m γ}.
```
In either form, this is the complex conjugate of the convention used
by versions of this package before 3.0, but is more standard in modern
references.  See the [conventions summary](@ref summary_wigner_D) for
the full list of properties, and the [comparison pages](@ref
"Comparisons") for how it relates to other sources.


## Evaluating the matrices

The simplest call takes a rotor and a maximum ``ℓ``:
```julia
using Quaternionic
using SphericalFunctions

R = randn(RotorF64)
ℓₘₐₓ = 8
𝔇 = D(R, ℓₘₐₓ)
```
A rotation is given as a `Rotor`, and only as a `Rotor`: that type is what says a
quaternion denotes a rotation.  A general `Quaternion` carries a magnitude
that the recurrence would simply divide out, and a `QuatVec` is a vector rather
than a rotation at all, so rather than reinterpret either one silently, both are
refused with a message naming what to write instead — `rotor(q)` for the
rotation a quaternion points at, `exp(v/2)` for the rotation a vector generates.

The result is indexed by ``ℓ`` and then by the two matrix indices, each
with its natural range:
```julia
𝔇[ℓ][m′, m]  # for ℓ ∈ 0:ℓₘₐₓ, m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
```
There is no index arithmetic to get wrong: for integer ``ℓ``, `𝔇[ℓ]` is an
`OffsetMatrix` whose axes really are `-ℓ:ℓ`.  (For half-integer ``ℓ`` it is a
[`WignerMatrix`](@ref) instead, indexed exactly the same way; see [Half-integer
indices](@ref half_integer_wigner) below.)  For the ``d`` matrices the interface is the
same, except that the argument is the angle ``β`` rather than a rotor, and the
values are real:
```julia
β = π * rand(Float64)
𝔡 = d(β, ℓₘₐₓ)
𝔡[ℓ][m′, m]
```
`d` also accepts the complex phase ``e^{iβ}`` directly, which is useful when
that is what you have, and a rotor, from which it takes the ``β`` angle.

Four keyword arguments — `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ` and `mₘᵢₙ` — restrict the
block of each matrix that is returned.  Restricting `m′` is the common case:
the spin-weighted spherical harmonics need only one column, and asking for
fewer columns makes the whole calculation cheaper as well as smaller.


## Reusing the workspace

`D` and `d` allocate a fresh workspace on every call, and copy their results out
of it.  When the matrices are needed at many rotors, allocate the workspace once
as a [`WignerDCalculator`](@ref) (or a [`WignerdCalculator`](@ref)) and iterate
over the values of ``ℓ``:
```julia
calculator = WignerDCalculator(R, ℓₘₐₓ)
for (ℓ, 𝔇ˡ) ∈ calculator
    𝔇ˡ[m′, m]                       # for m′, m ∈ -ℓ:ℓ
end
```
The rotor comes first, and it comes at construction: a calculator contains its
rotor from the moment it exists, which is also what fixes the element type,
because that type *is* the rotor's own.  A `Rotor{Float32}` gives a `Float32`
calculator, a `Rotor{Double64}` a `Double64` one, and
[`floattype`](@ref SphericalFunctions.floattype)`(calculator)` reports which it
is.  Nothing overrides that, since the type of the data is the
claim being made about the points, and tying the arithmetic to it keeps the two
from disagreeing; computing in another type therefore means converting the data
— `WignerDCalculator(Rotor{Double64}(R), ℓₘₐₓ)` — or, better, building the rotor
at that type in the first place.  A `WignerdCalculator` is built the same way
from ``β``, given as the angle, as the phase ``e^{iβ}``, or as a rotor of which
only ``β`` is used, and the same four keyword limits apply to both.

The loop yields `ℓ => block` pairs, because a calculator is keyed by ``ℓ``
rather than positioned by it.  That is the shape of a dictionary, and the usual
names follow it: `keys(calculator)` is the range of ``ℓ``,
`length(calculator)` is how many values of ``ℓ`` there are, and
`eltype(calculator)` is a `Pair` type.  Those last two are worth noting because
they mean something else one level down: on a block, `length` counts matrix
elements and `eltype` is the number type.

Each block is a view into the one buffer the calculator owns, and the next step
of the recurrence overwrites it.  Holding only one ``ℓ`` at a time is what makes
a large ``ℓₘₐₓ`` affordable, and it is why a whole pass over ``ℓ`` allocates
nothing at all.  Inside the loop body the view costs nothing; to keep a block,
`copy` it (which preserves the natural indices) or `collect` it (which gives an
ordinary 1-based `Matrix`).  To keep all of them, `collect` the calculator
itself:
```julia
blocks = collect(calculator)        # ℓ => block pairs, every block copied
```
Two things differ here from what `D` and `d` return.  The result is a `Vector`
indexed from 1 rather than by ``ℓ``; each pair has its own ``ℓ``, as
`first(blocks[i])`, and that is what identifies it.  And the copying is peculiar
to `collect`: the generic functions that merely store what they are handed, such
as `map`, `first(calculator, n)` and `Iterators.take`, give back that many
aliases of the same buffer, every one of them showing the last ``ℓ`` that was
computed rather than its own.

For the same reason, two loops over one calculator cannot be interleaved.  Each
step of either loop overwrites what the other is looking at, and nothing
complains; the answers are simply wrong.  `similar(calculator)` is the remedy —
a second workspace on the same rotor data, with nothing computed in it — and it
is why each thread should be given its own calculator.

[`eachℓ`](@ref) does the same iteration over a restricted range of ``ℓ``:
```julia
for (ℓ, 𝔇ˡ) ∈ eachℓ(calculator; ℓₘᵢₙ=2, ℓₘₐₓ=4)
    # ...
end
```
Starting above the calculator's own ``ℓₘᵢₙ`` costs nothing in accuracy — the
recurrence still runs through the intermediate values, and the blocks are
bit-for-bit those of a full pass — but for the same reason it does not save the
work of reaching ``ℓₘᵢₙ`` either.  `eachell` is the ASCII alias.

Underneath, each step of the iteration is a call to [`recurrence!`](@ref)
followed by an index, and that pair can still be written out by hand when a
`for` loop is not the shape wanted:
```julia
for ℓ ∈ 0:ℓₘₐₓ
    recurrence!(calculator, ℓ)
    𝔇ˡ = calculator[ℓ]              # 𝔇ˡ[m′, m], a view into the workspace
end
```
The two rules that iteration enforces are then the caller's: `recurrence!` comes
before the index, and the block must not be kept.  Asking for any ``ℓ`` other
than the one just computed is an error rather than a silently wrong answer.

Stepping through ``ℓ`` in increasing order is the cheap path, because each ``ℓ``
is computed from the one before.  Jumping to an arbitrary ``ℓ`` is allowed:
jumping forward simply runs through the intermediate values, while jumping
backward restarts the recursion from ``ℓₘᵢₙ``.  Either way the result is
bit-for-bit the same as a fresh calculator would give.

To move a calculator to another rotor, hand the new one to [`set_R!`](@ref) —
or, for a `WignerdCalculator` or a [`WignerHCalculator`](@ref), the new angle to
[`set_β!`](@ref).  The new data must be of the floating-point type the
calculator already works in, and a mismatch is an error rather than a silent
conversion: narrowing a `BigFloat` rotor into a `Float64` calculator would throw
away precision that nobody chose to lose.  Whatever the calculator was holding
is discarded, so the next pass over ``ℓ`` starts from ``ℓₘᵢₙ`` again:
```julia
calculator = WignerDCalculator(first(rotors), ℓₘₐₓ)
for R ∈ rotors
    set_R!(calculator, R)
    for (ℓ, 𝔇ˡ) ∈ calculator
        # ...
    end
end
```

That loop is not, however, the best way to handle many rotors.  Give the whole
collection to the constructor instead, and the calculator will run all
`length(rotors)` of them through the recurrence together, each block gaining a
rotor index in front:
```julia
calculator = WignerDCalculator(rotors, ℓₘₐₓ)
for (ℓ, 𝔇ˡ) ∈ calculator
    𝔇ˡ[iᵣ, m′, m]                   # for iᵣ ∈ 1:length(rotors)
end
```
The number of rotors is fixed when the calculator is built, and `set_R!` and
`set_β!` then expect exactly that many.

Batching this way is faster per rotor than looping over a single-rotor
calculator, and the reason is worth stating, because it is the one property of
these calculators that a caller can easily throw away.  For a single rotation
there is nothing here to do in parallel.  Every value the recurrence produces is
built from values it produced just before — the next ``ℓ`` from the previous one,
and each element of a matrix from its neighbours — so the work is inherently
sequential in ``ℓ``, ``m'`` and ``m`` alike.  The one independent dimension is
the rotation itself: the same arithmetic, on unrelated data.  The package is
therefore laid out so that the values belonging to different rotations sit next
to each other in memory, and the innermost loop of every step of the recurrence
runs across the batch.  That is the loop a processor can turn into vector
instructions, and it is the only one available.  Handing the batch to the
library is what lets that happen; a loop over rotors written outside the library
leaves that innermost loop one iteration long, and the vectorization is simply
lost.  For a batch of eight rotations, vectorizing that loop makes the dominant
step of the recurrence about 1.5 times faster, and the values it produces are
identical bit for bit.

Any floating-point element type will do: `Float32`, `Double64`, `BigFloat` and
`ForwardDiff` dual numbers all work, each of them inherited from the rotor data
the calculator was given.  Besides `similar(calculator)`, which keeps the rotor
data and computes nothing, `similar(calculator, R)` gives a fresh workspace on
new data — for the same number of rotors, since that is a property of the
storage, and in the same floating-point type, since a new workspace is no freer
than `set_R!` to change what the calculator works in.


## The ``H`` recursion underneath

Both calculators are built on a [`WignerHCalculator`](@ref), which computes the
real, symmetric ``H`` wedge that the Gumerov–Duraiswami recursion produces
before any phases are applied.  It is available directly for the rare cases
where the wedge itself is what is wanted:
```julia
h = WignerHCalculator(β, ℓₘₐₓ)
for ℓ ∈ 0:ℓₘₐₓ
    recurrence!(h, ℓ)
    h.Hˡ[iᵣ, m′, m]                 # always batched; only |m′| ≤ m ≤ ℓ is stored
end
```
This one is not iterable, and deliberately so: its wedge is a single mutable
object handed back as itself rather than a view that `copy` can preserve, so
there is nothing for a loop to yield that would behave like the blocks above.
Step it with `recurrence!` and read `h.Hˡ`, as here.

The wedge is stored as an [`HWedge`](@ref) (and, during the recursion, an
[`HAxis`](@ref)); the symmetries that relate the rest of the matrix to the
stored wedge are described in the notes on the [``H`` recursion](@ref
"Algorithm for computing ``H``").  Most users should prefer the ``𝔇`` and ``d``
calculators above, which apply those symmetries and the phases for you.


## [Half-integer indices](@id half_integer_wigner)

Wigner's matrices are defined for half-integer ``ℓ, m', m`` as well — they are the
irreducible representations of ``𝐒𝐩𝐢𝐧(3) ≅ 𝐒𝐔(2)`` that do not descend to ``𝐒𝐎(3)``
— and everything on this page works for them too.  Ask for them by giving a `Rational`
``ℓₘₐₓ`` whose denominator is exactly 2:
```julia
𝔇 = D(R, 7//2)           # ℓ = 1//2, 3//2, 5//2, 7//2
𝔡 = d(β, 7//2)
calculator = WignerDCalculator(rotors, 7//2)
```
`D(R, 3.5)` is a `MethodError`; `D(R, 3//1)` and `d(β, 4//2)` report that the index *"must
have denominator 2"*.  The four keyword limits must be half-integers too:
`D(R, 7//2; m′ₘₐₓ=1)` is an error, and `m′ₘₐₓ=1//2` is what was meant.  Both
``m' = ±1/2`` must lie in the requested `m′` range, because the recurrence seeds the whole
ladder from that pair of rows; that is the half-integer form of the integer requirement
that `m′ₘᵢₙ ≤ 0 ≤ m′ₘₐₓ`.

Internally the package stores these indices as [`HalfOddInteger`](@ref
SphericalFunctions.HalfOddInteger)s rather than as `Rational`s, which
is what lets the recurrences run entirely in integer arithmetic.  The
conversion happens at the boundary, so `7//2` remains the natural
thing to write; but the values that come *back* out — from [`ℓ`](@ref
SphericalFunctions.ℓ), [`m′ₘₐₓ`](@ref SphericalFunctions.m′ₘₐₓ),
`axes` and the rest — are of that type, but display as `7//2`.  Call
`Rational(x)` on one to get the actual type.

### What comes back

`OffsetArray` refuses non-integer axes, so half-integer results use a small family of
containers of this package's own instead.  They are indexed exactly as their integer
counterparts are:

| call | integer ``ℓ`` | half-integer ``ℓ`` |
|---|---|---|
| `D(R, ℓₘₐₓ)`, `d(β, ℓₘₐₓ)` | `OffsetVector` of blocks | [`WignerSeries`](@ref) of blocks |
| a block of `D` | `OffsetMatrix` | [`WignerDMatrix`](@ref) |
| a block of `d` | `OffsetMatrix` | [`WignerdMatrix`](@ref) |
| `calc[ℓ]`, `Nᵣ = 1` | `OffsetMatrix` | [`WignerMatrix`](@ref) |
| `calc[ℓ]`, `Nᵣ > 1` | 3-d `OffsetArray` | [`WignerMatrixBatch`](@ref) |

``ℓ`` runs over `1//2, 3//2, …, ℓₘₐₓ`, so `ℓₘᵢₙ` is `1//2` rather than `0`, and that is
where a loop over a calculator begins.  Indexing the series with something that is not one
of those values — `𝔇[3]`, `𝔇[1//1]` — is an error naming the values it does hold, not a
wrong block.

These containers support `[m′, m]` (or `[iᵣ, m′, m]`), `axes`, `size`, `size(w, d)`,
`length`, `ndims`, `eltype`, `parent`, `copy`, `similar`, `collect`, `Array`, `Matrix`,
`==`, iteration and `show`, and a `WignerMatrixBatch` additionally gives `w[iᵣ]` — a
[`WignerMatrix`](@ref) view of one rotor's block, which the integer path has no equivalent
of.  They are deliberately **not** `AbstractMatrix`es, because half-integer axes cannot
satisfy that interface, so linear algebra does not apply to them directly: `B'`, `B * C`,
`lu(B)` and the like are `MethodError`s, and broadcasting (`B .+ 1`) returns an ordinary
1-based `Matrix`, dropping the natural indices.  Call `Matrix(B)` first when you want to do
arithmetic on the block as a matrix.

### ``β``, the double cover, and the ``H`` wedge

Half-integer ``d`` has period ``4π`` in ``β`` rather than ``2π``:
``d^{(ℓ)}(β + 2π) = -d^{(ℓ)}(β)``.  An angle or a `Rotor` therefore pins the result down
unambiguously, but a bare phase ``e^{iβ}`` fixes ``β`` only modulo ``2π``; in that case the
branch ``β ∈ (-π, π]`` is used, so `d(cis(β), ℓₘₐₓ)` agrees with `d(β, ℓₘₐₓ)` bitwise inside
that branch and differs by ``(-1)^{2ℓ}`` outside it.

The double cover shows up in ``𝔇`` as ``𝔇^{(ℓ)}(-R) = -𝔇^{(ℓ)}(R)``, which holds
*exactly* — the two results are bitwise negatives of each other, with no rounding.  A
`Rotor` handed to `d` contributes only ``β ∈ [0, π]``, so `d(-R, ℓₘₐₓ) == d(R, ℓₘₐₓ)`: the
sign lives entirely in the ``α`` and ``γ`` phases.

One caveat applies to the ``H`` wedge underneath.  For half-integer indices ``H`` is
symmetric only up to the sign ``σ = \mathrm{sgn}(m)\,\mathrm{sgn}(m')``, so code that reads
`calc.Hˡ` and transposes it by hand will get the wrong sign for half of the elements; use
[`wedge_value`](@ref SphericalFunctions.wedge_value), which applies ``σ`` for you.  See the notes on the
[``H`` recursion](@ref "Algorithm for computing ``H`` (redesigned)").

### What is not supported

Half-integer indices are available for `D`, `d`, the three calculators, and
[`sYlmCalculator`](@ref) (see [``{}_{s}Y_{ℓ,m}`` functions](@ref interface_sYlm)).  They are
*not* available for the flat, mode-ordered interfaces — [`sYlm`](@ref), [`sYlm!`](@ref),
[`sYlm_matrix`](@ref), [`Ysize`](@ref), [`Yindex`](@ref), [`Yrange`](@ref),
[`ModeWeights`](@ref) — nor for the [transforms](@ref interface_transformations)
([`SSHT`](@ref), [`map2salm`](@ref)), all of which are built on the canonical integer
mode ordering.  Calling those with half-integer arguments is a `MethodError`.


## Docstrings

```@docs
D
d
WignerDCalculator
WignerdCalculator
WignerHCalculator
recurrence!
eachℓ
eachell
set_R!
set_β!
```


## Containers

The types that `D`, `d` and `calc[ℓ]` return for half-integer ``ℓ``, and the abstract type
they share with the workspaces below.

```@docs
AbstractWignerMatrix
WignerMatrix
WignerDMatrix
WignerdMatrix
WignerMatrixBatch
WignerSeries
WignerCalculator
```


## Workspaces

```@docs
HWedge
HAxis
```


## Methods of `Base` functions

Indexing a calculator or a container, copying one, emptying one, gathering every ``ℓ`` of
one, and converting one to an ordinary `Array` are all spelled with the usual `Base`
functions.  Every specialization the package defines is collected here, for the Wigner
calculators and containers and for the [`sYlmCalculator`](@ref) alike.

```@docs
Base.getindex
Base.similar
Base.fill!
Base.collect
Base.Matrix
```


## Accessors

The same handful of names reports the index ranges of every container and calculator in the
package.  Each has an ASCII alias, given in its docstring, for use where the subscripted
Unicode names are inconvenient.

```@docs
SphericalFunctions.ℓ
SphericalFunctions.ℓₘᵢₙ
SphericalFunctions.ℓₘₐₓ
SphericalFunctions.m′ₘₐₓ
SphericalFunctions.m′ₘᵢₙ
SphericalFunctions.mₘₐₓ
SphericalFunctions.mₘᵢₙ
SphericalFunctions.Nᵣ
SphericalFunctions.ishalfinteger
SphericalFunctions.isbatched
```
