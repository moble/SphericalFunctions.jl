# [``𝔇`` and ``d`` matrices, and ``{}_{s}Y_{ℓ,m}`` and ``Y_{ℓ,m}`` functions](@id interface_wigner_matrices)

```@meta
CurrentModule = SphericalFunctions
```

Wigner's ``𝔇`` matrices — and to a lesser extent, the related ``d``
matrices — are extremely important in the theory of rotations.  Each
element is, itself, a special function of the rotation group: in
particular, an eigenfunction of [the left- and right-Lie
derivatives](@ref "Differential operators"), and thus a spin-weighted
spherical function.  See the "Background" section, and particularly
[this page](@ref sYlm_and_Dlmpm), for details.  Collectively, they
describe how spin-weighted spherical functions transform under
rotation.  But their accurate and efficient computation is
surprisingly subtle.  This package implements the current
state-of-the-art techniques for their fast and accurate computation,
based on the [``H`` recursion](@ref "Algorithm for computing ``H``
(redesigned)") introduced by [Gumerov_2001](@citet).

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
In either form, this is the *complex conjugate of* the convention used
by versions of this package before 3.0, but is more standard in modern
references.  See the [conventions summary](@ref summary_wigner_D) for
the full list of properties, and the [comparison pages](@ref
"Comparisons") for how it relates to other sources.

The spin-weighted spherical harmonics are an important set of
[functions defined on](@cite Boyle_2016) that same group, particularly
useful in describing the angular dependence of polarized fields, like
the electromagnetic field and the gravitational-wave field.
Originally introduced by [Newman_1966](@citet), they are essentially
single columns of Wigner's ``𝔇`` matrices:
```math
{}_{s}Y_{ℓ,m}(𝐑)
  = (-1)^s \sqrt{\frac{2ℓ+1}{4π}} \, \overline{𝔇^{(ℓ)}_{m, -s}(𝐑)}.
```
(See the [conventions summary](@ref summary_swsh) for this and the
related definitions.)  They are therefore computed by the same
recursion, and because only the single column ``m' = -s`` is needed,
both the storage and the work are much smaller than for the full
matrices.  The standard (scalar) spherical harmonics are the special
case of spin weight ``s = 0``,
```math
Y_{ℓ,m}(𝐑) = {}_{0}Y_{ℓ,m}(𝐑),
```
with the usual spherical-coordinate form being defined as
```math
Y_{ℓ,m}(θ, ϕ) = {}_{0}Y_{ℓ,m}\left(
    \exp\left(\frac{ϕ}{2}𝐤\right) \exp\left(\frac{θ}{2}𝐣\right)
\right).
```

Everything below applies to all three sets of functions, and the few
differences among them are noted where they arise.

## [Evaluating the functions](@id interface_wigner_matrices_evaluation)

The simplest call takes a rotor and a maximum ``ℓ``:
```julia
using Quaternionic
using SphericalFunctions

R = randn(RotorF64)
ℓₘₐₓ = 8
𝔇 = D(R, ℓₘₐₓ)
```
A rotation *must* be given as a [`Quaternionic.Rotor`](@extref): that
type is what specifies that a quaternion denotes a rotation.  To make
this easier, the `Quaternionic` package includes several easy ways to
create a `Rotor`:
  - Explicitly from its components as `Rotor(w, x, y, z)`, or from a
    `Quaternion` as `Rotor(q)`.
  - [From Euler angles](@extref `Quaternionic.from_euler_angles`).
  - [From spherical coordinates](@extref
    `Quaternionic.from_spherical_coordinates`).
  - [From a rotation matrix](@extref
    `Quaternionic.from_rotation_matrix`).
  - From an angle `θ` and unit vector `v` as `exp(θ*QuatVec(v)/2)`,
    which is the right-handed rotation by `θ` about `v`.

The result is indexed by ``ℓ`` and then by the two matrix indices,
each with its natural range:
```julia
𝔇[ℓ][m′, m]  # for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ and m′,m ∈ -ℓ:ℓ
```
There is no index arithmetic to get wrong: for integer ``ℓ``, `𝔇[ℓ]`
is an
[`OffsetMatrix`](https://juliaarrays.github.io/OffsetArrays.jl/stable/)
whose axes are just `-ℓ:ℓ`.  (For half-integer ``ℓ`` it is a
[`WignerMatrix`](@ref) instead, which is designed to act just like an
`OffsetMatrix`, but can deal with half-integer indices; see
[Half-integer indices](@ref interface_half_integers).)  Note that
`ℓₘᵢₙ=0` for integers but `ℓₘᵢₙ=1//2` for half-integers.  To get the
underlying `Matrix`, call `parent(𝔇[ℓ])`.

For the ``d`` matrices the interface is the same, except that the
argument is the angle ``β`` rather than a rotor, and the values are
real:
```julia
β = π * rand(Float64)
𝔡 = d(β, ℓₘₐₓ)
𝔡[ℓ][m′, m]
```
`d` also accepts the complex phase ``e^{iβ}`` directly, which is
useful when that is what you have, and a rotor, from which it takes
the equivalent of the ``β`` angle.

Both `D` and `d` accept four keyword arguments — `m′ₘₐₓ`, `m′ₘᵢₙ`,
`mₘₐₓ`, and `mₘᵢₙ` — restricting the block of each matrix that is
returned.  Restricting `m′` is the common case: the spin-weighted
spherical harmonics need only one column, and asking for fewer columns
makes the whole calculation cheaper as well as smaller.

For those spin-weighted spherical harmonics, a more direct and
efficient method is provided by [`sYlm`](@ref), taking the spin weight
as a third argument:
```julia
s = -2
sY = sYlm(R, ℓₘₐₓ, s)
```
The `s` value can also be a range of spin weights, which is more
efficient than calling `sYlm` repeatedly for each spin weight:
```julia
sY = sYlm(R, ℓₘₐₓ, -2:2)
```

A harmonic has only one index besides ``ℓ``, so there is no matrix to
index into here, and the result is a single flat `Vector` of `Complex`
numbers ordered by ``ℓ`` and then by ``m``:
```julia
[sY[Yindex(ℓ, m, abs(s))] for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ] == sY
```
This is the canonical ordering of mode weights, described by
[`Ysize`](@ref), [`Yindex`](@ref) and [`Yrange`](@ref), and wrapped by
[`ModeWeights`](@ref) so that the index arithmetic need not be done by
hand.  Modes with ``ℓ < |s|`` do not exist (or are inherently zero),
so by default the vector starts at ``ℓ = |s|``.  Pass `ℓₘᵢₙ=0` to get
a vector that starts at ``ℓ = 0`` instead, with zeros in the
nonexistent modes; this is the layout that some downstream packages
use for every spin weight at once.  For ``s = 0`` these are the
ordinary scalar spherical harmonics ``Y_{ℓ,m}``, which [`Ylm`](@ref)
gives without the redundant argument: `Ylm(R, ℓₘₐₓ)` is exactly
`sYlm(R, ℓₘₐₓ, 0)`, and starts at ``ℓ = 0`` because no modes are
missing there.


## Iterating over ``ℓ`` and reusing the storage

The functions above allocate a fresh workspace on every call, and copy
the results *for all ``ℓ`` values* out of it.  When the total size is
too large, or the values are needed for many rotors, allocate the
workspace once as a [`WignerDCalculator`](@ref) (or a
[`WignerdCalculator`](@ref)) and iterate over the values of ``ℓ``:
```julia
calculator = WignerDCalculator(R, ℓₘₐₓ)
for (ℓ, 𝔇ˡ) ∈ calculator
    # 𝔇ˡ[m′, m] is available for m′, m ∈ -ℓ:ℓ
end
```
An [`sYlmCalculator`](@ref) is built and iterated in just the same
way, with the spin weight given as a third argument:
```julia
calculator = sYlmCalculator(R, ℓₘₐₓ, s)
for (ℓ, ₛYₗ) ∈ calculator
    # ₛYₗ[m] is available for m ∈ -ℓ:ℓ
end
```
The blocks for ``ℓ < |s|``, whose modes do not exist, come back full
of zeros rather than being skipped.

The `m′` keywords have their counterpart in the spin weight.  Because
``{}_{s}Y_{ℓ,m}`` is the ``m' = -s`` column of ``𝔇``, a range of spin
weights is a range of columns of one recursion, and an `sYlmCalculator`
will serve several of them at once, giving its blocks a spin axis
indexed by the spin weight itself:
```julia
calculator = sYlmCalculator(R, ℓₘₐₓ, -2:2)
for (ℓ, ₛYₗ) ∈ calculator
    # ₛYₗ[s, m] for s ∈ -2:2, m ∈ -ℓ:ℓ
end
```
What the recursion costs is governed by the *largest* ``|s|`` asked
for, so reading the rest of the range out of it is close to free;
naming a single spin weight is a saving in storage and in the final
assembly rather than in the recursion itself.  A single spin weight of
such a block is `ₛYₗ[s, :]`, and [`eachℓ`](@ref)`(calculator, s)`
iterates that row directly.  [`spins`](@ref SphericalFunctions.spins)
reports the range a calculator serves, and [`spin`](@ref) the one
value when there is only one.

Note that we begin with some `R`, which sets the underlying float type
of the calculator.  For example, an input `Rotor{Float32}` gives a
`Float32` calculator, which returns matrices of `Complex{Float32}`.
All subsequent calls *must* use that same type of `Rotor`; if you need
to use a different type, you have to create a new calculator.  For
example, when differentiating with `ForwardDiff`, the input rotor must
be a dual number.  A `WignerdCalculator` is built the same way from
``β``, ``e^{iβ}``, or a rotor.  The same four keywords are accepted by
both functions.  [`floattype`](@ref
SphericalFunctions.floattype)`(calculator)` reports the type in use.
Nothing but the input can set it — there is no type argument on any of
these constructors — so computing in a wider type means building the
rotor in that type, as in `sYlmCalculator(Rotor{BigFloat}(R), ℓₘₐₓ,
s)`.

!!! danger
    Each `𝔇ˡ` block is a *view* into the storage kept in the
    calculator.  The next step of the loop overwrites it, so you
    cannot keep a block between steps unless you `copy` it.

`copy` keeps the block's natural indices, while `collect` gives an
ordinary 1-based array; `collect` applied to the calculator itself
copies every block for you.

For the same reason, two loops over one calculator cannot be
interleaved.  Each step of either loop overwrites what the other is
looking at.  There is no way to warn about this behavior; the answers
will simply be wrong if you try this.  A simple way to get a second
calculator of the same type is `similar(calculator)`.

To reset the calculator to the beginning of the loop over ``ℓ``, and
change the `R` value, call [`set_R!`](@ref) on a `WignerDCalculator`
or an `sYlmCalculator`, or [`set_β!`](@ref) on a `WignerdCalculator`.
For example, given a collection of rotors, you can iterate over them
all like this:
```julia
calculator = WignerDCalculator(first(rotors), ℓₘₐₓ)
for R ∈ rotors
    set_R!(calculator, R)
    for (ℓ, 𝔇ˡ) ∈ calculator
        # 𝔇ˡ[m′, m] is available for m′, m ∈ -ℓ:ℓ with this value of R
    end
end
```

However, such a loop is not always the best way to handle many rotors.
You can also give a whole collection of rotors to the constructor, and
the calculator will use SIMD to compute the Wigner matrix elements for
all of them at once — which can be significantly faster than looping
over them one at a time.  This would not be accessible to the user by
external looping as above.  For example,
```julia
calculator = WignerDCalculator(rotors, ℓₘₐₓ)
for (ℓ, 𝔇ˡ) ∈ calculator
    # 𝔇ˡ[iᵣ, m′, m] is available for iᵣ ∈ 1:length(rotors) and m′, m ∈ -ℓ:ℓ
end
```
The type *and number* of rotors are fixed when the calculator is
built, but you can still use `set_R!` and `set_β!` to adjust their
values and restart the loop over ``ℓ``.  Each block of an
`sYlmCalculator` built this way gains the same leading rotor index, so
that it is `ₛYₗ[iᵣ, m]`, or `ₛYₗ[iᵣ, s, m]` for a range of spin
weights.

Finally, an `sYlmCalculator` also accepts real angles ``θ`` in place
of rotors, either at construction or later through [`set_θ!`](@ref),
and then evaluates the harmonics at ``(θ, ϕ=0)``:
```julia
calculator = sYlmCalculator(θ⃗, ℓₘₐₓ, -2)      # a vector of angles, or one θ
for (ℓ, ₛλₗ) ∈ calculator
    # ₛλₗ[iᵣ, m] for iᵣ ∈ 1:length(θ⃗), m ∈ -ℓ:ℓ
end
```
For integer spin weight the values there are real, though they are
still stored as complex numbers with zero imaginary part.  This is the
``{}_{s}λ_{ℓ,m}(θ)`` that ring-based transforms need — one ring of the
sphere for each angle, which is why a whole vector of them is the
natural input.  Angles fix the element type exactly as rotors do, and
`set_θ!` requires the same agreement as `set_R!`, so a `BigFloat`
calculator wants `big(θ)` rather than a bare literal.

## The underlying ``H`` recursion

All of these calculators are built on a [`WignerHCalculator`](@ref),
which computes the real, symmetric ``H`` wedge that the
[Gumerov–Duraiswami](@cite Gumerov_2006) recursion produces before any
phases are applied.  It is available directly for the rare cases where
the wedge itself is needed — probably as an optimization:
```julia
h = WignerHCalculator(β, ℓₘₐₓ)
for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ  # ℓₘᵢₙ=0 for integers or 1//2 for half-integers
    recurrence!(h, ℓ)
    # h.Hˡ[iᵣ, m′, m] is available; iᵣ is always present; only |m′| ≤ m ≤ ℓ is stored
end
```
This one is intentionally not iterable: it is a single mutable object
handed back.  It must also be stepped through manually with
`recurrence!` and read as `h.Hˡ`.

The wedge is stored as an [`HWedge`](@ref) (and, during the recursion,
an [`HAxis`](@ref)); the symmetries that relate the rest of the matrix
to the stored wedge are described in the notes on the [``H``
recursion](@ref "Algorithm for computing ``H``").  Most users should
prefer the ``𝔇``, ``d`` and ``{}_{s}Y_{ℓ,m}`` calculators above,
which apply the symmetries and the phases for you.


## Half-integer indices

Everything on this page works for half-integer ``ℓ, m', m`` and spin
weight ``s`` as well.  Ask for them by giving `Rational`s whose
denominator is exactly 2:
```julia
𝔇 = D(R, 7//2)                                # ℓ = 1//2, 3//2, 5//2, 7//2
𝔡 = d(β, 7//2)
calculator = sYlmCalculator(R, 7//2, -3//2:3//2)
```
The containers that come back in place of the `OffsetArray`s, the
index type behind them, and the handful of places where the
half-integer case genuinely differs are all described on the
[half-integer page](@ref interface_half_integers).


## Docstrings

```@docs
D
d
sYlm
sYlm!
Ylm
sYlm_matrix
WignerDCalculator
WignerdCalculator
WignerHCalculator
sYlmCalculator
recurrence!
eachℓ
eachell
set_R!
set_β!
set_θ!
```


## Containers

The types that `D`, `d` and `calc[ℓ]` return for half-integer ``ℓ``,
and the abstract type they share with the workspaces below.

```@docs
AbstractWignerMatrix
WignerMatrix
WignerDMatrix
WignerdMatrix
WignerMatrixBatch
WignerVector
WignerVectorBatch
SpinMatrix
SpinMatrixBatch
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
SphericalFunctions.spins
SphericalFunctions.sₘₐₓ
SphericalFunctions.sₘᵢₙ
SphericalFunctions.Nᵣ
SphericalFunctions.ishalfinteger
SphericalFunctions.isbatched
```
