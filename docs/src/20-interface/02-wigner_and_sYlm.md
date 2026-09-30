# [``𝔇`` and ``d`` matrices, and ``{}_{s}Y_{ℓ,m}`` and ``Y_{ℓ,m}`` functions](@id interface_wigner_matrices)

```@meta
CurrentModule = SphericalFunctions
```

Wigner's ``𝔇`` matrices — and to a lesser extent, the related ``d``
matrices — are extremely important in the theory of rotations.  Each
element is, itself, a special function of the rotation group: in
particular, an eigenfunction of [the left- and right-Lie
derivatives](@ref background_differential_operators), and thus a
spin-weighted spherical function.  See the "Background" section, and
particularly
[this page](@ref sYlm_and_Dlmpm), for details.  Collectively, they
describe how spin-weighted spherical functions transform under
rotation.  But their accurate and efficient computation is
surprisingly subtle.  This package implements the current
state-of-the-art techniques for their fast and accurate computation,
based on the [``H`` recursion](@ref "Algorithm for computing ``H``")
introduced by [Gumerov_2015](@citet).

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
recursion, and because only the single column ``m = -s`` is needed —
the second index of ``𝔇`` — both the storage and the work are much
smaller than for the full matrices.  The standard (scalar) spherical
harmonics are the special case of spin weight ``s = 0``,
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

For the most common of these, the conversion is built in: `D`,
`DCalculator`, `sYlm`, `Ylm`, `sYlmCalculator` and `YlmCalculator`
also accept Euler angles or spherical coordinates in place of a single
rotor.  So `D(α, β, γ, ℓₘₐₓ)` is `D(from_euler_angles(α, β, γ),
ℓₘₐₓ)`, and `sYlm(θ, ϕ, ℓₘₐₓ, s)` is
`sYlm(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ, s)`, with exactly the
same values.

The result is indexed by ``ℓ`` and then by the two matrix indices,
each with its natural range:
```julia
𝔇[ℓ][m′, m]  # for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ and m′,m ∈ -ℓ:ℓ
```
There is no index arithmetic to get wrong: `𝔇[ℓ]` is a
[`WignerMatrix`](@ref) whose axes are just `-ℓ:ℓ`, and the outer
container is a [`WignerSeries`](@ref).  Note that `ℓₘᵢₙ=0` for
integers but `ℓₘᵢₙ=1//2` for half-integers.

A block is deliberately **not** an `AbstractMatrix`, so linear algebra
does not apply to it directly; [`array_view`](@ref) gives a 1-based
`StridedArray` view of the same storage, which BLAS takes at full
speed, and [`relabel`](@ref) puts the natural indices back on the
result.  `Matrix(𝔇[ℓ])` gives an independent copy.  The reasons for
that arrangement are set out under [Containers](@ref
interface_containers) below.

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
`mₘₐₓ`, and `mₘᵢₙ`, which may also be passed as `mp_max`, `mp_min`,
`m_max` and `m_min` — restricting the block of each matrix that is
returned; each lower limit defaults to minus the corresponding upper
one, so that `D(R, ℓₘₐₓ; m′ₘₐₓ=2)` has the rows ``-2 ≤ m' ≤ 2``.
Restricting either `m′` (the rows) or `m` (the columns) to a narrow
band makes the calculation cheaper as well as smaller, because the
recurrence at each ``ℓ`` costs in proportion to the narrower of the
two ranges.  Each range must contain 0 — or both ``±1/2``, for
half-integer indices — because that is where the recurrence starts;
a single row or column elsewhere, such as ``m = 2``, is read from the
block afterwards.  (The spin-weighted spherical harmonics need a
single column, and are computed most cheaply by `sYlm`, below.)

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

The order of these arguments follows one rule throughout the package.
A function evaluated at a rotor takes the rotor first (or the angles
that stand in for it), then ``ℓₘₐₓ``, and then the spin weight when
there is one: `D(R, ℓₘₐₓ)`,
`sYlm(R, ℓₘₐₓ, s)`, `sYlm_matrix(R⃗, ℓₘₐₓ, s)`,
`sYlmCalculator(R, ℓₘₐₓ, s)` and the rest of that family, with ``ℓₘᵢₙ``
as a keyword.  An object labelled by a spin weight takes the spin
weight first, then ``ℓₘᵢₙ`` where it may be given, and then ``ℓₘₐₓ``:
`ModeWeights{T}(undef, s, ℓₘᵢₙ, ℓₘₐₓ)`, `SSHT(s, ℓₘₐₓ)`, the
differential operators `op(s, ℓₘᵢₙ, ℓₘₐₓ)`, and the pixelizations such
as `leja_rotors(s, ℓₘₐₓ)`.  Exchanging the arguments of one family
for those of the other asks for a spin weight larger than ``ℓₘₐₓ``,
unless the two are equal, and is usually refused for that reason.  The
exception is a spin weight larger by exactly one, for which the
operators and the `ModeWeights` constructors return an empty result,
since the range ``|s| ≤ ℓ ≤ ℓₘₐₓ`` is then empty rather than invalid.

A harmonic has only one index besides ``ℓ``, so a block is a vector
rather than a matrix, but the result is indexed the same way as
``𝔇``: by ``ℓ`` first, then naturally.  It is a
[`HarmonicValues`](@ref), and `sY[ℓ][m]` is one value; `sY[ℓ, :]` is
equivalent to the block `sY[ℓ]`.  As for ``𝔇``, iterating over the
result gives `ℓ => block` pairs, while `first`, `last` and `only` give
blocks.  A whole collection of rotors may be given instead of one, and
the spin weight may be a range, which between them give a block four
possible shapes:

| built for | `sY[ℓ]` is indexed |
|---|---|
| one rotor, one spin weight | `[m]` |
| many rotors, one spin weight | `[iᵣ, m]` |
| one rotor, a range of spin weights | `[s, m]` |
| many rotors, a range of spin weights | `[iᵣ, s, m]` |

Underneath, the values are held in one array whose *last* axis is the
modes in the canonical ordering described by [`Ysize`](@ref),
[`Yindex`](@ref) and [`Yrange`](@ref), and whose leading axes are the
rotors and spin weights.  [`array_view`](@ref) hands that array back:

```julia
array_view(sY)[Yindex(ℓ, m, abs(s))] == sY[ℓ][m]
```

That flat form is what a product with the plain vector of mode weights
takes, to synthesize a function at the rotors; [`sYlm_matrix`](@ref)
is the direct name for it, for those who want the bare array.  The
weights themselves are held in a [`ModeWeights`](@ref), which the
labelled `sY` multiplies directly, as described [below](@ref
mode_weight_operations).

Modes with ``ℓ < |s|`` do not exist (or are inherently zero), so by
default ``ℓ`` starts at ``|s|``.  Pass `ℓₘᵢₙ=0` (or `ell_min=0`) to
start at ``ℓ = 0`` instead, with zeros in the nonexistent modes; this
is the layout that some downstream packages use for every spin weight
at once.  For ``s = 0`` these are the ordinary scalar spherical
harmonics ``Y_{ℓ,m}``, which [`Ylm`](@ref) gives without the redundant
argument: `Ylm(R, ℓₘₐₓ)` is exactly `sYlm(R, ℓₘₐₓ, 0)`, and starts at
``ℓ = 0`` because no modes are missing there.


## Iterating over ``ℓ`` and reusing the storage

The functions above allocate a fresh workspace on every call, and copy
the results *for all ``ℓ`` values* out of it.  When the total size is
too large, or the values are needed for many rotors, allocate the
workspace once as a [`DCalculator`](@ref) (or a
[`dCalculator`](@ref)) and iterate over the values of ``ℓ``:
```julia
calculator = DCalculator(R, ℓₘₐₓ)
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
A calculator always starts at ``ℓ = 0`` (or ``1/2``), and takes no
`ℓₘᵢₙ` keyword; the blocks for ``ℓ < |s|``, whose modes do not exist,
come back full of zeros rather than being skipped.

The `m` keywords have their counterpart in the spin weight.  Because
``{}_{s}Y_{ℓ,m}`` is the ``m = -s`` column of ``𝔇``, a range of spin
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
such a block is `ₛYₗ[s, :]`, which is written the same way whichever
kind of index the calculator has.
[`spins`](@ref SphericalFunctions.spins)
reports the range a calculator serves, and [`spin`](@ref) the one
value when there is only one.

Note that we begin with some `R`, which sets the underlying float type
of the calculator.  For example, an input `Rotor{Float32}` gives a
`Float32` calculator, which returns matrices of `Complex{Float32}`.
All subsequent calls *must* use that same type of `Rotor`; if you need
to use a different type, you have to create a new calculator.  For
example, when differentiating with `ForwardDiff`, the input rotor must
be a dual number.  A `dCalculator` is built the same way from
``β``, ``e^{iβ}``, or a rotor.  The same four keywords are accepted by
both functions.  [`floattype`](@ref
SphericalFunctions.floattype)`(calculator)` reports the type in use.
Nothing but the input can set it — there is no type argument on any of
these constructors — so computing in a wider type means building the
rotor in that type, as in `sYlmCalculator(Rotor{BigFloat}(R), ℓₘₐₓ,
s)`.

!!! note "Derivatives at the poles"
    ``𝔇`` and the harmonics may be differentiated with respect to the
    rotor by automatic differentiation everywhere, including at rotors
    with ``β = 0`` or ``β = π``.  The same is not true of ``d`` and
    ``H`` *of a rotor*: some of their elements, such as
    ``d^{(1)}_{1,0}``, have no derivative at the poles as functions of
    the rotor, because ``β`` itself has none, and those derivatives
    are returned as `NaN`.  To differentiate ``d`` at the poles, give
    it the angle ``β`` or the phase ``e^{iβ}`` instead.  [This
    note](@ref automatic_differentiation) explains the details.

!!! danger
    Each `𝔇ˡ` block is a *view* into the storage kept in the
    calculator.  The next step of the loop overwrites it, so you
    cannot keep a block between steps unless you `copy` it.  The same
    applies to [`array_view`](@ref) of a block, which aliases that
    storage rather than copying it; `Matrix`, `Array` and `collect`
    are the forms that survive.

`copy` keeps the block's natural indices, while `collect` gives an
ordinary 1-based array; `collect` applied to the calculator itself
copies every block for you.

A calculator is not indexed, and there is no restricted form of the
iteration.  Both are the same `for` loop over [`recurrence!`](@ref),
which computes one ``ℓ`` and returns its block:
```julia
for ℓ ∈ 2:4
    𝔇ˡ = recurrence!(calculator, ℓ)
    # 𝔇ˡ[m′, m] for m′, m ∈ -ℓ:ℓ
end
```
Beginning above the calculator's own ``ℓₘᵢₙ`` costs nothing in
accuracy: the recursion runs through the values below either way, and
the result is bit-for-bit what a full sweep gives.  Values of ``ℓ``
taken in *decreasing* order are a different matter — each one restarts
the recursion from ``ℓₘᵢₙ``, so a loop that reads two neighboring
``ℓ`` together pays that restart at every step.

For the same reason, two loops over one calculator cannot be
interleaved.  Each step of either loop overwrites what the other is
looking at.  There is no way to warn about this behavior; the answers
will simply be wrong if you try this.  A simple way to get a second
calculator of the same type is `similar(calculator)`, or
`similar(calculator, R)` to give it other rotors.  The same holds for
tasks: a calculator is a mutable workspace, so one calculator must
never be used by two tasks at the same time, and each task that
computes in parallel with the others needs a calculator of its own.

To reset the calculator to the beginning of the loop over ``ℓ``, and
change the `R` value, call [`set_R!`](@ref) on a `DCalculator`
or an `sYlmCalculator`, or [`set_β!`](@ref) (also available as
`set_beta!`) on a `dCalculator`.
For example, given a collection of rotors, you can iterate over them
all like this:
```julia
calculator = DCalculator(first(R⃗), ℓₘₐₓ)
for R ∈ R⃗
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
calculator = DCalculator(R⃗, ℓₘₐₓ)
for (ℓ, 𝔇ˡ) ∈ calculator
    # 𝔇ˡ[iᵣ, m′, m] is available for iᵣ ∈ 1:length(R⃗) and m′, m ∈ -ℓ:ℓ
end
```
The type *and number* of rotors are fixed when the calculator is
built, but you can still use `set_R!` and `set_β!` to adjust their
values and restart the loop over ``ℓ``.  Each block of an
`sYlmCalculator` built this way gains the same leading rotor index, so
that it is `ₛYₗ[iᵣ, m]`, or `ₛYₗ[iᵣ, s, m]` for a range of spin
weights.  What decides this is that the rotors came as a vector, not
how many there are: a calculator built from a vector of one rotor is
still a batch, with a rotor index of length one, and
[`isbatched`](@ref SphericalFunctions.isbatched) says which kind a
calculator is.  The operations that need a single rotor — rotating a
`ModeWeights`, for example — refuse a batch of one, rather than
guessing which was meant.

## [The real harmonics ``{}_sλ_{ℓ,m}(θ)``](@id interface_real_harmonics)

A calculator also accepts real angles ``θ`` in place of rotors, either
at construction or later through [`set_θ!`](@ref) (also available as
`set_theta!`), and then evaluates the harmonics at ``(θ, ϕ=0)``.  That
is the ``{}_{s}λ_{ℓ,m}(θ)`` the ring-based transforms need — one ring of the
sphere for each angle, which is why a whole vector of them is the
natural input:
```julia
calculator = SphericalFunctions.sλlmCalculator(θ⃗, ℓₘₐₓ, -2)  # a vector of angles, or one θ
for (ℓ, ₛλₗ) ∈ calculator
    # ₛλₗ[iᵣ, m] for iᵣ ∈ 1:length(θ⃗), m ∈ -ℓ:ℓ
end
```
An [`sλlmCalculator`](@ref) stores its values as *real* numbers, and
[`sλlm`](@ref), [`sλlm!`](@ref) and [`sλlm_matrix`](@ref) are the flat
forms of it.  These names are public but not exported, and each has an
ASCII alias, `slambdalmCalculator`, `slambdalm`, `slambdalm!` and
`slambdalm_matrix`.  Everything else is as it is for the complex
family: the same blocks, the same iteration, the same half-integer
types, and the same containers, which are generic in the number type.

The two flavors share one struct, [`HarmonicCalculator`](@ref),
exactly as [`DCalculator`](@ref) and [`dCalculator`](@ref) do — and
for the same reason.  The underlying ``H`` recursion is real either
way; it is only the factor ``e^{-i(mα - sγ)}`` that ever makes a
result complex, and an angle sets ``α = γ = 0``.  So the real flavor
runs precisely the same recursion, allocates no phase tables at all,
and writes half as many numbers.

The definition is

```math
{}_sλ_{ℓ,m}(θ) = \begin{cases}
    {}_sY_{ℓ,m}(θ, 0), & s ∈ ℤ, \\
    {}_sY_{ℓ,m}(θ, 0) \big/ i^{2s}, & s ∈ ℤ + \tfrac{1}{2},
\end{cases}
```

which is real for both kinds of index.  For an integer spin weight the
prefactor ``(-1)^s`` in the definition of ``{}_sY_{ℓ,m}`` is ``\pm
1``, so ``{}_sY_{ℓ,m}(θ, 0)`` is already real, and ``{}_sλ_{ℓ,m}`` is
simply the harmonic at ``ϕ = 0``, as the literature writes it.  For a
half-odd spin weight that prefactor is ``i^{2s} = \pm i``, so
``{}_sY_{ℓ,m}(θ, 0)`` is imaginary rather than real, and dividing that
constant phase out is what leaves a real function behind.  (Dividing
by ``i^{2s}`` in both cases would give the wrong sign for odd integer
``s``, where ``i^{2s} = -1``.)

A `Rotor` is refused, by the constructor and by [`set_θ!`](@ref)
alike: it specifies the angles ``α`` and ``γ``, whose phases a real
calculator has nowhere to put.  Use an `sYlmCalculator` for that.
Conversely, [`set_R!`](@ref) takes only rotors, and refuses an angle
with a message naming `set_θ!`.  Angles fix the element type exactly
as rotors do, and `set_θ!` requires the same agreement as `set_R!`, so
a `BigFloat` calculator takes `big(θ)` rather than a bare literal.

## The underlying ``H`` recursion

All of these calculators are built on a [`HCalculator`](@ref),
which computes the real, symmetric ``H`` wedge that the
[Gumerov–Duraiswami](@cite Gumerov_2015) recursion produces before any
phases are applied.  It is available directly for the rare cases where
the wedge itself is needed — probably as an optimization:
```julia
h = HCalculator(β, ℓₘₐₓ)
for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ  # ℓₘᵢₙ=0 for integers or 1//2 for half-integers
    Hˡ = recurrence!(h, ℓ)
    # Hˡ[iᵣ, m′, m] is available; iᵣ is always present; only |m′| ≤ m ≤ ℓ is stored
end
```
This one is intentionally not iterable, and is the exception to
everything said above about blocks: the wedge is a single mutable
object, the calculator's own workspace, handed back by identity rather
than as a view.  The next step overwrites it, and its `ℓ` must not be
reassigned by hand.  Read the values out before stepping on, or keep
`copy(Hˡ)`, which is an independent wedge holding the same numbers.

The wedge is stored as an [`HWedge`](@ref) (and, during the recursion,
an [`HAxis`](@ref SphericalFunctions.HAxis), which is internal); the
symmetries that relate the rest of the matrix to the stored wedge are
described in the notes on the [``H`` recursion](@ref "Algorithm for
computing ``H``"), and [`wedge_value`](@ref) reads any element of the
matrix through them.  Most users should prefer the ``𝔇``, ``d`` and
``{}_{s}Y_{ℓ,m}`` calculators above, which apply the symmetries and
the phases for you.


## Half-integer indices

Everything on this page works for half-integer ``ℓ, m', m`` and spin
weight ``s`` as well.  Ask for them by giving `Rational{Int}`s whose
denominator is exactly 2, or [`HalfOddInteger`](@ref
SphericalFunctions.HalfOddInteger)s:
```julia
𝔇 = D(R, 7//2)                                # ℓ = 1//2, 3//2, 5//2, 7//2
𝔡 = d(β, 7//2)
calculator = sYlmCalculator(R, 7//2, -3//2:3//2)
```
The containers that come back are the same ones as for integer
indices.  The index type behind them, and the handful of places where
the half-integer case actually differs, are described on the
[half-integer page](@ref interface_half_integers).


## [Rotating and evaluating mode weights](@id mode_weight_operations)

Two products tie the containers together.  Multiplying a
[`ModeWeights`](@ref) by the Wigner matrices of a rotor rotates it,
and multiplying by the harmonics at some rotors evaluates it there:

```julia
w′ = D(R, ℓₘₐₓ) * w              # the weights of f′(𝐐) = f(𝐑⁻¹𝐐)
f  = sYlm(R, ℓₘₐₓ, s) * w         # the value of f at R
f⃗  = sYlm(R⃗, ℓₘₐₓ, s) * w        # ... and at each of many rotors
```

The calculator forms, `DCalculator(R, ℓₘₐₓ) * w` and
`sYlmCalculator(R, ℓₘₐₓ, s) * w`, compute the same things one ``ℓ`` at a
time rather than materializing every block.

Evaluation reads only the weights that belong to a function of spin
weight ``s``, those with ``ℓ ≥ \max(ℓₘᵢₙ(w), |s|)``, so the harmonics
need to cover only that range, and weights that hold nothing at or
above ``|s|`` evaluate to zero.  The flat [`sYlm_matrix`](@ref) is a
plain matrix, which has no labels to check, so its product with a
`ModeWeights` is refused rather than trusted.  The same synthesis is
written either with the labelled harmonics, as `sYlm(R⃗, ℓₘₐₓ, s) * f̃`,
which checks the spin weight and the range of ``ℓ``, or with the
matrix and the raw numbers, as `Y * array_view(f̃)`, for weights stored
from the matrix's own ``ℓₘᵢₙ``.  Weights over some other range of
``ℓ`` are copied into the one wanted with `ModeWeights(f̃; ℓₘᵢₙ,
ℓₘₐₓ)`, which fills the modes that `f̃` lacks with zeros.

Rotation is an ordinary matrix–vector product on each ``ℓ`` block, with
**no complex conjugate** — see [Rotation of mode
weights](@ref conv_rotation_of_modes) for the derivation.  Because
version 2 used the conjugate convention for ``𝔇``, code ported from it
must *drop* a `conj` rather than add one.  The matrices must cover the
weights' range of ``ℓ``, which is not the same as matching it, and
their blocks must be whole: a ``𝔇`` built with the `m′ₘₐₓ` or `mₘₐₓ`
restrictions cannot rotate anything, because every ``m`` mixes into
every ``m′``.

!!! warning
    Evaluation is written `*` and **not** `⋅`.  `⋅` is
    `LinearAlgebra.dot`, which conjugates its first argument;
    evaluation must not conjugate the harmonics.  Calling `dot` on
    these types raises an error saying so, rather than quietly
    returning an answer with the wrong phase.  For the conjugating
    inner product of two sets of mode weights, `dot(w₁, w₂)` is still
    what you want.

```@autodocs
Modules = [SphericalFunctions]
Pages = ["mode_weights/products.jl"]
```

## Docstrings

```@docs
D
d
sYlm
sYlm!
Ylm
sYlm_matrix
DCalculator
dCalculator
HCalculator
sYlmCalculator
YlmCalculator
sλlm
sλlm!
sλlm_matrix
sλlmCalculator
HarmonicCalculator
recurrence!
set_R!
set_β!
set_θ!
```


## [Containers](@id interface_containers)

The types that `D`, `d`, `sYlm` and [`recurrence!`](@ref) return, for
either kind of index, and the abstract types they share with the workspaces below.

These containers are deliberately **not** `AbstractArray`s.  Half-odd
indices cannot satisfy that interface at all — `axes` must be integer
ranges, and `-3//2:3//2` is not one — but the reason they are not
arrays on the integer path either is a sharper one.  There the natural
array would be an `OffsetArray`, and an `OffsetArray` with non-trivial
offsets *accepts* `*` and `mul!` and returns silently wrong answers:
a product of two blocks comes back as a 1-based `Matrix` of mostly
zeros, and an adjoint product comes back holding uninitialized memory.

!!! note "How to get results as `Array`s"
    The package provides [`array_view`](@ref) to get a 1-based
    `StridedArray` *view* of the storage, and [`relabel`](@ref) to put
    the natural indices back on the result.  `Matrix(𝔇[ℓ])` gives an
    independent *copy* of the storage.  The reason for that
    arrangement is discussed in the [Containers](@ref
    interface_containers) section.

```@docs
array_view
relabel
AbstractWignerMatrix
WignerMatrix
WignerDMatrix
WignerdMatrix
WignerMatrixBatch
DegreeBlock
DegreeBlockBatch
SpinMatrix
SpinMatrixBatch
WignerSeries
HarmonicValues
SphericalFunctions.AbstractModeContainer
WignerCalculator
```


## Workspaces

The wedge that an [`HCalculator`](@ref) returns is an `HWedge`, and
[`wedge_value`](@ref) reads any element of ``H^ℓ`` from it, applying
the symmetries described in the notes on the [``H`` recursion](@ref
"Algorithm for computing ``H``").  The axis that seeds the wedge
during the recursion is an internal type, described on the [internal
page](@ref "Internal functions").

```@docs
HWedge
wedge_value
```


## Methods of `Base` functions

Indexing a calculator or a container, copying one, emptying one, gathering every ``ℓ`` of
one, and converting one to an ordinary `Array` are all written with the usual `Base`
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

The same handful of names reports the index ranges of every container
and calculator in the package.  Each has an ASCII alias, given in its
docstring, for use where the subscripted Unicode names are
inconvenient: `ell`, `ell_min`, `ell_max`, `mp_max`, `mp_min`,
`m_max`, `m_min`, `s_max`, `s_min` and `Nr`.  The keyword arguments of
the same names are spelled the same way in ASCII, so that `D(R, ℓₘₐₓ;
mp_max=2)` is `D(R, ℓₘₐₓ; m′ₘₐₓ=2)` and `sYlm(R, ℓₘₐₓ, s; ell_min=0)`
is `sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=0)`.  The functions whose names are not
ASCII have aliases too — `set_beta!`, `set_theta!`, `slambdalm`,
`slambdalm!`, `slambdalm_matrix` and `slambdalmCalculator` here, and
those of the [differential operators](@ref
interface_differential_operators) — which are public but not exported,
and are mentioned in the docstrings of the functions they name.

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
