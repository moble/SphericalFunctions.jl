# Internal functions

```@meta
CurrentModule = SphericalFunctions
```

The functions on this page are not part of the public interface: their
names are not exported, they are not declared `public`, and their
signatures may change without a breaking release.  They are documented
because the public docstrings refer to them, and because anyone
reading the recurrence or the raw storage of a calculator needs them.

Everything here belongs to one machine.  An [`HCalculator`](@ref) holds
two kinds of buffer — an [`HWedge`](@ref) holding about a quarter of
``H^ℓ`` for a batch of rotors, and two [`HAxis`](@ref) buffers holding
successive orders of the ``m'=0``, ``m ≥ 0`` axis that seeds it — and
[`recurrence!`](@ref) walks them from one ``ℓ`` to the next through
the [``H`` recursion](@ref "Algorithm for computing ``H``") below.
The axis is always at integer order, labelled by `axis_ℓ`: for integer
``ℓ`` it is copied straight into the wedge's ``m'=0`` row, and for
half-integer ``ℓ``, where there is no such row, it seeds the two rows
``m' = ±1/2`` instead.  Everything the package computes — ``d``,
``𝔇``, ``{}_{s}Y_{ℓ,m}``, and hence the transforms — is that
recursion plus a phase.


## The ``H`` recursion

The recurrence is described by [Gumerov_2015](@citet), in the form
derived in the notes on the [``H`` recursion](@ref "Algorithm for
computing ``H``"), which is where the recurrence coefficients
themselves are written out.  Each step fills one region of ``H^ℓ``
from regions already filled: steps 1 and 2 run the ``ℓ`` ladder along
the ``m'=0`` axis, steps 3–5 run the ``m'`` ladders at fixed ``ℓ`` to
fill the wedge ``m ≥ |m'|``, and step 6 uses the symmetries to fill
the rest of the requested ``m'`` range.

The steps documented below are what an `HCalculator` runs.  Each takes an
[`HCalculator`](@ref) and works on its buffers for all `Nᵣ` rotors at
once, with the rotor index innermost: steps 1 and 2 fill its two axis
buffers, and steps 3–5 fill its batched quarter-wedge.  They handle
half-integer as well as integer ``ℓ``; for half-integer ``ℓ``,
[`recurrence_seed!`](@ref) takes the place of step 3.  Step 6 is
*never* applied to the wedge, because every element outside it is
supplied by the symmetries, including the sign ``σ``, when a block is
assembled from the wedge, or when an element is read through
[`wedge_value`](@ref).  The phases of step 7 — the ``ϵ`` signs and the
Euler phases ``e^{-im'α}`` and ``e^{-imγ}`` — enter in the
`materialize!` of each calculator, which writes every element of a
block through the one `materialize_element!` of
`src/calculators/engine.jl`.  The test suite, in
`test/wigner/recurrence.jl`, checks the `HCalculator` against a second
implementation of the same recurrence, which holds one whole ``H^ℓ``
of integer order for one rotor in a dense matrix and applies all seven
steps to it.

```@docs
SphericalFunctions.recurrence_step1!
SphericalFunctions.recurrence_step2!
SphericalFunctions.recurrence_step3!
SphericalFunctions.recurrence_seed!
SphericalFunctions.recurrence_step4!
SphericalFunctions.recurrence_step5!
```


## Rotor phases

The recurrence needs ``\cos(β/2)`` and ``\sin(β/2)`` and the two
spinor phases, never the Euler angles themselves: extracting ``α``,
``β``, ``γ`` from a rotor and re-forming ``e^{i(m'α+mγ)}`` would both
lose accuracy near the poles and throw away the rotor's double-cover
sign, which is exactly what half-integer ``ℓ`` needs.

```@docs
SphericalFunctions.spinor_phases
SphericalFunctions.half_angles
SphericalFunctions.axis_ℓ
```


## Reading the ``H`` wedge

Only about a quarter of each ``H^ℓ`` matrix is stored (see the notes
on the [``H`` recursion](@ref "Algorithm for computing ``H``")); every
other element is obtained from a stored one by symmetry, together with
a sign ``σ`` that is ``+1`` for integer indices and
``\mathrm{sgn}(m)\,\mathrm{sgn}(m')`` for half-integer ones.  The
public [`wedge_value`](@ref), documented with [`HWedge`](@ref), reads
one element that way, and the functions below define where it comes
from.  The calculators' assembly of their blocks applies the same
cases, a run of elements at a time rather than one by one, and is
tested element by element against `wedge_value`.

```@docs
SphericalFunctions.wedge_source
SphericalFunctions.wedge_offset
SphericalFunctions.transpose_sign
SphericalFunctions.HAxis
```

The recurrences are written in the natural indices ``ℓ``, ``m'`` and
``m`` throughout.  That costs nothing, because a sum or difference of
two [`HalfOddInteger`](@ref SphericalFunctions.HalfOddInteger)s is an
`Integer`: a coefficient such as ``(ℓ-m)(ℓ+m+1)`` is therefore integer
arithmetic for both index types, while being written exactly as the
references write it.  That type, and the `IntegerHalf` union it
belongs to, are described on the [half-integer page](@ref
interface_half_integers), together with [`@index_methods`](@ref
SphericalFunctions.@index_methods), the macro with which every public
function that takes an index is defined.  The methods it generates
convert a half-integer index passed as a `Rational{Int}` to a
`HalfOddInteger` before anything is computed, and refuse any other
index that is not an `Int` with the message built by
[`index_argument_error`](@ref
SphericalFunctions.index_argument_error).  Indexing a container is
not written with the macro, and neither are [`recurrence!`](@ref) and
[`wedge_value`](@ref), which take an index beside a calculator or a
wedge, because the kind of index they accept is fixed by the
container, the calculator, or the wedge rather than by the arguments
of the call.  Each of them converts an index as the macro does, and
refuses one that is not of the kind of the container, the calculator,
or the wedge with a message of the same form, so that it accepts
exactly the indices that the functions defined with the macro accept.

```@docs
SphericalFunctions.index_argument_error
SphericalFunctions.δ²
SphericalFunctions.sgn
SphericalFunctions.ϵ
```
