# Internal functions

```@meta
CurrentModule = SphericalFunctions
```

The functions on this page are not part of the public interface: their
names are not exported, they are not declared `public`, and their
signatures may change without a breaking release.  They are documented
because the public docstrings refer to them, and because anyone
reading the recurrence or the raw storage of a calculator needs them.

Everything here belongs to one machine.  A [`HCalculator`](@ref)
holds two buffers — an [`HWedge`](@ref) holding about a quarter of
``H^ℓ`` for a batch of rotors, and an [`HAxis`](@ref) holding the
``m'=0``, ``m ≥ 0`` axis that seeds it — and [`recurrence!`](@ref)
walks them from one ``ℓ`` to the next through the [``H``
recursion](@ref "Algorithm for computing ``H`` (redesigned)") below.
The axis is always at integer order, labelled by `axis_ℓ`: for integer
``ℓ`` it is copied straight into the wedge's ``m'=0`` row, and for
half-integer ``ℓ``, where there is no such row, it seeds the two rows
``m' = ±1/2`` instead.  Everything the package computes — ``d``,
``𝔇``, ``{}_{s}Y_{ℓ,m}``, and hence the transforms — is that
recursion plus a phase.


## The ``H`` recursion

The recurrence is described [Gumerov_2015](@citet), in the form
derived in the notes on the [``H`` recursion](@ref "Algorithm for
computing ``H`` (redesigned)"), which is where the recurrence
coefficients themselves are written out.  Each step fills one region
of ``H^ℓ`` from regions already filled: steps 1 and 2 run the ``ℓ``
ladder along the ``m'=0`` axis, steps 3–5 run the ``m'`` ladders at
fixed ``ℓ`` to fill the wedge ``m ≥ |m'|``, and step 6 uses the
symmetries to fill the rest of the requested ``m'`` range.

Each step name below has **two** methods, and they are not
interchangeable.

  - The methods taking a [`HCalculator`](@ref) are what the
    engine runs.  They work on the batched quarter-wedge for all `Nᵣ`
    rotors at once, with the rotor index innermost, and they handle
    half-integer as well as integer ``ℓ``.  There are only five of
    them: step 6 is *never* applied to the wedge, because every
    element outside it is read on demand through
    [`wedge_value`](@ref), which already contains the symmetry sign
    ``σ``.  They are internal to `src/wigner/wigner_H_calculator.jl`
    and have no docstrings.
  - The methods documented below take a single
    [`AbstractWignerMatrix`](@ref) holding one whole ``H^ℓ`` for one
    rotor, with ``\cos β`` and ``\sin β`` passed explicitly, and are
    restricted to integer indices.  They are the unbatched reference
    form of the same recurrence — a second implementation, useful for
    checking the engine, and tested against it by the "Dense H
    recurrence vs the batched engine" test item — but nothing in the
    package calls them.

The final two functions apply the phases of step 7, converting a
filled ``H^ℓ`` in place into ``d^ℓ`` or ``𝔇^ℓ``.  They are the
unbatched counterpart of the calculator's `materialize!`, which is
where the ``ϵ`` signs and the Euler phases ``e^{-im'α}``, ``e^{-imγ}``
actually enter for the engine (and, for ``{}_{s}Y_{ℓ,m}``, in
`src/sYlm/sYlm.jl`).

```@docs
SphericalFunctions.recurrence_step1!
SphericalFunctions.recurrence_step2!
SphericalFunctions.recurrence_step3!
SphericalFunctions.recurrence_step4!
SphericalFunctions.recurrence_step5!
SphericalFunctions.recurrence_step6!
SphericalFunctions.convert_H_to_d!
SphericalFunctions.convert_H_to_D!
```


## Rotor phases

The recurrence needs ``\cos(β/2)`` and ``\sin(β/2)`` and the two spinor phases, never the
Euler angles themselves: extracting ``α``, ``β``, ``γ`` from a rotor and re-forming
``e^{i(m'α+mγ)}`` would both lose accuracy near the poles and throw away the rotor's
double-cover sign, which is exactly what half-integer ``ℓ`` needs.

```@docs
SphericalFunctions.spinor_phases
SphericalFunctions.half_angles
SphericalFunctions.axis_ℓ
```


## Reading the ``H`` wedge

Only about a quarter of each ``H^ℓ`` matrix is stored (see the notes on the [``H``
recursion](@ref "Algorithm for computing ``H`` (redesigned)")); every other element is
obtained from a stored one by symmetry, together with a sign ``σ`` that is ``+1`` for integer
indices and ``\mathrm{sgn}(m)\,\mathrm{sgn}(m')`` for half-integer ones.  The functions below
are the only place in the package where those symmetries are encoded.  They are not part of
the public interface, but they are what [`HCalculator`](@ref)'s docstring refers to,
and anyone reading `calc.Hˡ` directly needs them.

```@docs
SphericalFunctions.wedge_value
SphericalFunctions.wedge_source
SphericalFunctions.wedge_offset
SphericalFunctions.transpose_sign
```

The recurrences are written in the natural indices ``ℓ``, ``m'`` and ``m`` throughout.
That costs nothing, because a sum or difference of two
[`HalfOddInteger`](@ref SphericalFunctions.HalfOddInteger)s is an `Integer`: a coefficient
such as ``(ℓ-m)(ℓ+m+1)`` is therefore integer arithmetic for both index types, while being
written exactly as the references write it.  That type, and the `IntegerHalf` union it
belongs to, are described on the [half-integer page](@ref interface_half_integers); the
functions below are the boundary machinery that produces them.  A public entry point accepts
a half-integer index passed as a `Rational`, and normalizes it with [`half_integer`](@ref
SphericalFunctions.half_integer) — or, where several indices arrive together and must be of
one kind, with [`half_integers`](@ref SphericalFunctions.half_integers) — before anything is
computed.

```@docs
SphericalFunctions.half_integer
SphericalFunctions.half_integers
SphericalFunctions.isindex
SphericalFunctions.δ²
SphericalFunctions.sgn
SphericalFunctions.ϵ
```
