# [``{}_{s}Y_{ℓ,m}`` functions](@id interface_sYlm)

```@meta
CurrentModule = SphericalFunctions
```

The spin-weighted spherical harmonics are an [important set of
functions defined on](@cite Boyle_2016) the rotation group
``𝐒𝐎(3)``, or more generally, the spin group ``𝐒𝐩𝐢𝐧(3)`` that
covers it.  They are eigenfunctions of [the left- and right-Lie
derivatives](@ref "Differential operators"), and are particularly
useful in describing the angular dependence of polarized fields, like
the electromagnetic field and gravitational-wave field.  Originally
introduced by [Newman_1966](@citet), they are essentially components
of Wigner's ``𝔇`` matrices:
```math
{}_{s}Y_{ℓ,m}(𝐑)
  = (-1)^s \sqrt{\frac{2ℓ+1}{4π}} \, \overline{𝔇^{(ℓ)}_{m, -s}(𝐑)}.
```
(See the [conventions summary](@ref summary_swsh) for this and the
related definitions.)  As such, they can be computed with the same
[``H`` recursion](@ref "Algorithm for computing ``H`` (redesigned)") algorithm as
the Wigner ``𝔇`` matrices.  But because only the single column
``m' = -s`` is needed, both the storage and the work are much smaller
than for the full matrices.


## Evaluating the harmonics

The basic call takes a rotor, a maximum ``ℓ`` and a spin weight, and
returns all the values for that spin weight:
```julia
using Quaternionic
using SphericalFunctions

R = randn(RotorF64)
ℓₘₐₓ = 8
s = -2
Y = sYlm(R, ℓₘₐₓ, s)
```
As for the ``𝔇`` matrices, the rotation must be a `Rotor` rather than any
quaternion that could be normalized into one; see [Evaluating the
matrices](@ref interface_wigner_matrices).

The result is a single `Vector` of `Complex` numbers whose base type
is `R`'s own — the input is the only thing that decides it — ordered
by ``ℓ`` and then by ``m``:
```julia
[Y[Yindex(ℓ, m, abs(s))] for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ] == Y
```
Modes with ``ℓ < |s|`` do not exist, so by default the vector starts
at ``ℓ = |s|``.  Pass `ℓₘᵢₙ=0` to get a vector that starts at
``ℓ = 0`` instead, with zeros in the nonexistent modes; this is the
layout that some downstream packages use for every spin weight at
once.  For ``s = 0`` these are the ordinary scalar spherical harmonics
``Y_{ℓ,m}``, which [`Ylm`](@ref) gives without the redundant
argument: `Ylm(R, ℓₘₐₓ)` is exactly `sYlm(R, ℓₘₐₓ, 0)`, and starts at
``ℓ = 0`` because no modes are missing there.  See [`Ysize`](@ref),
[`Yindex`](@ref) and [`Yrange`](@ref) for the indexing functions, and
[`ModeWeights`](@ref) for a wrapper that does the index arithmetic for
you.


## Reusing the workspace

`sYlm` allocates its workspace on every call.  When the harmonics are
needed at many rotors — the usual case — allocate the workspace once
as an [`sYlmCalculator`](@ref) and step through the values of ``ℓ``
yourself.  The calculator takes its rotor at construction, so it is
ready to use the moment it exists:
```julia
calculator = sYlmCalculator(first(rotors), ℓₘₐₓ, abs(s))
for R ∈ rotors
    set_R!(calculator, R)
    for (ℓ, ₛYₗ) ∈ eachℓ(calculator, s)
        # ₛYₗ[m] for m ∈ -ℓ:ℓ
    end
end
```
[`set_R!`](@ref) replaces the rotor of an existing calculator and
discards whatever it was holding, so the loop that follows starts from
the lowest ``ℓ`` again; [`eachℓ`](@ref) runs the recursion one ``ℓ`` at
a time, yielding `ℓ => block` pairs.  The blocks for ``ℓ < |s|``, whose
modes do not exist, come back full of zeros rather than being skipped.
The new rotor must have the same floating-point type as the one the
calculator was built from — that type is fixed at construction, as
described below — and a mismatch is an error rather than a conversion.

Unlike the Wigner calculators, which iterate directly, an
`sYlmCalculator` has to be told the spin weight.  One calculator serves
every spin weight up to the `sₘₐₓ` it was built with, at no extra cost
in the recursion, so a calculator built with `sₘₐₓ = 2` can produce the
``s = 0``, ``s = ±1`` and ``s = ±2`` harmonics from the same workspace;
there is therefore no single thing that bare iteration could yield, and
`eachℓ(calculator, s)` says which spin weight this loop is about.

When several spin weights are wanted at once, that shared recursion is
the whole point, and the way to use it is to read them all out of the
current ``ℓ`` with `calculator[ℓ, s]`:
```julia
calculator = sYlmCalculator(R, ℓₘₐₓ, 2)
for ℓ ∈ 0:ℓₘₐₓ
    recurrence!(calculator, ℓ)
    for s ∈ -2:2
        ₛYₗ = calculator[ℓ, s]  # ₛYₗ[m] for m ∈ -ℓ:ℓ
    end
end
```
This is the lower-level form of the same loop: [`recurrence!`](@ref)
advances the calculator to the next ``ℓ``, and indexing it reads one
spin weight of the result.  Either way the block is a view into the
calculator's storage, valid only until the next step; `copy` it (which
keeps the natural indices) or `collect` it (which gives an ordinary
1-based `Vector`) if you need to keep it.  Calling `collect` on the
iterator instead copies every block for you.

The element type is the rotor's own, so a `Rotor{Float32}` gives a
`Float32` calculator and blocks of `Complex{Float32}`, and
[`floattype`](@ref SphericalFunctions.floattype)`(calculator)` reports the type in use.
Nothing else can set it — there is no type argument on the
constructor — so to compute in a wider type, build the rotor in that
type or convert it: `sYlmCalculator(Rotor{BigFloat}(R), ℓₘₐₓ, sₘₐₓ)`.
That is also the honest way to say it, since the type of the data is
the claim being made about the points.

If the flat, mode-ordered vector is what the rest of the code wants,
[`sYlm!`](@ref) will fill one from a calculator, replacing the
calculator's rotor itself as it goes:
```julia
calculator = sYlmCalculator(first(rotors), ℓₘₐₓ, abs(s))
Y = Vector{ComplexF64}(undef, Ysize(abs(s), ℓₘₐₓ))
for R ∈ rotors
    sYlm!(Y, calculator, R, s)
    # Do something with Y before the next iteration overwrites it
end
```
The element type of `Y` must be `Complex` of the calculator's float
type — equivalently, of the rotor's — and a mismatch is an error
rather than a silent conversion, because the rotor decides the
arithmetic while `Y` is where the answer lands, so a disagreement
between them means one of the two is not what the caller thinks it is.
`ComplexF64` is right for the `Float64` rotors used throughout these
examples; for any other type, ask the calculator rather than guessing.
The same rule governs the calculator-free `sYlm!(Y, R, ℓₘₐₓ, s)`,
where the rotor alone decides.

Handing the constructor an `AbstractVector` of rotors in place of a
single one gives a calculator that evaluates the whole batch at once,
with each block gaining a leading rotor index:
```julia
calculator = sYlmCalculator(rotors, ℓₘₐₓ, 2)
for (ℓ, ₛYₗ) ∈ eachℓ(calculator, -2)
    # ₛYₗ[iᵣ, m] for iᵣ ∈ 1:length(rotors), m ∈ -ℓ:ℓ
end
```
The number of rotors is fixed when the calculator is built, so
`set_R!` then takes a vector of exactly that many.  Batching this way
is substantially faster per rotor than looping over a single-rotor
calculator, for a structural reason — the recursion is sequential in
every index that a single rotor has, which leaves the rotor index as
the only one a processor can vectorize over — explained in full under
"Reusing the workspace" on the [Wigner matrices](@ref
interface_wigner_matrices) page.

Finally, the calculator also accepts real angles ``θ`` in place of
rotors, either at construction or later through [`set_θ!`](@ref), and
then evaluates the harmonics at ``(θ, ϕ=0)``:
```julia
calculator = sYlmCalculator(θ⃗, ℓₘₐₓ, 2)      # a vector of angles, or one θ
for (ℓ, ₛλₗ) ∈ eachℓ(calculator, -2)
    # ₛλₗ[iᵣ, m] for iᵣ ∈ 1:length(θ⃗), m ∈ -ℓ:ℓ
end
```
For integer spin weight the values there are real (they are still
stored as complex numbers, with zero imaginary part).  This is the
``{}_{s}λ_{ℓ,m}(θ)`` that ring-based transforms need — one ring of the
sphere for each angle, which is why a whole vector of them is the
natural input.  Angles fix the element type exactly as rotors do, and
`set_θ!` requires the same agreement as `set_R!`: an angle of any other
floating-point type is an error, so a `BigFloat` calculator wants
`big(θ)` rather than a bare literal.  Note that it is the calculator
that accepts an angle: the flat `sYlm` and `sYlm!` take a rotor.


## The dense matrix of harmonics

Synthesis — evaluating a function from its mode weights — is a
matrix-vector product with the matrix of harmonics evaluated at the
sample points.  [`sYlm_matrix`](@ref) builds that matrix:
```julia
R⃗ = golden_ratio_spiral_rotors(s, ℓₘₐₓ)
𝐘 = sYlm_matrix(R⃗, ℓₘₐₓ, s)
f = 𝐘 * f̃          # synthesis
f̃ = lu(𝐘) \ f      # analysis, if there are as many points as modes
```
Its rows are exactly the vectors that `sYlm` returns:
`sYlm_matrix(R⃗, ℓₘₐₓ, s)[i, :] == sYlm(R⃗[i], ℓₘₐₓ, s)`.  Here too the
rotors decide the element type, so a matrix in a wider type comes from
rotors in that type.  The pixelizations do take a type argument,
because they create the points rather than being handed them, so
`sYlm_matrix(golden_ratio_spiral_rotors(s, ℓₘₐₓ, BigFloat), ℓₘₐₓ, s)`
gives a `Complex{BigFloat}` matrix.  This is the matrix behind the
`"Matrix"` method of [`SSHT`](@ref); see
[Transformations](@ref interface_transformations) for the transforms
built on it, which are usually what you want.


## Half-integer indices

An [`sYlmCalculator`](@ref) also handles half-integer ``ℓ``, ``m`` **and** spin weight
``s``, exactly as the Wigner calculators do (see [Half-integer indices](@ref
half_integer_wigner)).  Give it a half-integer `ℓₘₐₓ` and `sₘₐₓ`, spelled as `Rational`s with
denominator 2:
```julia
calculator = sYlmCalculator(R, 7//2, 3//2)
for ℓ ∈ 1//2:7//2
    recurrence!(calculator, ℓ)
    for s ∈ -3//2:3//2
        ₛYₗ = calculator[ℓ, s]   # ₛYₗ[m] for m ∈ -ℓ:ℓ
    end
end
```
`calc[ℓ, s]` then returns a [`WignerVector`](@ref) (or, for a
calculator built over several rotors, a [`WignerVectorBatch`](@ref))
in place of the `OffsetArray` the integer path gives, since
`OffsetArray` cannot use half-integer axes; both are indexed the same
way, and support the same `copy`, `collect`, `Matrix`/`Vector`,
`similar`, `==`, `size` and iteration as the Wigner containers.

Two things differ mathematically.  The prefactor ``(-1)^s`` in the definition above is
``\pm i`` for half-integer ``s``; this package uses the principal branch
``(-1)^s ≡ e^{iπs} = i^{2s}`` settled in the [conventions summary](@ref summary_swsh).
Consequently ``{}_{s}λ_{ℓ,m}(θ)`` is *not* real for half-integer ``s``, even though it is
for integer ``s``.

The flat interfaces — [`sYlm`](@ref), [`sYlm!`](@ref) and [`sYlm_matrix`](@ref) — remain
integer-only, as do [`Ysize`](@ref), [`Yindex`](@ref), [`Yrange`](@ref),
[`ModeWeights`](@ref) and the transforms, because the canonical mode ordering they share is
defined only for integer ``ℓ`` and ``m``.  [`Ylm`](@ref) has no half-integer analogue for a
second reason as well: half-integer ``ℓ`` goes with half-integer ``s``, and never with
``s = 0``.


## Docstrings

```@docs
sYlm
sYlm!
Ylm
sYlmCalculator
set_θ!
sYlm_matrix
SphericalFunctions.sₘₐₓ
```

(Indexing a calculator with `calc[ℓ, s]` is documented with the other `Base` methods, and
[`eachℓ`](@ref) and [`set_R!`](@ref) with the rest of the calculator machinery, on the
[Wigner matrices](@ref interface_wigner_matrices) page.)


## Containers

```@docs
WignerVector
WignerVectorBatch
```
