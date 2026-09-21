# Introduction

```@meta
CurrentModule = SphericalFunctions
DocTestSetup = quote
    using SphericalFunctions, Quaternionic
end
```

This is a Julia package for evaluating and transforming Wigner's 𝔇
matrices, and spin-weighted spherical harmonics ``{}_{s}Y_{ℓ,m}``
(which includes the ordinary scalar spherical harmonics).  Because
[*both* 𝔇 *and* the harmonics are most correctly considered](@cite
Boyle_2016) functions on the rotation group ``𝐒𝐎(3)`` — or more
generally, the spin group ``𝐒𝐩𝐢𝐧(3) \cong 𝐒𝐔(2)`` that covers it
— these functions are evaluated directly in terms of quaternions.
Concessions are also made for more standard forms of spherical
coordinates and Euler angles.[^1] Among other applications, those
functions permit "synthesis" (evaluation of the spin-weighted
spherical functions) of spin-weighted spherical harmonic coefficients
on regular or distorted grids.  This package also includes functions
enabling efficient "analysis" (decomposition into mode coefficients)
of functions evaluated on regular grids to high order and accuracy.


## Quick start

A handful of functions cover most first uses of the package.  Each
returns values for *every* ``ℓ`` up to a given ``ℓₘₐₓ``, rather than
for one ``ℓ`` at a time, because the recursion relations described in
the next section produce them that way.

The most direct of them is [`D`](@ref), which gives Wigner's
``𝔇^{(ℓ)}_{m',m}`` matrices for a single rotation.  Its result is
indexed by ``ℓ`` first, and then by the two matrix indices, each of
which runs over its own natural range — so there is no index
arithmetic to get wrong:

```jldoctest quickstart
julia> using SphericalFunctions, Quaternionic

julia> R = from_spherical_coordinates(π/3, π/4);  # the point (θ, ϕ), as a rotation

julia> ℓₘₐₓ = 8;

julia> 𝔇 = D(R, ℓₘₐₓ);

julia> size(𝔇[3])  # each block is (2ℓ+1)×(2ℓ+1)
(7, 7)

julia> axes(𝔇[3])  # and is indexed as 𝔇[ℓ][m′, m], with m′, m ∈ -ℓ:ℓ
(-3:1:3, -3:1:3)

julia> 𝔇[1][0, 0] ≈ cos(π/3)  # for m′ = m = 0 the phases drop out, leaving d = cos β
true
```

The [Wigner matrix interface](@ref interface_wigner_matrices)
describes what `D` returns in full, including the half-integer case.

Wigner's ``d^{(ℓ)}_{m',m}`` matrices are the factor of ``𝔇`` that
depends on the single Euler angle ``β``, and they are real.
Accordingly, [`d`](@ref) takes that angle in place of a rotation, and
is otherwise used exactly like `D`:

```jldoctest quickstart
julia> 𝔡 = d(π/3, ℓₘₐₓ);

julia> eltype(𝔡[3])
Float64

julia> round(𝔡[2][1, -1], digits=12)
0.5
```

The same page documents `d` beside `D`, along with the other ways of
giving ``β``.

The spin-weighted spherical harmonics are given by [`sYlm`](@ref),
which takes the spin weight ``s`` as a third argument.  A harmonic has
only one index besides ``ℓ``, so a block here is a vector rather than
a matrix, but the result is indexed the same way: by ``ℓ`` first, and
then naturally.  The lower limit defaults to ``|s|``, because every
harmonic below it vanishes.  The spin-weight-zero case is common
enough to have its own name, [`Ylm`](@ref):

```jldoctest quickstart
julia> Y = sYlm(R, ℓₘₐₓ, -2);  # spin weight -2, so ℓ starts at 2

julia> repr(Y)
"HarmonicValues{ComplexF64} for ℓ ∈ 2:8, s = -2"

julia> axes(Y[3])  # one ℓ, indexed by m ∈ -ℓ:ℓ
(-3:1:3,)

julia> sum(abs2, array_view(Y[3])) ≈ 7 / (4π)  # Σₘ |ₛYₗₘ|² = (2ℓ+1)/4π
true

julia> Ylm(R, ℓₘₐₓ)[0][0] ≈ 1 / √(4π)  # here ℓ starts at 0, and Y₀₀ = 1/√(4π)
true
```

Underneath, the values are held in one flat array in the canonical
ordering of mode weights, `[ₛYₗₘ for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]`,
and [`array_view`](@ref) hands that array back.  It is what a product
with a vector of mode weights takes, and [`Yindex`](@ref) gives the
position of any one mode in it:

```jldoctest quickstart
julia> length(array_view(Y)) == Ysize(2, ℓₘₐₓ)
true

julia> array_view(Y)[Yindex(3, -3, 2)] == Y[3][-3]
true
```

The [interface page](@ref interface_wigner_matrices) describes both
functions, alongside ``𝔇`` and ``d`` themselves.

Each of the functions above allocates its entire result and fills it,
which is convenient at moderate ``ℓₘₐₓ`` and wasteful at large
``ℓₘₐₓ``: storing every matrix up to ``ℓₘₐₓ`` takes ``O(ℓₘₐₓ^3)``
numbers, while the recursion that produces them needs only the current
``ℓ``, which is ``O(ℓₘₐₓ^2)``.  A *calculator* holds just that part of
the storage.  Iterating one walks through the values of ``ℓ`` in turn,
handing back the block for each, which is how to reach large ``ℓₘₐₓ``
without ever holding every matrix at once:

```jldoctest quickstart
julia> calc = DCalculator(R, ℓₘₐₓ);

julia> norms = Float64[];

julia> for (ℓ, 𝔇ˡ) ∈ calc
           push!(norms, sum(abs2, 𝔇ˡ))  # each 𝔇ˡ is unitary, so this is 2ℓ+1
       end

julia> norms ≈ [2ℓ + 1 for ℓ ∈ 0:ℓₘₐₓ]
true
```

The block returned by each iterations a view into the calculator,
which the next step overwrites, so `copy` it if it has to outlive the
iteration.  The calculators themselves, the `set_R!` family that
points an existing calculator at a new rotation, and the restricted
form of the iteration are all described on the [Wigner matrix
interface](@ref interface_wigner_matrices) page; an
[`sYlmCalculator`](@ref) does the same for the harmonics.

Finally, a calculator built from a *vector* of rotations evaluates all
of them together, adding the rotation index to the front of each
block.  Handing the whole batch to the library, rather than writing
the loop over rotations yourself, is what makes that worth doing: the
recursion is sequential in every index a single rotation has — each
``ℓ`` comes from the one before, and each element of a matrix from its
neighbours — so the rotation index is the only one along which the
same arithmetic can be done independently:

```jldoctest quickstart
julia> rotors = [from_spherical_coordinates(θ, π/4) for θ ∈ range(0, π, 8)];

julia> batch = DCalculator(rotors, ℓₘₐₓ);

julia> for (ℓ, 𝔇ˡ) ∈ batch
           @assert axes(𝔇ˡ) == (1:8, -ℓ:ℓ, -ℓ:ℓ)  # now indexed as 𝔇ˡ[iᵣ, m′, m]
       end

julia> size(recurrence!(batch, ℓₘₐₓ))
(8, 17, 17)
```

How much that arrangement is worth, and the structure of the recursion
it follows from, are described under [Reusing the workspace](@ref
interface_wigner_matrices).

Everything so far computes the harmonics themselves.  They are a
basis, so the other half of the story is the coefficients against
them.  A [`ModeWeights`](@ref) holds the ``f_{ℓ,m}`` of a
spin-weighted function ``f = \sum_{ℓ,m} f_{ℓ,m}\, {}_sY_{ℓ,m}``, in the
canonical ordering described above, together with the spin weight and
the range of ``ℓ`` they belong to.  Keeping those labels beside the
numbers is what lets the operations below know what they are acting
on, and refuse a combination that means nothing:

```jldoctest quickstart
julia> w = ModeWeights{ComplexF64}(undef, -2, 4);  # spin weight -2, so ℓ runs over 2:4

julia> w .= 0; w[2, -1] = 0.5; w[3, 2] = im;

julia> spin(w), length(modes(w)), w[2, -1]
(-2, 21, 0.5 + 0.0im)
```

Evaluating the function those weights describe is a call.  Writing the
sum out instead — as a product with the harmonics at the same point —
gives the same number, and is the form to prefer when the harmonics
are already in hand or are wanted for many sets of weights:

```jldoctest quickstart
julia> w(R) ≈ sYlm(R, 4, -2) * w
true

julia> w(R) ≈ sYlmCalculator(R, 4, -2) * w  # one ℓ at a time, for repeated use
true
```

Multiplying by Wigner matrices rotates the function rigidly.  The
result is a new `ModeWeights` with the same spin weight and the same
range of ``ℓ``, since a rotation changes neither, and it satisfies the
property that defines it:

```jldoctest quickstart
julia> Q = from_spherical_coordinates(π/5, π/6);

julia> w′ = D(R, 4) * w;

julia> w′(Q) ≈ w(inv(R) * Q)
true
```

The angular-momentum operators are applied the same way, as
[`ð`](@ref)`(w)` or `ð * w`.  Each gives a new `ModeWeights`, with the
spin weight adjusted where the operator changes it — the point of
labelling the weights in the first place, since ``ð`` maps a function
of spin weight ``s`` to one of spin weight ``s+1``:

```jldoctest quickstart
julia> spin(ð * w), spin(ð̄ * w)
(-1, -3)

julia> (L² * w)[3, 2] == 3 * (3 + 1) * w[3, 2]  # L² is diagonal, with eigenvalue ℓ(ℓ+1)
true
```

No matrix is built for any of this: the operator is applied by a loop,
so the only allocation is the result, and `mul!` into an existing
container allocates nothing at all.  The rest of what a `ModeWeights`
supports is described under [rotating and evaluating mode
weights](@ref mode_weight_operations), and the full list of operators
— along with the matrix forms of them, which act on a plain vector
instead — is on the [differential operators](@ref
interface_differential_operators) page.

These quantities are computed using recursion relations, which makes
it possible to compute to very high ℓ values.  Unlike direct
evaluation of individual elements, which would generally cause
overflow or underflow beyond ℓ≈30 when using double precision
(`Float64`), these recursion relations should be valid for far higher
ℓ values.  More precisely, when using *this* package, `Inf` values
appear starting at ℓ=128 for `Float16`, but I have not yet found any
for values up to at least ℓ=1024 with `Float32`, and presumably far
higher for `Float64`.  `BigFloat` also works, and presumably will not
overflow for any ℓ value that could reasonably fit into computer
memory — though it is far slower.  Also note that
[`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl) will
work, and achieve significantly greater accuracy (but no greater ℓ
range) than `Float64`.  In all cases, results are typically accurate
to roughly ℓ times the precision of the underlying float type.

Half-integer ``ℓ, m', m`` — the representations of ``𝐒𝐩𝐢𝐧(3)``
that do not descend to ``𝐒𝐎(3)`` — are supported throughout: by
``𝔇``, ``d`` and the spin-weighted harmonics, by the mode weights and
the operators on them, and by the transforms, at the same accuracy and
essentially the same speed as integer indices.  Pass a `Rational` with
denominator 2, as in `D(R, 7//2)` or `SSHT(1//2, 7//2)`.  See
[Half-integer indices](@ref interface_half_integers), and the
[half-integer section of the transforms page](@ref
transformations_half_integer) for what a function of half-integer spin
weight is a function *of*.

The conventions for this package diverge from its predecessors found
[here](https://moble.github.io/spherical_functions/) and
[here](https://moble.github.io/spherical/), but are described in
detail on [this page](@ref Summary) and the following pages, including
detailed comparisons to other sources that are tested automatically
with each change to this code.

Note that numerous other packages cover some of these use cases,
including
[`FastTransforms.jl`](https://JuliaApproximation.github.io/FastTransforms.jl/),
[`FastSphericalHarmonics.jl`](https://eschnett.github.io/FastSphericalHarmonics.jl/dev/),
[`WignerSymbols.jl`](https://github.com/Jutho/WignerSymbols.jl), and
[`WignerFamilies.jl`](https://github.com/xzackli/WignerFamilies.jl).
However, I need support for quaternions (via
[`Quaternionic.jl`](https://github.com/moble/Quaternionic.jl)) and for
higher-precision numbers — even at the cost of a very slight decrease
in speed in some cases — which are what this package provides.

[^1]: Euler angles are quite generally a very poor choice for
    computing with rotations.  (The only context in which they may be
    preferred is when *analytically* integrating some analytically
    known functions.)  Almost universally, it is best to use
    quaternions when computing with rotations.  All the computations
    done within this package use quaternions; the user interfaces
    involving Euler angles essentially convert to/from quaternions.
    While the calculations needed for those conversions would still
    need to be done if this package used Euler angles internally —
    meaning that this approach is as efficient as any — that work can
    be avoided entirely if you work with quaternions directly.
