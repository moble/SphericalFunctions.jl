# Changelog

## 3.0.0 (unreleased)

Version 3 is a rewrite.  The interface is new, the `Deprecated` module
that held version 2 is gone, and — most important for anyone upgrading
— **the convention for Wigner's ``𝔇`` matrices is the complex
conjugate of the one version 2 used.**

### Breaking

### Added

* **Rotating mode weights: `𝔇 * w`.**  Given the Wigner matrices of a rotor — as a
  `WignerSeries` from `D`, or a `DCalculator`, which streams one ℓ at a time — this
  gives the weights of the actively rotated function, ``f′(𝐐) = f(𝐑^{-1}𝐐)``.  There is
  **no complex conjugate**: version 2's ``𝔇`` was the conjugate of this one, so ported code
  must drop a `conj` rather than add one.  The derivation is now written out in the
  conventions section, under "Rotation of mode weights".  `mul!(w′, 𝔇, w)` writes into an
  existing container.
* **Evaluating a function: `Y * w`.**  Given the harmonics at one or more rotors — a
  `HarmonicValues` from `sYlm`, or an `sYlmCalculator` — this gives the function's values
  there.  It is written `*` and **not** `⋅`, because `⋅` is `LinearAlgebra.dot`, which
  conjugates its first argument, and evaluation must not; `dot` on these types raises an error
  saying so rather than answering with the wrong phase.  `w(R)` and `w(R⃗)` are the same
  product, written as a call.
* **The real harmonics: `sλlm`, `sλlm!`, `sλlm_matrix` and `sλlmCalculator`.**  These give
  ``{}_sλ_{ℓ,m}(θ) = {}_sY_{ℓ,m}(θ, 0) / i^{2s}``, which is real for integer and half-odd
  spin weights alike, in the same containers as `sYlm` and with the same iteration and
  indexing.  They share one struct, `HarmonicCalculator`, with the complex family, exactly as
  `dCalculator` shares `WignerCalculator` with `DCalculator`: the ``H`` recursion
  is real either way, and only the factor ``e^{-i(mα - sγ)}`` ever made a result complex.  A
  real calculator therefore allocates no phase tables at all.  It refuses a `Rotor`, at
  construction and through `set_R!`, because a rotor specifies the angles ``α`` and ``γ``
  whose phases it has nowhere to put; `set_θ!` is its setter.

  `SSHTRS` and `SSHTMinimal` now build their ``Λ`` tables with one of these instead of
  reading the real part out of a complex table at every access in their innermost loops.  The
  values are bit-for-bit what they were, so nothing downstream changes.

### Breaking

* **The differential operators are objects, not functions.**  `L²`, `Lz`, `L₊`, `L₋`, `Lx`,
  `Ly`, `R²`, `Rz`, `R₊`, `R₋`, `ð` and `ð̄` are now zero-size singleton instances of
  subtypes of `DifferentialOperator`, so an operator knows its own effect on the spin weight
  (`Δspin`) and its own band structure.  Every existing call still works — `ð(w)` for the
  labelled result and `ð(s, ℓₘᵢₙ, ℓₘₐₓ, [T])` for the matrix — and `ð * w` is added
  alongside `ð(w)`.  Applying one to mode weights no longer builds the matrix: a loop does
  it, so `ð * w` allocates only its result, and `mul!(w′, ð, w)` allocates nothing.  The
  results are bit-for-bit what the matrix product gave, because the matrix builders and the
  loops evaluate the same coefficient functions.
* **Underscore-prefixed internal names are gone.**  Where a worker could simply be another
  method of the public function — reached by dispatch, and undocumented rather than hidden —
  it now is; where it could not, because that would have added a public overload the design
  refuses, it has an ordinary descriptive name instead.
* **`OffsetArray`s are gone from every return value.**  Through
  version 2 the integer path returned `OffsetArray`s, so that a block
  could be indexed by its natural ``m'`` and ``m``; the half-integer
  path could not, and grew its own containers.  Both paths now return
  those containers.  The reason is not uniformity but safety: an
  `OffsetArray` with non-trivial offsets *accepts* `*` and `mul!` and
  returns silently wrong answers — a product of two blocks comes back
  as a 1-based `Matrix` of mostly zeros, and an adjoint product comes
  back holding uninitialized memory.  Written the obvious way, the
  composition law ``𝔇(R₁R₂) = 𝔇(R₁)𝔇(R₂)`` was wrong in the first
  digit with no error raised.  The containers are deliberately not
  `AbstractArray`s, so the same code is now a `MethodError`.
* **`strided` and `relabel` are the route to and from linear
  algebra.**  `strided(w)` gives a 1-based `StridedArray` **aliasing**
  the container's storage, which BLAS and LAPACK take at full speed —
  contiguity is not required, only a unit leading stride, which an
  unbatched block already has.  `relabel(w, A)` puts the natural
  indices back on a plain array.  `Matrix`, `Array` and `collect`
  remain the copying forms.  So `𝔇₁[ℓ] * 𝔇₂[ℓ]` becomes
  `strided(𝔇₁[ℓ]) * strided(𝔇₂[ℓ])`.
* **`D` and `d` return a `WignerSeries` for integer ``ℓ`` too**,
  rather than an `OffsetVector` of blocks.  It knows its own `ℓₘᵢₙ`
  and `ℓₘₐₓ`, and refuses an `ℓ` it does not hold with a sentence
  rather than a `BoundsError` about axes.
* **`sYlm` returns a `HarmonicValues`,** indexed by ``ℓ`` and then
  naturally, rather than a flat `Vector` reached through `Yindex`.
  `sYlm` also accepts a vector of rotors, so a block has four possible
  shapes: `[m]`, `[iᵣ, m]`, `[s, m]` and `[iᵣ, s, m]`.  The flat array
  is `strided(sY)`, in the same canonical ordering as before, and
  `sYlm_matrix` is unchanged — it remains the direct way to ask for
  the bare synthesis array.
* **`ModeWeights` is no longer an `AbstractVector`.**  It is an
  `AbstractModeContainer`, the supertype it shares with
  `HarmonicValues`.  Indexing, broadcasting (which still preserves the
  wrapper), `map`, `similar`, `copy`, reductions, `norm`, `dot`,
  `adjoint` and products with matrices all still work; `strided(w)` is
  the flat storage, and is what the transforms take.
* **``𝔇`` is conjugated.**  The convention is now
  ``𝔇^{(ℓ)}_{m',m}(α, β, γ) = e^{-i m' α} d^{(ℓ)}_{m',m}(β) e^{-i m
  γ}``, which agrees with LALSuite, Wikipedia, Sakurai, Shankar,
  Zettili and Varshalovich et al., and is the complex conjugate of
  what Wigner, Edmonds, Goldberg et al., Boyle (2016) and versions of
  this package before 3.0 used.  This closes issue #42.  Code that
  rotated mode weights with `conj(D_matrices(...))` should now use
  `D(...)` with no conjugation.  The ``d`` matrices and the
  spin-weighted spherical harmonics are numerically unchanged.
* **Calculators take their rotor data at construction, and take it
  first** — `DCalculator(R, ℓₘₐₓ)`, `dCalculator(β, ℓₘₐₓ)`,
  `HCalculator(β, ℓₘₐₓ)`, `sYlmCalculator(R, ℓₘₐₓ, s)`.  A
  calculator is therefore usable the moment it exists, and the
  `Nᵣ` keyword is gone: a vector argument is what makes one batched.
  The old argument order is a `MethodError` at the call site.
* **`sYlmCalculator` is built for the spin weights it will serve.**
  The third argument is a spin weight, or an ascending range of them:
  `sYlmCalculator(R, ℓₘₐₓ, 2)` serves ``s = 2`` alone, and
  `sYlmCalculator(R, ℓₘₐₓ, -2:2)` serves all five.  It follows that the
  calculator has something to yield, so it iterates like the Wigner
  ones — `for (ℓ, ₛYₗ) ∈ calc` — with the block indexed `ₛYₗ[m]` in the
  first case and `ₛYₗ[s, m]` in the second.  One spin weight of such a
  block is the slice `ₛYₗ[s, :]`, and `spins` and `spin` report what a
  calculator was built for.  The flat `sYlm`, `sYlm!` and `sYlm_matrix` take ranges
  too, laying the spin weights along a new axis of a plain array.
  Half-integer ranges are written the same way, `-3//2:3//2`, and their
  blocks are the new `SpinMatrix` and `SpinMatrixBatch` containers.
* **The element type is the input's, and there is no argument to
  override it.**  The positional element type is gone from every
  calculator constructor, and the `T` keyword from `sYlm`,
  `sYlm_matrix` and `Ylm`; `sYlm!` no longer takes its working type
  from the output buffer.  A result's type is decided by its input's
  type, if and only if — to compute in another type, convert the data,
  which is also the honest way to say it, since the type of the data is
  the claim being made about the points.  `SSHT`, the pixelizations,
  the weights and the operators keep their `T` arguments, because they
  build their own numbers rather than being handed any.
* **Type mismatches are errors rather than silent conversions.**
  `sYlm!` requires `eltype(Y)` to be `Complex` of its calculator's
  float type, and `set_R!`, `set_β!`, `set_θ!` and
  `similar(calc, data)` require the new data to match the calculator's.
* **Rotations must be `Rotor`s.**  A general `Quaternion` has a
  magnitude that these functions would divide out, and a `QuatVec` is a
  vector rather than a rotation; both are refused with a message naming
  `rotor(q)` or `exp(v/2)`.  Vectors of rotor data must also have a
  concrete element type, so a `Vector{Any}` is refused rather than
  guessed at.
* **The calculators are named for their functions.**
  `WignerDCalculator`, `WignerdCalculator` and `WignerHCalculator` are
  now `DCalculator`, `dCalculator` and `HCalculator`, so that every
  calculator is its function's name plus `Calculator`, as
  `sYlmCalculator` and `YlmCalculator` already were.
* **`recurrence!` returns the block** rather than the calculator, and
  **the calculators are no longer indexed.**  `calc[ℓ]` had to be given
  the ``ℓ`` just computed and threw for any other, so it asserted what
  the caller already knew rather than looking anything up; reading a
  result by hand is now one call, `𝔇ˡ = recurrence!(calc, ℓ)`.  This is
  also how an `HCalculator` hands back its wedge, in place of the field
  access `calc.Hˡ`, and how one spin weight is reached, as
  `recurrence!(calc, ℓ)[s, :]` in place of `calc[ℓ, s]`.
* **`eachℓ` and `eachell` are removed.**  `eachℓ(calc)` was bare
  iteration under another name; a restricted range and a single spin
  weight are both a `for` loop over `recurrence!`.
* **`SphericalFunctions.Deprecated` is removed**, and with it the
  whole version-2 API: `D_matrices`, `D_prep`, `D_iterator`,
  `d_matrices`, `d_prep`, `d_iterator`, `sYlm_values`, `sYlm_prep`,
  `sYlm_iterator`, `ₛ𝐘`, `H!`, `λ_iterator`, the
  `WignerHsize`/`WignerDindex` family, and the version-2 `SSHT` types.
  The table below gives the replacements.
* **`SSHTDirect` is now `SSHTMatrix`**, and the `SSHT` method name
  `"Direct"` is accepted with a deprecation warning.  The default
  method is now `"RS"` rather than the dense matrix.
* **`map2salm` output starts at ``ℓ = |s|``** rather than ``ℓ = 0`` (issue #59).
* The dependencies `ProgressMeter` and `LoopVectorization` are
  dropped, and `Quaternionic` 4 is now allowed.

### What replaces what

| Version 2 | Version 3 |
|---|---|
| `D_matrices(R, ℓₘₐₓ)` + `D_iterator` | `D(R, ℓₘₐₓ)`, indexed `𝔇[ℓ][m′, m]` (conjugated; see above) |
| `D_prep` + `D_matrices!` | `DCalculator(R, ℓₘₐₓ)`, iterated, or stepped with `recurrence!` |
| `d_matrices(β, ℓₘₐₓ)` | `d(β, ℓₘₐₓ)`, indexed `𝔡[ℓ][m′, m]` |
| `sYlm_values(R, ℓₘₐₓ, s)` | `sYlm(R, ℓₘₐₓ, s)` |
| `sYlm_prep(ℓₘₐₓ, sₘₐₓ)` + `sYlm_values!` | `sYlmCalculator(R, ℓₘₐₓ, s)` + `sYlm!(Y, calc, R)` |
| `ₛ𝐘(s, ℓₘₐₓ, T, R⃗)` | `sYlm_matrix(R⃗, ℓₘₐₓ, s)` |
| `SSHTDirect` | `SSHTMatrix` |
| `WignerHsize`, `WignerDindex`, … | not needed: matrices are indexed by ``(ℓ, m', m)`` through views |

### Added

* **Half-integer ``ℓ``.**  `d`, `D` and the calculators accept a
  `Rational` index type, so `D(R, 7//2)` and
  `DCalculator(R, 15//2)` work, and `sYlmCalculator` accepts
  half-integer spin weights.  Verified against two independent
  references to ``10^{-16}`` for ``J ≤ 31/2`` and by oracle-free
  identities to ``J = 101/2``.  (Issue #29.)  The same `Rational`
  indices are accepted by everything built on the canonical mode
  ordering: `sYlm`, `sYlm!` and `sYlm_matrix`; `Ysize`, `Yindex`,
  `Yrange` and `ModeWeights`; the angular-momentum operators; the
  pixelizations; and the `"RS"` and `"Matrix"` transforms, with
  `map2salm` and `salm2map`.  `Ylm` and the `"Minimal"` method remain
  integer-only, and say so.  A call that mixes integer and half-integer
  indices is refused with a message naming both kinds.
* **Batched calculators.**  Give a calculator a vector of rotors and it
  evaluates all of them at once, which is two to ten times faster per
  rotor than looping, and is what the transforms now use internally.
  Batching is internal because the recurrence is sequential in every
  index a single rotation has, leaving the rotation index as the only
  one that can be vectorized.  (Issue #32.)
* **Natural indexing.**  `D(R, ℓₘₐₓ)[ℓ][m′, m]` uses the real ranges
  `-ℓ:ℓ`; there is no index arithmetic to get wrong.  (Issues #41,
  #48.)
* `ModeWeights`, a vector of mode weights that knows its spin weight
  and its ``ℓ`` range, with `w[ℓ, m]`, `w[ℓ, :]`, `modes`, `spin` and
  evaluation `w(R)`.  Everything that takes mode weights still also
  accepts a plain vector in the canonical ordering.
* **The calculators iterate over ``ℓ``.**  `for (ℓ, 𝔇ˡ) ∈ calc` yields
  `ℓ => block` pairs, one ``ℓ`` at a time, which is how to reach large
  ``ℓₘₐₓ`` without holding every matrix at once; a whole pass allocates
  nothing.  With it come `keys`, `length`, `eltype`, `pairs`, and a
  `collect` that copies every block, since the blocks themselves are
  views that the next step overwrites.  A sweep over part of the range,
  or in some other order, is a `for` loop over `recurrence!`.
* `set_R!`, `set_β!` and `set_θ!` point an existing calculator at new
  data, each named for what its calculator actually holds.
* `Ylm(R, ℓₘₐₓ)`, the ordinary scalar spherical harmonics, which are
  the spin-weight-zero case of `sYlm`.
* `floattype(calc)` reports the floating-point type a calculator works
  in.
* `sYlm_matrix`, the dense synthesis matrix, as a public function.
* The angular-momentum operators `Lx` and `Ly`.  (There is
  deliberately no `Rx` or `Ry`; see the `Lx` docstring.)
* `SSHT` objects accept any number of trailing array dimensions and
  transform them independently.

### Fixed

* Size functions no longer return negative numbers for inverted ``ℓ``
  ranges; they throw.  (Issue #52.)
* Broadcasting a function over the axis of a half-integer container —
  `Rational.(axes(w[ℓ, :], 1))`, say — gives the axis values; it used
  to misread positions as values.
* `map2salm!` no longer raises a `BoundsError`.  (Issue #59.)
* Spinor phases are computed in the calculator's precision rather than
  the rotor's, so asking for `Double64` results from `Float64` rotors
  no longer silently loses half the digits.
* `complex_powers!` works for wrapper element types such as
  `ForwardDiff.Dual`, including at the phase exactly 1, where it used
  to produce `NaN` derivatives.
* The ring-based transform handles rings with different numbers of
  points, and uniform counts above ``2ℓₘₐₓ+1``; version 2 returned
  garbage for both.
