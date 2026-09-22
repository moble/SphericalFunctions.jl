# Changelog

## 3.0.0

Version 3 is a major rewrite of the code, with a new interface, new
capabilities, and a new convention for Wigner's ``𝔇`` matrices.  The
convention change deserves particular attention when porting, because
code that is updated only by renaming functions will run without
complaint, but will give the complex conjugate of the intended result.
Almost everything else in version 2's interface has been replaced, and
the table under "What replaces what" maps each old function to its
successor.  The most significant additions are support for
half-integer indices throughout the package, and the `ModeWeights`
container, which can be evaluated, rotated, and acted on by the
differential operators directly.  The entries below are relative to
version 2.2.9.

### Breaking

* **``𝔇`` is conjugated.**  The convention is now
  ``𝔇^{(ℓ)}_{m',m}(α, β, γ) = e^{-i m' α}\, d^{(ℓ)}_{m',m}(β)\, e^{-i m γ}``,
  which agrees with LALSuite, Wikipedia, Sakurai, Shankar, Zettili and
  Varshalovich et al.  It is the complex conjugate of what Wigner,
  Edmonds, Goldberg et al., Boyle (2016) and every earlier version of
  this package used.  Code that rotated mode weights with the conjugate
  of version 2's ``𝔇`` should now use `D` with no conjugation — or,
  more simply, `D(R, ℓₘₐₓ) * w`.  The ``d`` matrices and the
  spin-weighted spherical harmonics are numerically unchanged.  (Issue
  #42.)  The conventions are described in full in the documentation,
  along with comparisons to other sources that are tested
  automatically.
* **Julia 1.10 or later is required**, rather than 1.6.
* **The version-2 interface is removed.**  This includes `D_matrices`,
  `D_matrices!`, `D_prep`, `D_iterator`, `d_matrices`, `d_matrices!`,
  `d_prep`, `d_iterator`, `sYlm_values`, `sYlm_values!`, `sYlm_prep`,
  `sYlm_iterator`, `λ_iterator`, `ₛ𝐘`, `H!`, `H_recursion_coefficients`,
  the associated-Legendre functions `ALFRecursionCoefficients`,
  `ALFrecurse!`, `ALFcompute!` and `ALFcompute`, the index functions
  `WignerHsize`, `WignerHindex`, `WignerHrange`, `WignerDsize`,
  `WignerDindex` and `WignerDrange`, the helpers `deduce_limits`,
  `theta_phi` and `phi_theta`, and the legacy aliases `Diterator`,
  `diterator`, `Yiterator`, `λiterator`, `d!`, `D!`, `Y!`, `dprep`,
  `Dprep` and `Yprep`.  The table below gives the replacements.  Note
  that the names `D` and `d` survive with a new meaning: version 2's
  legacy `d(expiβ, ℓₘₐₓ)` returned a flat vector, while the new
  functions return the containers described next.
* **Results are indexed by their natural labels, rather than stored in
  flat vectors.**  `D` and `d` return a `WignerSeries`, indexed as
  `𝔇[ℓ][m′, m]` with ``m′`` and ``m`` running over `-ℓ:ℓ`, and `sYlm`
  returns a `HarmonicValues`, indexed as `Y[ℓ][m]`.  There is no longer
  any index arithmetic with `WignerDindex`.  The first index of each
  block is ``m′``; version 2's `D_iterator` and `d_iterator` returned
  the transpose of the blocks their documentation described, and warned
  so on every call.  (Issues #41 and #48.)  These containers are
  deliberately not `AbstractArray`s, since their indices may be
  half-integers; linear algebra on them goes through `array_view`
  (described under "Added"), so that `𝔇₁[ℓ] * 𝔇₂[ℓ]` is written
  `array_view(𝔇₁[ℓ]) * array_view(𝔇₂[ℓ])`.
* **Rotations must be given as `Rotor`s.**  The Euler-angle and
  spherical-coordinate forms, such as `D_matrices(α, β, γ, ℓₘₐₓ)` and
  `sYlm_values(θ, ϕ, ℓₘₐₓ, s)`, have no counterparts; convert with
  `from_euler_angles(α, β, γ)` or `from_spherical_coordinates(θ, ϕ)`
  from Quaternionic.  A general `Quaternion` or a `QuatVec` is refused,
  with a message suggesting `rotor(q)` or `exp(v/2)`.  The functions of
  ``β`` alone — `d`, `dCalculator` and `HCalculator` — still accept
  ``β``, ``e^{iβ}`` or a `Rotor`.
* **The element type of a result is that of its input.**  Version 2
  chose the type through arguments such as `D_prep(ℓₘₐₓ, T)`; now there
  is no such argument, and to compute in another type the rotor (or
  angle) must be converted.  A mismatch between the input and a
  preallocated output, as in `sYlm!`, is an error rather than a silent
  conversion.  `SSHT`, the pixelizations, the quadrature weights and the
  operator matrices keep their `T` arguments, since they construct
  their own numbers.
* **The differential operators are objects rather than functions.**
  `L²`, `Lz`, `L₊`, `L₋`, `R²`, `Rz`, `R₊`, `R₋`, `ð` and `ð̄` are now
  singleton instances of subtypes of `DifferentialOperator`.  Calling
  one as in version 2, `ð(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`, returns the same matrix
  as before; what is new is described under "Added".
* **The transforms have changed defaults and names.**  `SSHT` now uses
  the `"RS"` method by default, rather than the dense matrix.
  `SSHTDirect` is renamed `SSHTMatrix`, and the method name `"Direct"`
  is accepted with a deprecation warning in favor of `"Matrix"`.
  `SSHTMatrix` factorizes with `lu` rather than `qr` when the number of
  points equals the number of modes.  Analysis of a one-dimensional set
  of function values, `𝒯 \ f`, returns a `ModeWeights` rather than a
  `Vector`.
* **`map2salm` has a new signature and output.**  It is called as
  `map2salm(map, s, ℓₘₐₓ)` or `map2salm(map, 𝒯::SSHTRS)`; the
  `show_progress` argument is gone.  The output starts at ``ℓ = |s|``
  rather than ``ℓ = 0``, and is a `ModeWeights` for a single map.
  `map2salm!` and `plan_map2salm` are removed; the plan is now an
  `SSHTRS`, which the unexported `map2salm_plan` constructs.  Because
  `map2salm` is now implemented by the ring-based transform, the
  `BoundsError` that version 2's `map2salm!` raised under
  `--check-bounds=yes` is gone.  (Issue #59.)
* The dependencies `AbstractFFTs`, `DoubleFloats`, `Hwloc`,
  `LoopVectorization` and `ProgressMeter` are dropped, and
  `FixedSizeArrays` is added.

### What replaces what

| Version 2 | Version 3 |
|---|---|
| `D_matrices(R, ℓₘₐₓ)` + `D_iterator` | `D(R, ℓₘₐₓ)`, indexed `𝔇[ℓ][m′, m]` (conjugated; see above) |
| `D_matrices(α, β, γ, ℓₘₐₓ)` | `D(from_euler_angles(α, β, γ), ℓₘₐₓ)` |
| `D_prep` + `D_matrices!` | `DCalculator(R, ℓₘₐₓ)`, iterated or stepped with `recurrence!`, and reused with `set_R!` |
| `d_matrices(β, ℓₘₐₓ)` + `d_iterator` | `d(β, ℓₘₐₓ)`, indexed `𝔡[ℓ][m′, m]` |
| `d_prep` + `d_matrices!` | `dCalculator(β, ℓₘₐₓ)`, reused with `set_β!` |
| `H!`, `H_recursion_coefficients` | `HCalculator(β, ℓₘₐₓ)` |
| `sYlm_values(R, ℓₘₐₓ, s)` + `sYlm_iterator` | `sYlm(R, ℓₘₐₓ, s)`, indexed `Y[ℓ][m]` |
| `sYlm_values(θ, ϕ, ℓₘₐₓ, s)` | `sYlm(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ, s)` |
| `sYlm_prep` + `sYlm_values!` | `sYlmCalculator(R, ℓₘₐₓ, s)` + `sYlm!(Y, calc, R)` |
| `ₛ𝐘(s, ℓₘₐₓ, T, R⃗)` | `sYlm_matrix(R⃗, ℓₘₐₓ, s)` |
| `λ_iterator` | `SphericalFunctions.sλlm` or `SphericalFunctions.sλlmCalculator` |
| `ALFcompute` and relatives | no direct replacement; `sλlm(θ, ℓₘₐₓ, 0)` gives ``Y_{ℓ,m}(θ, 0)`` |
| `WignerDindex`, `WignerHsize`, … | not needed: blocks are indexed by ``(ℓ, m′, m)`` directly |
| `SSHTDirect` | `SSHTMatrix` |
| `map2salm(map, s, ℓₘₐₓ, show_progress)` | `map2salm(map, s, ℓₘₐₓ)` |
| `plan_map2salm` + `map2salm!` | `𝒯 = map2salm_plan(map, s, ℓₘₐₓ)` + `map2salm(map, 𝒯)` |

### Added

* **Half-integer indices.**  `D`, `d`, the calculators, `sYlm`,
  `sYlm!`, `sYlm_matrix`, `Ysize`, `Yindex`, `Yrange`, `ModeWeights`,
  the differential operators, the pixelizations, and the `"RS"` and
  `"Matrix"` transforms (with `map2salm` and `salm2map`) all accept
  half-integer ``ℓ``, ``m`` and ``s``, passed as `Rational`s with
  denominator 2 — as in `D(R, 7//2)` or `SSHT(1//2, 7//2)`.  The
  results are indexed exactly as in the integer case.  They have been
  verified against two independent references to ``10^{-16}`` for
  ``ℓ ≤ 31/2``, and by identities that need no reference to
  ``ℓ = 101/2``.  `Ylm` and the `"Minimal"` transform remain
  integer-only, and say so; a call that mixes integer and half-integer
  indices is refused.  (Issue #29.)
* **`ModeWeights`**, which holds the mode weights of a spin-weighted
  function in the canonical ordering, together with its spin weight and
  range of ``ℓ``.  It is indexed as `w[ℓ, m]` or `w[ℓ, :]`, and `modes`
  and `spin` report its labels.  It supports
  - evaluation at a rotor or a vector of rotors, `w(R)`, which is the
    same as the product `Y * w` with harmonics from `sYlm` or an
    `sYlmCalculator` (this is written `*` rather than `⋅`, since `dot`
    conjugates its first argument, and `dot` on these types raises an
    error saying so);
  - active rotation, `D(R, ℓₘₐₓ) * w` (or with a `DCalculator` in place
    of `D`), which gives the weights of ``f′(𝐐) = f(𝐑^{-1}𝐐)``; and
  - the differential operators, written `ð(w)` or `ð * w`, which return
    a `ModeWeights` with the spin weight adjusted appropriately.  These
    are applied by a loop rather than by building a matrix, so the only
    allocation is the result, and `mul!(w′, ð, w)` allocates nothing.

  Rotation also has an in-place form, `mul!(w′, 𝔇, w)`.  The
  transforms accept a `ModeWeights` directly, and check its range of
  ``ℓ`` against their own.
* **Calculators** — `DCalculator`, `dCalculator`, `HCalculator`,
  `sYlmCalculator` and `YlmCalculator` — take the same arguments as the
  corresponding functions, and compute one ``ℓ`` at a time, so that
  large ``ℓₘₐₓ`` can be reached without holding every matrix at once.
  Iterating one, as in `for (ℓ, 𝔇ˡ) ∈ calc`, yields each block as a
  view that the next step overwrites, and allocates nothing;
  `recurrence!(calc, ℓ)` computes a single step, and `collect` copies
  every block.  `set_R!`, `set_β!` and `set_θ!` point an existing
  calculator at new data.  (Issue #32.)
* **Batched evaluation.**  Given a vector of rotors, the calculators and
  `sYlm` evaluate all of them at once, with the rotor as the leading
  index of each block.  This is two to ten times faster per rotor than
  a loop, because the rotor index is the only one that the recursions
  allow to be vectorized.  The transforms use it internally.
* **Spinor phases.**  The phases of a rotor that the recursions need
  are computed directly from its components as
  ``e^{i(α±γ)/2}`` and ``\cos(β/2)``, ``\sin(β/2)``, rather than
  through Euler angles.  This is accurate near both poles, and gives
  ``𝔇(-R) = -𝔇(R)`` for half-integer indices.  (Issue #57.)
* **Restricted ranges.**  The keywords `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ` and
  `mₘᵢₙ` of `D`, `d` and their calculators limit the part of each matrix
  that is computed, and `sYlm` and its relatives accept an `ℓₘᵢₙ`
  keyword and an ascending range of spin weights such as `-2:2`.
* **`array_view` and `relabel`.**  `array_view(x)` gives the contents of
  a container as a 1-based `StridedArray` that aliases its storage, so
  that BLAS and LAPACK can operate on it without copying; for mode
  weights and harmonics this is the flat array in the canonical
  ordering.  `relabel(x, A)` puts the natural indices back on a plain
  array.  `Matrix`, `Array` and `collect` give copies.
* `Ylm`, the ordinary scalar spherical harmonics, and `sYlm_matrix`,
  the dense synthesis matrix.
* The real harmonics
  ``{}_sλ_{ℓ,m}(θ) = {}_sY_{ℓ,m}(θ, 0) / i^{2s}``, through the public
  but unexported `sλlm`, `sλlm!`, `sλlm_matrix` and `sλlmCalculator`.
* The angular-momentum operators `Lx` and `Ly`.  (There is deliberately
  no `Rx` or `Ry`; see the `Lx` docstring.)
* `salm2map`, the inverse of `map2salm`.
* The pixelizations `driscoll_healy_pixels`, `driscoll_healy_rotors`,
  `mcewen_wiaux_pixels` and `mcewen_wiaux_rotors`, which are public but
  unexported.
* `ComplexPowers`, an iterator over the powers of a unit complex
  number.
* The accessors `ℓₘᵢₙ`, `ℓₘₐₓ`, `spins`, `Nᵣ`, `floattype` and others,
  with ASCII aliases such as `ellmax`, are declared public but not
  exported.

### Fixed

* The size functions throw for an inverted range of ``ℓ``, rather than
  returning a negative number.  (Issue #52.)
* `complex_powers!` works for wrapper element types such as
  `ForwardDiff.Dual`, including at the phase 1, where it used to
  produce `NaN` derivatives.
* The ring-based transform handles rings with different numbers of
  points, and rings with more than ``2ℓₘₐₓ+1`` points; version 2 gave
  wrong results in both cases.
