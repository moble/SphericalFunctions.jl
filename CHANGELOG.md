# Changelog

## 3.0.0

Version 3 is a major rewrite of the code, with a new interface, new
capabilities, and a new convention for Wigner's ``𝔇`` matrices.  The
convention change deserves particular attention when porting, because
code that is updated only by renaming functions will run without
complaint, but will give the complex conjugate of the intended result.
Almost everything else in version 2's interface has been replaced, and
the table under "What replaces what" maps each old function to its
successor.  The most significant additions are a complete account of
the package's conventions, derived from first principles and compared
in detail with thirty other sources from Euler (1767) to current
software; support for half-integer indices throughout the package; and
the `ModeWeights` container, which can be evaluated, rotated, and
acted on by the differential operators directly.  The entries below
are relative to version 2.2.9.

### Breaking

* **``𝔇`` is conjugated.**  The convention is now
  ``𝔇^{(ℓ)}_{m',m}(α, β, γ) = e^{-i m' α}\, d^{(ℓ)}_{m',m}(β)\, e^{-i
  m γ}``, which agrees with LALSuite, Wikipedia, Sakurai, Shankar,
  Zettili and Varshalovich et al.  It is the complex conjugate of what
  Wigner, Edmonds, Goldberg et al., Boyle (2016) and every earlier
  version of this package used.  Code that rotated mode weights with
  the conjugate of version 2's ``𝔇`` should now use `D` with no
  conjugation — or, more simply, `D(R, ℓₘₐₓ) * w`.  The ``d`` matrices
  and the spin-weighted spherical harmonics are numerically unchanged.
  (Issue #42.)  The reasons for the change, and comparisons with the
  other sources, are in the new conventions documentation described
  under "Added".
* **Julia 1.10 or later is required**, rather than 1.6.
* **The version-2 interface is removed.**  This includes `D_matrices`,
  `D_matrices!`, `D_prep`, `D_iterator`, `d_matrices`, `d_matrices!`,
  `d_prep`, `d_iterator`, `sYlm_values`, `sYlm_values!`, `sYlm_prep`,
  `sYlm_iterator`, `λ_iterator`, `ₛ𝐘`, `H!`,
  `H_recursion_coefficients`, the associated-Legendre functions
  `ALFRecursionCoefficients`, `ALFrecurse!`, `ALFcompute!` and
  `ALFcompute`, the index functions `WignerHsize`, `WignerHindex`,
  `WignerHrange`, `WignerDsize`, `WignerDindex` and `WignerDrange`,
  the helpers `deduce_limits`, `theta_phi` and `phi_theta`, and the
  legacy aliases `Diterator`, `diterator`, `Yiterator`, `λiterator`,
  `d!`, `D!`, `Y!`, `dprep`, `Dprep` and `Yprep`.  The table below
  gives the replacements.  Note that the names `D` and `d` survive
  with a new meaning: version 2's legacy `d(expiβ, ℓₘₐₓ)` returned a
  flat vector, while the new functions return the containers described
  next.
* **Results are indexed by their natural labels, rather than stored in
  flat vectors.**  `D` and `d` return a `WignerSeries`, indexed as
  `𝔇[ℓ][m′, m]` with ``m′`` and ``m`` running over `-ℓ:ℓ`, and `sYlm`
  returns a `HarmonicValues`, indexed as `Y[ℓ][m]`.  Both iterate as
  `ℓ => block` pairs, exactly as the calculators described under
  "Added" do.  There is no longer any index arithmetic with
  `WignerDindex`.  The first index of each block is ``m′``; version
  2's `D_iterator` and `d_iterator` returned the transpose of the
  blocks their documentation described, and warned so on every call.
  (Issues #41 and #48.)  These containers are deliberately not
  `AbstractArray`s, since their indices may be half-integers; linear
  algebra on them goes through `array_view` (described under "Added"),
  so that `𝔇₁[ℓ] * 𝔇₂[ℓ]` is written `array_view(𝔇₁[ℓ]) *
  array_view(𝔇₂[ℓ])`.  The labels count wherever the containers meet:
  `==`, `isequal`, `hash` and `≈` compare the kind of block, its ``ℓ``
  and its ranges of indices as well as its numbers, and a broadcast
  over blocks, which returns an ordinary `Array`, refuses operands
  whose labels differ.  `first`, `last` and `only` of a series give
  blocks, as `Y[ℓ]` and `Y[ℓ, :]` do.
* **Rotations must be given as quaternions, or as angles.**  `D`,
  `DCalculator`, `sYlm`, `Ylm`, `sYlmCalculator` and `YlmCalculator`
  accept Euler angles or spherical coordinates in place of a single
  rotor, as in `D(α, β, γ, ℓₘₐₓ)` and `sYlm(θ, ϕ, ℓₘₐₓ, s)`;
  everything else takes a `Rotor`, converted with
  `from_euler_angles(α, β, γ)` or `from_spherical_coordinates(θ, ϕ)`
  from Quaternionic, or any other `Quaternion`, which denotes the
  rotation of its normalization.  A `QuatVec` is rejected, with a
  message suggesting `exp(v/2)`, by every function that takes a
  rotation — including `w(R)` and the `Rθϕ` keyword of `SSHTMatrix`.
  The functions of ``β`` alone — `d`, `dCalculator` and `HCalculator`
  — still accept ``β``, ``e^{iβ}`` or a quaternion.
* **An index can be an `Int` or a half-odd-integer.**  Every public
  function and constructor that takes an index, and every container
  indexed by one, accepts an `Int`, a `HalfOddInteger`, or a
  `Rational{Int}` with denominator 2, where version 2 took any
  `Integer`.  Any other value is refused with an `ArgumentError` that
  says how to write the index: a narrow integer such as an `Int8` or
  an `Int32`, an unsigned integer, a `Bool`, a `BigInt` or an
  `Int128`, a `Rational` of another type or with another denominator,
  a whole number written as a `Rational`, such as `2//1`, and a range
  of indices that is not a `UnitRange` written `a:b`, such as `0:2:4`,
  `0:1:4` or `Base.OneTo(4)`.  The index arithmetic is not closed
  under those types: ``ℓ^2`` overflows an `Int8` at ``ℓ = 12``, and
  ``-m`` wraps around in a `UInt`.  All the indices of one call must
  be of one kind, integer or half-odd.  The exceptions are
  `recurrence!` and `wedge_value`, which are given a calculator or a
  wedge whose kind of index is fixed, and which convert an index of
  any integer type, or a half-odd-integer `Rational` of any integer
  type, to that kind.
* **The element type of a result is that of its input.**  Version 2
  chose the type through arguments such as `D_prep(ℓₘₐₓ, T)`; now
  there is no such argument, and to compute in another type the rotor
  (or angle) must be converted.  A mismatch between the input and a
  preallocated output, as in `sYlm!`, is an error rather than a silent
  conversion.  The transforms, the pixelizations, the quadrature
  weights and the operator matrices take a positional `T` argument,
  since they construct their own numbers.  Sample points passed to a
  transform must be `Rotor`s or `Quaternion`s of its type, and
  colatitudes or weights must be of its type or integers, where
  version 2 converted them.
* **The differential operators are objects rather than functions.**
  `L²`, `Lz`, `L₊`, `L₋`, `R²`, `Rz`, `R₊`, `R₋`, `ð` and `ð̄` are now
  singleton instances of subtypes of `DifferentialOperator`.  Calling
  one as in version 2, `ð(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`, still returns a
  matrix; what is new is described under "Added".
* **`Rz`, `R₊` and `R₋` follow the sign of the spin weight.**  The
  eigenvalue of ``R_z`` on ``{}_sY_{ℓ,m}`` is now ``+s`` rather than
  ``-s``, so `R₊` raises the spin weight (it equals `ð`) and `R₋`
  lowers it, where version 2's `R₊` lowered it and its `R₋` raised it.
  The matrix of `Rz` is therefore the negative of version 2's, and the
  matrices of `R₊` and `R₋` are those of version 2's `R₋` and `R₊`;
  code ported from version 2 should replace `Rz` with `-Rz` and swap
  `R₊` with `R₋`.  The other seven operators return the same matrices
  as before.
* **The transforms have changed defaults and names.**  `SSHT` uses
  the `"RS"` method by default, rather than the dense matrix, and takes
  its element type as a positional argument, `SSHT(s, ℓₘₐₓ, T)`, where
  version 2 took the keyword `T=`; so do `SSHTRS`, `SSHTMatrix` and
  `SSHTMinimal`.  `SSHTDirect` is renamed `SSHTMatrix`, and the method
  name `"Direct"` is accepted with a deprecation warning in favor of
  `"Matrix"`.  `SSHTMatrix` factorizes with `qr` rather than `lu` when
  there are more points than modes, so that the analysis is a
  least-squares solve; version 2 used `lu` in every case.  Analysis of
  a one-dimensional set of function values, `𝒯 \ f`, returns a
  `ModeWeights` rather than a `Vector` — for the in-place `"Minimal"`
  and `"Matrix"` transforms, one that wraps the overwritten input —
  and in-place synthesis, which only `"Minimal"` performs, returns the
  overwritten storage as a plain array.  For ``s ≠ 0`` the `"Minimal"`
  transform samples on different rings (see "Fixed"), so its sample
  points, and the meaning of its `θ` keyword, have changed; for ``s =
  0`` the rings are those of version 2, taken in the order that
  `sorted_rings` gives (see below).  `SSHTMatrix` samples by default
  on the points of `leja_rotors` rather than on the golden-ratio
  spiral, which is far worse conditioned (see "Added"); pass
  `Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, T)` for the old points.
  Each transform runs on the thread of the task that calls it, where
  version 2 threaded the ring-based transform internally; transforms
  run in parallel when each task has its own `copy(𝒯)` (see "Added").
* **Two pixelizations give their points differently.**
  `golden_ratio_spiral_pixels` reduces each azimuth to ``[0, 2π)``,
  computing the reduction with enough extra precision that the result
  is correctly rounded, where version 2 returned ``kΔϕ`` unreduced.
  The rotors of `golden_ratio_spiral_rotors` are built from the
  reduced azimuths, so wherever the reduction subtracts an odd
  multiple of ``2π`` the rotor of a pixel is the negative of the one
  built from ``kΔϕ``, and the same holds for `leja_rotors`, whose
  points are drawn from the spiral.  The two rotors are the same
  rotation, but the sampled values of a function of half-integer spin
  weight differ in sign between them.  The transforms are unaffected,
  because the pixels of a transform and the rotors on which it
  evaluates the harmonics agree.
  `sorted_rings` orders its rings by exact comparisons of their
  positions, so that the order is the same in every floating-point
  type and, for ``s ≠ 0``, the orders for ``s`` and ``-s`` mirror each
  other exactly; for ``s = 0`` the larger ring of each pair placed
  symmetrically about the equator is the northern one, which
  conditions the `"Minimal"` transform better.  `sorted_ring_pixels` and
  `sorted_ring_rotors` follow the same order.
* **`map2salm` has a new signature and output.**  It is called as
  `map2salm(map, s, ℓₘₐₓ)` or `map2salm(map, 𝒯::SSHTRS)`; the
  `show_progress` argument is gone.  The output starts at ``ℓ = |s|``
  rather than ``ℓ = 0``, and is a `ModeWeights` for a single map.
  `map2salm!` and `plan_map2salm` are removed; the plan is now an
  `SSHTRS`, which `SphericalFunctions.map2salm_plan` constructs
  (public, but not exported).  Because `map2salm` is now implemented
  by the ring-based transform, the `BoundsError` that version 2's
  `map2salm!` raised under `--check-bounds=yes` is gone.  (Issue #59.)
* **`Yrange` returns a vector of pairs.**  `Yrange(ℓₘᵢₙ, ℓₘₐₓ)` is now
  a `Vector` of `(ℓ, m)` tuples, so that `Yrange(ℓₘᵢₙ, ℓₘₐₓ)[i]` is
  the pair stored at index `i`, where version 2 returned a matrix with
  one row per mode and the columns ``ℓ`` and ``m``.  Code that read
  `Yrange(…)[i, 1]` and `Yrange(…)[i, 2]` should read
  `Yrange(…)[i][1]` and `Yrange(…)[i][2]`; `stack(Yrange(…); dims=1)`
  gives the old matrix.
* The dependencies `AbstractFFTs`, `FastTransforms`, `Hwloc`,
  `LoopVectorization`, `OffsetArrays`, `ProgressMeter` and `Random`
  are dropped, and `FixedSizeArrays`, `GenericFFT`, `PrecompileTools`
  and `Serialization` are added.  `DoubleFloats` becomes a weak
  dependency: loading it adds exact conversions of a `HalfOddInteger`
  to its types.

### What replaces what

| Version 2 | Version 3 |
|---|---|
| `D_matrices(R, ℓₘₐₓ)` + `D_iterator` | `D(R, ℓₘₐₓ)`, indexed `𝔇[ℓ][m′, m]` (conjugated; see above) |
| `D_matrices(α, β, γ, ℓₘₐₓ)` | `D(α, β, γ, ℓₘₐₓ)` |
| `D_prep` + `D_matrices!` | `DCalculator(R, ℓₘₐₓ)`, iterated or stepped with `recurrence!`, and reused with `set_R!` |
| `d_matrices(β, ℓₘₐₓ)` + `d_iterator` | `d(β, ℓₘₐₓ)`, indexed `𝔡[ℓ][m′, m]` |
| `d_prep` + `d_matrices!` | `dCalculator(β, ℓₘₐₓ)`, reused with `set_β!` |
| `H!`, `H_recursion_coefficients` | `HCalculator(β, ℓₘₐₓ)` |
| `sYlm_values(R, ℓₘₐₓ, s)` + `sYlm_iterator` | `sYlm(R, ℓₘₐₓ, s)`, indexed `Y[ℓ][m]` |
| `sYlm_values(θ, ϕ, ℓₘₐₓ, s)` | `sYlm(θ, ϕ, ℓₘₐₓ, s)` |
| `sYlm_prep` + `sYlm_values!` | `sYlmCalculator(R, ℓₘₐₓ, s)` + `sYlm!(Y, calc, R; ℓₘᵢₙ=0)`, or for several spin weights `sYlmCalculator(R, ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)` + `sYlm!(Y, calc, R, s; ℓₘᵢₙ=0)` (see below) |
| `ₛ𝐘(s, ℓₘₐₓ, T, R⃗)` | `sYlm_matrix(R⃗, ℓₘₐₓ, s)` |
| `λ_iterator` | `SphericalFunctions.sλlm` or `SphericalFunctions.sλlmCalculator` |
| `ALFcompute` and relatives | no direct replacement; `sλlm(θ, ℓₘₐₓ, 0)` gives ``Y_{ℓ,m}(θ, 0)`` |
| `WignerDindex`, `WignerHsize`, … | not needed: blocks are indexed by ``(ℓ, m′, m)`` directly |
| `SSHT(s, ℓₘₐₓ; T=T, method)` | `SSHT(s, ℓₘₐₓ, T; method)` |
| `SSHTDirect` | `SSHTMatrix` |
| `map2salm(map, s, ℓₘₐₓ, show_progress)` | `map2salm(map, s, ℓₘₐₓ)` |
| `plan_map2salm` + `map2salm!` | `𝒯 = SphericalFunctions.map2salm_plan(map, s, ℓₘₐₓ)` + `map2salm(map, 𝒯)` |

The `ℓₘᵢₙ=0` in the `sYlm_prep` row matters.  Version 2's `sYlm_prep`
allocated storage starting at ``ℓ = 0``, while `sYlm!` starts at
``ℓₘᵢₙ = |s|`` unless told otherwise, and accepts a longer vector,
writing only its first `Ysize(ℓₘᵢₙ, ℓₘₐₓ)` elements.  A version-2
buffer, indexed with `Yindex(ℓ, m)` as before, would therefore be read
at the wrong modes without any error.  (The other `sYlm` functions of
version 2 already started at ``ℓ = |s|``.)

### Added

* **A complete account of the conventions.**  The documentation now
  derives every convention used in the package from first principles,
  starting from Cartesian coordinates and proceeding through
  quaternions, rotations, the angular-momentum operators, Wigner's
  matrices, and the spin-weighted spherical harmonics, for integer and
  half-integer indices alike.  Those conventions are then compared
  with thirty other sources, each on its own page.  They range from
  the founders — Euler (1767), Hamilton, Tait, Clifford, Gibbs and
  Wilson — through the standard texts and papers of quantum mechanics
  and relativity, including Whittaker, Condon and Shortley, Wigner,
  Edmonds, Newman and Penrose, Goldberg et al., Thorne, Varshalovich
  et al., Sakurai, Shankar, and Zettili, to current references and
  software, including the NIST DLMF, Wikipedia, LALSuite, Mathematica,
  SciPy and SymPy.  For the twenty-four sources whose formulas can be
  evaluated, those formulas are transcribed and checked numerically
  against this package.  The checks are test items that run with the
  test suite, as well as being rendered in the documentation, so a
  change that breaks agreement with any of these sources is caught.
* **Half-integer indices.**  `D`, `d`, the calculators, `sYlm`,
  `sYlm!`, `sYlm_matrix`, `Ysize`, `Yindex`, `Yrange`, `ModeWeights`,
  the differential operators, the golden-ratio, Leja and sorted-ring
  pixelizations, and the `"RS"` and `"Matrix"` transforms (with
  `map2salm` and `salm2map`) all accept half-integer ``ℓ``, ``m`` and
  ``s``, passed as `Rational{Int}`s with denominator 2 — as in `D(R,
  7//2)` or `SSHT(1//2, 7//2)` — or as `HalfOddInteger`s, the type to
  which such a `Rational` is converted.  The results are indexed
  exactly as in the integer case.  They have been verified against two
  independent references to ``10^{-16}`` for ``ℓ ≤ 31/2``, and by
  identities that need no reference to ``ℓ = 101/2``.  `Ylm` and
  `YlmCalculator`, the two equiangular grids (`driscoll_healy_pixels`
  and `mcewen_wiaux_pixels`, with their `_rotors` counterparts),
  `minimal_rings` and the `"Minimal"` transform remain integer-only,
  and say so; a call that mixes integer and half-integer indices is
  refused.  (Issue #29.)
* **The rules for indices are public.**  `HalfOddInteger` and
  `IntegerHalf` are public, as are the markers `IndexType`,
  `IndexRange` and `IndexOrRange` and the macro `@index_methods`, with
  which the public functions that take indices are defined.  A package
  building on this one can use them to define functions that accept
  exactly the indices these do, and refuse the others with the same
  messages.
* **`ModeWeights`**, which holds the mode weights of a spin-weighted
  function in the canonical ordering, together with its spin weight
  and range of ``ℓ``.  It is indexed as `w[ℓ, m]` or `w[ℓ, :]`, and
  `modes` and `spin` report its labels.  It supports
  - evaluation at a rotor or a vector of rotors, `w(R)`, which is the
    same as the product `Y * w` with harmonics from `sYlm` or an
    `sYlmCalculator` (this is written `*` rather than `⋅`, since `dot`
    conjugates its first argument, and `dot` on these types raises an
    error saying so);
  - active rotation, `D(R, ℓₘₐₓ) * w` (or with a `DCalculator` in
    place of `D`), which gives the weights of ``f′(𝐐) =
    f(𝐑^{-1}𝐐)``; and
  - the differential operators, written `ð(w)` or `ð * w`, which
    return a `ModeWeights` with the spin weight adjusted
    appropriately and the same range of ``ℓ`` as `w`.  These are
    applied by a loop rather than by building a matrix, so the only
    allocation is the result, and `mul!(w′, ð, w)` allocates nothing.
    A product of several operators, as in `ð̄ * ð * w`, is applied
    from right to left, and in a broadcast, as in `ð.(ws)`, an
    operator is a single value.

  Rotation also has an in-place form, `mul!(w′, 𝔇, w)`.  Evaluation
  reads only the weights with ``ℓ ≥ |s|``, which are those of a
  function of spin weight ``s``.  `ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` copies
  weights into another range of ``ℓ``, filling the modes it adds with
  zeros.  The transforms accept a `ModeWeights` directly and check its
  labels: synthesis takes weights of the transform's spin weight over
  any range of ``ℓ`` up to its ``ℓₘₐₓ``, while the outputs of analysis
  and the destinations of `mul!` must have exactly the transform's
  labels.  The linear algebra that keeps the labels is defined on a
  `ModeWeights` — `zero`, `fill!`, `rmul!`, `lmul!`, `copyto!`,
  `axpy!`, `axpby!` and `dot` — and `==` and `isequal` between two
  `ModeWeights` compare the labels as well as the numbers.  A plain
  matrix has no labels to check, so its product with a `ModeWeights`
  is refused with a message naming the labelled alternatives; the
  matrix forms of the operators, and `sYlm_matrix`, act on the plain
  vector `array_view(w)`.
* **Calculators** — `DCalculator`, `dCalculator`, `HCalculator`,
  `sYlmCalculator` and `YlmCalculator` — take the same positional
  arguments as the corresponding functions, and also accept a vector
  of rotors; a harmonic calculator always starts at ``ℓ = 0``, and
  takes no `ℓₘᵢₙ` keyword.  They compute one ``ℓ`` at a time, so that
  large ``ℓₘₐₓ`` can be reached without holding every matrix at once.
  Iterating one, as in `for (ℓ, 𝔇ˡ) ∈ calc`, yields each block as a
  view that the next step overwrites, and allocates nothing;
  `recurrence!(calc, ℓ)` computes a single step, and `collect` copies
  every block.  `set_R!`, `set_β!` and `set_θ!` point an existing
  calculator at new data.  A calculator is a mutable workspace, so
  each task that computes in parallel needs its own, which
  `similar(calc)` or `similar(calc, R)` provides.  (Issue #32.)
* **Batched evaluation.**  Given a vector of rotors, the calculators
  and `sYlm` evaluate all of them at once, with the rotor as the
  leading index of each block.  This is two to ten times faster per
  rotor than a loop, because the rotor index is the only one that the
  recursions allow to be vectorized.  The transforms use it
  internally.  What makes a batch is that the rotors come as a vector,
  so a vector of one rotor is a batch of one, and `isbatched` reports
  which kind a calculator is; the operations that need a single rotor,
  such as rotating a `ModeWeights`, refuse a batch.
* **Spinor phases.**  The phases of a rotor that the recursions need
  are computed directly from its components as ``e^{i(α±γ)/2}`` and
  ``\cos(β/2)``, ``\sin(β/2)``, rather than through Euler angles.
  This is accurate near both poles, and gives ``𝔇(-R) = -𝔇(R)`` for
  half-integer indices.  (Issue #57.)
* **Automatic differentiation with respect to the rotor.**  Package
  extensions for `ChainRulesCore`, `EnzymeCore`, `ForwardDiff`,
  `Mooncake`, and `ReverseDiff` supply rules that give the derivatives
  of ``𝔇`` and of the harmonics from the angular-momentum operators,
  as combinations of the values themselves, rather than by
  differentiating the recurrence, so the derivatives are as accurate
  as the values at every rotor, the poles included.  `ForwardDiff`,
  `Enzyme`, `Mooncake`, and `ReverseDiff` apply them to every step of
  a `DCalculator` or `sYlmCalculator`, batched or not, so that a loop
  over a calculator's blocks is differentiated one block at a time, as
  efficiently as it is evaluated; a calculator of dual numbers
  allocates nothing when stepped, and gives derivatives of every order
  exactly under nested differentiation.  `D`, `sYlm`, `Ylm`, and
  `sYlm_matrix` are computed by calculators, so the same rules serve
  them; `Zygote` uses rules for their arrays instead, and
  `ReverseDiff` uses those too.  The note on automatic differentiation
  in the documentation describes all of this.
* **Restricted ranges.**  The keywords `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ` and
  `mₘᵢₙ` of `D`, `d` and their calculators limit the part of each
  matrix that is computed, and restricting either index makes the
  recurrence itself cheaper.  Each lower limit defaults to minus the
  corresponding upper one, and each range must contain 0 (or both
  ``±1/2``), where the recurrence starts.  `sYlm` and its relatives
  accept an `ℓₘᵢₙ` keyword and a unit range of spin weights such as
  `-2:2`.
* **`array_view` and `relabel`.**  `array_view(x)` gives the contents
  of a container as a 1-based `StridedArray` that aliases its storage,
  so that BLAS and LAPACK can operate on it without copying; for mode
  weights and harmonics this is the flat array in the canonical
  ordering.  `relabel(x, A)` puts the natural indices back on a plain
  array.  `Matrix`, `Array` and `collect` give copies.  An array with
  offset axes, such as an `OffsetArray`, is refused by `array_view`,
  and so by the transforms, because the result is indexed from 1.
* `Ylm`, the ordinary scalar spherical harmonics, and `sYlm_matrix`,
  the dense synthesis matrix.
* `HWedge`, the part of ``H^ℓ`` that an `HCalculator` computes and
  stores, and `wedge_value`, which reads any element of ``H^ℓ`` from
  it through the symmetries, applying the sign they introduce for
  half-integer indices.
* The real harmonics ``{}_sλ_{ℓ,m}(θ)`` — ``{}_sY_{ℓ,m}(θ, 0)`` for
  integer ``s``, and ``{}_sY_{ℓ,m}(θ, 0) / i^{2s}`` for half-odd ``s``
  — through the public but unexported `sλlm`, `sλlm!`, `sλlm_matrix`
  and `sλlmCalculator`.
* The angular-momentum operators `Lx` and `Ly`.  (There is
  deliberately no `Rx` or `Ry`; see the `Lx` docstring.)
* `salm2map`, the inverse of `map2salm`, which also takes a
  `ModeWeights` directly, as `salm2map(w, Nϕ, Nθ)`.
* The pixelizations `driscoll_healy_pixels`, `driscoll_healy_rotors`,
  `mcewen_wiaux_pixels` and `mcewen_wiaux_rotors`, which are public
  but unexported.
* The pixelizations `leja_pixels` and `leja_rotors`, which choose
  exactly as many points as there are modes, as discrete Leja points
  drawn from a golden-ratio spiral, so that the matrix of harmonics on
  them is well conditioned.  The points are chosen in `Float64` and
  then evaluated in the requested type, so that every precision gets
  the same points.  They are the default sample points of
  `SSHTMatrix`; on the golden-ratio spiral, which version 2 used, a
  round trip loses all its digits by ``ℓₘₐₓ = 64``, and on these
  points it loses about 3.
* `ComplexPowers`, an iterator over the powers of a unit complex
  number, which runs the same recurrence as `complex_powers!`.
* `sqrtbinomial`, the square root of a binomial coefficient, computed
  through the logarithm of the beta function so that it stays finite
  and accurate where the coefficient itself overflows; it is public
  but not exported.
* **Transforms from several tasks.**  `copy(𝒯)` gives a transform
  that shares the read-only tables and FFT plans of `𝒯` with new
  workspace, for a small fraction of the cost of constructing one, and
  is the way to give each task its own.  An `SSHTMatrix` holds no
  workspace, so it may be shared, and its `copy` is itself.  A
  `deepcopy` is independent, and a transform may be serialized — to
  send it to another process, for example — and makes its FFT plans
  again where it arrives.
* The accessors `ℓₘᵢₙ`, `ℓₘₐₓ`, `spins`, `Nᵣ`, `floattype` and others
  are declared public but not exported, with the ASCII aliases `ell`,
  `ell_min`, `ell_max`, `mp_max`, `mp_min`, `m_max`, `m_min`, `s_max`,
  `s_min` and `Nr`.  The keyword arguments of the same names accept
  the same ASCII spellings, as in `D(R, ℓₘₐₓ; mp_max=2)` or `sYlm(R,
  ℓₘₐₓ, s; ell_min=0)`.
* ASCII aliases, public but not exported, for the other functions
  whose names are not ASCII: `L2`, `Lplus`, `Lminus`, `R2`, `Rplus`,
  `Rminus`, `eth` and `ethbar` for the operators, `Deltaspin` for
  `Δspin`, `set_beta!` and `set_theta!` for the setters, and
  `slambdalm`, `slambdalm!`, `slambdalm_matrix` and
  `slambdalmCalculator` for the real harmonics.

### Fixed

* The size functions throw for an inverted range of ``ℓ``, rather than
  returning a negative number.  (Issue #52.)
* `complex_powers!` works for wrapper element types such as
  `ForwardDiff.Dual`, including at the phase 1, where it used to
  produce `NaN` derivatives.
* The ring-based transform handles rings with different numbers of
  points, and rings with more than ``2ℓₘₐₓ+1`` points; version 2 gave
  wrong results in both cases.
* The `"Minimal"` transform is accurate for ``s ≠ 0``.  Its rings of
  ``2j+1`` points, which work well for ``s = 0``, are badly
  conditioned for any other spin weight: at ``s = 2`` a round trip
  lost half its digits by ``ℓₘₐₓ = 10`` and all of them by ``ℓₘₐₓ =
  14``, silently.  The rings are now arranged as `minimal_rings`
  describes, and the analysis solves small groups of ``m`` values
  together; at ``s = 2`` and ``ℓₘₐₓ = 16`` the error of a round trip
  falls from ``10^2`` to ``10^{-13}``.  The sample points still become
  badly conditioned at larger ``ℓₘₐₓ``, for every spin weight, and the
  constructor now warns when a round trip would lose more than half
  the digits.
