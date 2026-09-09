# SphericalFunctions.jl v3 — design memo

**Status:** authoritative design record for the v3 rewrite, written 2026-09-08.
Supersedes the earlier blueprint that lived in this file, `notes/docs/interface_design.md`,
`notes/redesign/*.md`, and the interface items in `TODO.md`.  Where any of those disagree
with this file, this file wins.  Companion decisions on *conventions* live in
`docs/src/conventions/{summary,details}.md`; comparisons against the literature live in
`docs/literate_input/conventions/comparisons/`.

Contents

1. Decisions in one screen
2. Where the code actually is
3. Bugs found while reviewing (not yet fixed)
4. Interface specification
5. Half-integer ℓ: the route, the math, the verification
6. Half-integer test layer
7. What is deprecated and what comes back
8. Reconciliation of the old notes
9. Performance constraints and benchmark results
10. Coordination with the conventions and comparisons work
11. Implementation checklist (ordered)
12. Provenance and dating of the notes


## 1. Decisions in one screen

| # | Decision |
|---|---|
| D1 | One engine: the rotor-batched, β-only `WignerHCalculator`.  Sizing and type at construction, rotor data at call time.  Convenience functions take the rotor first: `D(R, ℓₘₐₓ; …)`. |
| D2 | One stepping verb, `recurrence!(calc, ℓ)`.  Sequential ℓ is O(Nᵣℓ²); jumping restarts the axis (make it advance instead when possible).  `next!`/`compute!` from the old blueprint are dropped. |
| D3 | Internal storage is the rotor-first real `HWedge`.  Output is a per-ℓ view indexed by natural `[m′, m]` (half-integer allowed) that does **not** pretend to be an `AbstractMatrix` when its axes are not integers.  `Matrix(view)` materializes. |
| D4 | Half-integer ℓ is supported, with **no new ℓ-recursion**: reuse the integer axis, seed the rows m′ = ±1/2 by a Clebsch–Gordan step, then run the unchanged m′ ladder.  Verified numerically (§5). |
| D5 | sYlm is the same machinery with m′ limited to ±\|s\|, reading the column m = −s through one symmetry helper. |
| D6 | Mode weights: one canonical ordering (ℓ-major, m increasing) behind a small indexable type.  The three-ordering hierarchy is dropped. |
| D7 | Calculators *are* the workspace; add `similar(calc)` so pools of them fill a `Channel`.  No separate `workspace` function. |
| D8 | SSHT is rebuilt on top of this layer later; its requirements are recorded in §4.8. |
| D9 | Issue #41 is designed out by natural indexing.  Issue #42 is closed by the convention flip: v3 uses 𝔇_{m′m}(α,β,γ) = e^{−im′α} d_{m′m}(β) e^{−imγ}, which the new engine already implements. |
| D10 | `Deprecated` is a quarantine, not a verdict; §7 says what comes back. |


## 2. Where the code actually is

- `src/wigner/` holds the new engine: `AbstractWignerMatrix{IT,NT,ST}`, `WignerMatrix`
  (aliases `WignerDMatrix`, `WignerdMatrix`), `WignerRange`, `HWedge` (flat real storage
  `[iᵣ, m′, m]`, precomputed row offsets, mutable `ℓ`, not thread-safe), `HAxis`,
  `WignerCalculator` (single rotor, dense per ℓ; complete; agrees with `Deprecated.d_matrices`
  to 1e-16), and `WignerHCalculator` (Nᵣ rotors, β only; **incomplete** — step 6 and the
  phase conversion are commented out, so it produces only the H wedge, never D, d or sYlm).
- Half-integer indices: the containers, `WignerRange`, `validate_index_ranges` (requires
  denominator 2), `ℓₘᵢₙ(::Type{<:Rational}) = 1//2`, `isrational`, and printing are done and
  tested.  Every recursion function is still gated `where {IT<:Signed}`.
- Nothing outside `Deprecated` produces sYlm, indexes mode weights, or performs an SSHT.
  `src/ssht/` is an unreachable scaffold (its `include` paths resolve to files that do not
  exist).  The `mode_weights/` directory from the old blueprint was never created.
- Tests: hand-rolled `test/runtests.jl` plus TestItemRunner; there is no
  `JuliaTestItems.toml`; both `[targets] test` and `test/Project.toml` exist.  **No test
  exercises either calculator.**  Two `@testitem`s live inside `src/` (`ComplexPowers`,
  `WignerMatrix`).
- Earlier attempts at this redesign exist as branches: `simpler_API`, `matrix_order`,
  `requested_types`, `out_of_place`, `iterate_completely`, `threading`, `ssht`.  Diff them
  before re-solving a problem they may already have solved.


## 3. Bugs found while reviewing (not yet fixed)

| # | Where | What | Consequence |
|---|---|---|---|
| B1 | `src/wigner/wigner_H_calculator.jl:149-241`, `recurrence_step2!` | Uses the *source* ℓ where the *target* ℓ+1 is needed (coefficient `b̄ₗ` at ~177, loop bound at ~185, the peeled `m = ℓ−1` and `m = ℓ` blocks at ~206 and ~225). | Every batched axis is wrong for ℓ ≥ 2: from the ℓ=1 axis it yields H²₀₀ = cos²β instead of (3cos²β−1)/2 and never writes H²₀₂.  Max error vs `Deprecated.H!` ≈ 2.4 at ℓₘₐₓ = 7.  The single-rotor path in `recurrence.jl:37-83` is correct. |
| B2 | `wigner_H_calculator.jl:365`; `recurrence.jl:162, 210, 213, 232` | Steps 5, 6 and `convert_H_to_d!` assume `m′ₘᵢₙ = −m′ₘₐₓ`. | `WignerHCalculator(eⁱᵝ, 4, 4, -1)` then `recurrence!(w, 4)` throws `BoundsError` (caught only because `FixedSizeVector` bounds-checks under `@inbounds`). |
| B3 | `recurrence.jl:256` vs `src/deprecated/evaluate.jl:355-372` and `src/utilities/operators.jl` | New `convert_H_to_D!` applies e^{−i(m′α+mγ)}; deprecated `D!` applies e^{+i(m′α+mγ)}; `operators.jl` has `Rz → −s`, `R₊ = R_x − iR_y`. | Not a bug in the new code: the settled convention (2026-09-08) is the new one.  `Deprecated` and `operators.jl` carry the old convention and must be updated or retired.  Nothing tests this either way. |
| B4 | `src/wigner/wigner_matrix.jl:49, 81, 127` | `AbstractWignerMatrix{<:Rational} <: AbstractMatrix` with non-integer `axes`. | `Matrix(w)`, `copy(w)`, `w == w`, `sum(w)`, `w .+ 1`, `w[:,1]`, `collect(w)` all throw for `WignerDMatrix(rand(ComplexF64,2,2), 1//2)`; `show` is already special-cased to hide this. |
| B5 | `wigner_H_calculator.jl:243-282` (step 3 gate), `wigner_matrix.jl:156-256` (`validate_index_ranges`) | Step 3's gate `ℓ > 0 && m′ₘₐₓ ≥ 1` would fire for Rational ℓ; validation accepts `m′ₘᵢₙ = +1/2`. | Both block the half-integer route (§5 needs step 3 skipped and `m′ₘᵢₙ ≤ −1/2`). |
| B6 | `src/wigner/recurrence.jl:7` | `ϵ(m)` uses `isodd(m)`, which errors on `Rational`. | Blocks half-integer phase conversion. |
| B7 | `src/ssht/ssht.jl` | `include("ssht/direct.jl")` resolves to `src/ssht/ssht/direct.jl`. | Dead scaffold; `salm2map` is also defined twice in `huffenberger_wandelt.jl`. |


## 4. Interface specification

### 4.1 Construction and calls (D1)

```julia
calc = WignerDCalculator(ℓₘₐₓ, T; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ, Nᵣ=1)
calc = WignerdCalculator(ℓₘₐₓ, T; …)          # real output, β only
calc = sYlmCalculator(s, ℓₘₐₓ, T; Nᵣ=1)        # m′ limits fixed to ±|s| internally

recurrence!(calc, R, ℓ)      # R::Rotor or AbstractVector{<:Rotor} of length Nᵣ; ℓ may be any
                             # value in ℓₘᵢₙ:ℓₘₐₓ (sequential is cheap, a jump costs a restart)
calc[ℓ]                      # view for the current rotors: natural [m′, m] or [iᵣ, m′, m]
Matrix(calc[ℓ])              # materialize one rotor's ℓ block

D(R, ℓₘₐₓ; kwargs...)        # convenience: allocate, run every ℓ, return a container indexed by ℓ
d(β, ℓₘₐₓ; kwargs...)
sYlm(R, ℓₘₐₓ, s; kwargs...)
```

Rules that follow from the measured constraints (§9):

- Keyword arguments only at construction; none in any function called per ℓ, per m′ or per
  rotor.  No function-valued index arguments anywhere.
- `ℓₘₐₓ` and the four limits are all `IT` (an `Integer`, or a `Rational` with denominator 2);
  they fix `ℓₘᵢₙ ∈ {0, 1/2}`.  Mixed integer/half-integer arguments are an error.
- `WignerHCalculator` (rotor data and `Nᵣ` fixed at construction) becomes the engine.
  `recurrence!(calc, R, ℓ)` copies the rotor phases into the calculator's buffers and
  invalidates the axis state.  `WignerCalculator` either becomes a thin wrapper with `Nᵣ = 1`
  or is deleted; it duplicates all six recursion steps today.
- The rotor supplies three things per element: `eⁱᵝ` (for the axis recursion), the
  half-angle pair `(cos β/2, sin β/2)` computed as `(√(W²+Z²), √(X²+Y²))` (needed only for
  half-integer ℓ), and the two half-angle phases `zp`, `zm` of §5.3 (needed for D and sYlm).
  Add an accessor to Quaternionic.jl (`to_euler_phases` already computes all of these as
  locals) or compute them locally.

### 4.2 Stepping (D2)

`recurrence!(calc, ℓ)` (after rotors are set) computes the wedge for exactly one ℓ.  The
axis buffers hold ℓ and ℓ+1; advancing by one is O(Nᵣℓ²).  Today `_recurrence!`
(`wigner_H_calculator.jl:472-532`) restarts from `ℓₘᵢₙ` whenever the requested ℓ is not the
successor, even when the buffers already hold ℓ−1 or ℓ.  Change it to advance from the current
axis state whenever `ℓ(h⃗ˡ⁺¹) ≤ ℓ`, and restart only when going backwards or when the state
is invalid.  `Base.setproperty!(w, :ℓ, …)` relabels storage without recomputing and stays
internal.

### 4.3 Storage and views (D3)

- Internal: `HWedge`, real, rotor index innermost, wedge m ≥ |m′|.  This is what makes the
  inner loops vectorize (`@simd` over rotors) and what the vectorization notes argued for.
  It forecloses `Vector{Matrix{Complex}}` as *storage*; a nested container may still be the
  return type of the convenience function `D(R, ℓₘₐₓ)` for one rotor, built from views.
- Output view: natural indexing `view[m′, m]` (one rotor) or `view[iᵣ, m′, m]`.  Fix B4 by
  not subtyping `AbstractMatrix` when `IT <: Rational`; for integer `IT` an `OffsetArray`
  view is the idiomatic choice.  Provide `getindex`, `axes`-like accessors returning
  `WignerRange`, iteration, and `Matrix(view)`.
- Materialization: `fill_d!(Wˡ, Hˡ, iᵣ)` and `fill_D!(Wˡ, Hˡ, iᵣ, zp, zm)` apply ε, σ and
  phases while reading wedge entries through **one** helper `wedge_index_and_sign(H, m′, m)`
  that encodes the case analysis of `Deprecated.WignerHindex` (`indexing.jl:327-345`) plus the
  half-integer sign σ of §5.2.  No other code path may re-derive the symmetries.
- The old blueprint's lazy phase-on-the-fly wedges (`?WignerDWedge`, `?ₛYₗₘWedge`) are an
  optimization for SSHT inner loops, not the primary API.

### 4.4 sYlm (D5)

ₛYₗₘ(R) = (−1)^s √((2ℓ+1)/4π) conj(𝔇^{(ℓ)}_{m,−s}(R)) in the settled convention (see the
conventions pages for the equivalent form and the half-integer branch choice).  The first
index of 𝔇 is m, so the needed entries form the *column* m = −s of the wedge, reached from
the stored rows ±|s| by the transpose symmetry — which carries σ for half-integer s.
`sYlmCalculator(s, ℓₘₐₓ, T; Nᵣ)` fixes `m′ₘₐₓ = |s|`, `m′ₘᵢₙ = −|s|`, and applies the
prefactor and phases in its own `fill_Y!`.  Cost per rotor per ℓ is O(|s|·ℓ), i.e. O(ℓ²) for
all ℓ, as in v2's `Y!`.

### 4.5 Mode weights (D6)

One ordering: `[ f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ ]`, with `ℓₘᵢₙ` the smallest
representable ℓ for the index type (0 or 1/2) and a separate `ℓ₀ = |s|`-style lower cut
where an operator needs it.  A small type `ModeWeights{IT, T}` wraps a vector and offers
`w[ℓ, m]`, `size`, iteration over `(ℓ, m)`, and the `Yindex` closed form internally
(`Deprecated.indexing.jl` has the SymPy-derived formulas and their tests; reuse, do not
re-derive).  Operators (`L², Lz, L₊, L₋, R², Rz, R₊, R₋, ð, ð̄`, plus the planned
`Lx, Ly, Rx, Ry`) act on it.  Other orderings can be added later as new types; none exist now.

### 4.6 Parallelism (D7)

`HWedge` mutates its row-offset table when `ℓ` changes, so one wedge cannot be shared between
tasks.  Define `Base.similar(::WignerHCalculator)` reconstructing from the stored
`(RT, Nᵣ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)` (and deduplicate those fields against the wedge's `maxℓ`,
`maxm′ₘₐₓ`, `minm′ₘᵢₙ`).  `map2salm`'s existing `Channel` pool pattern
(`Deprecated/map2salm.jl`) then holds calculators.  There is no separate `workspace`
function; the TODO item proposing one is closed by this decision.

### 4.7 Convenience functions

`D(R, ℓₘₐₓ; …)`, `d(β, ℓₘₐₓ; …)`, `sYlm(R, ℓₘₐₓ, s; …)` allocate a calculator, run every ℓ,
and return a container indexed by ℓ whose elements are the views of §4.3.  They are documented
as the slow path for repeated use.  Argument order is rotor first (matches v2's
`D_matrices(R, ℓₘₐₓ)`; the `D(ℓₘₐₓ, R)` form in older notes is withdrawn).

### 4.8 What the SSHT layer will need (D8)

Recorded so the Wigner layer is not redesigned twice: a β-only H stage per ring with `Nᵣ` =
number of rings; batched rotors; `(zp, zm)` phases applied per pixel; the sYlm wedge with
m′ = ±|s|; rotor-first layout for the inner loops; `similar` for pools; and the
`FFTW`-style `*`/`\` semantics of the deprecated `SSHT` types, which are kept (§7).


## 5. Half-integer ℓ: the route, the math, the verification

### 5.1 Why no new ℓ-recursion is needed

Gumerov–Duraiswami's steps 4 and 5 (the m′ ladder at fixed ℓ) follow from equating ∂_β d^j
computed with J_y acting on the left index and on the right index:

    δ^{m′} d_{m′+1,m} = δ^{m′−1} d_{m′−1,m} + δ^{m−1} d_{m′,m−1} − δ^{m} d_{m′,m+1},
    δ^k ≡ √((j−k)(j+k+1)).

No parity assumption enters, so the identity and its coefficient magnitudes hold for
half-integer j.  Step 6 (symmetries) likewise follows from d_{m′m} = (−1)^{m′−m} d_{mm′} =
d_{−m,−m′}, valid for all j.  Only the *seed* is ℓ-dependent: for integer ℓ it is the ALF-type
axis recursion (step 2) plus step 3; for half-integer ℓ it is §5.3.  The stalled derivation in
`notes/redesign/redesign_half_integer.md`, which tried to build a self-contained half-integer
ℓ-recursion, is abandoned.  (Your `docs_outline.ipynb` also shows no half-step ladder
operator exists as a combination of L_x, L_y, so an operator route is closed.)

### 5.2 Sign bookkeeping that changes for half-integers

Write d_{m′m} = ε_{m′} ε_{−m} H_{m′m}.  For the H-form of the ladder to keep G&D's shape one
needs ε_{k+1}/ε_k = −sgn(k) (with sgn(0) = +1), and the unique extension with ε_k = 1 for
k ≤ 0 is

    ε(m) = (−1)^⌊m⌋  for m > 0,   ε(m) = 1  for m ≤ 0
    (Julia: ϵ(m) = ifelse(m > 0 && isodd(floor(Int, m)), -1, 1); reduces to (−1)^m on integers)

giving ε_{1/2} = +1, ε_{3/2} = −1, …  With it the recurrence reads

    sgn(m′) δ^{m′} H_{m′+1,m} = sgn(m′−1) δ^{m′−1} H_{m′−1,m} − sgn(m) δ^{m−1} H_{m′,m−1} + sgn(m+1) δ^{m} H_{m′,m+1}.

For integers this is G&D's Eq. (50).  Two things differ for half-integers:

1. **Step 4 is not sign-free.**  At m′ = 1/2 the coefficient of H_{m′−1,m} carries
   sgn(m′−1) = sgn(−1/2) = −1.  The code comment claiming "signs of m′ and m are always +1"
   is false there.  (The m-signs, sgn(m) and sgn(m+1), are +1 throughout the wedge ranges
   actually visited, m ≥ 3/2, so the code's `sgn(m−1)`, `sgn(m)` on those coefficients are
   harmless; only the m′ sign matters.)
2. **The transpose symmetry acquires a sign.**  For any real ε, H_{−m,−m′} = H_{m′,m}
   exactly, but H_{m,m′} = σ H_{m′,m} and H_{−m′,−m} = σ H_{m′,m} with

       σ(m′, m) = sgn(m)·sgn(m′)   for half-integers,   σ ≡ 1 for integers.

   This cannot be removed by any choice of ε: on the anti-diagonal m = −m′ the ε products
   are squares, and d^{1/2}_{1/2,−1/2} = −sin(β/2) = −d^{1/2}_{−1/2,1/2}.  Step 6 and every
   on-the-fly symmetry read (§4.3's helper, the sYlm column read) must carry σ.  The
   statement "H is symmetric" in `docs/src/notes/H_recurrence.md` needs the caveat.

Loop ranges become: step 4 starts at m′ = 1 − ℓₘᵢₙ (1 or 1/2); step 5 starts at
m′ = −ℓₘᵢₙ (0 or −1/2); step 6 runs m ∈ (1 − ℓₘᵢₙ):ℓ.  Step 3 is skipped for Rational ℓ.

### 5.3 The seed: two rows from the integer axis

Let J be half-integer, j = J − 1/2, and h_k = H^{j}_{0,k} = d^{j}_{0,k}(β) for k = 0…j the
integer axis that step 2 already produces (h_{j+1} ≡ 0).  With c = cos(β/2), s = sin(β/2),
Varshalovich Eqs. 4.8.2(14) and (15) — in the phase-consistent form already tested in
`test/conventions/varshalovich.jl` (the transcriptions in the notes carry a sign error) —
applied at M′ = ∓1/2 and transposed into this package's (m′, m) order give, for m ∈ 1/2 : J,

    H^J_{+1/2, m} = [ √(J+m) · c · h_{m−1/2}  −  √(J−m) · s · h_{m+1/2} ] / √(J + 1/2)
    H^J_{−1/2, m} = [ √(J+m) · s · h_{m−1/2}  +  √(J−m) · c · h_{m+1/2} ] / √(J + 1/2)

(and on these two rows H = d, since ε_{±1/2} ε_{−m} = 1 for m > 0).  The denominator is
j + 1 ≥ 1, so nothing is singular; the full wedge range m = 1/2…J is covered from
h_0…h_j alone; cost is O(Nᵣ J), cheaper than integer step 3 (which needs the ℓ+1 axis).
Both rows are required: the corner H_{−1/2,1/2} cannot be reached by step 5 from the +1/2 row
without leaving the wedge.  Hence `validate_index_ranges` must require m′ₘᵢₙ ≤ −1/2 for
Rational index types (B5).  Then step 4 climbs from m′ = 1/2 using the −1/2 row as
"m′ − 1", and step 5 descends from m′ = −1/2 using the +1/2 row as "m′ + 1".

Consequences for the calculator: the axis buffers `h⃗ᵃ, h⃗ᵇ` stay **integer-indexed**
(`HAxis{Int,RT}` regardless of `IT`; either add a type parameter or make `HAxis`
integer-only), `consistent_ℓ` becomes `ℓ(h⃗ˡ) == ℓ(Hˡ) − ℓₘᵢₙ(IT)`, and the driver's
`ℓ == ℓₘᵢₙ` branch for Rational is "step 1 on the j = 0 axis, then the seed".  The
calculator also needs `(c, s)` per rotor, best taken from the quaternion as
(√(W²+Z²), √(X²+Y²)) rather than from cos β (which loses relative accuracy near β = 0).

### 5.4 Phases without square roots

For half-integer m′, m the factors e^{−im′α}, e^{−imγ} are not integer powers of eⁱᵅ, eⁱᵞ.
But m′ ± m are integers, and with

    zp = (W + iZ)/√(W²+Z²) = e^{i(α+γ)/2},     zm = (Y − iX)/√(X²+Y²) = e^{i(α−γ)/2}

one has e^{i(m′α + mγ)} = zp^{m′+m} · zm^{m′−m} with integer exponents.  These are exactly the
`zp`, `zm` locals of `Quaternionic.to_euler_phases`.  The same formula holds for integer
indices, so `convert_H_to_D!` can be written once, driven by two `ComplexPowers` iterators
over zp and zm (powers up to 2ℓ instead of ℓ; negligible).  For half-integer indices one of
the two exponents is odd, so D(−R) = −D(R) is automatic — the double cover is respected
without any branch choice.  The singular fallbacks zp = 1 (β = π) and zm = 1 (β = 0) are
harmless because the corresponding d entries vanish.

The (−1)^s prefactor on sYlm is ±i for half-integer s.  **Resolved by the conventions pages
(details.md, "Half-integer indices", 2026-09-08):** the principal branch (−1)^s ≡ e^{iπs} = i^{2s}
is chosen, the conjugate form ₛYₗₘ = (−1)^s √((2ℓ+1)/4π) conj(𝔇_{m,−s}) is *the* definition,
and the 𝔇_{−m,s} form is an integer-index corollary (it is off by (−1)^{2s} otherwise).
`fill_Y!` must implement the conjugate form with `i^(2s)` computed from the integer `2s`.

### 5.5 Verification

Prototype: `notes/half-integer/prototype_2026-09-08.jl` (integer axis from `Deprecated.d_matrices`,
seed rows of §5.3, steps 4/5 of §5.2, symmetry fill with σ, phases of §5.4), run 2026-09-08 on
Julia 1.12.7.  Worst-case deviations over β ∈ {0.3, 1.1, 1.5, 2.0, 2.9}:

| Check | Range | Max deviation |
|---|---|---|
| vs Boyle (2016) `WignerDElement` (`test/conventions/boyle2016.jl`), d and full D with §5.4 phases, random rotors | J ≤ 15/2 | 2.5e-15 |
| vs Varshalovich Eq. 4.3.1(2) closed form (`test/conventions/varshalovich.jl`) | J ≤ 15/2 | 2.9e-15 |
| vs Varshalovich Tables 4.3–4.12 (720 transcribed entries) | J ≤ 9/2 | 3.6e-14 |
| character identity Σₘ d^J_{mm}(β) = sin((2J+1)β/2)/sin(β/2) (oracle-free) | J ≤ 61/2 | 1.8e-14 |
| Float64 vs BigFloat, absolute / relative on \|d\| > 1e-3 | J ∈ {21/2, 41/2, 61/2} | 8.9e-16 / 1.7e-13 |
| D(−R) + D(R) | J ≤ 15/2 | 0 (exact) |

An independent BigFloat prototype during review reproduced Varshalovich's closed forms to
1e-76 for J ≤ 25/2.  Stability is no worse than the integer G&D scheme: the seed has one
subtraction (the +1/2 row) with coefficients bounded by √2.


## 6. Half-integer test layer

Today there is **no test of package half-integer values** (nothing computes them).  What
exists: `test/conventions/varshalovich.jl` checks Varshalovich's own D against his index
relations for j ≤ 7 and his d against his tables for J ≤ 3/2 only (tables are transcribed to
9/2 but return `nothing` for the rows the book omits, and the negative-M reflection then
fails); `test/conventions/boyle2016.jl` checks `WignerDElement` against `Deprecated` for
integer ℓ only, so its half-integer branch is unverified against anything; the container
tests cover indexing only.

Stage 1 — reference vs reference (can be done now, in the comparisons work): put
`WignerDElement` and `d_½_explicit` in a shared test module; extend the table check to
J = 9/2 over all transcribed entries, guarding the untranscribed rows; cross-check the two
independent references for J ≤ 15/2 over β and random rotors, including D(−R) = −D(R).
This becomes the half-integer oracle.

Stage 2 — engine vs oracle (once `IT <: Rational` runs): `calc[ℓ][m′, m]` against the
stage-1 oracle; plus oracle-free metamorphic tests valid at large j: the character identity,
D(−R) = −D(R), the representation property D(R₁R₂) = D(R₁)D(R₂), unitarity, and Float64 vs
BigFloat at J ∈ {21/2, 41/2, 61/2}.  The prototype of §5.5 is the first version of these.
The Literate Varshalovich page then gains a "verified against the package" section.


## 7. What is deprecated and what comes back

`SphericalFunctions.Deprecated` was created by moving *everything* out of the way to start
clean, not as a judgment on each piece.  Code outside it is more authoritative; pieces come
back as follows.

| `src/deprecated/…` | Verdict | Why |
|---|---|---|
| `indexing.jl` (`Yindex`, `WignerDindex`, `WignerHindex`, sizes, ranges) | **Reuse the closed forms** inside the new indexable types; do not re-export the bare functions | Correct and tested (SymPy derivations in comments).  Flat-index call sites are what produced #41. |
| `iterators.jl` (`D_iterator`, `d_iterator`, `sYlm_iterator`, `λ_iterator`) | **Retire** | Superseded by `calc[ℓ]` views; `D_iterator` *is* #41. |
| `Hrecursion.jl` (`H!`, `H_recursion_coefficients`) | **Keep as the benchmark/regression baseline** until the new engine is verified against it; then retire | It is the whole-array reference the 2022 measurement favoured. |
| `evaluate.jl` (`d!`, `D!`, `Y!`, `sYlm_values`, `ₛ𝐘`, …) | **Retire**; port the `spin_factor`/ε logic into `fill_Y!` | Carries the old D convention (B3). |
| `associated_legendre.jl` (`ALFcompute!`, …) | **Undeprecate as a utility** if the SSHTRS λ-recursion still needs it; else retire | Independent of the D redesign. |
| `map2salm.jl` (Channel pool, `plan_map2salm`) | **Undeprecate the pattern**, re-target onto calculators (§4.6) | Already the intended parallelism model. |
| `ssht.jl`, `ssht/{direct,minimal,rs}.jl` | **Undeprecate**, re-target onto the new Wigner layer; rename `SSHTDirect → SSHTMatrix`; `SSHTOperator` (LinearOperators + GMRES/BiCGSTAB) stays a later idea | Working implementations; only their Wigner dependency changes. |
| `rotate.jl` (already commented out) | **Retire** or rewrite on the new views | Trivial once §4.3 exists. |
| utilities (`pixelizations`, `weights`, `complex_powers`, `operators`) | Already live outside `Deprecated` | `operators.jl` still needs the convention flip (B3). |
| `src/ssht/` scaffold (`huffenberger_wandelt.jl`, empty modules) | **Delete**; nothing in it is reachable | The deprecated SSHT code is the real starting point. |


## 8. Reconciliation of the old notes

| Topic | Variants found | Resolution |
|---|---|---|
| Argument order | `D(ℓₘₐₓ, R)` (`TODO.md`, `docs/src/index.md` TODO, notes `docs/todo.md`) vs `D(R, ℓₘₐₓ, …)` (old blueprint here, `interface_design.md`) | **`D(R, ℓₘₐₓ; …)`** |
| Limit arguments | `m′ₘₐₓ, m′ₘᵢₙ` positional (blueprint, `WignerHCalculator`) vs four keywords (`WignerCalculator`, `redesign2.jl`) vs `m′ₘₐₓ, mₘₐₓ` (docs TODO) | Four keywords, constructor-time only |
| Method names | `next!`/`compute!`/`calc[ℓ]` (blueprint) vs `recurrence!`/`w(ℓ)` (shipped) | `recurrence!` and `calc[ℓ]`; `next!`, `compute!` dropped |
| Return type | `Vector{Matrix{Complex}}` / `Matrix{Matrix}` (blueprint) vs rotor-first flat (vectorization notes, `HWedge`) | Rotor-first storage; views for access; nested containers only as convenience-function output |
| Type hierarchy | (a) `AbstractD/DMatrix/DOperator, Abstractd/…, AbstractP/…, AbstractH/HMatrix/HOperator/HWedge/HStorage, AbstractY/YVector` (original `redesign_notes.md`); (b) `WignerMatrix{NT,IT}` + `AbstractWignerMatrices{NT,IT,MT}` (same file, second sketch); (c) shipped `AbstractWignerMatrix{IT,NT,ST}` | (c) is canonical; (a) and (b) are dead.  `notes/docs/trash.md` describes (a) as "partially done" — it was abandoned, and the *shipped* scheme is what is half-built. |
| Rotor data | at construction (`WignerHCalculator`) vs at call (`WignerCalculator`) | One engine, rotors at call (§4.1) |
| Streaming per ℓ | headline feature (blueprint, `redesign_notes.md`) vs "4–5× slower at large ℓ, I give up" (`IterativeHrecursion.ipynb`, 2022) | Reproduced for Nᵣ = 1; **reversed** for Nᵣ ≥ 8 (§9).  Memory decides for per-ℓ regardless. |
| Workspace / Channel | active TODO (`TODO.md`, 2026-01) vs "speculative" (`trash.md`) | Calculators are the workspace (§4.6) |
| Half-integer status | "deferred" (`half_integer.md`, `trash.md`) vs half-built containers in `src/` | Route verified (§5); implement after the integer engine is complete |
| Half-integer relation | `redesign_half_integer.md`: Eq. (14) "wrong way", Eq. (15) "right way", valid for M′ ≠ J; `half_integer.md`: valid for M′ ≠ −J | Both are used, one per seed row; (15) is singular at M′ = −J.  Both transcriptions carry a sign error; the tested form in `test/conventions/varshalovich.jl` is correct. |
| Mode-weight orderings | three types (blueprint) vs none | One ordering, extensible type |
| Docs layout | `TODO.md` (intro/background/tutorials/…) vs Divio (`docs_outline.ipynb`) vs shipped `make.jl` | Keep the shipped skeleton; add a Tutorials section; low priority |
| SSHT names/files | `SSHTRS`, `SSHTDirect → SSHTMatrix`, `SSHTOperator`, an unplanned `huffenberger_wandelt.jl` | See §7 |
| Calculator prototypes | `HCalculator(T, ℓₘₐₓ)` (2022 notebooks), `DˡComputer` (`redesign.jl`), `WignerCalculator` (`redesign2.jl`), `WignerHCalculator` (`redesign3.jl`) | Lineage only; `WignerHCalculator` is the survivor |


## 9. Performance constraints and benchmark results

Constraints, all measured in the notes and still valid:

- No keyword arguments and no function-valued index arguments in hot loops: an unspecialized
  keyword index function cost 60–90× and up to 6 GiB of allocation
  (`notes/Hrecursor_allocations/`, Julia issue #45162).  The current hot path allocates 0 bytes.
- Innermost loop over rotors ⇒ rotor-first real storage.
- `Rational` loop indices cost a gcd per arithmetic operation, O(ℓ²) per ℓ — negligible for
  Nᵣ ≳ 8, visible at Nᵣ = 1.  Convert to `Int` "twice-index" offsets at loop entry.
- `sqrt` on the fly is not the bottleneck; steps 4/5 are (2022 profile).  Step 5 currently
  divides inside the rotor loop without `@fastmath`; hoist `inv(d̄ₗᵐ′⁻¹)`.
- LoopVectorization is effectively unmaintained on Julia ≥ 1.11; design so `@simd`/`@inbounds`
  suffice and treat `@turbo` as optional.
- Never form matrix inverses (5% error vs a few ε for LU); call `BLAS.set_num_threads`
  (43× on 56 cores); use `logbinomial`/`sqrtbinomial` from `utils.jl`, never `binomial`
  (overflows `Int` at ℓ ≈ 66 and Float64 at ℓ ≈ 1026).
- Float16 overflows at ℓ = 128; Float32 survives ℓ = 1024 on random rotors.

Benchmark (`notes/redesign/benchmark_per_ell_2026-09-08.jl` (output alongside), 2026-09-08, Julia 1.12.7, one thread, Float64,
β = 1.1; **preliminary** because B1 makes the batched axis slightly cheaper than it should be):
ns per (H element · rotor), `Deprecated.H!` whole-array vs `WignerHCalculator` looped over ℓ.

| ℓₘₐₓ | m′ₘₐₓ | Nᵣ = 1 | Nᵣ = 8 | Nᵣ = 64 |
|---|---|---|---|---|
| 8 | 8 | 3.65 → 4.39 (1.20×) | 3.56 → 0.95 (0.27×) | 3.58 → 0.38 (0.11×) |
| 64 | 64 | 0.85 → 2.24 (2.63×) | 0.86 → 0.68 (0.79×) | 0.92 → 0.31 (0.33×) |
| 64 | 2 | 1.03 → 3.01 (2.91×) | 1.12 → 0.84 (0.75×) | 1.01 → 0.42 (0.42×) |
| 256 | 256 | 0.50 → 2.25 (4.47×) | 0.51 → 0.69 (1.36×) | 0.51 → 0.30 (0.58×) |
| 256 | 2 | 0.85 → 2.68 (3.14×) | 0.85 → 0.81 (0.95×) | 0.84 → 0.42 (0.50×) |

Reading: the 2022 finding ("per-ℓ is 4–5× slower") is real for a single rotor at large ℓ,
and inverts for batches — parity near Nᵣ = 8, 2–3× faster per rotor at Nᵣ = 64.  The SSHT
use case (Nᵣ = number of rings ~ ℓₘₐₓ, m′ₘₐₓ = |s|) lands on the fast side.  The whole-array
approach needs ~2.7 GB at ℓₘₐₓ = 1000, which settles the design regardless.  Re-run the full
grid (ℓₘₐₓ ∈ {8, 64, 256, 1024}, m′ₘₐₓ ∈ {ℓₘₐₓ, 2}, Nᵣ ∈ {1, 8, 64, 512}) after B1 is fixed;
if Nᵣ = 1 stays > 2× slower, add coefficient precomputation to the checklist.


## 10. Coordination with the conventions and comparisons work

Three parallel efforts touch disjoint files: this memo, `TODO.md`, and the private notes repo
(interface); `docs/src/conventions/{summary,details,outline}.md` and
`docs/literate_input/conventions/calculations/` (conventions); `docs/literate_input/
conventions/comparisons/`, `comparisons.md`, `references.bib`, `CondaPkg.toml`, a scheduled
CI workflow, and the deletion of `test/conventions/` (comparisons).  The comparisons work
depends on the conventions work (its `@ref` targets and oracle) and on this memo (it deletes
`test/conventions/{boyle2016,varshalovich}.jl`, cited above; their contents move to
`boyle_2016.jl` and `varshalovich_1988.jl` Literate pages).

Hand-offs from this memo:

- To conventions: the two forms of the sYlm relation differ by (−1)^{2s} for half-integer s
  (§5.4); name one as the definition and fix the branch of (−1)^s.  **Done** — see §5.4.
- To comparisons: implement §6 Stage 1 in the Varshalovich and Boyle 2016 pages; keep
  `WignerDElement` half-integer capable; the (14)/(15) forms in `varshalovich.jl` are the
  correct ones, the notes' transcriptions are not.
- From conventions (received): the D convention is e^{−im′α} d e^{−imγ}; the R sign is
  `R_𝐮 f = −i d/dε f(𝐑 e^{−ε𝐮/2})` so `R_z = i∂_γ` with eigenvalue +s and `ð = R_x + iR_y`;
  `operators.jl` and `Deprecated` are stale on both counts (B3).


## 11. Implementation checklist (ordered)

1. Tests for both calculators against `Deprecated.H!` / `d_matrices` (there are none).
2. Fix B1 (step 2 uses ℓ for ℓ+1).
3. Fix B2 (`m′ₘᵢₙ` asymmetry in steps 5/6 and `convert_H_to_d!`).
4. Flip `operators.jl` to the settled R convention; mark `Deprecated.D!` as old-convention
   in its docstring (it is retired anyway).
5. `wedge_index_and_sign` helper; `fill_d!`, `fill_D!`; wire step 6 and phase conversion into
   `WignerHCalculator` so it produces d and D.
6. View type that does not violate the `AbstractArray` contract (B4); `Matrix(view)`.
7. `Base.similar(::WignerHCalculator)`; deduplicate size fields.
8. Rewrite `convert_H_to_D!` on `(zp, zm)` powers (§5.4); add the Quaternionic accessor.
9. Half-integer: generalized `ϵ` (B6), σ in step 6 and the helper, loop starts by `ℓₘᵢₙ`,
   the seed step of §5.3, integer-indexed axes, `validate_index_ranges` (B5), skip step 3;
   then §6 Stage 2 tests.
10. `sYlmCalculator` and `fill_Y!`.
11. `ModeWeights` type; move operators onto it; add `Lx, Ly, Rx, Ry`.
12. Un-deprecate `map2salm` and the SSHT types onto the new layer (§7); delete `src/ssht/`.
13. Re-run the §9 benchmark grid; decide on coefficient precomputation.
14. Retire the remaining `Deprecated` code; single test-project setup (`JuliaTestItems.toml`,
    drop `[targets]`).


## 12. Provenance and dating of the notes

`notes/` in this repository is a symlink into the private repo
`~/Research/Code/Notes/SphericalFunctions.jl`, whose history is a single bulk import
(2026-05-11) followed by a reorganization the same day that reduced `redesign/redesign_notes.md`,
`redesign/redesign2.md`, `operators/raising_lowering_K.md` and the notes `TODO.md` to one-line
stubs (originals recoverable with `git show 241438a:SphericalFunctions.jl/<path>`).  Git dates
are therefore meaningless for the notes; the usable stratigraphy is internal:

| Evidence | Date | Notes |
|---|---|---|
| Julia 1.7.2 kernels, Pluto 0.19.4, Julia issue #45162, discourse 80354 | spring 2022 | `Hrecursor_allocations/`, `IterativeHrecursion.ipynb` (the 4–5× result), `SpecializedArrays.ipynb`, `Untitled.ipynb` |
| Julia 1.8 kernel | 2022 | `Operators.ipynb` |
| Julia 1.9.3 kernels / Manifests | 2023 | `docs_outline.ipynb`, `find_overflow*.ipynb`, `operators/explicit_definition.jl` |
| Julia 1.10.x | 2024 | `golden_ratio_redux.ipynb`, `super_fibonacci.ipynb`, `issue_40/` |
| Julia 1.11.2 Manifests | late 2024 / early 2025 | notes root project, `SymPy/` |
| "Keefe sent to me on January 12, 2025" | 2025-01-12 | `lorentz/paper-6.txt` |
| `src/redesign/README.md` added / moved | 2025-12-11 / 2026-01-06 | package repo commits `5d29114`, `1722efa`; `Deprecated` created in `572c4c8` |
| Julia 1.12.x Manifests, `redesign{,2,3}.jl` naming lineage | late 2025 – early 2026 | `notes/redesign/`, `vectorization/` |
| `TODO.md` mtime | 2026-01-22 | untracked file; its conventions filenames (`euler_1776`, `whittaker_1927`, `wilson_1929`) are misremembered — the real files are `euler_1767`, `whittaker_1947`, `wilson_1921` |
| `lorentz/Jan_16_26-1.pdf` | 2026-01-16 | filename |
| `docs/*.md` consolidations | 2026-05-11 | derivative of the above; several statuses in them are corrected by this memo |
