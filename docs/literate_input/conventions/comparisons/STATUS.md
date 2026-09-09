# Status of the literature-comparison pages

Working notes for converting the comparisons in `docs/src/conventions/comparisons.md` and the
test modules in `test/conventions/` into Literate pages in this directory.  Each page renders
as a documentation page *and* runs as a `@testitem`.  Update this table as work proceeds so a
later session can resume.

Ground truth: `docs/src/conventions/{summary,details}.md`.  Oracle: the `D`, `d`, `Y`
functions in `ConventionsUtilities.jl`, which implement the settled conventions (currently via
`conj(Deprecated.D)`; see the TODO there).

Column meanings:

- **Page**: Literate `.jl` page written in this directory.
- **Tests**: the page's `@testitem`s pass locally.
- **Old .md**: the files under `docs/src/conventions/comparisons/` are *generated* by
  `docs/make_literate.jl` and gitignored; there was never a hand-written one to delete, so this
  column is always "n/a".
- **comparisons.md**: the reference's section in `docs/src/conventions/comparisons.md` has been
  removed and a row added to the summary table there.
- **test/conventions**: the corresponding `test/conventions/<ref>.jl` file has been deleted
  (or never existed).

## Already finished before this effort (12 pages)

`blanchet_2024`, `clifford_1878`, `cohen_tannoudji_1991`, `condon_shortley_1935`,
`euler_1767`, `gibbs_1881`, `hamilton_1844`, `lalsuite_2025`, `ninja_2011`, `tait_1868`,
`whittaker_1947`, `wilson_1921`.  (`lalsuite_2025` was updated in this effort to compare
against the settled-convention oracle instead of `conj(Deprecated.D)`.)

## References being converted

| # | Reference | Page | Tests | Old .md | comparisons.md | test/conventions | Notes |
|---|-----------|------|-------|---------|----------------|------------------|-------|
| 1 | Newman & Penrose (1966) | ☑ `newman_penrose_1966.jl` | ☑ | n/a | ☑ | never existed | only Eq. (3.8) numbered; other NP equation numbers need confirming against the paper (`#src` TODO in page) |
| 2 | Goldberg et al. (1967) | ☑ `goldberg_et_al_1967.jl` | ☑ | n/a | ☑ | ☑ deleted | ``ₛY`` = ``(-1)^m`` ours; ``D`` = conj with α↔γ |
| 3 | Thorne (1980) | ☑ `thorne_1980.jl` | ☑ | n/a | ☑ | ☑ deleted | agrees; MTW out of scope |
| 4 | Wikipedia | ☑ `wikipedia_2026.jl` | ☑ | n/a | ☑ | ☑ deleted | agrees; accessed 2026-09-08; 𝒫 = −R |
| 5 | Sakurai (1994) | ☑ `sakurai_1994.jl` | ☑ | n/a | ☑ | ☑ deleted | agrees |
| 6 | Shankar (1994) | ☑ `shankar_1994.jl` | ☑ | n/a | ☑ | ☑ deleted | agrees; D via matrix exponentials of J |
| 7 | Zettili (2009) | ☑ `zettili_2009.jl` | ☑ | n/a | ☑ | never existed | agrees, incl. rotation law (7.70) |
| 8 | Edmonds (1960) | ☑ `edmonds_1960.jl` | ☑ | n/a | ☑ | ☑ deleted | 𝒟 = 𝔇(R⁻¹) = conj transpose; d transposed; bib key `Edmonds_2016` |
| 9 | Varshalovich et al. (1988) | ☑ `varshalovich_1988.jl` | ☑ | n/a | ☑ | ☑ deleted | agrees; rendered `@testmodule Varshalovich`; Stage 1 half-integer oracle (tables to 9/2, recursion, cross-check vs Boyle 2016 to 15/2); Ĵ′ = (−R_x, R_y, −R_z) tested |
| 10 | Wigner (1959) | ☑ `wigner_1959.jl` | ☑ | n/a | ☑ | ☑ deleted | (−1)^{μ′−μ} conj(𝔇) |
| 11 | Torres del Castillo (2003) | ☑ `torres_del_castillo_2003.jl` | ☑ | n/a | ☑ | ☑ deleted | agrees |
| 12 | Griffiths (1995) | ☑ `griffiths_1995.jl` | ☑ | n/a | ☑ | never existed | agrees |
| 13 | Le Bellac (2006) | ☑ `le_bellac_2006.jl` | ☑ | n/a | ☑ | never existed | agrees; D via matrix exponentials of J |
| 14 | Mathematica | ☑ `mathematica_2026.jl` | ☑ | n/a | ☑ | ☑ deleted | docs accessed 2026-09-08; WignerD = 𝔇_{−m₁,−m₂} (identities only) |
| 15 | SymPy | ☑ `sympy_2026.jl` | ☑ (incl. `:python` cross-check vs SymPy 1.14.0, run locally 2026-09-09) | n/a | ☑ | never existed | wigner_d = 𝔇_{−m′,−m} |
| 16 | SciPy | ☑ `scipy_2026.jl` | ☑ (incl. `:python` cross-check vs SciPy 1.18.0, run locally 2026-09-09) | n/a | ☑ | never existed | agrees; docs quoted from 1.18.1 |
| 17 | NIST DLMF | ☑ `nist_dlmf_2026.jl` | ☑ | n/a | ☑ (no section existed) | ☑ deleted | agrees; Release 1.2.7 (2026-06-15) |
| 18 | Boyle (2016) | ☑ `boyle_2016.jl` | ☑ | n/a | ☑ (no section existed) | ☑ deleted | conj of ours; rendered `@testmodule Boyle2016`; Stage 1 cross-check lives on Varshalovich page |

## Infrastructure

| Item | Done |
|------|------|
| Oracle `D`, `d`, `Y` in `ConventionsUtilities.jl` + pinning testitem | ☑ (tests pass) |
| `lalsuite_2025.jl` uses the oracle (no `conj`, no `@test_broken`) | ☑ (tests pass; NINJA, Condon-Shortley, Cohen-Tannoudji re-run and pass) |
| `test/CondaPkg.toml` (sympy, scipy for `:python` items) | ☑ |
| `.github/workflows/scheduled.yml` (quarterly: full suite, `:python` items, docs; LTS + current) | ☑ written (not yet exercised) |
| `references.bib` entries for Wikipedia, Mathematica, SymPy, SciPy | ☑ (plus DLMF release bumped to 1.2.7) |
| `comparisons.md`: sections removed, summary table added | ☑ all per-reference sections removed; table has 31 rows |
| `test/conventions/` directory removed | ☑ |
| Docs build clean | ☑ no warnings from the comparison pages (2026-09-09); remaining warnings are pre-existing (`local_notes`, duplicate operator docstrings, unresolved "Differential operators" refs in `interface/`) |

## Items for the author to confirm against the sources

These are places where the page text goes beyond what the repository or the fetched
documentation recorded; each is also marked with a non-rendered `#src` TODO line in the page.

- **Newman & Penrose (1966)**: only the ð/ð̄ definition carries an equation number, (3.8).
  Needed: the equation numbers (or section) for the stereographic coordinate ζ, the ``m^μ``
  quotation and the spin-weight law ``η' = e^{isψ}η``, and the "∝" expression for ``ₛY_{ℓm}``
  in terms of ζ.  The page states the relation ``ₛY_ζ = (-1)^ℓ e^{-isϕ}
  √(4π/[(2ℓ+1)(ℓ+m)!(ℓ-m)!])\, ₛY_{ℓm}`` as found numerically from the transcribed formula.
- **Goldberg et al. (1967)**: the raising relation ``ð ₛY = √… ₛ₊₁Y`` is attributed only
  collectively to "Eqs. (2.7)" (the lowering one is (2.7b) per details.md); confirm (2.7a) and
  the equation number of their ð definition.
- **Le Bellac (2006)**: the page does not transcribe a ``d`` formula and tests his definition
  via matrix exponentials; confirm whether he gives one worth transcribing.
- **Varshalovich et al. (1988)**: confirm the section/page of the quoted definition of
  ``Ĵ`` and of "Eq. (12)" for the commutators (both quoted from the old comparisons.md
  without a location).
- **Wigner (1959)**: Eq. (A.11) is quoted in the schematic form from the old notes
  (``Y = c (-1)^m 𝔇(ϕ,θ,0)_{m0}``); confirm the exact statement and normalization.
- **Mathematica**: the `WignerD` relation ``𝔇_{-m_1,-m_2}`` is inferred from documented
  identities only.  Evaluating, e.g., `WignerD[{2, 1, -1}, ψ, θ, ϕ]` and `WignerD[{1, 1, 0},
  ψ, θ, ϕ]` in Mathematica and comparing with ``𝔇_{-m_1,-m_2}`` would pin the general
  element; also confirm `SphericalHarmonicY[1, -1, θ, ϕ]` is ``+√(3/8π) sin θ e^{-iϕ}``
  (the documentation's formula was fetched with ``|m|`` in one rendering and without in the
  old notes).
- **SciPy**: the legacy `sph_harm` page is gone from the current docs; the deprecation text is
  quoted from the 1.15.0 documentation.

## Out-of-scope items noticed (owned elsewhere)

- `src/utilities/operators.jl` (`Rz`, `R₊`, `R₋`, `ð`, `ð̄`) and `Deprecated.D!` still carry
  the pre-3.0 conventions; flip is item 4 of the checklist in `src/redesign/README.md`.
- `docs/src/background/sYlm_and_Dlmpm.md` says `L_z 𝔇 = m′ 𝔇` (settled: `−m′`).
- `docs/src/notes/normalization.md` relates `ₛY` to `𝔇_{m,s}` (settled: `(−1)^s √… conj(𝔇_{m,−s})`).
- Root `TODO.md` lists conventions file names that do not exist (`euler_1776`, `whittaker_1927`,
  `wilson_1929`).
