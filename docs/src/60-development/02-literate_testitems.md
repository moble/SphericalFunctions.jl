# Literate TestItems

This package computes things like Wigner's 𝔇 matrices.  If you've
ever dealt with these things — or even just rotations generally —
you'll know that the literature is an absolute quagmire of subtly
differing conventions wrapped up in terminology and notation from
hundreds of years ago.  In an effort to sort some of this out, I
decided to carefully compare the conventions of as many significant
sources as I could — everything from Wikipedia and Mathematica, to
current quantum-mechanics textbooks, all the way back to the original
books and papers that introduced some of these concepts as early as
1767!  So I went through each reference, and wrote a documentation
page for each one that carefully laid out the conventions used in that
source.  Those pages are collected in the "Comparisons" part of the
[Conventions](../30-conventions/01-summary.md) section of these docs.

But I didn't just want to *document* these conventions; I wanted to
*test* how they compared to the implementations in this package.  I
wanted to feed actual numbers into the actual expressions written down
by all these sources, and see if they agreed (or all too commonly, how
they disagreed) with this package's output.  This page describes how
that was done: each comparison page is a
[Literate.jl](https://github.com/fredrikekre/Literate.jl) script that
is simultaneously the source of a documentation page and a
[TestItems.jl](https://github.com/julia-vscode/TestItems.jl) test
item.  Getting the two tools to cooperate took a few tricks, which
are laid out below in the hope that they will be useful to anyone
else who wants documentation that tests itself this thoroughly.


## Why do it this way?

The pages and the tests have to be the same thing, because the tests
determine what the pages say.  While working through any given
reference, I would transcribe a formula, compare it numerically
against this package, and only then write down whether it agreed —
and if not, by what factor or sign.  Each page's summary box is the
conclusion of its tests.  Keeping the code and the prose in one file
also means that the conclusions are re-checked whenever the package
changes, so a page cannot quietly become false.

Documenter's doctests may seem like the obvious tool for code in
documentation, but they are not up to this job, for two reasons.  The
first is the sheer complexity of the code and tests.  A page typically
defines a whole module of functions transcribed from its reference,
and then loops over many angles and indices comparing those functions
to this package's, to within floating-point tolerances.  Doctests
compare *printed output*, which is the wrong notion of agreement for
floating-point results, and they are awkward for code that runs to a
couple hundred lines.  The second reason is granularity.  While
writing a page, I needed to run just the tests for the reference I was
working on — not the entire test suite — so that I could iterate
quickly; some of these literal transcriptions are also slow enough
that they should not hold anything else up.  Test items give exactly
this: each page is one item that can be run on its own, from the
editor or the command line, in a worker process that stays warm
between runs.

Those two requirements lead to Literate and TestItems.  Test items are
found by statically scanning Julia source files, so the tests have to
live in `.jl` files; the critical feature of Literate is that it takes
actual runnable Julia scripts as input — rather than markdown
containing fenced Julia code — and turns them into markdown for
Documenter.  The same file can therefore be scanned by the test
runners and converted into a page.  This does not appear to be
possible with Quarto.  It looks like it *would* work with Weave.jl; I
just happen not to use it, because I'm not sure how easily it
integrates with Documenter.jl and
[DocumenterCitations.jl](https://juliadocs.org/DocumenterCitations.jl/stable/),
on which the comparison pages rely heavily for their citations.


## What a page looks like

The sources of the comparison pages live in
`docs/literate_input/30-conventions/10-comparisons/`, and
`docs/make_literate.jl` converts them into markdown under `docs/src/`
at the start of each documentation build.  Stripped of its content, a
page looks like this:

```julia
md"""
# Some Author (1999)

!!! info "Summary"
    Some Author's spherical harmonics agree with this package's.

Several paragraphs of prose, with math and citations like
[Some Author's book](@cite SomeAuthor_1999)...
"""

using TestItems: @testitem  #hide
@testitem "Some Author conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module SomeAuthor
#+

# Prose quoting equation (12) of the reference, followed by its transcription:
import ..ConventionsUtilities: 𝒾, ❗
function Y(ℓ, m, θ, ϕ)
    # ...
end
#+

end  # module SomeAuthor
#+

# ## Tests
#
# Prose explaining what is compared, and why:
for (θ, ϕ) ∈ θϕrange()
    for (ℓ, m) ∈ ℓmrange(4)
        @test SomeAuthor.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

end  #hide
```

The long introduction is written as an `md"""..."""` string, which
Literate accepts when it is given the `mdstrings=true` option; this is
much easier to write than hundreds of lines each beginning with `# `.
Everything after it is a single test item, whose body alternates
between code and markdown comments.  The body is not indented, so that
the rendered code blocks read as ordinary top-level code rather than
as the inside of a macro call.  The
[Condon-Shortley](../30-conventions/10-comparisons/condon_shortley_1935.md)
page is a short, complete example, and its [source
file](https://github.com/moble/SphericalFunctions.jl/blob/main/docs/literate_input/30-conventions/10-comparisons/condon_shortley_1935.jl)
shows what it looks like before Literate converts it.

The rest of this page explains the odd-looking pieces of that
skeleton.


## The tricks

### The test-item wrapper is hidden, not removed

`make_literate.jl` calls `Literate.markdown` with `documenter=true`
and `execute=false`.  Literate itself therefore runs nothing; it
writes each run of code into an `@example` block, and it is
*Documenter* that evaluates those blocks when it builds the docs.
That would ordinarily mean that every test ran again in the
documentation build, in a sandbox module where the setup modules do
not exist, and failed.

The way out is that `@testitem`, `@testmodule`, and `@testsnippet` are
macros from TestItems.jl that expand to `nothing`.  They only exist so
that the code parses; the runners find test items by reading the
source, not by evaluating these macros.  So if Documenter evaluates
the `@testitem` line along with the rest of the page, the entire body
is swallowed by the macro and discarded, and nothing in the test item
runs during the build.  (This is also why `TestItems` is a dependency
of the docs environment.)

Of course, readers should not have to see that wrapper.  So the
`@testitem ... begin` line and its closing `end` go on lines of their
own, each ending in `#hide`.  Documenter evaluates lines marked
`#hide` but leaves them out of the rendered page.  (Literate keeps
those lines in its output only because `execute=false`; when Literate
executes the code itself, it strips them.)  The tempting alternative,
`#src`, is exactly wrong: Literate deletes `#src` lines entirely, so
Documenter would never see the macro, and would try to evaluate the
body as ordinary code.

### Every code chunk continues into the next

Literate splits its output into a new code block every time a markdown
comment interrupts the code, and Documenter evaluates each `@example`
block on its own.  The first block of a page would then contain an
opening `@testitem ... begin` with no matching `end`, which is a parse
error.

The fix is to put `#+` on a line by itself immediately after every run
of code that is followed by markdown.  Literate then marks that block
`continued = true`, which tells Documenter not to evaluate it yet, but
to hold onto its code and evaluate it together with the next
`@example` block of the same name.  The generated markdown looks like
this:

`````markdown
````@example condon_shortley_1935; continued = true
using TestItems: @testitem  #hide
@testitem "Condon-Shortley conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module CondonShortley
````

We'll also use some predefined utilities to make the code look more like the equations.

````@example condon_shortley_1935; continued = true
import ..ConventionsUtilities: 𝒾, ❗, dʲsin²ᵏθdcosθʲ
````
`````

Only the final block, containing the hidden `end`, is not continued,
so the whole test item is assembled and evaluated in one piece — at
which point `@testitem` discards it.  A forgotten `#+` shows up only
in the documentation build, as a parse error there; the test runners
never notice, because the source file is perfectly valid Julia.

### One module per reference

The expressions in the literature are, naturally, given very similar
or identical names — nearly every source has its own `Y` or `D` or `d`
— so each reference's formulas are defined in a module of their own,
which also lets the tests read as a direct comparison like
`CondonShortley.𝜙(...) ≈ ConventionsUtilities.Y(...)`.  Of course, I
couldn't have all these extra inefficient and inaccurate functions
cluttering up the actual package, so these modules are defined inside
the test items.  Because such a module is nested inside the module
that the test runner creates for the test item, it reaches the setup
modules with two dots, as in `import ..ConventionsUtilities: 𝒾`.

Occasionally one reference's formulas are needed by another page;
Varshalovich's tests, for example, also use the formulas from my own
2016 paper.  In that case the page defines its module with
`@testmodule` (hidden in the same way as `@testitem`), and any test
item that needs it lists it in its `setup`.

### Shared machinery lives beside the pages

Two files in the comparisons directory are not pages at all, and
`make_literate.jl` lists them in `skip_input_files` so that Literate
leaves them alone:

- `ConventionsSetup.jl` is a `@testsnippet` that seeds the random
  number generator.
- `ConventionsUtilities.jl` is a `@testmodule` holding everything the
  pages share.

The most important contents of `ConventionsUtilities` are the
functions `D`, `d`, and `Y`, which are the standard every page is
compared against.  They call this package, but they have their own
test item pinning them to the explicit formulas on the conventions
[Summary](../30-conventions/01-summary.md) page, so the pages are
guaranteed to compare the literature against the *documented*
conventions rather than against whatever the code happens to compute
today.  The module also defines a few notational conveniences that let
transcribed formulas look like the originals: `𝒾` is the imaginary
unit, and `❗` is a singleton object whose multiplication method
computes a factorial in `BigInt` arithmetic, so that `(ℓ+m)❗` — which
Julia parses as juxtaposition, and hence multiplication — reads just
like ``(ℓ+m)!``.  Finally, the generic angle and index ranges
(`θϕrange`, `ℓmrange`, and so on) come from the `Utilities` snippet in
`test/utilities/utilities.jl`, which is shared with the rest of the
test suite.

The pages import these names explicitly, even though TestItems.jl
already brings setup modules into scope with `using`.  The explicit
imports show the reader where each symbol came from; for the same
reason, the setup modules export nothing, so that nothing they define
can silently shadow a name defined in a page.

### Derivatives are computed, not rewritten

Many of these sources define their functions in terms of derivatives.
Condon and Shortley, for example, give the ``θ`` dependence of the
spherical harmonics in terms of ``d^{ℓ-m} \sin^{2ℓ}θ / d(\cos
θ)^{ℓ-m}``, and others use Rodrigues' formula for the Legendre
polynomials.  Rewriting such expressions in terms of standard
functions would subject each comparison to another round of convention
ambiguity, which is exactly what the pages are trying to eliminate.
Instead, `ConventionsUtilities` uses
[FastDifferentiation.jl](https://github.com/brianguenter/FastDifferentiation.jl),
which differentiates the expression symbolically and compiles the
result into an ordinary function.  It does a very good job of turning
these general formulas into useful implementations; the helpers
`dʲsin²ᵏθdcosθʲ` and `∂ⁿ` wrap it, and they have test items of their
own.

### Formulas are transcribed literally

The formulas on each page are copied verbatim from the reference,
including its notation, and changed only as far as Julia demands — for
example, by making arguments explicit where the reference leaves them
implicit.  When a comparison fails, the formula is not massaged until
it agrees; the disagreement is the finding, and it is what the page
reports.  The tests are what get adjusted, to express the relation
between the reference's conventions and this package's (a factor of
``(-1)^m``, say, or a complex conjugate).

### Finding and running the tests

`JuliaTestItems.toml` in the package root tells the runners where to
look for test items, and it lists `docs/literate_input/` alongside
`src/` and `test/`.  It deliberately does not list the generated
copies under `docs/src/`, which would otherwise be found twice.
(Those generated files are also added to `.gitignore` by
`make_literate.jl`, so that the Literate script is the only version
under version control.)  TestItemRunner.jl honors the same file, so
`Pkg.test` and continuous integration run the comparison pages along
with the rest of the suite.

While writing a page, the page alone can be run from VS Code's Testing
panel, or from the command line with, for example,

```bash
juliati --filter 'Condon-Shortley'
```

A few pages also cross-check against Python libraries (SciPy and
SymPy).  Those comparisons are separate test items tagged `:python`
and `:skipci`, because continuous integration does not have the Python
environment they need.


## Loose ends

- The documentation build does not currently run these tests; it
  relies on the test suite having done so.  TestItemRunner.jl
  [supports
  filtering](https://www.julia-vscode.org/docs/stable/userguide/testitems/#Filtering-support-in-TestItemRunner.jl),
  so the docs build could run just the comparison items before
  building the pages.
- Every page begins with `using TestItems: @testitem  #hide`, which
  only Documenter's evaluation needs.  Literate's `preprocess` option
  could insert that line during conversion instead, at the cost of
  making each source file slightly less self-explanatory.
