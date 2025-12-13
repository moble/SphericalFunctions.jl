# Literate TestItems

TestItemRunner.jl + Literate.jl = 💪

(With honorable mention going to the excellent
[DocumenterCitations.jl](https://juliadocs.org/DocumenterCitations.jl/stable/)
and
[FastDifferentiation.jl](https://github.com/brianguenter/FastDifferentiation.jl)
packages.)

I have a package called
[SphericalFunctions.jl](https://moble.github.io/SphericalFunctions.jl/dev/)
that computes things like Wigner D matrices and spherical harmonics.
If you've ever dealt with these things — or even just rotations
generally — you'll know that the literature is an absolute quagmire of
subtly differing conventions wrapped up in terminology and notation
from hundreds of years ago.  In an effort to sort some of this out, I
decided to carefully compare the conventions of as many significant
sources as I could — everything from Wikipedia and Mathematica, to
current quantum-mechanics textbooks, all the way back to the original
books and papers that introduced some of these concepts.  So I went
through each reference, and wrote a documentation page for each one
that carefully laid out the conventions used in that source.  The
documentation for my package has a whole section collecting all of
these different pages.

But I didn't just want to *document* these conventions; I wanted to
*test* how they compared to the implementations in my package.  I
wanted to feed actual numbers into the actual expressions written down
by all these sources, and see if they agreed (or all too commonly, how
they disagreed) with my package's output.  But I wanted thorough tests
— too much for simple doctests.

All these expressions are, naturally, given very similar or identical
names in the literature, which could lead to collisions; we have to
define them in a separate module for each reference.  Of course, I
couldn't have all these extra inefficient and inaccurate functions
cluttering up my actual package, so these modules would have to be
defined elsewhere.

I wanted to be able to run the tests automatically, as part of the
documentation-build process and on their own.  Of course, I also
wanted to be able to test what I was writing *as I was going through*
and writing the page for any given reference, because the results of
those tests would determine what I should write in the documentation.
This meant that I needed the test suite to be as modular and fast as
possible — preferably testing just the reference I was working on,
rather than my entire test suite, so that I could iterate quickly.

* The critical feature of Literate.jl is that it takes actual runnable
  Julia scripts as input — rather than markdown that has fenced Julia
  code, for example — which allows TestItemRunner.jl to find the tests
  and run them properly.  This does not appear to be possible with
  Quarto.  It looks like it *would* work with Weave.jl; I just happen
  to not use it in my workflow (because I'm not sure how easy it is to
  integrate with Documenter.jl and DocumenterCitations.jl).
* TODO: Use TestItemRunner as part of doc building process, with
  [filtering](https://www.julia-vscode.org/docs/stable/userguide/testitems/#Filtering-support-in-TestItemRunner.jl)
  to select tests in the docs.
* TODO: I'm adding `using TestItems: @testmodule, @testitem  #hide` at
  the top of each page; could this be put somewhere else?
* `@testmodule` and `@testitem` and their corresponding `end`
  statements go on their own lines, followed by `#hide` to make the
  output look nicer, while still functioning properly.  Note that
  `#src` is *not* the right thing to do, as it will fail to provide
  those lines to Literate's execution, which needs them in order to
  skip them — which we presumably want to do so that the tests won't
  run and possibly pollute the namespace.
* Inside those blocks, all the code output has to end up in the same
  `@example` section, because they're all part of the same `@testitem`
  or whatever, which means we have to put `#+` on its own line
  immediately after any code and before any non-code (markdown).
* Mention how FastDifferentiation makes it easy to test expressions
  that explicitly use derivatives.
* Mention DocumenterCitations.jl for managing citations in the docs.
