# Common development tasks

## Running tests

The suite is made of [TestItems.jl](https://github.com/julia-vscode/TestItems.jl)
test items, so any of the usual runners will do.  Day to day, the quickest is
the `juliati` command-line runner, from the package root.  Its `--filter` option
takes a Julia expression in the variables `name`, `tags`, `filename` and
`package_name`, and runs the items for which it is true:

```bash
juliati                                             # every item, even the `:python` ones
juliati --filter '!(:python in tags)'               # every item but those, for everyday runs
juliati --filter 'occursin("ComplexPowers", name)'  # the items whose names mention ComplexPowers
juliati --filter ':slow in tags'
juliati --filter '!(:slow in tags)'
```

Unlike `Pkg.test`, described below, `juliati` knows nothing of the tags this
package gives special meaning, so a run without a filter includes the `:python`
items, and builds their Conda environment the first time.

VS Code's Julia extension lists the same items in its Testing panel
and runs them individually, and the `julia` MCP server exposes them to
editors and agents.  All three keep worker processes warm between
runs, so repeated runs skip Julia's startup and most compilation.

For coverage, and for the whole suite in one go, there is a script:

```bash
julia -t auto scripts/test.jl            # the whole suite
julia -t auto scripts/test.jl --coverage # ... and write lcov.info
```

It runs `Pkg.test`, so it skips the `:python` items described below
unless `:python` is passed as an argument, and it exits with a nonzero
status when any test fails — after writing `lcov.info`, if coverage
was requested.

Finally, `Pkg.test` works, because `test/runtests.jl` is a thin shim
over `@run_package_tests`:

```julia
using Pkg
Pkg.test("SphericalFunctions")
Pkg.test("SphericalFunctions"; test_args=[":python"])  # only the `:python` items
```

That shim supports filtering by tag only, written as `:sometag`, and
ignores any other argument with a warning; richer filtering belongs to
`juliati` rather than to a hand-written argument parser.  Under the
shim, items tagged `:python` run only when that tag is asked for,
because they build a Conda environment on first use; the scheduled
workflow asks for them.  Items tagged `:skipci` are skipped
automatically when the environment variable `CI` is `"true"`, unless
that tag is what was asked for — they need something continuous
integration does not have.  These rules belong to the shim alone;
other runners apply only the filters they are given.

Which files are searched for test items is set by
`JuliaTestItems.toml` in the package root: `src/`, `test/`, and the
Literate sources under `docs/literate_input/`.  The list is explicit
so that nothing else in the checkout is searched — in particular not
`notes/` or `docs/src/local_notes`, which are symbolic links to a
separate repository of notes — and `test/.CondaPkg/`, the Conda
environment that the `:python` items build, is excluded from `test/`.
The pages generated from the Literate sources are Markdown files, which
no runner searches for test items.


## Writing tests and coverage

Tags can be added to individual test items, which can then be used
either in the VS Code interface or the command line to include or
exclude certain tests.

```julia
@testitem "My testitem" tags=[:skipci, :slow] begin
    @test my_function() == expected_value
end
```

It's a well hidden fact that you can turn coverage on and off by
adding certain comments around the code you don't want to check:

```julia
# COV_EXCL_START
untested_code_that_wont_show_up_in_coverage()
# COV_EXCL_STOP
```


## Precompiling during development

The package has a precompile workload, written with
[PrecompileTools](https://github.com/JuliaLang/PrecompileTools.jl) in
`src/precompile.jl`.  It makes the calls with which most sessions
begin — `D`, `d`, `sYlm`, the operations on `ModeWeights`, and the
transforms — while the package is being precompiled, so that the code
they need is saved with the package, and their first calls in a
session take milliseconds rather than seconds.  The cost is paid at
each precompilation, which the workload makes more than twice as
long.  While `src/` is being edited, the workload can be turned off
with a preference, set in a file `LocalPreferences.toml` in the
package root:

```toml
[SphericalFunctions]
precompile_workload = false
```

That file is listed in `.gitignore`, so it applies to one checkout
only, and installations of the registered package still run the
workload.  From Julia 1.12 on, the test, docs and benchmark projects
are members of the package's workspace, and read the same file.
Julia 1.10 and 1.11 have no workspaces, so for the test project under
those versions the same two lines must also be in
`test/LocalPreferences.toml`, which `.gitignore` covers as well.


## Building the documentation

To build the documentation locally, run the following command from the
package root:

```bash
julia --project=. scripts/docs.jl
```

By default, this will build the documentation, run the doctests, and
launch a local server to view the docs in your web browser.
