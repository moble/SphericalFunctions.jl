# Common development tasks

## Running tests

The suite is made of [TestItems.jl](https://github.com/julia-vscode/TestItems.jl)
test items, so any of the usual runners will do.  Day to day, the quickest is
the `juliati` command-line runner, from the package root:

```bash
juliati                        # everything
juliati --filter 'ComplexPowers'
juliati --filter ':slow'
juliati --filter '!:slow'
```

VS Code's Julia extension lists the same items in its Testing panel
and runs them individually, and the `julia` MCP server exposes them to
editors and agents.  All three keep worker processes warm between
runs, so repeated runs skip Julia's startup and most compilation.

For coverage, and for the whole suite in one go, there is a script:

```bash
julia -t auto scripts/test.jl            # the whole suite
julia -t auto scripts/test.jl --coverage # ... and write lcov.info
```

Finally, `Pkg.test` works, because `test/runtests.jl` is a thin shim
over `@run_package_tests`:

```julia
using Pkg
Pkg.test("SphericalFunctions")
Pkg.test("SphericalFunctions"; test_args=[":python"])  # only the `:python` items
```

That shim supports filtering by tag only, written as `:sometag`;
richer filtering belongs to `juliati` rather than to a hand-written
argument parser.  Items tagged `:skipci` are skipped automatically
when the environment variable `CI` is `"true"`, unless that tag is
what was asked for — they need something continuous integration does
not have, such as a Python environment.

Which files are searched for test items is set by
`JuliaTestItems.toml` in the package root: `src/`, `test/`, and the
Literate sources under `docs/literate_input/`.  The generated copies
under `docs/src/` and `docs/build/` are deliberately not searched, and
neither is `notes/`, which is a symbolic link to a separate
repository.


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


## Building the documentation

To build the documentation locally, run the following command from the
package root:

```bash
julia --project=. scripts/docs.jl
```

By default, this will build the documentation, run the doctests, and
launch a local server to view the docs in your web browser.
