# Common development tasks

## Running tests

You can run all tests with coverage and process that coverage *from
the package root* with:

```bash
julia --project=. scripts/test.jl
```

Optionally, at the end of this command, specify names of individual
tests to run (in quotes if there are spaces), tags of tests to run
(which must start with a colon), or files to run all tests in (which
must end with `.jl`).  If any are specified, only matching tests or
files will be run.  By default, all tests in all files will be run.

Optionally, either with or without any of the above specifications,
add `--skip` followed by one or more tests, tags, or files to skip.
These override any inclusion criteria specified earlier in the
command.

The names of individual tests or files can be given as regex patterns
(probably in quotes), and all such matches will be via `occursin`
matching, so that partial matches will work.  Tags must be given
exactly as they appear in the code (including the colon).

Note that the tests with the `:skipci` tag will be skipped whenever
the environment variable `CI` is set to "true" (which is the case on
GitHub Actions), even if it is explicitly included in the command
line.

Here are some example invocations:

```bash
# Run all tests in all files
julia --project=. scripts/test.jl

# Run everything in complex_powers.jl
julia --project=. scripts/test.jl complex_powers.jl

# Run only the ComplexPowers test
julia --project=. scripts/test.jl ComplexPowers

# Run everything in complex_powers.jl except ComplexPowers
julia --project=. scripts/test.jl complex_powers.jl --skip ComplexPowers

# Run everything in every file except ComplexPowers
julia --project=. scripts/test.jl --skip ComplexPowers

# Run only tests tagged :fast
julia --project=. scripts/test.jl :fast

# Run everything except tests tagged :slow
julia --project=. scripts/test.jl --skip :slow
```


## Writing tests and coverage

Tags can be added to individual test items, which can then be used either in the VS Code interface or the command line to include or exclude certain tests.

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
