#! /usr/bin/env bash

# This runs the test suite (or any tests matching the args to this command)
# with coverage enabled, then processes that coverage into `lcov.info`, and
# finally deletes all the *.cov files generated during the tests.

time julia --project=test --code-coverage=user test/runtests.jl "$@"
julia --project=test test/coverage.jl
find . -name \*.cov -exec /bin/rm {} \+
