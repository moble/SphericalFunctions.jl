# Run this script from the top-level directory as
#
#   julia -t 4 --project=. scripts/docs.jl
#
# The docs will build and the browser should open automatically.  `LiveServer`
# will monitor the docs for any changes, then rebuild them and refresh the browser
# until this script is stopped.

import Revise
Revise.revise()

import Dates
println("Building docs starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

import Pkg
cd((@__DIR__) * "/..")
Pkg.activate("docs")

# Keep the `LiveServer` loop alive through half-written pages: broken doctests, dead
# cross-references and missing docstrings become warnings rather than errors.  A plain
# `docs/make.jl` build — and CI — still fails on all of them.
ENV["SPHERICALFUNCTIONS_DOCS_WARNONLY"] = "true"

import LiveServer: servedocs
literate_input = joinpath(pwd(), "docs", "literate_input")
@info "Using input for Literate.jl from $literate_input"
servedocs(
    include_dirs=["src/"],  # So that docstring changes are picked up
    include_files=["docs/make_literate.jl"],
    skip_files=["docs/src/30-conventions/10-comparisons/lalsuite_SphericalHarmonics.md"],
    literate_dir=literate_input,
    launch_browser=true,
)
