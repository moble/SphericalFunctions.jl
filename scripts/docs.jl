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

import LiveServer: servedocs
literate_input = joinpath(pwd(), "docs", "literate_input")
literate_output = joinpath(pwd(), "docs", "src", "literate_output")
@info "Using input for Literate.jl from $literate_input"
servedocs(
    include_dirs=["src/"],  # So that docstring changes are picked up
    include_files=["docs/make_literate.jl"],
    skip_files=["docs/src/conventions/comparisons/lalsuite_SphericalHarmonics.md"],
    literate_dir=literate_input,
    launch_browser=true,
)
