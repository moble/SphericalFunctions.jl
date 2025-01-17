# Run this script from the top-level directory as
#
#   julia -t 4 --project=. scripts/docs.jl
#
# The docs will build and the browser should open automatically.  `LiveServer`
# will monitor the docs for any changes, then rebuild them and refresh the browser
# until this script is stopped.

using Revise

import Dates
println("Building docs starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd((@__DIR__) * "/..")
Pkg.activate("docs")

using LiveServer
literate_input = joinpath(pwd(), "docs", "literate_input")
literate_output = joinpath(pwd(), "docs", "src", "literate_output")
@info "Using input for Literate.jl from $literate_input"
servedocs(
    literate_dir = literate_input,
    skip_dir = literate_output,
    launch_browser=true
)
