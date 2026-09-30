# Run this script from the top-level directory as
#
#   julia -t 4 --project=. scripts/docs.jl
#
# The docs will build and the browser should open automatically.  `LiveServer`
# will monitor the docs for any changes, then rebuild them and refresh the browser
# until this script is stopped.
#
# A change to a file in `src/` also starts a rebuild, but the rebuild shows the new
# docstrings only if `Revise` has updated the package that is already loaded.  No project
# in this repository lists `Revise`, so it is used when it can be loaded from the global
# environment; without it, the docs still build, and edits to docstrings appear once this
# script is restarted.

revise_available = try
    import Revise
    true
catch
    @info "Revise is not available; edits to docstrings will appear only after a restart."
    false
end

import Dates
println("Building docs starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

import Pkg
cd((@__DIR__) * "/..")
Pkg.activate("docs")

# Keep the `LiveServer` loop alive through half-written pages: broken doctests, dead
# cross-references and missing docstrings become warnings rather than errors.  A plain
# `docs/make.jl` build — and CI — still fails on all of them.
ENV["SPHERICALFUNCTIONS_DOCS_WARNONLY"] = "true"

# Outside the REPL, `Revise` applies an edit only when `Revise.revise()` is called, and
# nothing in a rebuild calls it.  This task calls it as soon as `Revise` sees a change to a
# file it tracks, which is well before the rebuild that `LiveServer` starts for the same
# change reaches the docstrings.  An edit that does not parse is reported and skipped, and
# the loop continues.
if revise_available
    @async while true
        wait(Revise.revision_event)
        Revise.revise()
    end
end

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
