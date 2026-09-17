# Run with
#   julia -t 4 --project=. scripts/docs.jl
# assuming you are in this top-level directory

# Pretty-print the current time
using Dates
println("\n")
@info """Building docs starting at $(Dates.format(Dates.now(), "HH:MM:SS"))."""

start = time()  # We'll display the total after everything has finished

using Documenter
using Literate
using DocumenterCitations


docs_src_dir = joinpath(@__DIR__, "src")
package_root = dirname(@__DIR__)

# Run `make_literate.jl` to generate the literate files
include(joinpath(@__DIR__, "make_literate.jl"))
include("local_notes.jl")
(notes_pages, notes_remotes) = local_notes()


# Documenter treats a broken doctest, a dead cross-reference or a missing docstring as an
# error, by default.  That is what we want from a normal build and from CI, so the default
# here is to fail.  While drafting, though, the `LiveServer` loop in `scripts/docs.jl` has
# to survive half-written pages, so that script sets the environment variable below.  A
# one-off lenient build can also be had with
#
#   julia --project=docs docs/make.jl warnonly
#
const warnonly = ("warnonly" in ARGS) || get(ENV, "SPHERICALFUNCTIONS_DOCS_WARNONLY", "false") == "true"
@info "Building with warnonly=$warnonly"


bib = CitationBibliography(
    joinpath(docs_src_dir, "references.bib");
    #style=:authoryear,
)

using SphericalFunctions

DocMeta.setdocmeta!(
    SphericalFunctions,
    :DocTestSetup,
    :(using SphericalFunctions);
    recursive=true,
    warn=false,
)

makedocs(
    plugins=[bib],
    sitename="SphericalFunctions.jl",
    modules = [SphericalFunctions],
    remotes=notes_remotes,
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),  # Use clean URLs, unless built as a "local" build
        edit_link = "main",  # Link out to "main" branch on github
        canonical = "https://moble.github.io/SphericalFunctions.jl/stable/",
        assets = String["assets/citations.css", "assets/extras.css"],
    ),
    pages = [
        "index.md",
        "Background" => [
            "background/domain.md",
            "background/operators.md",
            "background/sYlm_and_Dlmpm.md",
            "background/mode_weights.md",
        ],
        "Interface" => [
            "interface/wigner_matrices.md",
            "interface/sYlm.md",
            "interface/transformations.md",
            "interface/operators.md",
            "interface/utilities.md",
        ],
        "Conventions" => [
            "conventions/summary.md",
            "conventions/details.md",
            "Comparisons" => map(
                s -> joinpath("conventions", "comparisons", s),
                sort(
                    filter(
                        s -> s != "lalsuite_SphericalHarmonics.md",
                        readdir(joinpath(docs_src_dir, "conventions", "comparisons"))
                    )
                )
            ),
            "Calculations" => map(
                s -> joinpath("conventions", "calculations", s),
                sort(readdir(joinpath(docs_src_dir, "conventions", "calculations")))
            ),
        ],
        "API" => [
            "api/internal.md",
            "api/functions.md",
        ],
        "Notes" => map(
            s -> joinpath("notes", s),
            sort(readdir(joinpath(docs_src_dir, "notes")))
        ),
        "Development" => [
            "development/index.md",
            "development/literate_testitems.md",
        ],
        "index_of_docstrings.md",
        "References" => "references.md",
        notes_pages...,
    ],
    warnonly=warnonly,
    #doctest = false,
    #draft=true,  # Skips running code in the docs for speed
)

deploydocs(
    repo="github.com/moble/SphericalFunctions.jl",
    devbranch="main",
    push_preview=true
)

println("Docs built in ", time() - start, " seconds.\n")
