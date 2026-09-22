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
using DocumenterInterLinks
using DocumenterCodeBlocks

docs_src_dir = joinpath(@__DIR__, "src")
package_root = dirname(@__DIR__)

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

# Plugin configuration
links = InterLinks(
    "Quaternionic" => "https://moble.github.io/Quaternionic.jl/stable/",
    "Julia" => "https://docs.julialang.org/en/v1/",
)
bib = CitationBibliography(
    joinpath(docs_src_dir, "references.bib");
    #style=:authoryear,
)
codeblocks = CodeBlocks(;line_counter=:continue)

# Run `make_literate.jl` to generate the literate files
include(joinpath(@__DIR__, "make_literate.jl"))
include("local_notes.jl")
(notes_pages, notes_remotes) = local_notes()

using SphericalFunctions

DocMeta.setdocmeta!(
    SphericalFunctions,
    :DocTestSetup,
    :(using SphericalFunctions);
    recursive=true,
    warn=false,
)

makedocs(
    plugins=[links, bib, codeblocks],
    sitename="SphericalFunctions.jl",
    modules = [SphericalFunctions],
    remotes=notes_remotes,
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),  # Use clean URLs, unless built as a "local" build
        edit_link = "main",  # Link out to "main" branch on github
        canonical = "https://moble.github.io/SphericalFunctions.jl/stable/",
        assets = String["assets/citations.css", "assets/extras.css"],
        size_threshold = 300 * 2^10, # 300 KiB
        size_threshold_warn = 200 * 2^10, # 200 KiB
        search_size_threshold_warn = 1 * 2^20, # 1 MiB
    ),
    pages = [
        "index.md",
        "Background" => [
            "10-background/01-domain.md",
            "10-background/02-operators.md",
            "10-background/03-sYlm_and_Dlmpm.md",
            "10-background/04-mode_weights.md",
        ],
        "Interface" => [
            "20-interface/01-half_integers.md",
            "20-interface/02-wigner_and_sYlm.md",
            "20-interface/03-transformations.md",
            "20-interface/04-operators.md",
            "20-interface/05-utilities.md",
        ],
        "Conventions" => [
            "30-conventions/01-summary.md",
            "30-conventions/02-details.md",
            "Comparisons" => map(
                s -> joinpath("30-conventions", "10-comparisons", s),
                sort(
                    filter(
                        s -> s != "lalsuite_SphericalHarmonics.md",
                        readdir(joinpath(docs_src_dir, "30-conventions", "10-comparisons"))
                    )
                )
            ),
            "Calculations" => map(
                s -> joinpath("30-conventions", "20-calculations", s),
                sort(readdir(joinpath(docs_src_dir, "30-conventions", "20-calculations")))
            ),
        ],
        "API" => [
            "40-api/01-internal.md",
            "40-api/02-functions.md",
        ],
        "Notes" => map(
            s -> joinpath("50-notes", s),
            sort(readdir(joinpath(docs_src_dir, "50-notes")))
        ),
        "Development" => [
            "60-development/01-index.md",
            "60-development/02-literate_testitems.md",
        ],
        "index_of_docstrings.md",
        "References" => "references.md",
        notes_pages...,
    ],
    warnonly=warnonly,
    #doctest = false,
    #draft=true,  # Skips running code in the docs for speed
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo="github.com/moble/SphericalFunctions.jl",
        devbranch="main",
        push_preview=true
    )
else
    @info "Skipping docs deployment because CI environment variable is not set."
end

println("Docs built in ", time() - start, " seconds.\n")
