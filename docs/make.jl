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


bib = CitationBibliography(
    joinpath(docs_src_dir, "references.bib");
    #style=:authoryear,
)

using SphericalFunctions
DocMeta.setdocmeta!(SphericalFunctions, :DocTestSetup, :(using SphericalFunctions); recursive=true)

makedocs(
    plugins=[bib],
    sitename="SphericalFunctions.jl",
    modules = [SphericalFunctions],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),  # Use clean URLs, unless built as a "local" build
        edit_link = "main",  # Link out to "main" branch on github
        canonical = "https://moble.github.io/SphericalFunctions.jl/stable/",
        assets = String["assets/citations.css", "assets/extras.css"],
    ),
    pages = [
        "index.md",
        "transformations.md",
        "wigner_matrices.md",
        "sYlm.md",
        "operators.md",
        "utilities.md",
        "API" => [
            "internal.md",
            "functions.md",
        ],
        "Conventions" => [
            "conventions/summary.md",
            "conventions/details.md",
            "conventions/comparisons.md",
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
        "Notes" => map(
            s -> joinpath("notes", s),
            sort(readdir(joinpath(docs_src_dir, "notes")))
        ),
        "References" => "references.md",
    ],
    #warnonly=true,
    #doctest = false,
    #draft=true,  # Skips running code in the docs for speed
)

deploydocs(
    repo="github.com/moble/SphericalFunctions.jl",
    devbranch="main",
    push_preview=true
)

println("Docs built in ", time() - start, " seconds.\n")
