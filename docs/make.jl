# Run with
#   julia -t 4 --project=. scripts/docs.jl
# assuming you are in this top-level directory

using SphericalFunctions
using Documenter
using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    #style=:authoryear,
)

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
            "conventions/conventions.md",
            "conventions/comparisons.md",
        ],
        "Notes" => map(
            s -> "notes/$(s)",
            sort(readdir(joinpath(@__DIR__, "src/notes")))
        ),
        "References" => "references.md",
    ],
    #warnonly=true,
    #doctest = false
)

deploydocs(
    repo="github.com/moble/SphericalFunctions.jl",
    devbranch="main",
    push_preview=true
)
