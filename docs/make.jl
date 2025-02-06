# Run with
#   julia -t 4 --project=. scripts/docs.jl
# assuming you are in this top-level directory

# Pretty-print the current time
using Dates
@info """Building docs starting at $(Dates.format(Dates.now(), "HH:MM:SS"))."""

start = time()  # We'll display the total after everything has finished

using Documenter
using Literate

docs_src_dir = joinpath(@__DIR__, "src")

# See LiveServer.jl docs for this: https://juliadocs.org/LiveServer.jl/dev/man/ls+lit/
literate_input = joinpath(@__DIR__, "literate_input")
literate_output = joinpath(docs_src_dir, "literate_output")
for (root, _, files) ∈ walkdir(literate_input), file ∈ files
    # ignore non julia files
    splitext(file)[2] == ".jl" || continue
    # full path to a literate script
    input_path = joinpath(root, file)
    # generated output path
    output_path = splitdir(replace(input_path, literate_input=>literate_output))[1]
    # generate the markdown file calling Literate
    Literate.markdown(input_path, output_path, documenter=true, mdstrings=true)
end
relative_literate_output = relpath(literate_output, docs_src_dir)

using DocumenterCitations

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
            "Calculations" => [
                joinpath(relative_literate_output, "euler_angular_momentum.md"),
            ],
        ],
        "Notes" => map(
            s -> "notes/$(s)",
            sort(readdir(joinpath(docs_src_dir, "notes")))
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

println("Docs built in ", time() - start, " seconds.")
