# Run with
#   time julia --project=. make.jl && julia --project=. -e 'using LiveServer; serve(dir="build")'
# assuming you are in this `docs` directory (otherwise point the project argument here)

using SphericalFunctions
using Documenter

DocMeta.setdocmeta!(SphericalFunctions, :DocTestSetup, :(using SphericalFunctions); recursive=true)

makedocs(
    sitename="SphericalFunctions.jl",
    modules = [SphericalFunctions],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),  # Use clean URLs, unless built as a "local" build
        edit_link = "main",  # Link out to "main" branch on github
        canonical = "https://moble.github.io/SphericalFunctions.jl/stable/",
    ),
    pages = [
        "Introduction" => "index.md",
        "Primary functions" => "manual.md",
        "Differential operators" => "operators.md",
        "Utilities" => "utilities.md",
        "Notes" => map(
            s -> "notes/$(s)",
            sort(readdir(joinpath(@__DIR__, "src/notes")))
        ),
    ],
    #doctest = false
)

deploydocs(
    repo="github.com/moble/SphericalFunctions.jl",
    devbranch="main",
    push_preview=true
)
