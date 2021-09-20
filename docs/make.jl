using Documenter, SphericalFunctions

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
        "Basics" => "manual.md",
        "Utilities" => "utilities.md",
    ],
    # doctest = false
)

deploydocs(
    repo="github.com/moble/SphericalFunctions.jl.git",
    devbranch="main"
)
