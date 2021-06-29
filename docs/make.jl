using Documenter, Spherical

DocMeta.setdocmeta!(Spherical, :DocTestSetup, :(using Spherical); recursive=true)

makedocs(
    sitename="Spherical.jl",
    modules = [Spherical],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),  # Use clean URLs, unless built as a "local" build
        edit_link = "main",  # Link out to "main" branch on github
        #= UNCOMMENT ONCE STABLE:
        canonical = "https://moble.github.io/Spherical.jl/stable/",
        warn_outdated = true,
        :UNCOMMENT ONCE STABLE =#
    ),
    pages = [
        "Introduction" => "index.md",
        "Basics" => "manual.md",
    ],
    # doctest = false
)

deploydocs(
    repo="github.com/moble/Spherical.jl.git",
    devbranch="main"
)
