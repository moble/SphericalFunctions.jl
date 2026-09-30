# This file is intended to be included in the `docs/make.jl` file, and is responsible for
# converting the Literate-formatted scripts in `docs/literate_input` into
# Documenter-friendly markdown files in `docs/src`.

# Set up which files will be converted, and which will be skipped
skip_input_files = (  # Non-.jl files will be skipped anyway
    "ConventionsUtilities.jl",  # Used for TestItemRunners.jl
    "ConventionsSetup.jl",  # Used for TestItemRunners.jl
)
literate_input = joinpath(@__DIR__, "literate_input")

# The directories of `docs/src` that hold only generated pages.  The `.gitignore` file
# ignores both of them as a whole, so nothing written here is ever tracked.
generated_dirs = (
    joinpath(docs_src_dir, "30-conventions", "10-comparisons"),
    joinpath(docs_src_dir, "30-conventions", "20-calculations"),
)

# A page left behind by a Literate script that has since been renamed or deleted would
# still be listed in the navigation, because `make.jl` lists every Markdown file in these
# directories.  So each Markdown file without a source in `literate_input` is removed
# before the pages are generated.  The LALSuite source page is generated from a `.c` file,
# below.
for dir ∈ generated_dirs
    isdir(dir) || continue
    sourcedir = joinpath(literate_input, relpath(dir, docs_src_dir))
    for file ∈ readdir(dir)
        endswith(file, ".md") || continue
        file == "lalsuite_SphericalHarmonics.md" && continue
        isfile(joinpath(sourcedir, splitext(file)[1] * ".jl")) && continue
        @info "Removing $(joinpath(relpath(dir, docs_src_dir), file)), which has no Literate source"
        rm(joinpath(dir, file))
    end
end

# Generate markdown file for Documenter.jl from a Literate script
function generate_markdown(inputfile)
    # I've written the docs specifically to be consumed by Documenter; setting this option
    # enables lots of nice conversions.
    documenter=true
    # To support markdown strings, as in md""" ... """, we need to set this option.
    mdstrings=true
    # We *don't* want to execute the code in the literate script, because they are meant to
    # be used with TestItems.jl, and we don't want to run the tests here.
    execute=false
    # Output will be generated here.  The path is built from the part below
    # `literate_input`, so that the checkout's own path — often ending in
    # `SphericalFunctions.jl` — is never rewritten.
    outputdir = joinpath(docs_src_dir, dirname(relpath(inputfile, literate_input)))
    # Generate the markdown file calling Literate
    Literate.markdown(inputfile, outputdir; documenter, mdstrings, execute)
end

# Now, just walk through the literate_input directory and generate the markdown files for
# each literate script.
for (root, _, files) ∈ walkdir(literate_input), file ∈ files
    # Skip some files
    if splitext(file)[2] != ".jl" || file ∈ skip_input_files
        continue
    end
    # Full path to the literate script
    inputfile = joinpath(root, file)
    # Run the conversion
    generate_markdown(inputfile)
end

# Make "lalsuite_SphericalHarmonics.c" available in the docs
let
    inputfile = joinpath(literate_input, "30-conventions", "10-comparisons", "lalsuite_SphericalHarmonics.c")
    outputfile = joinpath(docs_src_dir, "30-conventions", "10-comparisons", "lalsuite_SphericalHarmonics.md")
    lalsource = read(inputfile, String)
    write(
        outputfile,
        "# LALSuite: Spherical Harmonics original source code\n"
        * "The official repository is [here]("
        * "https://git.ligo.org/lscsoft/lalsuite/-/blob/22e4cd8fff0487c7b42a2c26772ae9204c995637/lal/lib/utilities/SphericalHarmonics.c"
        * ")\n"
        * "```c\n$lalsource\n```\n"
    )
end
