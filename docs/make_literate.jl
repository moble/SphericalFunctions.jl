# This file is intended to be included in the `docs/make.jl` file, and is responsible for
# converting the Literate-formatted scripts in `docs/literate_input` into
# Documenter-friendly markdown files in `docs/src`.

# Set up which files will be converted, and which will be skipped
skip_input_files = (  # Non-.jl files will be skipped anyway
    "ConventionsUtilities.jl",  # Used for TestItemRunners.jl
    "ConventionsSetup.jl",  # Used for TestItemRunners.jl
)
literate_input = joinpath(@__DIR__, "literate_input")

# Ensure a file is listed in the .gitignore file
function ensure_in_gitignore(file_path)
    gitignore_path = joinpath(package_root, ".gitignore")
    if isfile(gitignore_path)
        existing_entries = readlines(gitignore_path)
        if file_path in existing_entries
            return  # File is already listed
        end
    end
    open(gitignore_path, "a") do io
        write(io, file_path * "\n")
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
    # Output will be generated here:
    outputfile = replace(inputfile, "literate_input"=>"src")
    outputdir = dirname(outputfile)
    # Ensure the output path is in .gitignore
    ensure_in_gitignore(relpath(replace(outputfile, ".jl"=>".md"), package_root))
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
    inputfile = joinpath(literate_input, "conventions", "comparisons", "lalsuite_SphericalHarmonics.c")
    outputfile = joinpath(docs_src_dir, "conventions", "comparisons", "lalsuite_SphericalHarmonics.c")
    ensure_in_gitignore(relpath(replace(outputfile, ".c"=>".md"), package_root))
    lalsource = read(
        joinpath(literate_input, "conventions", "comparisons", "lalsuite_SphericalHarmonics.c"),
        String
    )
    write(
        joinpath(docs_src_dir, "conventions", "comparisons", "lalsuite_SphericalHarmonics.md"),
        "# LALSuite: Spherical Harmonics original source code\n"
        * "The official repository is [here]("
        * "https://git.ligo.org/lscsoft/lalsuite/-/blob/22e4cd8fff0487c7b42a2c26772ae9204c995637/lal/lib/utilities/SphericalHarmonics.c"
        * ")\n"
        * "```c\n$lalsource\n```\n"
    )
end
