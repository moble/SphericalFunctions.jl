### Currently, all the generated output goes into the `docs/src/literate_output` directory.
### This is nice just because it lets me add just that directory to the `.gitignore` file;
### I can't add `docs/src` to the `.gitignore` file because it would ignore all the
### non-generated files in that directory.  However, it would be nice to have the generated
### files in directories that are more consistent with the documentation structure.  For
### example, the `docs/src/literate_output/conventions/comparisons` directory contains files
### that really should be in the `docs/src/conventions/comparisons` directory.  I intend to
### reorganize the functionality in this file so that the outputs are in directories that
### are more consistent with the documentation structure.
###
### Instead of just plain for loops doing all the work, I will create functions that
### encapsulate things, and then call those functions in the for loops.  This will make it
### easier to add new functionality in the future, and make it easier to read the code.
###
### To deal with the gitignore issue, I will add a step that ensures the output file is
### listed in the `.gitignore` file.
###
### The files in the `docs/src/literate_input` directory will be rearranged so that they
### are in directories that are more consistent with the documentation structure.  For
### example, the `docs/literate_input/conventions/comparisons` directory will be
### moved to `docs/literate_input/conventions/comparisons`, and the output for every file in
### that directory will be sent to `docs/src/conventions/comparisons`.  There will no longer
### be a `docs/src/literate_output` directory; all the output will be in the same
### directory as the non-generated files.

literate_input = joinpath(@__DIR__, "literate_input")
skip_input_files = (  # Non-.jl files will be skipped anyway
    "ConventionsUtilities.jl",  # Used for TestItemRunners.jl
    "ConventionsSetup.jl",  # Used for TestItemRunners.jl
    "conventions_install_lalsuite.jl",  # lalsuite_2025.jl
)

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

for (root, _, files) ∈ walkdir(literate_input), file ∈ files
    # Skip some files
    if splitext(file)[2] != ".jl" || file ∈ skip_input_files
        continue
    end
    # full path to a literate script
    inputfile = joinpath(root, file)
    generate_markdown(inputfile)
end



# # See LiveServer.jl docs for this: https://juliadocs.org/LiveServer.jl/dev/man/ls+lit/
# literate_input = joinpath(@__DIR__, "literate_input")
# literate_output = joinpath(docs_src_dir, "literate_output")
# relative_literate_output = relpath(literate_output, docs_src_dir)
# relative_convention_comparisons = joinpath(relative_literate_output, "conventions", "comparisons")
# rm(literate_output; force=true, recursive=true)
# skip_input_files = (  # Non-.jl files will be skipped anyway
#     "ConventionsUtilities.jl",  # Used for TestItemRunners.jl
#     "ConventionsSetup.jl",  # Used for TestItemRunners.jl
# )
# for (root, _, files) ∈ walkdir(literate_input), file ∈ files
#     # Skip some files
#     if splitext(file)[2] != ".jl" || file ∈ skip_input_files
#         continue
#     end
#     # full path to a literate script
#     inputfile = joinpath(root, file)
#     # generated output path
#     output_path = splitdir(replace(inputfile, literate_input=>literate_output))[1]
#     # Ensure the output path is in .gitignore
#     ensure_in_gitignore(relpath(output_path, @__DIR__))
#     # generate the markdown file calling Literate
#     Literate.markdown(inputfile, output_path; documenter, mdstrings, execute)
# end

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
