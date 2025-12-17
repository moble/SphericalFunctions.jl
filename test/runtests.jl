# See docs/src/developments/index.md for details of how to run tests with this script.

using TestItemRunner
using ArgParse

function parse_commandline()
    settings = ArgParseSettings()
    @add_arg_table! settings begin
        # Collect everything before the optional `--skip`
        "run"
            nargs    = '*'
        # Collect everything after the optional `--skip`
        "--skip"
            nargs    = '*'
            default  = String[]
    end
    parsed_args = parse_args(settings)
    run_files = Tuple(Regex(s) for s ∈ parsed_args["run"] if endswith(s, ".jl"))
    run_tags = Tuple(Symbol(s[2:end]) for s ∈ parsed_args["run"] if startswith(s, ":"))
    run_tests = Tuple(Regex(s) for s ∈ parsed_args["run"] if !endswith(s, ".jl") && !startswith(s, ":"))
    skip_files = Tuple(Regex(s) for s ∈ parsed_args["skip"] if endswith(s, ".jl"))
    skip_tags = Tuple(Symbol(s[2:end]) for s ∈ parsed_args["skip"] if startswith(s, ":"))
    skip_tests = Tuple(Regex(s) for s ∈ parsed_args["skip"] if !endswith(s, ".jl") && !startswith(s, ":"))
    return run_files, run_tags, run_tests, skip_files, skip_tags, skip_tests
end

const run_files, run_tags, run_tests, skip_files, skip_tags, skip_tests = parse_commandline()

# Get the `CI` environment variable, defaulting to "false" if not set
const CI = get(ENV, "CI", "false") == "true"

# Create the function that will filter which tests to run
function filter(testitem)
    # Destructure the input NamedTuple.  Note that `filename` is the full path to the file
    # containing the test item, and `name` is the full string used to name the test item,
    # while `tags` is a vector of `Symbol`s that can be used to tag test items.
    (; filename, name, tags) = testitem

    for skip ∈ skip_files
        if occursin(skip, filename)
            @info "Skipping test '$name' in file '$(relpath(filename))' " *
                "due to skip filter '$skip'."
            return false
        end
    end

    for skip ∈ skip_tags
        if tag ∈ tags
            @info "Skipping test '$name' tagged '$tag' due to skip filter."
            return false
        end
    end

    for skip ∈ skip_tests
        if occursin(skip, name)
            @info "Skipping test '$name' due to skip filter '$skip'."
            return false
        end
    end

    if !isempty(run_files) || !isempty(run_tags) || !isempty(run_tests)
        for run ∈ run_files
            if occursin(run, filename)
                # @warn "Dry run: including test '$name' in file '$(relpath(filename))' " *
                #     "due to run filter '$run'."
                # return false
                return true
            end
        end
        for run ∈ run_tags
            if run ∈ tags
                # @warn "Dry run: including test '$name' tagged '$run' due to run filter."
                # return false
                return true
            end
        end
        for run ∈ run_tests
            if occursin(run, name)
                # @warn "Dry run: including test '$name' due to run filter '$run'."
                # return false
                return true
            end
        end
        # @info "Excluding test '$name' in file '$(relpath(filename))' because " *
        #     "it does not match any requested tests."
        return false
    end

    if CI && :skipci ∈ tags && :skipci ∉ run_tags
        @info "Skipping test '$name' tagged ':skipci' because `CI` is true."
        return false
    end

    return true
end

@info "Filtering tests with" run_files run_tags run_tests skip_files skip_tags skip_tests CI

@run_package_tests verbose=true filter=filter
