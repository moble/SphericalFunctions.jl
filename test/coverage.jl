# To see coverage, first run tests as
#     time julia --project=test --code-coverage=user test/runtests.jl
# Then run this script as
#     julia --project=test test/coverage.jl
# You may want to then delete all other coverage files with
#     find . -name \*.cov -exec /bin/rm {} \+
# Also note that the resulting `lcov.info` file can be used by VS Code's "Coverage Gutters"
# extension to view coverage in the files themselves.


using Coverage

coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
coverage = merge_coverage_counts(coverage, filter!(
    let prefixes = (joinpath(pwd(), "src", ""),)
        c -> any(p -> startswith(c.filename, p), prefixes)
    end,
    LCOV.readfolder("test")))

# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)
println("$(100covered_lines/total_lines)% coverage")
println("Use lcov.info for full coverage info")

LCOV.writefile("lcov.info", coverage)
