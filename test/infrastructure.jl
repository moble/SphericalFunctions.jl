# Tests of the files that run the test suite, rather than of the package.

@testitem "runtests.jl parses" begin
    import SphericalFunctions

    # `test/runtests.jl` is executed only by `Pkg.test` — which is what CI runs, through
    # julia-runtest — and the test-item runners never read it, so a statement in it that does
    # not parse would not show up in an ordinary run of the suite.  Under `Pkg.test` the file
    # is parsed one top-level statement at a time, so every test runs first, and the whole
    # run then fails with a `LoadError`.
    path = joinpath(pkgdir(SphericalFunctions), "test", "runtests.jl")
    parsed = Meta.parseall(read(path, String); filename=path)
    unparsed = Any[]
    function find_unparsed!(unparsed, ex)
        if ex isa Expr
            ex.head ∈ (:error, :incomplete) && push!(unparsed, ex)
            foreach(a -> find_unparsed!(unparsed, a), ex.args)
        end
        unparsed
    end
    @test isempty(find_unparsed!(unparsed, parsed))
end
