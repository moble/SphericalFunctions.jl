# Tests that every docstring in `src` documents something.
#
# A docstring documents only the expression that immediately follows its closing quotes.  With
# anything in between — a comment, a blank line — the parser leaves the string as a statement
# of its own, which evaluates to itself and is discarded, so that its text is registered
# nowhere; and `@doc raw"""…"""` followed by a comment becomes an `@doc` call with a single
# argument, which looks the text up rather than documenting anything with it.  Neither is an
# error, and Documenter's `checkdocs` does not notice either, because the binding the text was
# written for usually has another docstring of its own.  Both are visible in the parsed source,
# however, which is what these items examine.

@testitem "Docstrings: every docstring in src is attached to a definition" begin
    import SphericalFunctions

    # The statements that are evaluated when the file is, walked through the constructs that
    # hold further such statements: modules, `begin` and `let` blocks, loops, conditionals,
    # and `@static`.  Function bodies and quoted code are not entered, because a string there
    # is an ordinary value.
    is_doc_macro(x) = x === Symbol("@doc") || x == GlobalRef(Core, Symbol("@doc")) ||
        (x isa Expr && x.head === :. && x.args[end] == QuoteNode(Symbol("@doc")))
    is_string_literal(x) = x isa AbstractString ||
        (x isa Expr && x.head === :string) ||
        (x isa Expr && x.head === :macrocall && x.args[1] isa Symbol &&
            endswith(String(x.args[1]), "_str"))
    function detached_docstrings!(found, ex, file, line)
        if ex isa Expr && ex.head ∈ (:toplevel, :block)
            for a ∈ ex.args
                if a isa LineNumberNode
                    line = a.line
                else
                    detached_docstrings!(found, a, file, line)
                end
            end
        elseif ex isa Expr && ex.head === :module
            detached_docstrings!(found, ex.args[3], file, line)
        elseif ex isa Expr && ex.head ∈ (:let, :for, :while)
            detached_docstrings!(found, ex.args[2], file, line)
        elseif ex isa Expr && ex.head ∈ (:if, :elseif)
            for a ∈ ex.args[2:end]
                detached_docstrings!(found, a, file, line)
            end
        elseif ex isa Expr && ex.head === :macrocall
            arguments = filter(a -> !(a isa LineNumberNode), ex.args[2:end])
            if is_doc_macro(ex.args[1])
                length(arguments) < 2 && push!(found, (file, line, "@doc documents nothing"))
            elseif ex.args[1] === Symbol("@static")
                for a ∈ arguments
                    detached_docstrings!(found, a, file, line)
                end
            end
        elseif is_string_literal(ex)
            push!(found, (file, line, "a string that documents nothing"))
        end
        found
    end
    detached_docstrings(source, file="none") =
        detached_docstrings!(Any[], Meta.parseall(source; filename=file), file, 0)

    # The check itself: what it must find ...
    q = "\"\"\""
    @test length(detached_docstrings("$q\nDocs.\n$q\n# a comment\nfunction f end\n")) == 1
    @test length(detached_docstrings("$q\nDocs.\n$q\n\nfunction f end\n")) == 1
    @test length(detached_docstrings("@doc raw$q\nDocs.\n$q\n# a comment\nf(x) = x\n")) == 1
    @test length(detached_docstrings("\"Docs \$x.\"\n# a comment\nf(x) = x\n")) == 1
    @test length(detached_docstrings("module M\n$q\nDocs.\n$q\n# c\nf() = 1\nend\n")) == 1
    @test length(detached_docstrings("begin\n$q\nDocs.\n$q\n# c\nf() = 1\nend\n")) == 1
    # Inside a conditional, a `let` block or a loop the parser attaches no docstring at all,
    # so that there only an explicit `@doc` documents anything
    @test length(detached_docstrings("@static if true\n$q\nDocs.\n$q\nf() = 1\nend\n")) == 1
    @test length(detached_docstrings("let\n$q\nDocs.\n$q\nf() = 1\nend\n")) == 1
    @test length(detached_docstrings("for T ∈ (1,)\n$q\nDocs.\n$q\nf() = 1\nend\n")) == 1
    # ... and what it must not
    @test isempty(detached_docstrings("$q\nDocs.\n$q\nfunction f end\n"))
    @test isempty(detached_docstrings("@doc raw$q\nDocs.\n$q\nf(x) = x\n"))
    @test isempty(detached_docstrings("# a comment\n$q\nDocs.\n$q\nf(x) = x\n"))
    @test isempty(detached_docstrings("module M\n$q\nDocs.\n$q\nf() = 1\nend\n"))
    @test isempty(detached_docstrings("begin\n$q\nDocs.\n$q\nf() = 1\nend\n"))
    @test isempty(detached_docstrings("@static if true\n@doc raw$q\nDocs.\n$q\nf() = 1\nend\n"))
    @test isempty(detached_docstrings("f(x) = \"a value\"\nconst s = \"another\"\n"))

    # The package's source
    root = joinpath(pkgdir(SphericalFunctions), "src")
    found = Any[]
    for (directory, _, files) ∈ walkdir(root), file ∈ sort(files)
        endswith(file, ".jl") || continue
        path = joinpath(directory, file)
        detached_docstrings!(
            found, Meta.parseall(read(path, String); filename=path), relpath(path, root), 0
        )
    end
    @test isempty(found)
end

@testitem "Docstrings: the SSHT constructor and w[ℓ, :] are documented" begin
    import SphericalFunctions
    import SphericalFunctions: SSHT

    # `?SSHT` shows the constructor's docstring beside the abstract type's, and the
    # constructor's is the one that describes its keywords.  The docstrings are read from
    # the package's metadata, because `@doc SSHT` returns only the abstract type's docstring
    # when the REPL is not loaded, as under `Pkg.test` from Julia 1.12 on.
    SSHT_docs = Base.Docs.meta(SphericalFunctions)[Base.Docs.Binding(SphericalFunctions, :SSHT)]
    @test any(d -> occursin(r"method\s*=", join(string.(d.text))), values(SSHT_docs.docs))

    # The block accessor of `ModeWeights` is one of the `getindex` methods documented by the
    # package, each under a header naming its call
    getindex_docs = Base.Docs.meta(SphericalFunctions)[Base.Docs.Binding(Base, :getindex)]
    headers = [
        first(split(strip(join(string.(d.text))), '\n')) for d ∈ values(getindex_docs.docs)
    ]
    @test "w[ℓ, :]" ∈ headers
    @test "w[ℓ, m]" ∈ headers
end
