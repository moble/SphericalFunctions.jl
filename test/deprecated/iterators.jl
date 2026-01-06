@testitem "D iterators" begin
    import SphericalFunctions: Deprecated
    import Logging: with_logger, NullLogger

    ℓₘₐₓ = 20
    Drange = mapslices(
        ℓm′m -> tuple(ℓm′m...),
        Deprecated.WignerDrange(ℓₘₐₓ);
        dims=[2]
    )[:, 1]
    let
        𝔇 = with_logger(NullLogger()) do
            Deprecated.D_iterator(Drange, ℓₘₐₓ)
        end
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ -= 1
            Dˡ = [
                (ℓ, m′, m)
                for m in -ℓ:ℓ, m′ in -ℓ:ℓ
            ]
            @test 𝔇ˡ == Dˡ
        end
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ -= 1
            Dˡ = [
                (ℓ, m′, m)
                for m′ in -ℓ:ℓ, m in -ℓ:ℓ
            ]
            @test 𝔇ˡ == Dˡ broken=(ℓ ≠ 0)
        end
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ -= 1
            Dˡ = Matrix{Any}(undef, 2ℓ+1, 2ℓ+1)
            for m in -ℓ:ℓ
                for m′ in -ℓ:ℓ
                    Dˡ[m′+ℓ+1, m+ℓ+1] = (ℓ, m′, m)
                end
            end
            @test 𝔇ˡ == Dˡ broken=(ℓ ≠ 0)
        end
        @test Base.IteratorSize(typeof(𝔇)) == Base.HasShape{1}()
        @test Base.IteratorEltype(typeof(𝔇)) == Base.HasEltype()
        collection = collect(𝔇)
        @test length(collection) == length(𝔇)
        @test size(collection, 1) == size(𝔇, 1)
        @test size(collection) == size(𝔇)
        @test eltype(collection) == eltype(𝔇)
    end

    for ℓₘᵢₙ in 0:ℓₘₐₓ
        𝔇 = with_logger(NullLogger()) do
            Deprecated.D_iterator(Drange, ℓₘₐₓ, ℓₘᵢₙ)
        end
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ += ℓₘᵢₙ - 1
            Dˡ = [
                (ℓ, m′, m)
                for m in -ℓ:ℓ, m′ in -ℓ:ℓ
            ]
            @test 𝔇ˡ == Dˡ
        end
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ += ℓₘᵢₙ - 1
            Dˡ = [
                (ℓ, m′, m)
                for m in -ℓ:ℓ, m′ in -ℓ:ℓ
            ]
        end
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ += ℓₘᵢₙ - 1
            Dˡ = Matrix{Any}(undef, 2ℓ+1, 2ℓ+1)
            for m in -ℓ:ℓ
                for m′ in -ℓ:ℓ
                    Dˡ[m′+ℓ+1, m+ℓ+1] = (ℓ, m′, m)
                end
            end
            @test 𝔇ˡ == Dˡ broken=(ℓ ≠ 0)
        end
    end
end

@testitem "d iterators" begin
    import SphericalFunctions: Deprecated
    import Logging: with_logger, NullLogger

    ℓₘₐₓ = 20
    drange = mapslices(
        ℓm′m -> tuple(ℓm′m...),
        Deprecated.WignerDrange(ℓₘₐₓ);
        dims=[2]
    )[:, 1]
    let
        𝔡 = with_logger(NullLogger()) do
            Deprecated.d_iterator(drange, ℓₘₐₓ)
        end
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ -= 1
            dˡ = [
                (ℓ, m′, m)
                for m′ in -ℓ:ℓ, m in -ℓ:ℓ
            ]
            @test 𝔡ˡ == dˡ broken=(ℓ ≠ 0)
        end
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ -= 1
            dˡ = [
                (ℓ, m′, m)
                for m′ in -ℓ:ℓ
                for m in -ℓ:ℓ
            ]
            @test 𝔡ˡ == dˡ broken=true
        end
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ -= 1
            dˡ = Matrix{Any}(undef, 2ℓ+1, 2ℓ+1)
            for m′ in -ℓ:ℓ
                for m in -ℓ:ℓ
                    dˡ[m′+ℓ+1, m+ℓ+1] = (ℓ, m′, m)
                end
            end
            @test 𝔡ˡ == dˡ broken=(ℓ ≠ 0)
        end
        @test Base.IteratorSize(typeof(𝔡)) == Base.HasShape{1}()
        @test Base.IteratorEltype(typeof(𝔡)) == Base.HasEltype()
        collection = collect(𝔡)
        @test length(collection) == length(𝔡)
        @test size(collection, 1) == size(𝔡, 1)
        @test size(collection) == size(𝔡)
        @test eltype(collection) == eltype(𝔡)
    end

    for ℓₘᵢₙ in 0:ℓₘₐₓ
        𝔡 = with_logger(NullLogger()) do
            Deprecated.d_iterator(drange, ℓₘₐₓ, ℓₘᵢₙ)
        end
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ += ℓₘᵢₙ - 1
            dˡ = [
                (ℓ, m′, m)
                for m′ in -ℓ:ℓ, m in -ℓ:ℓ
            ]
            @test 𝔡ˡ == dˡ broken=(ℓ ≠ 0)
        end
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ += ℓₘᵢₙ - 1
            dˡ = [
                (ℓ, m′, m)
                for m′ in -ℓ:ℓ
                for m in -ℓ:ℓ
            ]
            @test 𝔡ˡ == dˡ broken=true
        end
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ += ℓₘᵢₙ - 1
            dˡ = Matrix{Any}(undef, 2ℓ+1, 2ℓ+1)
            for m′ in -ℓ:ℓ
                for m in -ℓ:ℓ
                    dˡ[m′+ℓ+1, m+ℓ+1] = (ℓ, m′, m)
                end
            end
            @test 𝔡ˡ == dˡ broken=(ℓ ≠ 0)
        end
    end
end

@testitem "Y iterators" begin
    import SphericalFunctions: Deprecated
    ℓₘₐₓ = 20
    Yrange = mapslices(
        ℓm -> tuple(ℓm...),
        Deprecated.Yrange(ℓₘₐₓ);
        dims=[2]
    )[:, 1]
    𝔜 = Deprecated.sYlm_iterator(Yrange, ℓₘₐₓ)
    for (ℓ, 𝔜ˡ) in enumerate(𝔜)
        ℓ -= 1
        Yˡ = [
            (ℓ, m)
            for m in -ℓ:ℓ
        ]
        @test 𝔜ˡ == Yˡ
    end
    @test Base.IteratorSize(typeof(𝔜)) == Base.HasShape{1}()
    @test Base.IteratorEltype(typeof(𝔜)) == Base.HasEltype()
    collection = collect(𝔜)
    @test length(collection) == length(𝔜)
    @test size(collection, 1) == size(𝔜, 1)
    @test size(collection) == size(𝔜)
    @test eltype(collection) == eltype(𝔜)

    for ℓₘᵢₙ in 0:ℓₘₐₓ
        local 𝔜 = Deprecated.sYlm_iterator(Yrange, ℓₘₐₓ, ℓₘᵢₙ)
        for (ℓ, 𝔜ˡ) in enumerate(𝔜)
            ℓ += ℓₘᵢₙ - 1
            Yˡ = [
                (ℓ, m)
                for m in -ℓ:ℓ
            ]
            @test 𝔜ˡ == Yˡ
        end
        iₘᵢₙ = Deprecated.Ysize(ℓₘᵢₙ-1)+1
        @test Yrange[iₘᵢₙ] == (ℓₘᵢₙ, -ℓₘᵢₙ)
        𝔜 = Deprecated.sYlm_iterator(Yrange[iₘᵢₙ:end], ℓₘₐₓ, ℓₘᵢₙ, 1)
        for (ℓ, 𝔜ˡ) in enumerate(𝔜)
            ℓ += ℓₘᵢₙ - 1
            Yˡ = [
                (ℓ, m)
                for m in -ℓ:ℓ
            ]
            @test 𝔜ˡ == Yˡ
        end
    end
end
