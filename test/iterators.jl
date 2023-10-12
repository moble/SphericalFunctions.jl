@testset verbose=true "Iterators" begin
    ℓₘₐₓ = 20

    @testset "D" begin
        Drange = mapslices(
            ℓm′m -> tuple(ℓm′m...),
            SphericalFunctions.WignerDrange(ℓₘₐₓ);
            dims=[2]
        )[:, 1]
        𝔇 = D_iterator(Drange, ℓₘₐₓ)
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ -= 1
            Dˡ = [
                (ℓ, m′, m)
                for m in -ℓ:ℓ, m′ in -ℓ:ℓ
            ]
            @test 𝔇ˡ == Dˡ
        end
        @test Base.IteratorSize(typeof(𝔇)) == Base.HasShape{1}()
        @test Base.IteratorEltype(typeof(𝔇)) == Base.HasEltype()
        collection = collect(𝔇)
        @test length(collection) == length(𝔇)
        @test size(collection, 1) == size(𝔇, 1)
        @test size(collection) == size(𝔇)
        @test eltype(collection) == eltype(𝔇)

        for ℓₘᵢₙ in 0:ℓₘₐₓ
            𝔇 = D_iterator(Drange, ℓₘₐₓ, ℓₘᵢₙ)
            for (ℓ, 𝔇ˡ) in enumerate(𝔇)
                ℓ += ℓₘᵢₙ - 1
                Dˡ = [
                    (ℓ, m′, m)
                    for m in -ℓ:ℓ, m′ in -ℓ:ℓ
                ]
                @test 𝔇ˡ == Dˡ
            end
        end
    end

    @testset "d" begin
        drange = mapslices(
            ℓm′m -> tuple(ℓm′m...),
            SphericalFunctions.WignerDrange(ℓₘₐₓ);
            dims=[2]
        )[:, 1]
        𝔡 = d_iterator(drange, ℓₘₐₓ)
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ -= 1
            dˡ = [
                (ℓ, m′, m)
                for m in -ℓ:ℓ, m′ in -ℓ:ℓ
            ]
            @test 𝔡ˡ == dˡ
        end
        @test Base.IteratorSize(typeof(𝔡)) == Base.HasShape{1}()
        @test Base.IteratorEltype(typeof(𝔡)) == Base.HasEltype()
        collection = collect(𝔡)
        @test length(collection) == length(𝔡)
        @test size(collection, 1) == size(𝔡, 1)
        @test size(collection) == size(𝔡)
        @test eltype(collection) == eltype(𝔡)

        for ℓₘᵢₙ in 0:ℓₘₐₓ
            𝔡 = d_iterator(drange, ℓₘₐₓ, ℓₘᵢₙ)
            for (ℓ, 𝔡ˡ) in enumerate(𝔡)
                ℓ += ℓₘᵢₙ - 1
                dˡ = [
                    (ℓ, m′, m)
                    for m in -ℓ:ℓ, m′ in -ℓ:ℓ
                ]
                @test 𝔡ˡ == dˡ
            end
        end
    end

    @testset "Y" begin
        Yrange = mapslices(
            ℓm -> tuple(ℓm...),
            SphericalFunctions.Yrange(ℓₘₐₓ);
            dims=[2]
        )[:, 1]
        𝔜 = sYlm_iterator(Yrange, ℓₘₐₓ)
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
            𝔜 = sYlm_iterator(Yrange, ℓₘₐₓ, ℓₘᵢₙ)
            for (ℓ, 𝔜ˡ) in enumerate(𝔜)
                ℓ += ℓₘᵢₙ - 1
                Yˡ = [
                    (ℓ, m)
                    for m in -ℓ:ℓ
                ]
                @test 𝔜ˡ == Yˡ
            end
            iₘᵢₙ = Ysize(ℓₘᵢₙ-1)+1
            @test Yrange[iₘᵢₙ] == (ℓₘᵢₙ, -ℓₘᵢₙ)
            𝔜 = sYlm_iterator(Yrange[iₘᵢₙ:end], ℓₘₐₓ, ℓₘᵢₙ, 1)
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

end
