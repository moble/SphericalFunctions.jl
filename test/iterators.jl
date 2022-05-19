@testset verbose=true "Iterators" begin
    ℓₘₐₓ = 20

    @testset "D" begin
        Drange = mapslices(
            ℓm′m -> tuple(ℓm′m...),
            SphericalFunctions.WignerDrange(ℓₘₐₓ);
            dims=[2]
        )[:, 1]
        𝔇 = Diterator(Drange, ℓₘₐₓ)
        for (ℓ, 𝔇ˡ) in enumerate(𝔇)
            ℓ -= 1
            Dˡ = [
                (ℓ, m′, m)
                for m in -ℓ:ℓ, m′ in -ℓ:ℓ
            ]
            @test 𝔇ˡ == Dˡ
        end
        for ℓₘᵢₙ in 0:ℓₘₐₓ
            𝔇 = Diterator(Drange, ℓₘₐₓ, ℓₘᵢₙ)
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
        𝔡 = diterator(drange, ℓₘₐₓ)
        for (ℓ, 𝔡ˡ) in enumerate(𝔡)
            ℓ -= 1
            dˡ = [
                (ℓ, m′, m)
                for m in -ℓ:ℓ, m′ in -ℓ:ℓ
            ]
            @test 𝔡ˡ == dˡ
        end
        for ℓₘᵢₙ in 0:ℓₘₐₓ
            𝔡 = diterator(drange, ℓₘₐₓ, ℓₘᵢₙ)
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
        𝔜 = Yiterator(Yrange, ℓₘₐₓ)
        for (ℓ, 𝔜ˡ) in enumerate(𝔜)
            ℓ -= 1
            Yˡ = [
                (ℓ, m)
                for m in -ℓ:ℓ
            ]
            @test 𝔜ˡ == Yˡ
        end
        for ℓₘᵢₙ in 0:ℓₘₐₓ
            𝔜 = Yiterator(Yrange, ℓₘₐₓ, ℓₘᵢₙ)
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
            𝔜 = Yiterator(Yrange[iₘᵢₙ:end], ℓₘₐₓ, ℓₘᵢₙ, 1)
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
