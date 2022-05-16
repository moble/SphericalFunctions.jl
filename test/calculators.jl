@testset verbose=true "calculators" begin

    @testset "io" begin
        io = IOBuffer()

        for T in [Float64, Float32, BigFloat]
            for nₘₐₓ in 0:12
                Base.show(io, MIME("text/plain"), DCalculator(T, nₘₐₓ))
                @test String(take!(io)) == "DCalculator($T, $nₘₐₓ)"
                for m′ₘₐₓ in 0:nₘₐₓ
                    Base.show(io, MIME("text/plain"), HCalculator(T, nₘₐₓ, m′ₘₐₓ))
                    @test String(take!(io)) == "HCalculator($T, $nₘₐₓ, $m′ₘₐₓ)"
                end
            end
        end
    end

    @testset "offsets (HCalculator)" begin
        nₘₐₓ = 8

        for n in 1:nₘₐₓ
            for m′ₘₐₓ in 0:n
                HC = HCalculator(Float64, n, m′ₘₐₓ)
                indices = [
                    (n, m′, m)
                    for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                    for m in abs(m′):n
                ]
                i = 1
                for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                    for m in abs(m′):n
                        @test offset(HC, n, m′, m) == i
                        if m ≥ abs(m′+1) && abs(m′+1) ≤ m′ₘₐₓ
                            @test (n, m′+1, m) == indices[i+m′offset₊(HC, n, m′, m)]
                        end
                        if m ≥ abs(m′-1) && abs(m′-1) ≤ m′ₘₐₓ
                            @test (n, m′-1, m) == indices[i-m′offset₋(HC, n, m′, m)]
                        end
                        i += 1
                    end
                end
            end
        end
    end

    @testset "offsets (DCalculator)" begin
        nₘₐₓ = 8
        DC = DCalculator(Float64, nₘₐₓ)

        for n in 1:nₘₐₓ
            indices = [
                (n, m′, m)
                for m′ in -n:n
                for m in -n:n
            ]
            i = 1
            for m′ in -n:n
                for m in -n:n
                    @test offset(DC, n, m′, m) == i
                    if abs(m′+1) ≤ n
                        @test (n, m′+1, m) == indices[i+m′offset₊(DC, n, m′, m)]
                    end
                    if abs(m′-1) ≤ n
                        @test (n, m′-1, m) == indices[i-m′offset₋(DC, n, m′, m)]
                    end
                    i += 1
                end
            end
        end
    end

end
