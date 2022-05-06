@testset verbose=true "indexing" begin
    import SphericalFunctions: WignerHsize, WignerHrange, WignerHindex
    import SphericalFunctions: WignerDsize, WignerDrange, WignerDindex
    import SphericalFunctions: Ysize, Yrange, Yindex, deduce_limits
    import SphericalFunctions: theta_phi

    ell_max = 16
    ell_max_slow = ell_max ÷ 2

    @testset "WignerHrange" begin
        r1(mp_max, ell_max) = hcat([
            [ell, mp, m] for ell in 0:ell_max
            for mp in -min(ell, mp_max):min(ell, mp_max)
            for m in abs(mp):ell
        ]...)'
        for ell_max in 0:ell_max
            a = WignerHrange(ell_max)  # Implicitly, mp_max=ell_max
            b = r1(ell_max, ell_max)
            @test a == b
            for mp_max in 0:ell_max
                a = WignerHrange(ell_max, mp_max)
                b = r1(mp_max, ell_max)
                @test a == b
            end
        end
    end

    @testset "WignerHsize" begin
        for i in -5:-1
            @test WignerHsize(i) == 0
            for j in -5:5
                @test WignerHsize(i, j) == 0
            end
        end
        for ell_max in 0:ell_max
            @test WignerHsize(ell_max) == size(WignerHrange(ell_max, ell_max), 1)
            for mp_max in ell_max
                @test WignerHsize(ell_max, mp_max) == size(WignerHrange(ell_max, mp_max), 1)
            end
        end
    end

    @testset "WignerHindex" begin
        function fold_H_indices(ell, mp, m)
            if m < -mp
                if m < mp
                    return [ell, -mp, -m]
                else
                    return [ell, -m, -mp]
                end
            else
                if m < mp
                    return [ell, m, mp]
                else
                    return [ell, mp, m]
                end
            end
        end

        for ell_max_i in 0:ell_max_slow
            r = WignerHrange(ell_max_i)
            for ell in 0:ell_max_i
                for mp in -ell:ell
                    for m in -ell:ell
                        i = WignerHindex(ell, mp, m)
                        @test r[i, :] == fold_H_indices(ell, mp, m)
                    end
                end
            end
            for mp_max in 0:ell_max_i
                r = WignerHrange(ell_max_i, mp_max)
                for ell in 0:ell_max_i
                    for mp in -ell:ell
                        for m in -ell:ell
                            if abs(mp) > mp_max && abs(m) > mp_max
                                continue
                            end
                            i = WignerHindex(ell, mp, m, mp_max)
                            @test r[i, :] == fold_H_indices(ell, mp, m)
                        end
                    end
                end
            end
        end
    end

    @testset "WignerDrange" begin
        function r2(ell_min, mp_max, ell_max)
            a = [
                [ℓ, mp, m]
                for ℓ in ell_min:ell_max
                for mp in -min(ℓ, mp_max):min(ℓ, mp_max)
                for m in -ℓ:ℓ
            ]
            collect(transpose(reshape(collect(Iterators.flatten(a)), (3, length(a)))))
        end

        for ell_max in 0:ell_max_slow÷2
            let ell_min = 0
                a = WignerDrange(ell_max)  # Implicitly, mp_max=ell_max
                b = r2(ell_min, ell_max, ell_max)
                @test a == b
                for mp_max in 0:ell_max
                    a = WignerDrange(ell_max, mp_max)
                    b = r2(ell_min, mp_max, ell_max)
                    @test a == b
                end
            end
        end
    end

    @testset "WignerDsize" begin
        for ell_max in 0:ell_max
            let ell_min = 0
                for mp_max in 0:ell_max
                    a = WignerDsize(ell_max, mp_max)
                    b = size(WignerDrange(ell_max, mp_max), 1)
                    @test a == b
                end
            end
        end
        for ell_max in 0:ell_max
            let ell_min = 0
                for mp_max in 0:ell_max
                    a = WignerDsize(ell_max, mp_max)
                    b = size(WignerDrange(ell_max, mp_max), 1)
                    @test a == b
                end
            end
        end
        for ell_max in 0:ell_max
            let ell_min = 0
                for mp_max in [ell_max]
                    a = WignerDsize(ell_max, mp_max)
                    # a = WignerDsize_ellmin(ell_min, ell_max)
                    b = size(WignerDrange(ell_max, mp_max), 1)
                    @test a == b
                end
            end
        end
        for ell_max in 0:ell_max
            let ell_min = 0
                for mp_max in [ell_max]
                    a = WignerDsize(ell_max, mp_max)
                    # a = WignerDsize(ell_max)
                    b = size(WignerDrange(ell_max, mp_max), 1)
                    @test a == b
                end
            end
        end
    end

    @testset "WignerDindex" begin
        for ellmax in 0:ell_max_slow
            r = WignerDrange(ellmax)
            for ell in 0:ellmax
                for mp in -ell:ell
                    for m in -ell:ell
                        i = WignerDindex(ell, mp, m)
                        @test r[i, :] == [ell, mp, m]
                    end
                end
            end
            let ell_min = 0
                r = WignerDrange(ellmax)
                for ell in ell_min:ellmax
                    for mp in -ell:ell
                        for m in -ell:ell
                            i = WignerDindex(ell, mp, m)
                            @test r[i, :] == [ell, mp, m]
                        end
                    end
                end
                for mp_max in 0:ellmax
                    r = WignerDrange(ellmax, mp_max)
                    for ell in ell_min:ellmax
                        for mp in -min(ell, mp_max):min(ell, mp_max)
                            for m in -ell:ell
                                i = WignerDindex(ell, mp, m, mp_max)
                                @test r[i, :] == [ell, mp, m]
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Yrange" begin
        function r3(ell_min, ell_max)
            a = [
                [ℓ, m]
                for ℓ in ell_min:ell_max
                for m in -ℓ:ℓ
            ]
            collect(transpose(reshape(collect(Iterators.flatten(a)), (2, length(a)))))
        end
        for ell_max in 0:ell_max
            for ell_min in 0:ell_max
                a = Yrange(ell_min, ell_max)
                b = r3(ell_min, ell_max)
                @test a == b
            end
        end
    end

    @testset "Ysize" begin
        for ell_max in 0:ell_max
            for ell_min in 0:ell_max
                a = Ysize(ell_min, ell_max)
                b = size(Yrange(ell_min, ell_max))[1]
                @test a == b
            end
        end
    end

    @testset "deduce_limits" begin
        for ℓmax in 0:4096
            for ℓmin in 0:min(2, ℓmax)
                deduced = deduce_limits(Ysize(ℓmin, ℓmax), nothing)
                @test deduced == (ℓmin, ℓmax)
            end
        end
        for ℓmax in 0:ell_max
            for ℓmin in [0]
                deduced = deduce_limits(Ysize(ℓmin, ℓmax))
                @test deduced == (ℓmin, ℓmax)
            end
        end
        for ℓmax in 0:ell_max
            for ℓmin in 0:ℓmax
                deduced = deduce_limits(Ysize(ℓmin, ℓmax), ℓmin)
                @test deduced == (ℓmin, ℓmax)
            end
        end
        for ℓmax in 1:ell_max
            prev_size = Ysize(ℓmax-1)
            this_size = Ysize(ℓmax)
            if abs(prev_size - this_size) > 1
                mid_size = (prev_size + this_size) ÷ 2
                @test_throws ErrorException deduce_limits(mid_size)
            end
        end
    end

    @testset "Yindex" begin
        for ell_max in 0:ell_max
            let ell_min = 0
                r = Yrange(ell_min, ell_max)
                for ell in ell_min:ell_max
                    for m in -ell:ell
                        i = Yindex(ell, m)
                        @test r[i, :] == [ell, m]
                    end
                end
                s = Yrange(ell_max)
                @test array_equal(r, s)
            end
            for ell_min in 0:ell_max
                r = Yrange(ell_min, ell_max)
                for ell in ell_min:ell_max
                    for m in -ell:ell
                        i = Yindex(ell, m, ell_min)
                        @test r[i, :] == [ell, m]
                    end
                end
            end
        end
    end

    @testset "phi_theta <-> theta_phi" for T in [Float64, Float32, BigFloat]
        for nθ in 2:(2ell_max+1)
            for nϕ in 2:(2ell_max+1)
                θϕ = theta_phi(nθ, nϕ, T)
                ϕθ = permutedims(phi_theta(nϕ, nθ, T), [2, 1, 3])[:, :, 2:-1:1]
                @test array_equal(θϕ, ϕθ)
            end
        end
    end

end
