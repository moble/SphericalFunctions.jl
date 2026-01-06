@testitem "WignerHrange" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    r1(mp_max, ell_max) = hcat([
        [ell, mp, m] for ell in 0:ell_max
        for mp in -min(ell, mp_max):min(ell, mp_max)
        for m in abs(mp):ell
    ]...)'
    for ell_max in 0:ell_max
        a = Deprecated.WignerHrange(ell_max)  # Implicitly, mp_max=ell_max
        b = r1(ell_max, ell_max)
        @test a == b
        for mp_max in 0:ell_max
            a = Deprecated.WignerHrange(ell_max, mp_max)
            b = r1(mp_max, ell_max)
            @test a == b
        end
    end
end

@testitem "WignerHsize" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    for i in -5:-1
        @test Deprecated.WignerHsize(i) == 0
        for j in -5:5
            @test Deprecated.WignerHsize(i, j) == 0
        end
    end
    for ell_max in 0:ell_max
        @test Deprecated.WignerHsize(ell_max) == size(Deprecated.WignerHrange(ell_max, ell_max), 1)
        for mp_max in ell_max
            @test Deprecated.WignerHsize(ell_max, mp_max) == size(Deprecated.WignerHrange(ell_max, mp_max), 1)
        end
    end
end

@testitem "WignerHindex" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    ell_max_slow = ell_max ÷ 2
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
        r = Deprecated.WignerHrange(ell_max_i)
        for ell in 0:ell_max_i
            for mp in -ell:ell
                for m in -ell:ell
                    i = Deprecated.WignerHindex(ell, mp, m)
                    @test r[i, :] == fold_H_indices(ell, mp, m)
                end
            end
        end
        for mp_max in 0:ell_max_i
            r = Deprecated.WignerHrange(ell_max_i, mp_max)
            for ell in 0:ell_max_i
                for mp in -ell:ell
                    for m in -ell:ell
                        if abs(mp) > mp_max && abs(m) > mp_max
                            continue
                        end
                        i = Deprecated.WignerHindex(ell, mp, m, mp_max)
                        @test r[i, :] == fold_H_indices(ell, mp, m)
                    end
                end
            end
        end
    end
end

@testitem "WignerDrange" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    ell_max_slow = ell_max ÷ 2
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
            a = Deprecated.WignerDrange(ell_max)  # Implicitly, mp_max=ell_max
            b = r2(ell_min, ell_max, ell_max)
            @test a == b
            for mp_max in 0:ell_max
                a = Deprecated.WignerDrange(ell_max, mp_max)
                b = r2(ell_min, mp_max, ell_max)
                @test a == b
            end
        end
    end
end

@testitem "WignerDsize" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    for ell_max in 0:ell_max
        let ell_min = 0
            for mp_max in 0:ell_max
                a = Deprecated.WignerDsize(ell_max, mp_max)
                b = size(Deprecated.WignerDrange(ell_max, mp_max), 1)
                @test a == b
            end
        end
    end
    for ell_max in 0:ell_max
        let ell_min = 0
            for mp_max in 0:ell_max
                a = Deprecated.WignerDsize(ell_max, mp_max)
                b = size(Deprecated.WignerDrange(ell_max, mp_max), 1)
                @test a == b
            end
        end
    end
    for ell_max in 0:ell_max
        let ell_min = 0
            for mp_max in [ell_max]
                a = Deprecated.WignerDsize(ell_max, mp_max)
                # a = WignerDsize_ellmin(ell_min, ell_max)
                b = size(Deprecated.WignerDrange(ell_max, mp_max), 1)
                @test a == b
            end
        end
    end
    for ell_max in 0:ell_max
        let ell_min = 0
            for mp_max in [ell_max]
                a = Deprecated.WignerDsize(ell_max, mp_max)
                # a = WignerDsize(ell_max)
                b = size(Deprecated.WignerDrange(ell_max, mp_max), 1)
                @test a == b
            end
        end
    end
end

@testitem "WignerDindex" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    ell_max_slow = ell_max ÷ 2
    for ellmax in 0:ell_max_slow
        r = Deprecated.WignerDrange(ellmax)
        for ell in 0:ellmax
            for mp in -ell:ell
                for m in -ell:ell
                    i = Deprecated.WignerDindex(ell, mp, m)
                    @test r[i, :] == [ell, mp, m]
                end
            end
        end
        let ell_min = 0
            r = Deprecated.WignerDrange(ellmax)
            for ell in ell_min:ellmax
                for mp in -ell:ell
                    for m in -ell:ell
                        i = Deprecated.WignerDindex(ell, mp, m)
                        @test r[i, :] == [ell, mp, m]
                    end
                end
            end
            for mp_max in 0:ellmax
                r = Deprecated.WignerDrange(ellmax, mp_max)
                for ell in ell_min:ellmax
                    for mp in -min(ell, mp_max):min(ell, mp_max)
                        for m in -ell:ell
                            i = Deprecated.WignerDindex(ell, mp, m, mp_max)
                            @test r[i, :] == [ell, mp, m]
                        end
                    end
                end
            end
        end
    end
end

@testitem "Yrange" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
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
            a = Deprecated.Yrange(ell_min, ell_max)
            b = r3(ell_min, ell_max)
            @test a == b
        end
    end
end

@testitem "Ysize" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    for ell_max in 0:ell_max
        for ell_min in 0:ell_max
            a = Deprecated.Ysize(ell_min, ell_max)
            b = size(Deprecated.Yrange(ell_min, ell_max))[1]
            @test a == b
        end
    end
end

@testitem "deduce_limits" begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    for ℓmax in 0:4096
        for ℓmin in 0:min(2, ℓmax)
            deduced = Deprecated.deduce_limits(Deprecated.Ysize(ℓmin, ℓmax), nothing)
            @test deduced == (ℓmin, ℓmax)
        end
    end
    for ℓmax in 0:ell_max
        for ℓmin in [0]
            deduced = Deprecated.deduce_limits(Deprecated.Ysize(ℓmin, ℓmax))
            @test deduced == (ℓmin, ℓmax)
        end
    end
    for ℓmax in 0:ell_max
        for ℓmin in 0:ℓmax
            deduced = Deprecated.deduce_limits(Deprecated.Ysize(ℓmin, ℓmax), ℓmin)
            @test deduced == (ℓmin, ℓmax)
        end
    end
    for ℓmax in 1:ell_max
        prev_size = Deprecated.Ysize(ℓmax-1)
        this_size = Deprecated.Ysize(ℓmax)
        if abs(prev_size - this_size) > 1
            mid_size = (prev_size + this_size) ÷ 2
            @test_throws ErrorException Deprecated.deduce_limits(mid_size)
        end
    end
end

@testitem "Yindex" setup=[Utilities] begin
    import SphericalFunctions: Deprecated
    ell_max = 16
    for ell_max in 0:ell_max
        let ell_min = 0
            r = Deprecated.Yrange(ell_min, ell_max)
            for ell in ell_min:ell_max
                for m in -ell:ell
                    i = Deprecated.Yindex(ell, m)
                    @test r[i, :] == [ell, m]
                end
            end
            s = Deprecated.Yrange(ell_max)
            @test array_equal(r, s)
        end
        for ell_min in 0:ell_max
            r = Deprecated.Yrange(ell_min, ell_max)
            for ell in ell_min:ell_max
                for m in -ell:ell
                    i = Deprecated.Yindex(ell, m, ell_min)
                    @test r[i, :] == [ell, m]
                end
            end
        end
    end
end

@testitem "phi_theta <-> theta_phi" setup=[Utilities] begin
    for T in [Float64, Float32, BigFloat]
        import SphericalFunctions.Deprecated: theta_phi, phi_theta
        ell_max = 16
        for nθ in 2:(2ell_max+1)
            for nϕ in 2:(2ell_max+1)
                θϕ = theta_phi(nθ, nϕ, T)
                ϕθ = permutedims(phi_theta(nϕ, nθ, T), [2, 1, 3])[:, :, 2:-1:1]
                @test array_equal(θϕ, ϕθ)
            end
        end
    end
end
