@testset verbose=true "indexing" begin
    ell_max = 16

    @testset "WignerHrange" begin
        r(mp_max, ell_max) = hcat([
            [ell, mp, m] for ell in 0:ell_max
                for mp in -min(ell, mp_max):min(ell, mp_max)
                    for m in abs(mp):ell
                        ]...)'
        for ell_max in 0:ell_max
            a = WignerHrange(ell_max)  # Implicitly, mp_max=ell_max
            b = r(ell_max, ell_max)
            @test a == b
            for mp_max in 0:ell_max
                a = WignerHrange(mp_max, ell_max)
                b = r(mp_max, ell_max)
                @test a == b
            end
        end
    end

    @testset "WignerHsize" begin
        for ell_max in 0:ell_max
            @test WignerHsize(ell_max) == size(WignerHrange(ell_max, ell_max))[1]
            for mp_max in ell_max
                @test WignerHsize(mp_max, ell_max) == size(WignerHrange(mp_max, ell_max))[1]
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

        ell_max_slow = ell_max ÷ 2
        for ell_max_i in 0:ell_max_slow
            r = WignerHrange(ell_max_i)
            for ell in 0:ell_max_i
                for mp in -ell:ell
                    for m in -ell:ell
                        i = WignerHindex(ell, mp, m)
                        # println(i)
                        # println(fold_H_indices(ell, mp, m))
                        # println(r)
                        # println()
                        @test r[i, :] == fold_H_indices(ell, mp, m)
                    end
                end
            end
            for mp_max in 0:ell_max_i
                r = WignerHrange(mp_max, ell_max_i)
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

end
