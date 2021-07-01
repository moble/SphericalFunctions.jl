@doc raw"""
    map2salm(map_values, spin, ℓmax, [ℓmin])

Transform values sampled on the sphere to ``{}_sa_{\ell, m}`` modes.

The `map_values` array should have size Nφ along its first dimension and Nϑ
along its second; any number of dimensions may follow.

"""
function map2salm(map_values::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
    Nφ, Nϑ, Nextra... = size(map_values)
    m′max = abs(spin)
    wigner = WignerMatrixCalculator(ℓmin, ℓmax, m′max, T)
    weight = clenshaw_curtis(Nϑ, T)
    expiθ = complex_powers(exp(im * (π / T(Nϑ-1))), Nϑ-1)
    ϵs = Spherical.WignerMatrices.ϵ(-spin)

    G = Array{complex(T)}(undef, (Nφ, Nϑ, Nextra...));

    for extra = Base.Iterators.product((1:e for e in Nextra)...)
        for ϑ = 1:Nϑ
            G[:, ϑ, extra...] = weight[ϑ] * fft(map_values[:, ϑ, extra...])
        end
    end

    modes = zeros(complex(T), (Ysize(ℓmin, ℓmax), Nextra...))

    for ϑ = 1:Nϑ
        H!(wigner, expiθ[ϑ])  # Not thread safe
        for extra = Base.Iterators.product((1:e for e in Nextra)...)
            Base.Threads.@threads for ℓ ∈ abs(spin):ℓmax
                sqrt_factor = √((2ℓ+1)*T(π)) / (2ℓmax+1)

                let m=0
                    λ_factor = ϵs * sqrt_factor

                    i_H = Spherical.WignerMatrices.WignerHindex(ℓ, spin, -m, m′max)
                    sλlm = λ_factor * wigner.Hwedge[i_H]
                    modes[Yindex(ℓ, m, ℓmin), extra...] += G[m+1, ϑ, extra...] * sλlm
                end

                for m ∈ 1:ℓ
                    # λ_factor = m%2==0 ? ϵs * sqrt_factor : -ϵs * sqrt_factor
                    λ_factor = ϵs * sqrt_factor

                    i_H = Spherical.WignerMatrices.WignerHindex(ℓ, spin, -m, m′max)
                    sλlm = Spherical.WignerMatrices.ϵ(m) * λ_factor * wigner.Hwedge[i_H]
                    modes[Yindex(ℓ, m, ℓmin), extra...] += G[m+1, ϑ, extra...] * sλlm

                    i_H = Spherical.WignerMatrices.WignerHindex(ℓ, spin, m, m′max)
                    sλlm = Spherical.WignerMatrices.ϵ(-m) * λ_factor * wigner.Hwedge[i_H]
                    modes[Yindex(ℓ, -m, ℓmin), extra...] += G[Nφ-m+1, ϑ, extra...] * sλlm
                end

            end
        end
    end

    return modes
end
