# struct PlanMap2Salm
#     T
#     map_values
#     spin
#     ℓmax
#     ℓmin
#     Nφ
#     Nϑ
#     Nextra
#     G
#     modes
#     m′max
#     wigner
#     weight
#     expiθ
#     ϵs
#     extras
# end

# macro unpackPlan(q)  # Stolen from https://stackoverflow.com/a/67138749/1194883
#     code =  Expr(:block, [ :($field = $q.$field) for field in fieldnames(PlanMap2Salm) ]...)
#     esc(code)
# end

# map2salm(map_values::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
# plan_map2salm(map_values::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
# map2salm!(modes::AbstractArray{Complex{T}}, map_values::AbstractArray{Complex{T}}, plan::PlanMap2Salm) where {T<:Real}
# map2salm(map_values::AbstractArray{Complex{T}}, plan::PlanMap2Salm) where {T<:Real}
# _map2salm!(modes::AbstractArray{Complex{T}}, map_values::AbstractArray{Complex{T}}, plan::PlanMap2Salm) where {T<:Real}



@doc raw"""
    map2salm(map_values, spin, ℓmax, [ℓmin])

Transform values sampled on the sphere to ``{}_sa_{\ell, m}`` modes.

The `map_values` array should have size Nφ along its first dimension and Nϑ
along its second; any number of dimensions may follow.

"""
function map2salm(map_values::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
    goodies = map2salm_goodies(map_values, spin, ℓmax, ℓmin)
    _map2salm(goodies...)
end

function map2salm_goodies(map_values::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
    Nφ, Nϑ, Nextra... = size(map_values)
    G = Array{complex(T)}(undef, (Nφ, Nϑ, Nextra...));
    modes = zeros(complex(T), (Ysize(ℓmin, ℓmax), Nextra...))
    m′max = abs(spin)
    wigner = WignerMatrixCalculator(ℓmin, ℓmax, m′max, T)
    weight = clenshaw_curtis(Nϑ, T)
    expiθ = complex_powers(exp(im * (π / T(Nϑ-1))), Nϑ-1)
    ϵs = Spherical.WignerMatrices.ϵ(-spin)
    extras = Base.Iterators.product((1:e for e in Nextra)...)

    (T, map_values, spin, ℓmax, ℓmin, Nφ, Nϑ, Nextra, G, modes, m′max, wigner, weight, expiθ, ϵs, extras)
end

function _map2salm(T, map_values, spin, ℓmax, ℓmin, Nφ, Nϑ, Nextra, G, modes, m′max, wigner, weight, expiθ, ϵs, extras)
    # @inbounds for extra ∈ extras
    #     for ϑ ∈ 1:Nϑ
    #         @views G[:, ϑ, extra...] = weight[ϑ] * fft(map_values[:, ϑ, extra...])
    #     end
    # end

    # plan = plan_fft(map_values[:, 1, first(extras)...])
    # @inbounds for extra ∈ extras
    #     for ϑ ∈ 1:Nϑ
    #         @views G[:, ϑ, extra...] = weight[ϑ] * (plan * map_values[:, ϑ, extra...])
    #     end
    # end

    plan = plan_fft(map_values[:, 1, first(extras)...])
    @inbounds for extra ∈ extras
        for ϑ ∈ 1:Nϑ
            @views mul!(G[:, ϑ, extra...], plan, map_values[:, ϑ, extra...])
            @views G[:, ϑ, extra...] *= weight[ϑ]
        end
    end

    # G = fft(map_values)
    # @inbounds for extra ∈ extras
    #     for ϑ ∈ 1:Nϑ
    #         G[:, ϑ, extra...] *= weight[ϑ]
    #     end
    # end

    # G = fft(map_values)
    # @inbounds for ϑ ∈ 1:Nϑ
    #     w = weight[ϑ]
    #     for extra ∈ extras
    #         @views G[:, ϑ, extra...] *= w
    #     end
    # end

    @inbounds for ϑ ∈ 1:Nϑ
        H!(wigner, expiθ[ϑ])  # Not thread safe
        for extra ∈ extras
            for ℓ ∈ abs(spin):ℓmax
                sqrt_factor = √((2ℓ+1)*T(π)) / (2ℓmax+1)

                let m=0
                    λ_factor = ϵs * sqrt_factor

                    modes[Yindex(ℓ, m, ℓmin), extra...] +=
                        G[m+1, ϑ, extra...] * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, -m, m′max)]
                end

                for m ∈ 1:ℓ
                    λ_factor = ϵs * sqrt_factor

                    modes[Yindex(ℓ, m, ℓmin), extra...] +=
                        G[m+1, ϑ, extra...] * ϵ(m) * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, -m, m′max)]

                    modes[Yindex(ℓ, -m, ℓmin), extra...] +=
                        G[Nφ-m+1, ϑ, extra...] * ϵ(-m) * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, m, m′max)]
                end

            end
        end
    end

    return modes
end
