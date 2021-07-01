struct PlanMap2Salm{T<:Real}
    spin
    ℓmax
    ℓmin
    Nφ
    Nϑ
    Nextra
    G
    m′max
    wigner
    weight
    expiθ
    ϵs
    extra_dims
end

macro unpackPlan(q)  # Stolen from https://stackoverflow.com/a/67138749/1194883
    code =  Expr(:block, [ :($field = $q.$field) for field in fieldnames(PlanMap2Salm) ]...)
    esc(code)
end

@doc raw"""
    map2salm(map, spin, ℓmax, [ℓmin])
    map2salm(map, plan)

Transform `map` values sampled on the sphere to ``{}_sa_{\ell, m}`` modes.

The `map` array should have size Nφ along its first dimension and Nϑ along its
second; any number of dimensions may follow.  The `spin` must be entered
explicitly, and `ℓmax` is the highest ℓ value you want in the output.  The
`ℓmin` represents the smallest ℓ value in the output and defaults to
`abs(spin)`; it generally should never be larger than this.

For repeated applications of this function with different values of `map`, it
is more efficient to pre-compute `plan` using [`plan_map2salm`](@ref).  These
functions will create a new `salm` array on each call.  To operate in place on
a pre-allocated `salm` array, use [`map2salm!`](@ref).

"""
function map2salm(
    map::AbstractArray{Complex{T}},
    spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)
) where {T<:Real}
    Nφ, Nϑ, Nextra... = size(map)
    salm = zeros(complex(T), (Ysize(ℓmin, ℓmax), Nextra...))
    map2salm!(salm, map, spin, ℓmax, ℓmin)
    return salm
end

@doc raw"""
    map2salm!(salm, map, spin, ℓmax, [ℓmin])
    map2salm!(salm, map, plan)

Transform `map` values sampled on the sphere to ``{}_sa_{\ell, m}`` modes in
place.

For details, see [`map2salm`](@ref).

"""
function map2salm!(
    salm::AbstractArray{Complex{T}},
    map::AbstractArray{Complex{T}},
    spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)
) where {T<:Real}
    plan = plan_map2salm(map, spin, ℓmax, ℓmin)
    map2salm!(salm, map, plan)
end

"""
    plan_map2salm(map, spin, ℓmax, [ℓmin])

Precompute values to use in executing [`map2salm`](@ref) or
[`map2salm!`](@ref).

The arguments to this function exactly mirror those of the first form
[`map2salm`](@ref).

"""
function plan_map2salm(map::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
    Nφ, Nϑ, Nextra... = size(map)
    G = Array{complex(T)}(undef, (Nφ, Nϑ, Nextra...));
    m′max = abs(spin)
    wigner = WignerMatrixCalculator(ℓmin, ℓmax, m′max, T)
    weight = clenshaw_curtis(Nϑ, T)
    expiθ = complex_powers(exp(im * (π / T(Nϑ-1))), Nϑ-1)
    ϵs = Spherical.WignerMatrices.ϵ(-spin)
    extra_dims = Base.Iterators.product((1:e for e in Nextra)...)

    return (spin, ℓmax, ℓmin, Nφ, Nϑ, Nextra, G, m′max, wigner, weight, expiθ, ϵs, extra_dims)
end

function map2salm!(
    salm::AbstractArray{Complex{T}},
    map::AbstractArray{Complex{T}},
    (spin, ℓmax, ℓmin, Nφ, Nϑ, Nextra, G, m′max, wigner, weight, expiθ, ϵs, extra_dims)
) where {T<:Real}
    @assert size(salm) == (Ysize(ℓmin, ℓmax), Nextra...)

    fftplan = plan_fft(map[:, 1, first(extra_dims)...])
    @inbounds for extra ∈ extra_dims
        for ϑ ∈ 1:Nϑ
            @views mul!(G[:, ϑ, extra...], fftplan, map[:, ϑ, extra...])
            @views G[:, ϑ, extra...] *= weight[ϑ]
        end
    end

    @inbounds for ϑ ∈ 1:Nϑ
        H!(wigner, expiθ[ϑ])  # Not thread safe
        for extra ∈ extra_dims
            for ℓ ∈ abs(spin):ℓmax
                sqrt_factor = √((2ℓ+1)*T(π)) / (2ℓmax+1)

                let m=0
                    λ_factor = ϵs * sqrt_factor

                    salm[Yindex(ℓ, m, ℓmin), extra...] +=
                        G[m+1, ϑ, extra...] * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, -m, m′max)]
                end

                for m ∈ 1:ℓ
                    λ_factor = ϵs * sqrt_factor

                    salm[Yindex(ℓ, m, ℓmin), extra...] +=
                        G[m+1, ϑ, extra...] * ϵ(m) * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, -m, m′max)]

                    salm[Yindex(ℓ, -m, ℓmin), extra...] +=
                        G[Nφ-m+1, ϑ, extra...] * ϵ(-m) * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, m, m′max)]
                end

            end
        end
    end
end

function map2salm(map::AbstractArray{Complex{T}}, plan) where {T<:Real}
    salm = zeros(complex(T), (Ysize(plan[3], plan[2]), plan[6]...))
    map2salm!(salm, map, plan)
    return salm
end


# # function map2salm(map::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
# #     goodies = map2salm_goodies(map, spin, ℓmax, ℓmin)
# #     _map2salm(goodies...)
# # end

# function map2salm_goodies(map::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, ℓmin::Int=abs(spin)) where {T<:Real}
#     Nφ, Nϑ, Nextra... = size(map)
#     G = Array{complex(T)}(undef, (Nφ, Nϑ, Nextra...));
#     modes = zeros(complex(T), (Ysize(ℓmin, ℓmax), Nextra...))
#     m′max = abs(spin)
#     wigner = WignerMatrixCalculator(ℓmin, ℓmax, m′max, T)
#     weight = clenshaw_curtis(Nϑ, T)
#     expiθ = complex_powers(exp(im * (π / T(Nϑ-1))), Nϑ-1)
#     ϵs = Spherical.WignerMatrices.ϵ(-spin)
#     extra_dims = Base.Iterators.product((1:e for e in Nextra)...)

#     (T, map, spin, ℓmax, ℓmin, Nφ, Nϑ, Nextra, G, modes, m′max, wigner, weight, expiθ, ϵs, extra_dims)
# end

# function _map2salm(T, map, spin, ℓmax, ℓmin, Nφ, Nϑ, Nextra, G, modes, m′max, wigner, weight, expiθ, ϵs, extra_dims)
#     plan = plan_fft(map[:, 1, first(extra_dims)...])
#     @inbounds for extra ∈ extra_dims
#         for ϑ ∈ 1:Nϑ
#             @views mul!(G[:, ϑ, extra...], plan, map[:, ϑ, extra...])
#             @views G[:, ϑ, extra...] *= weight[ϑ]
#         end
#     end

#     @inbounds for ϑ ∈ 1:Nϑ
#         H!(wigner, expiθ[ϑ])  # Not thread safe
#         for extra ∈ extra_dims
#             for ℓ ∈ abs(spin):ℓmax
#                 sqrt_factor = √((2ℓ+1)*T(π)) / (2ℓmax+1)

#                 let m=0
#                     λ_factor = ϵs * sqrt_factor

#                     modes[Yindex(ℓ, m, ℓmin), extra...] +=
#                         G[m+1, ϑ, extra...] * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, -m, m′max)]
#                 end

#                 for m ∈ 1:ℓ
#                     λ_factor = ϵs * sqrt_factor

#                     modes[Yindex(ℓ, m, ℓmin), extra...] +=
#                         G[m+1, ϑ, extra...] * ϵ(m) * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, -m, m′max)]

#                     modes[Yindex(ℓ, -m, ℓmin), extra...] +=
#                         G[Nφ-m+1, ϑ, extra...] * ϵ(-m) * λ_factor * wigner.Hwedge[WignerHindex(ℓ, spin, m, m′max)]
#                 end

#             end
#         end
#     end

#     return modes
# end
