@doc raw"""
    map2salm(map, spin, ℓmax)
    map2salm(map, plan)

Transform `map` values sampled on the sphere to ``{}_sa_{\ell, m}`` modes.

The `map` array should have size Nφ along its first dimension and Nϑ along its
second; any number of dimensions may follow.  The `spin` must be entered
explicitly, and `ℓmax` is the highest ℓ value you want in the output.

For repeated applications of this function with different values of `map`, it
is more efficient to pre-compute `plan` using [`plan_map2salm`](@ref).  These
functions will create a new `salm` array on each call.  To operate in place on
a pre-allocated `salm` array, use [`map2salm!`](@ref).

The core of this function follows the method described by [Reinecke and
Seljebotn](https://dx.doi.org/10.1051/0004-6361/201321494).

"""
function map2salm(map::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int, show_progress=false) where {T<:Real}
    Nφ, Nϑ, Nextra... = size(map)
    salm = zeros(complex(T), (Ysize(ℓmax), Nextra...))
    map2salm!(salm, map, spin, ℓmax, show_progress)
    return salm
end

@doc raw"""
    map2salm!(salm, map, spin, ℓmax)
    map2salm!(salm, map, plan)

Transform `map` values sampled on the sphere to ``{}_sa_{\ell, m}`` modes in
place.

For details, see [`map2salm`](@ref).

"""
function map2salm!(
    salm::AbstractArray{Complex{T}},
    map::AbstractArray{Complex{T}},
    spin::Int, ℓmax::Int, show_progress=false
) where {T<:Real}
    plan = plan_map2salm(map, spin, ℓmax)
    map2salm!(salm, map, plan, show_progress)
end

"""
    plan_map2salm(map, spin, ℓmax)

Precompute values to use in executing [`map2salm`](@ref) or
[`map2salm!`](@ref).

The arguments to this function exactly mirror those of the first form of
[`map2salm`](@ref), and all but the first argument in the first form of
[`map2salm!`](@ref).  The `plan` returned by this function can be passed to the
second forms of those functions to avoid some computation and allocation costs.

Note that the `plan` object is not thread safe; a separate `plan` should be
created for each thread that will use one, or locks should be used to ensure
that a single `plan` is not used at the same time on different threads.

"""
function plan_map2salm(map_data::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int) where {T<:Real}
    Nφ, Nϑ, Nextra... = size(map_data)
    Gs = [Array{complex(T)}(undef, (Nφ,)) for i = 1:nthreads()]
    m′max = abs(spin)
    Hwedge = Array{T}(undef, WignerHsize(ℓmax, m′max))
    H_rec_coeffs = H_recursion_coefficients(ℓmax, T)
    weight = clenshaw_curtis(Nϑ, T)
    expiθ = complex_powers(cis(π / T(Nϑ-1)), Nϑ-1)
    ϵs = SphericalFunctions.ϵ(-spin)
    extra_dims = Base.Iterators.product((1:e for e in Nextra)...)
    fftplan = T<:MachineFloat ? plan_fft(map_data[:, 1, first(extra_dims)...]) : nothing

    return (spin, ℓmax, Nφ, Nϑ, Nextra, Gs, m′max, Hwedge, H_rec_coeffs, weight, expiθ, ϵs, extra_dims, fftplan)
end


function computeG!(
    G::AbstractArray{Complex{T}},
    map::AbstractArray{Complex{T}},
    weight::T, fftplan
) where {T<:Real}
    G[:] = weight * fft(map)
end


function computeG!(
    G::AbstractArray{Complex{T}},
    map::AbstractArray{Complex{T}},
    weight::T, fftplan
) where {T<:MachineFloat}
    @views mul!(G[:], fftplan, map)
    @views G[:] *= weight
end


function map2salm!(
    salm::AbstractArray{Complex{T}}, map::AbstractArray{Complex{T}},
    (spin, ℓmax, Nφ, Nϑ, Nextra, Gs, m′max, Hwedge, H_rec_coeffs, weight, expiθ, ϵs, extra_dims, fftplan),
    show_progress=false
) where {T<:Real}
    s1 = size(salm)
    s2 = (Ysize(ℓmax), Nextra...)
    @assert s1==s2 "size(salm)=$s1  !=  (Ysize(ℓmax), Nextra...)=$s2"

    absspin = abs(spin)
    progress = Progress(Nϑ * prod(Nextra); showspeed=true, enabled=show_progress)

    @inbounds for ϑ ∈ 1:Nϑ
        H!(Hwedge, expiθ[ϑ], ℓmax, m′max, H_rec_coeffs, WignerHindex)
        @threads for extra ∈ collect(extra_dims)
            # NOTE: We can't thread at a higher level because each thread could access the
            # same element of `salm` simultaneously below; by threading at this level, we
            # are assured that the `extra...` index used below is unique to each thread.
            iₜ = Threads.threadid()
            G = Gs[iₜ]
            computeG!(G, map[:, ϑ, extra...], weight[ϑ], fftplan)
            for ℓ ∈ absspin:ℓmax
                λ_factor = ϵs * √((2ℓ+1)*T(π)) / Nφ

                i₀ = WignerHindex(ℓ, spin, 0, m′max)

                let m=0
                    salm[Yindex(ℓ, m), extra...] +=
                        G[m+1] * λ_factor * Hwedge[i₀]
                end

                i₊ = i₀
                i₋ = i₀
                if !signbit(spin)
                    for m ∈ 1:min(ℓ, absspin)
                        i₊ -= ℓ-m+2
                        i₋ += ℓ-m+1
                        salm[Yindex(ℓ, m), extra...] +=
                            G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                        salm[Yindex(ℓ, -m), extra...] +=
                            G[Nφ-m+1] * λ_factor * Hwedge[i₋]
                    end
                else
                    for m ∈ 1:min(ℓ, absspin)
                        i₊ += ℓ-m+1
                        i₋ -= ℓ-m+2
                        salm[Yindex(ℓ, m), extra...] +=
                            G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                        salm[Yindex(ℓ, -m), extra...] +=
                            G[Nφ-m+1] * λ_factor * Hwedge[i₋]
                    end
                end
                for m ∈ absspin+1:ℓ
                    i₊ += 1
                    i₋ += 1
                    salm[Yindex(ℓ, m), extra...] +=
                        G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                    salm[Yindex(ℓ, -m), extra...] +=
                        G[Nφ-m+1] * λ_factor * Hwedge[i₋]
                end
            end
            next!(progress)
        end
    end
end

function map2salm(map::AbstractArray{Complex{T}}, plan, show_progress=false) where {T<:Real}
    salm = zeros(complex(T), (Ysize(plan[2]), plan[5]...))
    map2salm!(salm, map, plan, show_progress)
    return salm
end
