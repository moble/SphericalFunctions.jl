@doc raw"""
    map2salm(map, spin, ℓmax)
    map2salm(map, plan)

Transform `map` values sampled on the sphere to ``{}_sa_{\ell, m}`` modes.

The `map` array should have size Nφ along its first dimension and Nϑ along its second; any
number of dimensions may follow.  The `spin` must be entered explicitly, and `ℓmax` is the
highest ℓ value you want in the output.

For repeated applications of this function with different values of `map`, it is more
efficient to pre-compute `plan` using [`plan_map2salm`](@ref).  These functions will create
a new `salm` array on each call.  To operate in place on a pre-allocated `salm` array, use
[`map2salm!`](@ref).

The core of this function follows the method described by [Reinecke and Seljebotn](@cite
Reinecke_2013).

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

Transform `map` values sampled on the sphere to ``{}_sa_{\ell, m}`` modes in place.

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

# Per-worker workspace of preallocated memory to modify during map2salm computations.
struct Workspace_map2salm{T<:Real,P}
    G::Vector{Complex{T}}
    fftplan::P
end

@inline function with_workspace(f, pool::Channel)
    ws = take!(pool)
    try
        return f(ws)
    finally
        put!(pool, ws)
    end
end

"""
    plan_map2salm(map, spin, ℓmax)

Precompute values to use in executing [`map2salm`](@ref) or [`map2salm!`](@ref).

The arguments to this function exactly mirror those of the first form of [`map2salm`](@ref),
and all but the first argument in the first form of [`map2salm!`](@ref).  The `plan`
returned by this function can be passed to the second forms of those functions to avoid some
computation and allocation costs.

Note that the `plan` object is not thread safe; a separate `plan` should be created for each
thread that will use one, or locks/Channels should be used to ensure that a single `plan` is
not used at the same time on different threads.

"""
function plan_map2salm(map_data::AbstractArray{Complex{T}}, spin::Int, ℓmax::Int) where {T<:Real}
    Nφ, Nϑ, Nextra... = size(map_data)
    m′max = abs(spin)
    Hwedge = Array{T}(undef, WignerHsize(ℓmax, m′max))
    H_rec_coeffs = H_recursion_coefficients(ℓmax, T)
    weight = clenshaw_curtis(Nϑ, T)
    expiθ = complex_powers(cis(π / T(Nϑ-1)), Nϑ-1)
    ϵs = ϵ(-spin)

    # Indexable without allocation (unlike collect(product(...))).
    extra_dims = CartesianIndices(Nextra)

    # Number of workers seen by `@thread`
    Nworkers = threadpoolsize(:default)

    workspace_pool = if T <: MachineFloat
        # Build FFT plan from a representative slice
        proto_idx = extra_dims[1]
        proto = @views map_data[:, 1, proto_idx.I...]
        fftplan₁ = plan_fft(proto)
        P = typeof(fftplan₁)
        workspace_pool = Channel{Workspace_map2salm{T,P}}(Nworkers)
        for i ∈ 1:Nworkers
            Gᵢ = Vector{Complex{T}}(undef, Nφ)
            planᵢ = plan_fft(proto)
            put!(workspace_pool, Workspace_map2salm{T,P}(Gᵢ, planᵢ))
        end
        workspace_pool
    else
        P = Nothing
        workspace_pool = Channel{Workspace_map2salm{T,Nothing}}(Nworkers)
        for i ∈ 1:Nworkers
            Gᵢ = Vector{Complex{T}}(undef, Nφ)
            planᵢ = nothing
            put!(workspace_pool, Workspace_map2salm{T,P}(Gᵢ, planᵢ))
        end
        workspace_pool
    end

    return (
        spin, ℓmax, Nφ, Nϑ, Nextra, workspace_pool, m′max,
        Hwedge, H_rec_coeffs, weight, expiθ, ϵs, extra_dims
    )
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
    (
        spin, ℓmax, Nφ, Nϑ, Nextra, pool, m′max, Hwedge,
        H_rec_coeffs, weight, expiθ, ϵs, extra_dims
    ),
    show_progress=false
) where {T<:Real}
    s1 = size(salm)
    s2 = (Ysize(ℓmax), Nextra...)
    @assert s1==s2 "size(salm)=$s1  !=  (Ysize(ℓmax), Nextra...)=$s2"

    absspin = abs(spin)
    progress = Progress(Nϑ * prod(Nextra); showspeed=true, enabled=show_progress)
    proglock = ReentrantLock()

    @inbounds for ϑ ∈ 1:Nϑ
        H!(Hwedge, expiθ[ϑ], ℓmax, m′max, H_rec_coeffs, WignerHindex)

        # NOTE: We can't thread at a higher level because each thread could access the
        # same element of `salm` simultaneously below; by threading at this level, we
        # are assured that the `extra.I...` index used below is unique to each thread.
        @threads for extra ∈ extra_dims
            # extra is a CartesianIndex; indices are in extra.I
            with_workspace(pool) do ws
                G = ws.G
                @views computeG!(G, map[:, ϑ, extra.I...], weight[ϑ], ws.fftplan)
                for ℓ ∈ absspin:ℓmax
                    λ_factor = ϵs * √((2ℓ+1)*T(π)) / Nφ

                    i₀ = WignerHindex(ℓ, spin, 0, m′max)

                    let m=0
                        salm[Yindex(ℓ, m), extra.I...] +=
                            G[m+1] * λ_factor * Hwedge[i₀]
                    end

                    i₊ = i₀
                    i₋ = i₀
                    if !signbit(spin)
                        for m ∈ 1:min(ℓ, absspin)
                            i₊ -= ℓ-m+2
                            i₋ += ℓ-m+1
                            salm[Yindex(ℓ, m), extra.I...] +=
                                G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                            salm[Yindex(ℓ, -m), extra.I...] +=
                                G[Nφ-m+1] * λ_factor * Hwedge[i₋]
                        end
                    else
                        for m ∈ 1:min(ℓ, absspin)
                            i₊ += ℓ-m+1
                            i₋ -= ℓ-m+2
                            salm[Yindex(ℓ, m), extra.I...] +=
                                G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                            salm[Yindex(ℓ, -m), extra.I...] +=
                                G[Nφ-m+1] * λ_factor * Hwedge[i₋]
                        end
                    end
                    for m ∈ absspin+1:ℓ
                        i₊ += 1
                        i₋ += 1
                        salm[Yindex(ℓ, m), extra.I...] +=
                            G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                        salm[Yindex(ℓ, -m), extra.I...] +=
                            G[Nφ-m+1] * λ_factor * Hwedge[i₋]
                    end
                end
            end

            if show_progress
                lock(proglock) do
                    next!(progress)
                end
            end
        end
    end
end

function map2salm(map::AbstractArray{Complex{T}}, plan, show_progress=false) where {T<:Real}
    salm = zeros(complex(T), (Ysize(plan[2]), plan[5]...))
    map2salm!(salm, map, plan, show_progress)
    return salm
end
