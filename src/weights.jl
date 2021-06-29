"""
    fejer1(n, [T])

Compute weights for Fejér's first rule, wₖᶠ¹ for k ∈ 0:n

This function uses [Waldvogel's method](https://doi.org/10.1007/s10543-006-0045-4).

The type `T` may be any `AbstractFloat`, but defaults to `Float64`.
"""
function fejer1(n, ::Type{T}=Float64) where {T<:AbstractFloat}
    v = Vector{complex(T)}(undef, n)
    complex_powers!(@view(v[1:(n-1)÷2+1]), exp(im*(π/T(n))))
    @inbounds begin
        v[1] = 2
        for k = 1 : (n-1)÷2
            v[k+1] *= 2 / T(1 - 4k^2)
        end
        if iseven(n)
            let k = n ÷ 2
                v[k+1] = 0
            end
        end
        for k = 1 : (n-1)÷2
            v[n-k+1] = conj(v[k+1])
        end
        w = real(ifft(v))
    end
    [w; w[1]]
end


"""
    fejer2(n, [T])

Compute weights for Fejér's second rule, wₖᶠ² for k ∈ 1:n-1.

This function uses [Waldvogel's method](https://doi.org/10.1007/s10543-006-0045-4).
However, contrary to Waldvogel's notation, this routine *does not* return the
weight corresponding to the ϑ=0 or π nodes, which both have weight 0.

The type `T` may be any `AbstractFloat`, but defaults to `Float64`.

"""
function fejer2(n, ::Type{T}) where {T<:AbstractFloat}
    # General function for any type of float
    v = Vector{complex(T)}(undef, n)
    @inbounds begin
        v[1] = 2
        for k = 1 : n÷2-1
            v[k+1] = 2 / T(1 - 4k^2)
        end
        let k = n ÷ 2
            v[k+1] = ((n-2k-2) / T(2k-1))
        end
        for k = 1 : (n-1)÷2
            v[n-k+1] = v[k+1]
        end
        w = real(ifft(v))
    end
    w[2:end]
end

function fejer2(n, ::Type{T}=Float64) where {T<:MachineFloat}
    # Specialized to "machine" floats; significant reduction in memory and increase in speed
    v = Vector{T}(undef, n÷2 + 1)
    @inbounds begin
        v[1] = 2
        for k = 1 : n÷2-1
            v[k+1] = 2 / T(1 - 4k^2)
        end
        let k = n ÷ 2
            v[k+1] = ((n-2k-2) / T(2k-1))
        end
    end
    w = irfft(v, n)
    w[2:end]
end


"""
    clenshaw_curtis(n, [T])

Compute weights for the Clenshaw-Curtis rule, wₖᶜᶜ for k ∈ 0:n.

This function uses [Waldvogel's method](https://doi.org/10.1007/s10543-006-0045-4).

The type `T` may be any `AbstractFloat`, but defaults to `Float64`.

"""
function clenshaw_curtis(n, ::Type{T}) where {T<:AbstractFloat}
    # General function for any type of float
    nmod2 = n % 2
    w₀ᶜᶜ = inv(T(n^2 - 1 + nmod2))
    v = Vector{complex(T)}(undef, n)
    @inbounds begin
        v[1] = 2 - w₀ᶜᶜ
        for k = 1 : n÷2-1
            v[k+1] = 2 / T(1 - 4k^2) - w₀ᶜᶜ
        end
        let k = n ÷ 2
            v[k+1] = ((n-2k-2) / T(2k-1)) + w₀ᶜᶜ * ((2-nmod2)*n-1)
        end
        for k = 1 : (n-1)÷2
            v[n-k+1] = v[k+1]
        end
        w = real(ifft(v))
        w[1] = w₀ᶜᶜ
    end
    [w; w[1]]
end

function clenshaw_curtis(n, ::Type{T}=Float64) where {T<:MachineFloat}
    # Specialized to "machine" floats; significant reduction in memory and increase in speed
    nmod2 = n % 2
    w₀ᶜᶜ = inv(T(n^2 - 1 + nmod2))
    v = Vector{T}(undef, n÷2 + 1)
    @inbounds begin
        v[1] = 2 - w₀ᶜᶜ
        for k = 1 : n÷2-1
            v[k+1] = 2 / T(1 - 4k^2) - w₀ᶜᶜ
        end
        let k = n ÷ 2
            v[k+1] = ((n-2k-2) / T(2k-1)) + w₀ᶜᶜ * ((2-nmod2)*n-1)
        end
        w = irfft(v, n)
        w[1] = w₀ᶜᶜ
    end
    [w; w[1]]
end
