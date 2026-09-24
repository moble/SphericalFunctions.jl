# Each rule below fills a buffer whose length is set by the number of nodes `n`, starting with
# its first element and under `@inbounds`, so `n` is checked at the top of every method, before
# anything is allocated.

@doc raw"""
    fejer1(n, [T])

Compute `n` weights for Fejér's first rule, corresponding to `n` evenly spaced
nodes between 0 and π exclusive, offset by half a spacing from each end.  That
is, the nodes are located at
```math
θ_k = (2k+1) \frac{π}{2n} \quad k=0, \ldots, n-1,
```
as given by [`fejer1_rings`](@ref).

This function uses [Waldvogel's method](@cite Waldvogel_2006).

The number of nodes `n` must be at least 1.  The type `T` may be any `AbstractFloat`, but
defaults to `Float64`.
"""
function fejer1(n, ::Type{T}=Float64) where {T<:AbstractFloat}
    if n < 1
        throw(ArgumentError("`fejer1` needs at least one node; got n=$n."))
    end
    v = Vector{complex(T)}(undef, n)
    complex_powers!(@view(v[1:(n-1)÷2+1]), cis(π/T(n)))
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
            v[n+1-k] = conj(v[k+1])
        end
        w = real(ifft(v))
    end
    w
end


@doc raw"""
    fejer2(n, [T])

Compute `n` weights for Fejér's second rule, corresponding to `n` evenly spaced
nodes between 0 and π exclusive.  That is, the nodes are located at
```math
θ_k = k \frac{π}{n+1} \quad k=1, \ldots, n.
```

This function uses [Waldvogel's method](@cite Waldvogel_2006).  However,
contrary to Waldvogel's notation, this routine *does not* include the weight
corresponding to the ϑ=0 or π nodes, which both have weight 0.

The number of nodes `n` must be at least 1.  The type `T` may be any `AbstractFloat`, but
defaults to `Float64`.

"""
function fejer2(n, ::Type{T}) where {T<:AbstractFloat}
    # General function for any type of float
    if n < 1
        throw(ArgumentError("`fejer2` needs at least one node; got n=$n."))
    end
    v = Vector{complex(T)}(undef, n+1)
    @inbounds begin
        v[1] = 2
        for k = 1 : (n+1)÷2-1
            v[k+1] = 2 / T(1 - 4k^2)
        end
        let k = (n+1) ÷ 2
            v[k+1] = ((n-2k-1) / T(2k-1))
        end
        for k = 1 : n÷2
            v[n-k+2] = v[k+1]
        end
        w = real(ifft(v))
    end
    w[2:end]
end

function fejer2(n, ::Type{T}=Float64) where {T<:IEEEFloat}
    # Specialized to "machine" floats; significant reduction in memory and increase in speed
    if n < 1
        throw(ArgumentError("`fejer2` needs at least one node; got n=$n."))
    end
    v = Vector{T}(undef, (n+1)÷2 + 1)
    @inbounds begin
        v[1] = 2
        for k = 1 : (n+1)÷2-1
            v[k+1] = 2 / T(1 - 4k^2)
        end
        let k = (n+1) ÷ 2
            v[k+1] = ((n-2k-1) / T(2k-1))
        end
    end
    w = irfft(v, n+1)
    w[2:end]
end


@doc raw"""
    clenshaw_curtis(n, [T])


Compute `n` weights for the Clenshaw-Curtis rule, corresponding to `n` evenly
spaced nodes from 0 to π inclusive.  That is, the nodes are located at
```math
θ_k = k \frac{π}{n-1} \quad k=0, \ldots, n-1.
```

This function uses [Waldvogel's method](@cite Waldvogel_2006).

The number of nodes `n` must be at least 2, since the nodes include both poles.  The type
`T` may be any `AbstractFloat`, but defaults to `Float64`.

"""
function clenshaw_curtis(n, ::Type{T}) where {T<:AbstractFloat}
    # General function for any type of float
    if n < 2
        throw(ArgumentError(
            "`clenshaw_curtis` needs at least two nodes, one at each pole; got n=$n."
        ))
    end
    nmod2 = (n-1) % 2
    w₀ᶜᶜ = inv(T((n-1)^2 - 1 + nmod2))
    v = Vector{complex(T)}(undef, n-1)
    @inbounds begin
        v[1] = 2 - w₀ᶜᶜ
        for k = 1 : (n-1)÷2-1
            v[k+1] = 2 / T(1 - 4k^2) - w₀ᶜᶜ
        end
        let k = (n-1) ÷ 2
            v[k+1] = ((n-2k-3) / T(2k-1)) + w₀ᶜᶜ * ((2-nmod2)*(n-1)-1)
        end
        for k = 1 : (n-2)÷2
            v[n-k] = v[k+1]
        end
        w = real(ifft(v))
        w[1] = w₀ᶜᶜ
    end
    [w; w[1]]
end

function clenshaw_curtis(n, ::Type{T}=Float64) where {T<:IEEEFloat}
    # Specialized to "machine" floats; significant reduction in memory and increase in speed
    if n < 2
        throw(ArgumentError(
            "`clenshaw_curtis` needs at least two nodes, one at each pole; got n=$n."
        ))
    end
    nmod2 = mod(n-1, 2)
    w₀ᶜᶜ = inv(T((n-1)^2 - 1 + nmod2))
    v = Vector{T}(undef, (n-1)÷2 + 1)
    @inbounds begin
        v[1] = 2 - w₀ᶜᶜ
        for k = 1 : (n-1)÷2-1
            v[k+1] = 2 / T(1 - 4k^2) - w₀ᶜᶜ
        end
        let k = (n-1) ÷ 2
            v[k+1] = ((n-2k-3) / T(2k-1)) + w₀ᶜᶜ * ((2-nmod2)*(n-1)-1)
        end
        w = irfft(v, n-1)
        w[1] = w₀ᶜᶜ
    end
    [w; w[1]]
end
