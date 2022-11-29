@doc raw"""
    golden_ratio_spiral(s, ℓₘₐₓ, [T=Float64])

Cover the sphere 𝕊² with pixels generated by the golden-ratio spiral.
Successive pixels are separated by the azimuthal angle ``Δϕ =
2π(2-φ)``, and are uniformly distributed in ``\cos θ``.

Visually, this is a very reasonable-looking pixelization, with fairly uniform
distance between neighbors, and approximate isotropy.  No two pixels will share
the same values of either ``θ`` or ``ϕ``.  Also note that no point is present
on either the North or South poles.

The returned quantity is a vector of `Rotor`s providing each pixel.

"""
function golden_ratio_spiral(s, ℓₘₐₓ, T=Float64)
    let π = T(π), φ = T(MathConstants.φ)
        N = (ℓₘₐₓ+1)^2 - s^2
        Δϕ = 2π * (2 - φ)
        ϕ = (0:N-1) * Δϕ
        cosθ = range(1, -1, length=N+1)[begin:end-1] .- 1/T(N)
        [from_spherical_coordinates(acos(cosθ), ϕ) for (cosθ, ϕ) in zip(cosθ, ϕ)]
    end
end

@doc raw"""
    sorted_rings(s, ℓₘₐₓ, [T=Float64])

Compute locations of a series of rings labelled by ``j ∈ |s|:ℓₘₐₓ`` (analogous
to ``ℓ``), where each ring will contain ``k = 2j+1`` (analogous to ``m``)
pixels distributed evenly around the ring.  These rings are then sorted, so
that the ring with the most pixels (``j = ℓₘₐₓ``) is closest to the equator,
and the next-largest ring is placed just above or below the equator (depending
on the sign of ``s``), the next just below or above, and so on.  This is
generally a fairly good first guess when minimizing the condition number of
matrices used to solve for mode weughts from function values.  In particular, I
use this to initialize the modified EKKM algorithm, which is then fed into an
optimizer to fine-tune the positions of the rings.

This function does not provide the individual pixels; it just provides the
colatitude values of the rings on which the pixels will be placed.  The pixels
themselves are provided by [`sorted_ring_pixels`](@ref).

"""
function sorted_rings(s, ℓₘₐₓ, T=Float64)
    let πo2 = prevfloat(T(π)/2, s)
        sort(
            collect(LinRange(T(0), T(π), 2+ℓₘₐₓ-abs(s)+1))[begin+1:end-1],
            lt=(x,y)->(abs(x-πo2)<abs(y-πo2)),
            rev=true
        )
    end
end


"""
    sorted_ring_pixels(s, ℓₘₐₓ, [T=Float64])

Cover the sphere 𝕊² with ``(ℓₘₐₓ+1)²-s²`` pixels distributed in rings provided
by [`sorted_rings`](@ref); see that function's documentation for more
description.

The returned quantity is a vector of `Rotor`s providing each pixel.

"""
function sorted_ring_pixels(s, ℓₘₐₓ, T=Float64)
    θrings = sorted_rings(s, ℓₘₐₓ, T)
    [
        from_spherical_coordinates(θ, ϕ)
        for (j,θ) ∈ zip(abs(s):ℓₘₐₓ, θrings)
        for ϕ ∈ LinRange(T(0), 2T(π), 2j+2)[begin:end-1]
    ]
end

"""
    fejer1_rings(s, ℓₘₐₓ, [T=Float64])

Values of the colatitude coordinate (``θ``) appropriate for quadrature by
Fejér's first rule, using weights provided by [`fejer1`](@ref).
"""
function fejer1_rings(s, ℓₘₐₓ, T=Float64)
    # Eq. (12) of Reinecke and Seljebotn
    N = 2ℓₘₐₓ + 1
    let π = T(π)
        [(2n+1)*π/2N for n ∈ 0:N-1]
    end
end

"""
    fejer2_rings(s, ℓₘₐₓ, [T=Float64])

Values of the colatitude coordinate (``θ``) appropriate for quadrature by
Fejér's second rule, using weights provided by [`fejer2`](@ref).
"""
function fejer2_rings(s, ℓₘₐₓ, T=Float64)
    # Eq. (13) of Reinecke and Seljebotn
    N = 2ℓₘₐₓ + 1
    let π = T(π)
        [n*π/N for n ∈ 1:N-1]
    end
end

"""
    clenshaw_curtis_rings(s, ℓₘₐₓ, [T=Float64])

Values of the colatitude coordinate (``θ``) appropriate for quadrature by the
Clenshaw-Curtis rule, using weights provided by [`clenshaw_curtis`](@ref).
"""
function clenshaw_curtis_rings(s, ℓₘₐₓ, T=Float64)
    # Eq. (14) of Reinecke and Seljebotn
    N = 2ℓₘₐₓ + 1
    let π = T(π)
        [n*π/N for n ∈ 0:N]
    end
end
