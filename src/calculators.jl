export HCalculator, DCalculator
export m′offset₊, m′offset₋, offset

abstract type WignerCalculator end

struct HCalculator{T<:Real} <: WignerCalculator
    ℓₘₐₓ::Int
    m′ₘₐₓ::Int
    Hₙ₊₁⁰::Vector{T}
    Hₙ::Vector{T}
    sqrt3::T
    invsqrt2::T
end

Base.show(io::IO, hc::HCalculator{T}) where T =
    print(io, "HCalculator($(T), $(hc.ℓₘₐₓ), $(hc.m′ₘₐₓ))")

"""
    HCalculator(T, ℓₘₐₓ, [m′ₘₐₓ=ℓₘₐₓ])

Construct an HCalculator object.  This object pre-allocates some memory for
computing the "H wedge" for a given ℓ value, and stores some useful constants.

"""
function HCalculator(T, ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)
    @assert ℓₘₐₓ ≥ 0
    m′ₘₐₓ = min(abs(m′ₘₐₓ), ℓₘₐₓ)
    Hₙsize = m′ₘₐₓ * (2*ℓₘₐₓ - m′ₘₐₓ + 1) + ℓₘₐₓ + 1
    Hₙ₊₁⁰ = Vector{T}(undef, ℓₘₐₓ+2)
    Hₙ = Vector{T}(undef, Hₙsize)
    HCalculator{T}(ℓₘₐₓ, m′ₘₐₓ, Hₙ₊₁⁰, Hₙ, √T(3), inv(√T(2)))
end

struct DCalculator{T<:Real} <: WignerCalculator
    ℓₘₐₓ::Int
    Hₙ₊₁⁰::Vector{T}
    Hₙ::Vector{Complex{T}}
    expimα::Vector{Complex{T}}
    expimγ::Vector{Complex{T}}
    sqrt3::T
    invsqrt2::T
end

Base.show(io::IO, dc::DCalculator{T}) where T =
    print(io, "DCalculator($(T), $(dc.ℓₘₐₓ))")

"""
    DCalculator(T, ℓₘₐₓ)

Construct a DCalculator object.  This object pre-allocates some memory for
computing Wigner's 𝔇 matrix for a given ℓ value, and stores some useful
constants.

"""
function DCalculator(T, ℓₘₐₓ)
    @assert ℓₘₐₓ ≥ 0
    Hₙsize = (2ℓₘₐₓ + 1)^2
    Hₙ₊₁⁰ = Vector{T}(undef, ℓₘₐₓ+2)
    Hₙ = Vector{complex(T)}(undef, Hₙsize)
    expimα = Vector{complex(T)}(undef, ℓₘₐₓ+1)
    expimγ = Vector{complex(T)}(undef, ℓₘₐₓ+1)
    DCalculator{T}(ℓₘₐₓ, Hₙ₊₁⁰, Hₙ, expimα, expimγ, √T(3), inv(√T(2)))
end


"""
    m′offset₊(WC, ℓ, m′, m)

Find the number of elements between (ℓ, m′, m) and (ℓ, m′+1, m).

If `Hₙ[i]` refers to the (ℓ, m′, m) element, then `Hₙ[i+m′offset₊(...)]`
refers to the (ℓ, m′+1, m) element.

Note that no testing is done to ensure that any of these elements actually
exist, or that the resulting index will be inbounds.

"""
function m′offset₊(::HCalculator, ℓ, m′, m)
    ifelse(m′ ≥ 0, ℓ-m′, ℓ+m′+2)
end
function m′offset₊(::DCalculator, ℓ, m′, m)
    2ℓ+1
end

"""
    m′offset₋(WC, ℓ, m′, m)

Find the number of elements between (ℓ, m′, m) and (ℓ, m′-1, m).

If `Hₙ[i]` refers to the (ℓ, m′, m) element, then `Hₙ[i-m′offset₋(...)]`
refers to the (ℓ, m′-1, m) element.

Note that no testing is done to ensure that any of these elements actually
exist, or that the resulting index will be inbounds.

"""
function m′offset₋(::HCalculator, ℓ, m′, m)
    ℓ-abs(m′)+1
end
function m′offset₋(::DCalculator, ℓ, m′, m)
    2ℓ+1
end

"""
    offset(WC, ℓ, m′, m)

Find the linear index of element (ℓ, m′, m).
"""
function offset(HC::HCalculator, ℓ, m′, m)
    m′ₘₐₓ = min(HC.m′ₘₐₓ, ℓ)
    if m′<1
        (
            (m′ₘₐₓ + m′) * (2ℓ - m′ₘₐₓ + m′ + 1)
            + 2*(m + m′)
        ) ÷ 2 + 1
    else
        (
            (m′ₘₐₓ + 1) * (2ℓ - m′ₘₐₓ + 2)
            + (m′ - 1) * (2ℓ - m′ + 2)
            + 2*(m - m′)
        ) ÷ 2 + 1
    end
end
function offset(::DCalculator, ℓ, m′, m)
    (ℓ+m′) * (2ℓ+1) + ℓ + m + 1
end
