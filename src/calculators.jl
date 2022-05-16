export HCalculator
#export H_recursion_coefficients, HCalculator, DCalculator

abstract type WignerCalculator end

struct HCalculator{T<:Real} <: WignerCalculator
    ℓₘₐₓ::Int
    m′ₘₐₓ::Int
    Hₙ₊₁⁰::Vector{T}
    Hₙ::Vector{T}
    sqrt3::T
    invsqrt2::T
end

function HCalculator(T, ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)
    @assert ℓₘₐₓ ≥ 0
    m′ₘₐₓ = min(abs(m′ₘₐₓ), ℓₘₐₓ)
    Hₙsize = m′ₘₐₓ * (2*ℓₘₐₓ - m′ₘₐₓ + 1) + ℓₘₐₓ + 1
    Hₙ₊₁⁰ = Vector{T}(undef, ℓₘₐₓ+2)
    Hₙ = Vector{T}(undef, Hₙsize)
    HCalculator{T}(ℓₘₐₓ, m′ₘₐₓ, Hₙ₊₁⁰, Hₙ, √T(3), inv(√T(2)))
end

Base.show(io::IO, hc::HCalculator{T}) where T =
    print(io, "HCalculator($(T), $(hc.ℓₘₐₓ), $(hc.m′ₘₐₓ)) ")

"""
    m′offset(HCalculator, ℓ, m′, m)

Find the number of elements between (ℓ, m′, m) and (ℓ, m′+1, m).

"""
function m′offset(HCalculator, ℓ, m′, m)
    m′ ≥ 0 ? ℓ-m′ : ℓ+m′+2
end

"""
    offset(HCalculator, ℓ, m′, m)

Find the linear index of element (ℓ, m′, m).
"""
function offset(HCalculator, ℓ, m′, m)
    m′ₘₐₓ = min(HCalculator.m′ₘₐₓ, ℓ)
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

# struct DCalculator{T<:Real} <: WignerCalculator
#     ℓₘₐₓ::Int
#     Hc::HCalculator{T}
#     𝔇ˡ::OffsetMatrix{Complex{T}}
#     expimα::Vector{Complex{T}}
#     expimγ::Vector{Complex{T}}
# end
# Base.show(io::IO, dc::DCalculator{T}) where T = print(io, "DCalculator($(T), $(dc.ℓₘₐₓ)) ")

# function DCalculator(T, ℓₘₐₓ)
#     Hc = HCalculator(T, ℓₘₐₓ)
#     A = Array{complex(T)}(undef, 2ℓₘₐₓ+1, 2ℓₘₐₓ+1)
#     𝔇ˡ = OffsetMatrix(transpose(A), -ℓₘₐₓ-1, -ℓₘₐₓ-1)
#     expimα = Vector{complex(T)}(undef, ℓₘₐₓ+1)
#     expimγ = Vector{complex(T)}(undef, ℓₘₐₓ+1)
#     DCalculator{T}(ℓₘₐₓ, Hc, 𝔇ˡ, expimα, expimγ)
# end
