export HCalculator
#export H_recursion_coefficients, HCalculator, DCalculator

# """
#     H_recursion_coefficients(ℓₘₐₓ, T)
#
# Pre-compute constants used in Wigner H recursion.
#
# """
# function H_recursion_coefficients(ℓₘₐₓ, ::Type{T}) where {T<:Real}
#     aₙᵐ = T[√((n+1+m)*(n+1-m)/T((2n+1)*(2n+3))) for n in 0:ℓₘₐₓ+1 for m in 0:n]
#     bₙᵐ = T[(m<0 ? -1 : 1) * √((n-m-1)*(n-m)/T((2n-1)*(2n+1))) for n in 0:ℓₘₐₓ+1 for m in -n:n]
#     dₙᵐ = T[(m<0 ? -1 : 1) * (√T((n-m)*(n+m+1))) / 2 for n in 0:ℓₘₐₓ+1 for m in -n:n]
#     (aₙᵐ, bₙᵐ, dₙᵐ)
# end

struct HCalculator{T<:Real}
    ℓₘₐₓ::Int
    m′ₘₐₓ::Int
    Hₙ₊₁⁰::Vector{T}
    Hₙ::Vector{T}
end
Base.show(io::IO, hc::HCalculator{T}) where T = print(io, "HCalculator($(T), $(hc.ℓₘₐₓ), $(hc.m′ₘₐₓ)) ")

function HCalculator(T, ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)
    @assert ℓₘₐₓ ≥ 0
    m′ₘₐₓ = min(abs(m′ₘₐₓ), ℓₘₐₓ)
    Hₙsize = m′ₘₐₓ * (2*ℓₘₐₓ - m′ₘₐₓ + 1) + ℓₘₐₓ + 1
    Hₙ₊₁⁰ = Vector{T}(undef, ℓₘₐₓ+2)
    Hₙ = Vector{T}(undef, Hₙsize)
    HCalculator{T}(ℓₘₐₓ, m′ₘₐₓ, Hₙ₊₁⁰, Hₙ)
end

# struct DCalculator{T<:Real}
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
