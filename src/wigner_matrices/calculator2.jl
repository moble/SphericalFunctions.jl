
struct WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real}
    Hwedge::Vector{T}
    Hv::Vector{T}
    Hextra::Vector{T}
    zₐpowers::Vector{Complex{T}}
    zᵧpowers::Vector{Complex{T}}
    z::Vector{Complex{T}}
    function WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T}() where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real}
        if ℓₘᵢₙ < 0
            throw(DomainError("ℓₘᵢₙ = $ℓₘᵢₙ",
                    "This function is only valid for non-negative values of ℓ."
            ))
        end
        if ℓₘₐₓ < 0
            throw(DomainError("ℓₘₐₓ = $ℓₘₐₓ",
                    "This function is only valid for non-negative values of ℓ."
            ))
        end
        if ℓₘᵢₙ > ℓₘₐₓ
            throw(DomainError("ℓₘᵢₙ = $ℓₘᵢₙ",
                    "The input ℓₘᵢₙ must be no greater than the input ℓₘₐₓ=$ℓₘₐₓ."
            ))
        end
        new{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T}(
            zeros(T, WignerHsize(m′ₘₐₓ, ℓₘₐₓ)),
            zeros(T, (ℓₘₐₓ + 1)^2),
            zeros(T, ℓₘₐₓ + 2),
            zeros(Complex{T}, ℓₘₐₓ + 1),
            zeros(Complex{T}, ℓₘₐₓ + 1),
            zeros(Complex{T}, 3)
        )
    end
end
WignerMatrixCalculator(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ=typemax(Int), T=Float64) = WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, min(ℓₘₐₓ, abs(m′ₘₐₓ)), T}()

# Access type parameters
@generated ℓₘᵢₙ(::Type{<:WignerMatrixCalculator{ℓₘᵢₙW, ℓₘₐₓ, m′ₘₐₓ, T}}) where {ℓₘᵢₙW, ℓₘₐₓ, m′ₘₐₓ, T<:Real} = ℓₘᵢₙW
@generated ℓₘₐₓ(::Type{<:WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓW, m′ₘₐₓ, T}}) where {ℓₘᵢₙ, ℓₘₐₓW, m′ₘₐₓ, T<:Real} = ℓₘₐₓW
@generated m′ₘₐₓ(::Type{<:WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓW, T}}) where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓW, T<:Real} = m′ₘₐₓW
@generated T(::Type{<:WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, TW}}) where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, TW<:Real} = TW
ℓₘᵢₙ(w) = ℓₘᵢₙ(typeof(w))
ℓₘₐₓ(w) = ℓₘₐₓ(typeof(w))
m′ₘₐₓ(w) = m′ₘₐₓ(typeof(w))
T(w) = T(typeof(w))

# Access pre-computed array sizes
@generated WignerHsize(::Type{W}) where {W<:WignerMatrixCalculator} = WignerHsize(m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated Wignerdsize(::Type{W}) where {W<:WignerMatrixCalculator} = Wignerdsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated WignerDsize(::Type{W}) where {W<:WignerMatrixCalculator} = WignerDsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated Ysize(::Type{W}) where {W<:WignerMatrixCalculator} = Ysize(ℓₘᵢₙ(W), ℓₘₐₓ(W))
WignerHsize(w::WignerMatrixCalculator) = WignerHsize(typeof(w))
Wignerdsize(w::WignerMatrixCalculator) = Wignerdsize(typeof(w))
WignerDsize(w::WignerMatrixCalculator) = WignerDsize(typeof(w))
Ysize(w::WignerMatrixCalculator) = Ysize(typeof(w))

# Access pre-computed helper data
@inline a(n, m, T) = √((n+1+m) * (n+1-m) / T((2n+1)*(2n+3)))
@inline b(n, m, T) = (m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T((2n-1)*(2n+1)))
@inline d(n, m, T) = (m<0 ? -1 : 1) * (√T((n-m) * (n+m+1))) / 2

function Base.show(io::IO, w::WignerMatrixCalculator)
    print(io, typeof(w))
end
