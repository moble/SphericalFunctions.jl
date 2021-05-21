### TODO: Merge the Wigner and WignerWorkspace types; there's no real reason for the singleton


struct Wigner{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real} end
function Wigner(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ=typemax(Int), T::Type{<:Real}=Float64)
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
    Wigner{ℓₘᵢₙ, ℓₘₐₓ, min(ℓₘₐₓ, abs(m′ₘₐₓ)), T}
end

# Access type parameters
@generated ℓₘᵢₙ(::Type{W}) where {W<:Wigner} = W.parameters[1]
@generated ℓₘₐₓ(::Type{W}) where {W<:Wigner} = W.parameters[2]
@generated m′ₘₐₓ(::Type{W}) where {W<:Wigner} = W.parameters[3]
@generated T(::Type{W}) where {W<:Wigner} = W.parameters[4]

# Access pre-computed array sizes
@generated WignerHsize(::Type{W}) where {W<:Wigner} = WignerHsize(m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated Wignerdsize(::Type{W}) where {W<:Wigner} = Wignerdsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated WignerDsize(::Type{W}) where {W<:Wigner} = WignerDsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated Ysize(::Type{W}) where {W<:Wigner} = Ysize(ℓₘᵢₙ(W), ℓₘₐₓ(W))

# Access pre-computed helper data
@generated function a(::Type{W}) where {W<:Wigner}
    [√((n+1+m) * (n+1-m) / T(W)((2n+1)*(2n+3))) for n in 0:ℓₘₐₓ(W)+1 for m in 0:n]
end
@generated function b(::Type{W}) where {W<:Wigner}
    [(m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T(W)((2n-1)*(2n+1))) for n in 0:ℓₘₐₓ(W)+1 for m in -n:n]
end
@generated function d(::Type{W}) where {W<:Wigner}
    [(m<0 ? -1 : 1) * (√T(W)((n-m) * (n+m+1))) / 2 for n in 0:ℓₘₐₓ(W)+1 for m in -n:n]
end

struct WignerWorkspace{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real}
    W::Type{Wigner{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T}}
    Hwedge::Vector{T}
    Hv::Vector{T}
    Hextra::Vector{T}
    zₐpowers::Vector{Complex{T}}
    zᵧpowers::Vector{Complex{T}}
    z::Vector{Complex{T}}
    function WignerWorkspace(W::Type{Wigner{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T}}) where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real}
        new{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T}(
            W,
            zeros(T, WignerHsize(W)),
            zeros(T, (ℓₘₐₓ + 1)^2),
            zeros(T, ℓₘₐₓ + 2),
            zeros(Complex{T}, ℓₘₐₓ + 1),
            zeros(Complex{T}, ℓₘₐₓ + 1),
            zeros(Complex{T}, 3)
        )
    end
end

function Base.show(io::IO, workspace::WignerWorkspace)
    print(io, typeof(workspace))
end

include("wignerHrecursions.jl")
include("wignerfunctions.jl")
