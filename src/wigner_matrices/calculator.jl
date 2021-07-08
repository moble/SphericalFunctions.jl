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

"""
   WignerMatrixCalculator(ℓmin, ℓmax, [m′max], [T])

Construct a `WignerMatrixCalculator` to use in calculating ``𝔇``, ``d``, or
``{}_{s}Y_{ℓ,m}``.

If `m′max` is given, values may be calculated with `ℓ ∈ ℓmin:ℓmax` and `m ∈
-ℓ:ℓ` but only `m′ ∈ ℓmin:m′max.  This will reduce the storage and computation
costs.  Note that if your goal is just to compute spin-weighted spherical
harmonics up to some spin `s`, you should use `m′max=abs(s)`.

If `T<:Real` is given, all calculations will be done using that type; the
default is `Float64`.

See also:
  * [`D!`](@ref) to compute the full ``𝔇`` matrix
  * [`d!`](@ref) to compute the small ``d`` matrix
  * [`Y!`](@ref) to compute spin-weighted spherical harmonic values

"""
WignerMatrixCalculator(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ=typemax(Int), T=Float64) = WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, min(ℓₘₐₓ, abs(m′ₘₐₓ)), T}()

# Access type parameters
ℓₘᵢₙ(::Type{<:WignerMatrixCalculator{ℓₘᵢₙW, ℓₘₐₓ, m′ₘₐₓ, T}}) where {ℓₘᵢₙW, ℓₘₐₓ, m′ₘₐₓ, T<:Real} = ℓₘᵢₙW
ℓₘₐₓ(::Type{<:WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓW, m′ₘₐₓ, T}}) where {ℓₘᵢₙ, ℓₘₐₓW, m′ₘₐₓ, T<:Real} = ℓₘₐₓW
m′ₘₐₓ(::Type{<:WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓW, T}}) where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓW, T<:Real} = m′ₘₐₓW
T(::Type{<:WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, TW}}) where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, TW<:Real} = TW
ℓₘᵢₙ(::WignerMatrixCalculator{ℓₘᵢₙW, ℓₘₐₓ, m′ₘₐₓ, T}) where {ℓₘᵢₙW, ℓₘₐₓ, m′ₘₐₓ, T<:Real} = ℓₘᵢₙW
ℓₘₐₓ(::WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓW, m′ₘₐₓ, T}) where {ℓₘᵢₙ, ℓₘₐₓW, m′ₘₐₓ, T<:Real} = ℓₘₐₓW
m′ₘₐₓ(::WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓW, T}) where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓW, T<:Real} = m′ₘₐₓW
T(::WignerMatrixCalculator{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, TW}) where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, TW<:Real} = TW
# ℓₘᵢₙ(w) = ℓₘᵢₙ(typeof(w))
# ℓₘₐₓ(w) = ℓₘₐₓ(typeof(w))
# m′ₘₐₓ(w) = m′ₘₐₓ(typeof(w))
# T(w) = T(typeof(w))

# Access pre-computed array sizes
WignerHsize(::Type{W}) where {W<:WignerMatrixCalculator} = WignerHsize(m′ₘₐₓ(W), ℓₘₐₓ(W))
Wignerdsize(::Type{W}) where {W<:WignerMatrixCalculator} = Wignerdsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
WignerDsize(::Type{W}) where {W<:WignerMatrixCalculator} = WignerDsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
Ysize(::Type{W}) where {W<:WignerMatrixCalculator} = Ysize(ℓₘᵢₙ(W), ℓₘₐₓ(W))
WignerHsize(w::WignerMatrixCalculator) = WignerHsize(typeof(w))
Wignerdsize(w::WignerMatrixCalculator) = Wignerdsize(typeof(w))
WignerDsize(w::WignerMatrixCalculator) = WignerDsize(typeof(w))
Ysize(w::WignerMatrixCalculator) = Ysize(typeof(w))

# Access pre-computed helper data
function a(::Type{W}) where {W<:WignerMatrixCalculator}
    [√((n+1+m) * (n+1-m) / T(W)((2n+1)*(2n+3))) for n in 0:ℓₘₐₓ(W)+1 for m in 0:n]
end
a(w) = a(typeof(w))
function b(::Type{W}) where {W<:WignerMatrixCalculator}
    [(m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T(W)((2n-1)*(2n+1))) for n in 0:ℓₘₐₓ(W)+1 for m in -n:n]
end
b(w) = b(typeof(w))
function d(::Type{W}) where {W<:WignerMatrixCalculator}
    [(m<0 ? -1 : 1) * (√T(W)((n-m) * (n+m+1))) / 2 for n in 0:ℓₘₐₓ(W)+1 for m in -n:n]
end
d(w) = d(typeof(w))

function Base.show(io::IO, w::WignerMatrixCalculator)
    print(io, typeof(w))
end
