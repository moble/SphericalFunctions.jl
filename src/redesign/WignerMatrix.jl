abstract type WignerMatrix{NT, IT} end

### General methods for all WignerMatrix types

data(w::WignerMatrix{NT, IT}) where {NT, IT} = w.data
ℓ(w::WignerMatrix{NT, IT}) where {NT, IT} = w.ℓ
m′ₘₐₓ(w::WignerMatrix{NT, IT}) where {NT, IT} = w.m′ₘₐₓ
mₘₐₓ(w::WignerMatrix{NT, IT}) where {NT, IT} = w.mₘₐₓ

ℓₘᵢₙ(::WignerMatrix{NT, IT}) where {NT, IT<:Integer} = zero(IT)
ℓₘᵢₙ(::WignerMatrix{NT, IT}) where {NT, IT<:Rational} = IT(1//2)

is_rational(::WignerMatrix{NT, IT}) where {NT, IT<:Integer} = false
is_rational(::WignerMatrix{NT, IT}) where {NT, IT<:Rational} = true

Base.eltype(::WignerMatrix{NT, IT}) where {NT, IT} = NT
Base.size(w::WignerMatrix{NT, IT}) where {NT, IT} = size(data(w))

function Base.getindex(w::WignerMatrix{NT, IT}, i::Int) where {NT, IT}
    data(w)[i]
end
function Base.getindex(w::WignerMatrix{NT, IT}, m′::IT, m::IT) where {NT, IT}
    data(w)[Int(m′+m′ₘₐₓ(w))+1, Int(m+mₘₐₓ(w))+1]
end

function Base.setindex!(w::WignerMatrix{NT, IT}, v::NT, i::Int) where {NT, IT}
    data(w)[i] = v
end
function Base.setindex!(w::WignerMatrix{NT, IT}, v::NT, m′::IT, m::IT) where {NT, IT}
    data(w)[Int(m′+m′ₘₐₓ(w))+1, Int(m+mₘₐₓ(w))+1] = v
end

function Base.axes(w::WignerMatrix{NT, IT}) where {NT, IT}
    (-m′ₘₐₓ(w):m′ₘₐₓ(w), -mₘₐₓ(w):mₘₐₓ(w))
end


### Specialize to D and d matrices

struct WignerDMatrix{NT, IT} <: WignerMatrix{NT, IT}
    data::Matrix{NT}
    ℓ::IT
    m′ₘₐₓ::IT
    mₘₐₓ::IT
    function WignerDMatrix{NT, IT}(data::Matrix{NT}, ℓ::IT, m′ₘₐₓ::IT=ℓ, mₘₐₓ::IT=ℓ) where {NT, IT}
        if !(NT <: NTuple{3, IT}) && complex(NT) ≢ NT
            throw(ErrorException(
                "WignerDMatrix only supports complex types; the input type is $NT.\n"
                * "Perhaps you meant to use WignerdMatrix?"
            ))
        end
        if IT <: Rational
            if denominator(ℓ) ≠ 2 || denominator(m′ₘₐₓ) ≠ 2 || denominator(mₘₐₓ) ≠ 2
                throw(ErrorException(
                    "Index limits must be either integers or half-integer Rationals; "
                    * "the inputs are Rationals $ℓ, $m′ₘₐₓ, $mₘₐₓ."
                ))
            end
        end
        new(data, ℓ, abs(m′ₘₐₓ), abs(mₘₐₓ))
    end
end
function WignerDMatrix(
    data::Matrix{NT}, ℓ::IT;
    mpmax::IT=ℓ, mmax::IT=ℓ,
    m′ₘₐₓ::IT=mpmax, mₘₐₓ::IT=mmax
) where {NT, IT}
    WignerDMatrix{NT, IT}(data, ℓ, m′ₘₐₓ, mₘₐₓ)
end


struct WignerdMatrix{NT, IT} <: WignerMatrix{NT, IT}
    data::Matrix{NT}
    ℓ::IT
    m′ₘₐₓ::IT
    mₘₐₓ::IT
    function WignerdMatrix{NT, IT}(data::Matrix{NT}, ℓ::IT, m′ₘₐₓ::IT, mₘₐₓ::IT) where {NT, IT}
        if !(NT <: NTuple{3, IT}) && real(NT) ≢ NT
            throw(ErrorException(
                "WignerdMatrix only supports real types; the input type is $NT.\n"
                * "Perhaps you meant to use WignerDMatrix?"
            ))
        end
        if IT <: Rational
            if denominator(ℓ) ≠ 2 || denominator(m′ₘₐₓ) ≠ 2 || denominator(mₘₐₓ) ≠ 2
                throw(ErrorException(
                    "Index limits must be either integers or half-integer Rationals; "
                    * "the inputs are Rationals $ℓ, $m′ₘₐₓ, $mₘₐₓ."
                ))
            end
        end
        new(data, ℓ, abs(m′ₘₐₓ), abs(mₘₐₓ))
    end
end
function WignerdMatrix(
    data::Matrix{NT}, ℓ::IT;
    mpmax::IT=ℓ, mmax::IT=ℓ,
    m′ₘₐₓ::IT=mpmax, mₘₐₓ::IT=mmax
) where {NT, IT}
    WignerdMatrix{NT, IT}(data, ℓ, m′ₘₐₓ, mₘₐₓ)
end


@testitem "WignerMatrix" begin
    import SphericalFunctions.Redesign: WignerDMatrix, WignerdMatrix

    for ℓ ∈ 0:8
        for mₘₐₓ ∈ 0:ℓ
            for m′ₘₐₓ ∈ 0:ℓ
                data = [
                    (ℓ, m′, m)
                    for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ, m ∈ -mₘₐₓ:mₘₐₓ
                ]
                for D ∈ (WignerDMatrix(data, ℓ; m′ₘₐₓ, mₘₐₓ), WignerdMatrix(data, ℓ; m′ₘₐₓ, mₘₐₓ))
                    @test D.data == data
                    @test D.ℓ == ℓ
                    @test D.m′ₘₐₓ == m′ₘₐₓ
                    @test D.mₘₐₓ == mₘₐₓ
                    for m ∈ -mₘₐₓ:mₘₐₓ
                        for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ
                            @test D[m′, m] == (ℓ, m′, m)
                        end
                    end
                end
            end
        end
    end
end
