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

@propagate_inbounds function Base.getindex(w::WignerMatrix{NT, IT}, i::Int) where {NT, IT}
    data(w)[i]
end
@propagate_inbounds function Base.getindex(w::WignerMatrix{NT, IT}, m′::IT, m::IT) where {NT, IT}
    @boundscheck begin
        if abs(m′) > m′ₘₐₓ(w)
            throw(BoundsError("m′=$m′ out of bounds for WignerMatrix with m′ₘₐₓ=$(m′ₘₐₓ(w))."))
        end
        if abs(m) > mₘₐₓ(w)
            throw(BoundsError("m=$m out of bounds for WignerMatrix with mₘₐₓ=$(mₘₐₓ(w))."))
        end
    end
    @inbounds data(w)[Int(m′+m′ₘₐₓ(w))+1, Int(m+mₘₐₓ(w))+1]
end

@propagate_inbounds function Base.setindex!(w::WignerMatrix{NT, IT}, v::NT, i::Int) where {NT, IT}
    data(w)[i] = v
end
@propagate_inbounds function Base.setindex!(w::WignerMatrix{NT, IT}, v::NT, m′::IT, m::IT) where {NT, IT}
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
    function WignerDMatrix{NT, IT}(data::Matrix{NT}, ℓ::IT) where {NT, IT}
        if ℓ < 0
            throw(ErrorException("ℓ=$ℓ should be non-negative."))
        end
        if size(data, 1) == 0
            throw(ErrorException("Input data has 0 extent along first dimension."))
        end
        if size(data, 2) == 0
            throw(ErrorException("Input data has 0 extent along second dimension."))
        end
        m′ₘₐₓ = IT((size(data, 1) - 1) // 2)
        mₘₐₓ = IT((size(data, 2) - 1) // 2)
        if ℓ < max(m′ₘₐₓ, mₘₐₓ)
            throw(ErrorException(
                "ℓ=$ℓ should be greater than or equal to both m′ₘₐₓ=$m′ₘₐₓ and mₘₐₓ=$mₘₐₓ."
            ))
        end
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
                    * "the inputs are Rationals: $ℓ, $m′ₘₐₓ, $mₘₐₓ."
                ))
            end
        end
        new(data, ℓ, abs(m′ₘₐₓ), abs(mₘₐₓ))
    end
end
function WignerDMatrix(data::Matrix{NT}, ℓ::IT) where {NT, IT}
    WignerDMatrix{NT, IT}(data, ℓ)
end


struct WignerdMatrix{NT, IT} <: WignerMatrix{NT, IT}
    data::Matrix{NT}
    ℓ::IT
    m′ₘₐₓ::IT
    mₘₐₓ::IT
    function WignerdMatrix{NT, IT}(data::Matrix{NT}, ℓ::IT) where {NT, IT}
        if ℓ < 0
            throw(ErrorException("ℓ=$ℓ should be non-negative."))
        end
        if size(data, 1) == 0
            throw(ErrorException("Input data has 0 extent along first dimension."))
        end
        if size(data, 2) == 0
            throw(ErrorException("Input data has 0 extent along second dimension."))
        end
        m′ₘₐₓ = IT((size(data, 1) - 1) // 2)
        mₘₐₓ = IT((size(data, 2) - 1) // 2)
        if ℓ < max(m′ₘₐₓ, mₘₐₓ)
            throw(ErrorException(
                "ℓ=$ℓ should be greater than or equal to both m′ₘₐₓ=$m′ₘₐₓ and mₘₐₓ=$mₘₐₓ."
            ))
        end
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
                    * "the inputs are Rationals: $ℓ, $m′ₘₐₓ, $mₘₐₓ."
                ))
            end
        end
        new(data, ℓ, m′ₘₐₓ, mₘₐₓ)
    end
end
function WignerdMatrix(data::Matrix{NT}, ℓ::IT) where {NT, IT}
    WignerdMatrix{NT, IT}(data, ℓ)
end


@testitem "WignerMatrix" begin
    import SphericalFunctions.Redesign: WignerDMatrix, WignerdMatrix

    # Check that a negative ℓ value throws an error
    @test_throws "should be non-negative." WignerDMatrix(rand(ComplexF64, 3, 3), -1)
    @test_throws "should be non-negative." WignerdMatrix(rand(Float64, 3, 3), -1)
    @test_throws "should be non-negative." WignerDMatrix(rand(ComplexF64, 2, 2), -1//2)
    @test_throws "should be non-negative." WignerdMatrix(rand(Float64, 2, 2), -1//2)

    for ℓ ∈ Any[collect(0:8); collect(1//2:15//2)]
        # Check that ℓ < m′ₘₐₓ and ℓ < mₘₐₓ throw errors
        @test_throws "both m′ₘₐₓ=" WignerDMatrix(Array{Float64}(undef, Int(2ℓ)+3, Int(2ℓ)+1), ℓ)
        @test_throws "both m′ₘₐₓ=" WignerDMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+3), ℓ)
        @test_throws "both m′ₘₐₓ=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+3, Int(2ℓ)+1), ℓ)
        @test_throws "both m′ₘₐₓ=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+3), ℓ)

        # Check that a mismatch between integer/half-integer throws an error
        if ℓ>0 && ℓ isa Int
            @test_throws "InexactError: Int64(1//2)" WignerDMatrix(rand(ComplexF64, 2, 2), ℓ)
            @test_throws "InexactError: Int64(1//2)" WignerdMatrix(rand(Float64, 2, 2), ℓ)
        elseif ℓ isa Rational
            @test_throws "either integers or half-integer" WignerDMatrix(rand(ComplexF64, 1, 1), ℓ)
            @test_throws "either integers or half-integer" WignerdMatrix(rand(Float64, 1, 1), ℓ)
        end

        for mₘₐₓ ∈ (ℓ isa Rational ? (1//2:ℓ) : (0:ℓ))
            # Check a data array with a dimension of 0 extent throws an error.
            # (Note that we're pretending mₘₐₓ is m′ₘₐₓ for two cases, just for efficiency.)
            @test_throws "along second dim" WignerDMatrix(Array{Float64}(undef, Int(2mₘₐₓ)+1, 0), ℓ)
            @test_throws "along first dim" WignerDMatrix(Array{Float64}(undef, 0, Int(2mₘₐₓ)+1), ℓ)
            @test_throws "along second dim" WignerdMatrix(Array{Float64}(undef, Int(2mₘₐₓ)+1, 0), ℓ)
            @test_throws "along first dim" WignerdMatrix(Array{Float64}(undef, 0, Int(2mₘₐₓ)+1), ℓ)

            for m′ₘₐₓ ∈ (ℓ isa Rational ? (1//2:ℓ) : (0:ℓ))
                # Make a big, dumb array full of the explicit indices to check that indexing
                # works as expected.
                data = [
                    (ℓ, m′, m)
                    for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ, m ∈ -mₘₐₓ:mₘₐₓ
                ]
                for w ∈ (WignerDMatrix(data, ℓ), WignerdMatrix(data, ℓ))
                    @test w.data == data
                    @test w.ℓ == ℓ
                    @test w.m′ₘₐₓ == m′ₘₐₓ
                    @test w.mₘₐₓ == mₘₐₓ
                    for m ∈ -mₘₐₓ:mₘₐₓ
                        for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ
                            @test w[m′, m] == (ℓ, m′, m)
                        end
                    end
                end
            end
        end
    end
end
