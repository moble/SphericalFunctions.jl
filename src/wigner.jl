struct Wugner{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real} end
function Wugner(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ=typemax(Int), T::Type{<:Real}=Float64)
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
    Wugner{ℓₘᵢₙ, ℓₘₐₓ, min(ℓₘₐₓ, abs(m′ₘₐₓ)), T}
end

# Access type parameters
@generated ℓₘᵢₙ(::Type{W}) where {W<:Wugner} = W.parameters[1]
@generated ℓₘₐₓ(::Type{W}) where {W<:Wugner} = W.parameters[2]
@generated m′ₘₐₓ(::Type{W}) where {W<:Wugner} = W.parameters[3]
@generated T(::Type{W}) where {W<:Wugner} = W.parameters[4]

# Access pre-computed array sizes
@generated WignerHsize(::Type{W}) where {W<:Wugner} = WignerHsize(m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated Wignerdsize(::Type{W}) where {W<:Wugner} = Wignerdsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated WignerDsize(::Type{W}) where {W<:Wugner} = WignerDsize(ℓₘᵢₙ(W), m′ₘₐₓ(W), ℓₘₐₓ(W))
@generated Ysize(::Type{W}) where {W<:Wugner} = Ysize(ℓₘᵢₙ(W), ℓₘₐₓ(W))

# Access pre-computed helper data
@generated function a(::Type{W}) where {W<:Wugner}
    [√((n+1+m) * (n+1-m) / T(W)((2*n+1)*(2*n+3))) for n in 0:ℓₘₐₓ(W)+1 for m in 0:n]
end
@generated function b(::Type{W}) where {W<:Wugner}
    [(m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T(W)((2*n-1)*(2*n+1))) for n in 0:ℓₘₐₓ+1 for m in -n:n]
end
@generated function d(::Type{W}) where {W<:Wugner}
    [(m<0 ? -0.5 : 0.5) * √T(W)((n-m) * (n+m+1)) for n in 0:ℓₘₐₓ+1 for m in -n:n]
end
@generated function g(::Type{W}) where {W<:Wugner}
    [2(m+1) / √T(W)((n-m)*(n+m+1)) for n in 0:ℓₘₐₓ+1 for m in -n:n]
end
@generated function h(::Type{W}) where {W<:Wugner}
    [n==m ? NaN : √((n+m+2)*(n-m-1) / T(W)((n-m)*(n+m+1))) for n in 0:ℓₘₐₓ+1 for m in -n:n]
end






struct Wagner{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real}
    Hsize::Int
    dsize::Int
    Dsize::Int
    Ysize::Int
    a::Vector{T}
    b::Vector{T}
    d::Vector{T}
    g::Vector{T}
    h::Vector{T}
    function Wagner{ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T}() where {ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T<:Real}
        # Some generators to use below
        nm = ((n, m) for n in 0:ℓₘₐₓ+1 for m in -n:n)
        absnm = ((n, m) for n in 0:ℓₘₐₓ+1 for m in 0:n)

        # Pre-compute these useful arrays
        a = [√((n+1+m) * (n+1-m) / T((2*n+1)*(2*n+3))) for (n, m) in absnm]
        b = [(m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T((2*n-1)*(2*n+1))) for (n, m) in nm]
        d = [(m<0 ? -0.5 : 0.5) * √T((n-m) * (n+m+1)) for (n, m) in nm]
        g = [2(m+1) / √T((n-m)*(n+m+1)) for (n, m) in nm]
        h = [n==m ? NaN : √((n+m+2)*(n-m-1) / T((n-m)*(n+m+1))) for (n, m) in nm]

        new(
            WignerHsize(m′ₘₐₓ, ℓₘₐₓ),
            Wignerdsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ),
            WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ),
            Ysize(ℓₘᵢₙ, ℓₘₐₓ),
            a, b, d, g, h
        )
    end
end

function Wagner(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ=typemax(Int), T::Type{<:Real}=Float64)
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
    Wagner{ℓₘᵢₙ, ℓₘₐₓ, min(ℓₘₐₓ, abs(m′ₘₐₓ)), T}()
end

function Base.show(io::IO, w::Wagner)
    print(io, typeof(w))
end



struct Wogner{Lmin, Lmax, Mmax, T<:Real, Hsize, dsize, Dsize, Ysize, a, b, d, g, h}
end

function Wogner(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ=typemax(Int), T::Type{<:Real}=Float64)
    m′ₘₐₓ = min(ℓₘₐₓ, abs(m′ₘₐₓ))

    # Some generators to use below
    nm = ((n, m) for n in 0:ℓₘₐₓ+1 for m in -n:n)
    absnm = ((n, m) for n in 0:ℓₘₐₓ+1 for m in 0:n)

    Wogner{
        ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, T,
        WignerHsize(m′ₘₐₓ, ℓₘₐₓ),
        WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ),
        WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ),
        Ysize(ℓₘᵢₙ, ℓₘₐₓ),
        tuple([√((n+1+m) * (n+1-m) / T((2*n+1)*(2*n+3))) for (n, m) in absnm]...),
        tuple([(m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T((2*n-1)*(2*n+1))) for (n, m) in nm]...),
        tuple([(m<0 ? -0.5 : 0.5) * √T((n-m) * (n+m+1)) for (n, m) in nm]...),
        tuple([2(m+1) / √T((n-m)*(n+m+1)) for (n, m) in nm]...),
        tuple([n==m ? NaN : √((n+m+2)*(n-m-1) / T((n-m)*(n+m+1))) for (n, m) in nm]...)
    }()
end

function Base.show(io::IO, w::Wogner)
    wtype = typeof(w)
    print(io, "Wogner($(wtype.parameters[1]), $(wtype.parameters[2]), $(wtype.parameters[3]), $(wtype.parameters[4]))")
end

struct Wigner{T<:Real}
    ℓₘᵢₙ::Int
    ℓₘₐₓ::Int
    m′ₘₐₓ::Int
    Hsize::Int
    dsize::Int
    Dsize::Int
    Ysize::Int
    a::Vector{T}
    b::Vector{T}
    d::Vector{T}
    g::Vector{T}
    h::Vector{T}
    function Wigner{T}(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ=typemax(Int)) where {T<:Real}
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

        # Some generators to use below
        nm = ((n, m) for n in 0:ℓₘₐₓ+1 for m in -n:n)
        absnm = ((n, m) for n in 0:ℓₘₐₓ+1 for m in 0:n)

        # Pre-compute these useful arrays
        a = [√((n+1+m) * (n+1-m) / T((2*n+1)*(2*n+3))) for (n, m) in absnm]
        b = [(m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T((2*n-1)*(2*n+1))) for (n, m) in nm]
        d = [(m<0 ? -0.5 : 0.5) * √T((n-m) * (n+m+1)) for (n, m) in nm]
        g = [2(m+1) / √T((n-m)*(n+m+1)) for (n, m) in nm]
        h = [n==m ? NaN : √((n+m+2)*(n-m-1) / T((n-m)*(n+m+1))) for (n, m) in nm]

        new(
            ℓₘᵢₙ, ℓₘₐₓ, min(ℓₘₐₓ, abs(m′ₘₐₓ)),
            WignerHsize(m′ₘₐₓ, ℓₘₐₓ),
            WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ),
            WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ),
            Ysize(ℓₘᵢₙ, ℓₘₐₓ),
            a, b, d, g, h
        )
    end
end


function Base.show(io::IO, w::Wigner{T}) where {T<:Real}
    print(io, "Wigner{$T}($(w.ℓₘᵢₙ), $(w.ℓₘₐₓ), $(w.m′ₘₐₓ))")
end
