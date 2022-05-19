"""
    Diterator(D, ℓₘₐₓ, [ℓₘᵢₙ])

Construct an Iterator that returns sub-matrices of `D`, each of which consists
of elements ``(ℓ,-ℓ,-ℓ)`` through ``(ℓ,ℓ,ℓ)``, for ``ℓ`` from `ℓₘᵢₙ` through
`ℓₘₐₓ`.

By default, `ℓₘᵢₙ` is 0.

"""
struct Diterator{VT}
    D::VT
    ℓₘₐₓ::Int
    ℓₘᵢₙ::Int
    function Diterator{VT}(D, ℓₘₐₓ, ℓₘᵢₙ=0) where VT
        #@assert ℓₘₐₓ ≥ ℓₘᵢₙ ≥ 0
        new{VT}(D, ℓₘₐₓ, ℓₘᵢₙ)
    end
end
Diterator(D::VT, ℓₘₐₓ, ℓₘᵢₙ=0) where VT = Diterator{VT}(D, ℓₘₐₓ, ℓₘᵢₙ)

"""
    diterator(d, ℓₘₐₓ, [ℓₘᵢₙ])

Construct an Iterator that returns sub-matrices of `d`, each of which consists
of elements ``(ℓ,-ℓ,-ℓ)`` through ``(ℓ,ℓ,ℓ)``, for ``ℓ`` from `ℓₘᵢₙ` through
`ℓₘₐₓ`.

By default, `ℓₘᵢₙ` is 0.

"""
struct diterator{VT}
    d::VT
    ℓₘₐₓ::Int
    ℓₘᵢₙ::Int
    function diterator{VT}(d, ℓₘₐₓ, ℓₘᵢₙ=0) where VT
        #@assert ℓₘₐₓ ≥ ℓₘᵢₙ ≥ 0
        new{VT}(d, ℓₘₐₓ, ℓₘᵢₙ)
    end
end
diterator(d::VT, ℓₘₐₓ, ℓₘᵢₙ=0) where VT = diterator{VT}(d, ℓₘₐₓ, ℓₘᵢₙ)

"""
    Yiterator(Y, ℓₘₐₓ, [ℓₘᵢₙ, [iₘᵢₙ]])

Construct an Iterator that returns sub-vectors of `Y`, each of which consists
of elements ``(ℓ,-ℓ)`` through ``(ℓ,ℓ)``, for ``ℓ`` from `ℓₘᵢₙ` through `ℓₘₐₓ`.

By default, `Y` is assumed to contain all possible values, beginning with
`(0,0)`.  However, if `ℓₘᵢₙ` is not 0, this can be ambiguous: do we mean that
`Y` really starts with the `(0,0)` element and we are just asking to begin the
iteration higher?  Or do we mean that `Y` doesn't even contain data for lower
`ℓ` values?  We can resolve this using `iₘᵢₙ`, which gives the index of `ℓₘᵢₙ`
in `Y`.  By default, we assume the first case, and set `iₘᵢₙ=Ysize(ℓₘᵢₙ-1)+1`.
However, if `Y` doesn't contain data below `ℓₘᵢₙ`, we could use `iₘᵢₙ=1` to
indicate the index in `Y` at which to find ``(ℓₘᵢₙ,-ℓₘᵢₙ)``.

"""
struct Yiterator{VT}
    Y::VT
    ℓₘₐₓ::Int
    ℓₘᵢₙ::Int
    iₘᵢₙ::Int
    function Yiterator{VT}(Y, ℓₘₐₓ, ℓₘᵢₙ=0, iₘᵢₙ=Ysize(ℓₘᵢₙ-1)+1) where VT
        #@assert ℓₘₐₓ ≥ ℓₘᵢₙ ≥ 0
        new{VT}(Y, ℓₘₐₓ, ℓₘᵢₙ, iₘᵢₙ)
    end
end
Yiterator(Y::VT, ℓₘₐₓ, ℓₘᵢₙ=0) where VT = Yiterator{VT}(Y, ℓₘₐₓ, ℓₘᵢₙ)
Yiterator(Y::VT, ℓₘₐₓ, ℓₘᵢₙ, iₘᵢₙ) where VT = Yiterator{VT}(Y, ℓₘₐₓ, ℓₘᵢₙ, iₘᵢₙ)

function Base.iterate(
    Di::Diterator{VT},
    state=(Di.ℓₘᵢₙ,WignerDsize(Di.ℓₘᵢₙ-1)+1)
) where VT
    if state[1] ≥ Di.ℓₘₐₓ
        nothing
    else
        ℓ = state[1]
        i1 = state[2]
        i2 = i1 + (2ℓ+1)^2 - 1
        Dˡ = reshape(@view(Di.D[i1:i2]), 2ℓ+1, 2ℓ+1)
        (Dˡ, (ℓ+1, i2+1))
    end
end

function Base.iterate(
    di::diterator{VT},
    state=(di.ℓₘᵢₙ,WignerDsize(di.ℓₘᵢₙ-1)+1)
) where VT
    if state[1] ≥ di.ℓₘₐₓ
        nothing
    else
        ℓ = state[1]
        i1 = state[2]
        i2 = i1 + (2ℓ+1)^2 - 1
        dˡ = reshape(@view(di.d[i1:i2]), 2ℓ+1, 2ℓ+1)
        (dˡ, (ℓ+1, i2+1))
    end
end

function Base.iterate(
    Yi::Yiterator{VT},
    state=(Yi.ℓₘᵢₙ,Yi.iₘᵢₙ)
) where VT
    if state[1] ≥ Yi.ℓₘₐₓ
        nothing
    else
        ℓ = state[1]
        i1 = state[2]
        i2 = i1 + (2ℓ+1) - 1
        Yˡ = @view(Yi.Y[i1:i2])
        (Yˡ, (ℓ+1, i2+1))
    end
end
