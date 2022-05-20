"""
    Diterator(D, ℓₘₐₓ, [ℓₘᵢₙ])

Construct an Iterator that returns sub-matrices of `D`, each of which consists
of elements ``(ℓ,-ℓ,-ℓ)`` through ``(ℓ,ℓ,ℓ)``, for ``ℓ`` from `ℓₘᵢₙ` through
`ℓₘₐₓ`.  By default, `ℓₘᵢₙ` is 0.

Note that the returned objects are *views* into the original `D` data — meaning
that you may alter their values.

Because the result is a matrix restricted to a particular ``ℓ`` value, you can
index the ``(ℓ, m′, m)`` element as `[ℓ+m′+1, ℓ+m+1]`.  For example, you might
use this as something like

    for (ℓ, Dˡ) in zip(ℓₘᵢₙ:ℓₘₐₓ, Diterator(D, ℓₘₐₓ))
        for m′ in -ℓ:ℓ
            for m in -ℓ:ℓ
                Dˡ[ℓ+m′+1, ℓ+m+1] = <...>
            end
        end
    end

Also note that no bounds checking is done, either at instantiation time or
during iteration.  You are responsible for ensuring that the size of `D` and
the values of `ℓₘₐₓ` and `ℓₘᵢₙ` make sense.

"""
struct Diterator{VT<:Vector}
    D::VT
    ℓₘₐₓ::Int
    ℓₘᵢₙ::Int
    function Diterator{VT}(D, ℓₘₐₓ, ℓₘᵢₙ=0) where VT
        #@assert ℓₘₐₓ ≥ ℓₘᵢₙ ≥ 0
        new{VT}(D, ℓₘₐₓ, ℓₘᵢₙ)
    end
end
Diterator(D::VT, ℓₘₐₓ, ℓₘᵢₙ=0) where VT = Diterator{VT}(D, ℓₘₐₓ, ℓₘᵢₙ)

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
Base.IteratorSize(::Type{<:Diterator}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:Diterator}) = Base.HasEltype()
Base.eltype(::Type{Diterator{VT}}) where VT = Base.ReshapedArray{eltype(VT), 2, SubArray{eltype(VT), 1, VT, Tuple{UnitRange{Int64}}, true}, Tuple{}}
Base.length(Di::Diterator) = Di.ℓₘₐₓ - Di.ℓₘᵢₙ + 1
Base.size(Di::Diterator) = (length(Di),)
Base.size(Di::Diterator, dim) = dim > 1 ? 1 : length(Di)


"""
    diterator(d, ℓₘₐₓ, [ℓₘᵢₙ])

Construct an Iterator that returns sub-matrices of `d`, each of which consists
of elements ``(ℓ,-ℓ,-ℓ)`` through ``(ℓ,ℓ,ℓ)``, for ``ℓ`` from `ℓₘᵢₙ` through
`ℓₘₐₓ`.  By default, `ℓₘᵢₙ` is 0.

Note that the returned objects are *views* into the original `d` data — meaning
that you may alter their values.

Because the result is a matrix restricted to a particular ``ℓ`` value, you can
index the ``(ℓ, m′, m)`` element as `[ℓ+m′+1, ℓ+m+1]`.  For example, you might
use this as something like

    for (ℓ, dˡ) in zip(ℓₘᵢₙ:ℓₘₐₓ, diterator(d, ℓₘₐₓ))
        for m′ in -ℓ:ℓ
            for m in -ℓ:ℓ
                dˡ[ℓ+m′+1, ℓ+m+1] = <...>
            end
        end
    end

Also note that no bounds checking is done, either at instantiation time or
during iteration.  You are responsible for ensuring that the size of `d` and
the values of `ℓₘₐₓ` and `ℓₘᵢₙ` make sense.

"""
struct diterator{VT<:Vector}
    d::VT
    ℓₘₐₓ::Int
    ℓₘᵢₙ::Int
    function diterator{VT}(d, ℓₘₐₓ, ℓₘᵢₙ=0) where VT
        #@assert ℓₘₐₓ ≥ ℓₘᵢₙ ≥ 0
        new{VT}(d, ℓₘₐₓ, ℓₘᵢₙ)
    end
end
diterator(d::VT, ℓₘₐₓ, ℓₘᵢₙ=0) where VT = diterator{VT}(d, ℓₘₐₓ, ℓₘᵢₙ)

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
Base.IteratorSize(::Type{<:diterator}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:diterator}) = Base.HasEltype()
Base.eltype(::Type{diterator{VT}}) where VT = Base.ReshapedArray{eltype(VT), 2, SubArray{eltype(VT), 1, VT, Tuple{UnitRange{Int64}}, true}, Tuple{}}
Base.length(di::diterator) = di.ℓₘₐₓ - di.ℓₘᵢₙ + 1
Base.size(di::diterator) = (length(di),)
Base.size(di::diterator, dim) = dim > 1 ? 1 : length(di)


"""
    Yiterator(Y, ℓₘₐₓ, [ℓₘᵢₙ, [iₘᵢₙ]])

Construct an Iterator that returns sub-vectors of `Y`, each of which consists
of elements ``(ℓ,-ℓ)`` through ``(ℓ,ℓ)``, for ``ℓ`` from `ℓₘᵢₙ` through `ℓₘₐₓ`.

Note that the returned objects are *views* into the original `Y` data — meaning
that you may alter their values.

Because the result is a vector restricted to a particular ``ℓ`` value, you can
index the ``(ℓ, m)`` element as `[ℓ+m+1]`.  For example, you might
use this as something like

    for (ℓ, Yˡ) in zip(ℓₘᵢₙ:ℓₘₐₓ, Yiterator(Y, ℓₘₐₓ))
        for m in -ℓ:ℓ
            Yˡ[ℓ+m+1] = <...>
        end
    end

By default, `Y` is assumed to contain all possible values, beginning with
`(0,0)`.  However, if `ℓₘᵢₙ` is not 0, this can be ambiguous: do we mean that
`Y` really starts with the `(0,0)` element and we are just asking to begin the
iteration higher?  Or do we mean that `Y` doesn't even contain data for lower
`ℓ` values?  We can resolve this using `iₘᵢₙ`, which gives the index of `ℓₘᵢₙ`
in `Y`.  By default, we assume the first case, and set `iₘᵢₙ=Ysize(ℓₘᵢₙ-1)+1`.
However, if `Y` doesn't contain data below `ℓₘᵢₙ`, we could use `iₘᵢₙ=1` to
indicate the index in `Y` at which to find ``(ℓₘᵢₙ,-ℓₘᵢₙ)``.

Also note that no bounds checking is done, either at instantiation time or
during iteration.  You are responsible for ensuring that the size of `Y` and
the values of `ℓₘₐₓ`, `ℓₘᵢₙ`, and `iₘᵢₙ` make sense.

"""
struct Yiterator{VT<:Vector}
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
Base.IteratorSize(::Type{<:Yiterator}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:Yiterator}) = Base.HasEltype()
Base.eltype(::Type{Yiterator{VT}}) where VT = SubArray{eltype(VT), 1, VT, Tuple{UnitRange{Int64}}, true}
Base.length(Yi::Yiterator) = Yi.ℓₘₐₓ - Yi.ℓₘᵢₙ + 1
Base.size(Yi::Yiterator) = (length(Yi),)
Base.size(Yi::Yiterator, dim) = dim > 1 ? 1 : length(Yi)
