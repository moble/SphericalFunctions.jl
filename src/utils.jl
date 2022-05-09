using Base: @propagate_inbounds


mutable struct OffsetVec{T,VV<:AbstractVector{T}} <: AbstractVector{T}
    parent::VV
    offset::Int
end
offset(V::OffsetVec) = V.offset
Base.parent(V::OffsetVec) = V.parent
@propagate_inbounds Base.getindex(V::OffsetVec, i::Int) = parent(V)[i-offset(V)]
@propagate_inbounds function Base.setindex!(V::OffsetVec, val, i::Int)
    parent(V)[i-offset(V)] = val
    V
end


mutable struct OffsetMat{T,MM<:AbstractMatrix{T}} <: AbstractMatrix{T}
    parent::MM
    offset1::Int
    offset2::Int
end
offset1(M::OffsetMat) = M.offset1
offset2(M::OffsetMat) = M.offset2
Base.parent(M::OffsetMat) = M.parent
@propagate_inbounds Base.getindex(M::OffsetMat, i::Int, j::Int) = parent(M)[i-offset1(M), j-offset2(M)]
#@propagate_inbounds Base.getindex(M::OffsetMat, ::Colon, j::Int) = parent(M)[:, j-offset2(M)]
