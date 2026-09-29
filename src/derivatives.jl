### Derivatives with respect to the rotor
#
# Wigner's 𝔇 matrices and the spin-weighted harmonics are differentiated here by the
# angular-momentum operators, rather than through the recurrence that computes them.  The
# derivative of a block of degree ℓ along any direction is a combination of values of that
# same block, which are exact wherever the values are, so the derivatives are as accurate as
# the values at every rotor, the poles included.  The kernels below act on one block at a
# time, which is what lets a calculator produce the derivatives of each block as it produces
# the block itself.  The calculators of dual numbers use them to lift each block of values
# into a block of duals (see `lift!` below), and the extensions for Enzyme and Mooncake use
# them in rules for the calculators' steps, while the extensions for ChainRulesCore and
# ReverseDiff use them in rules for `D_array`, `sYlm_array`, and `sYlm_matrix`.
#
# The derivative of 𝔇 along a rotation is given by its generators.  For R(t) = exp(t𝐮/2) R,
# with 𝐮 a vector,
#
#     d𝔇/dt = -i (𝐮 ⋅ 𝐉) 𝔇,        at t = 0,
#
# where 𝐉 is the angular momentum acting on the index m′, with ⟨m′|J_z|m′⟩ = m′ and
# ⟨m′±1|J_±|m′⟩ = √((ℓ∓m′)(ℓ±m′+1)); `test/wigner/poles.jl` checks exactly this.  A general
# tangent Ṙ at R — any quaternion, because 𝔇 is taken to depend on R only through R/‖R‖ —
# is of that form to first order: writing Ṙ = q R, with q = Ṙ R̄ / ‖R‖², the scalar part of
# q changes only the norm of R and drops out, and the vector part v gives 𝐮 = 2v.  So, with
# w = v_x + i v_y,
#
#     d𝔇^ℓ_{m′,m} = -i [2 v_z m′ 𝔇^ℓ_{m′,m} + w̄ a_{m′} 𝔇^ℓ_{m′-1,m} + w b_{m′} 𝔇^ℓ_{m′+1,m}],
#
#     a_m = √((ℓ-m+1)(ℓ+m)),        b_m = √((ℓ+m+1)(ℓ-m)),
#
# in which the terms beyond ±ℓ vanish with their coefficients.  This is the derivative from
# the left, and it couples each element to its neighbors in the same column.  Writing
# instead Ṙ = R q′, with q′ = R̄ Ṙ / ‖R‖² and v′ its vector part, gives the derivative from
# the right, which couples each element to its neighbors in the same row,
#
#     d𝔇^ℓ_{m′,m} = -i [2 v′_z m 𝔇^ℓ_{m′,m} + w′ a_m 𝔇^ℓ_{m′,m-1} + w̄′ b_m 𝔇^ℓ_{m′,m+1}].
#
# A block of 𝔇 restricted in m′ but not in m is differentiated from the right, and one
# restricted in m but not in m′ from the left, so that every neighbor needed is in the block
# already; only a block restricted in both needs values beyond its limits, and a calculator
# computes one row more on each side of such a block.  The harmonics are ₛY_{ℓ,m} = (-1)^s
# √((2ℓ+1)/4π) conj(𝔇^ℓ_{m,-s}), a conjugated row of 𝔇 at a fixed column, so they are
# differentiated from the left, which keeps the spin weight fixed,
#
#     dₛY_{ℓ,m} = i [2 v_z m ₛY_{ℓ,m} + w a_m ₛY_{ℓ,m-1} + w̄ b_m ₛY_{ℓ,m+1}];
#
# from the right, the derivative would couple the harmonics of weight s to those of weights
# s ± 1.
#
# The reverse-mode rules need the adjoint of these linear maps.  A complex cotangent z̄ is
# taken to be ∂L/∂(Re z) + i ∂L/∂(Im z), for a real function L, which is the convention of
# ChainRules, Enzyme, Mooncake, and ReverseDiff alike, so that dL = Σ Re(conj(z̄) dz).  With
# K = Σ conj(z̄) 2m z, P = Σ conj(z̄) a z₋, and Q = Σ conj(z̄) b z₊, where z₋ and z₊ are the
# neighbors that the derivative couples, and m, a, and b are those of the index that it
# steps, that sum is g ⋅ v for
#
#     g = (Im(P + Q), Re(Q - P), Im K)        for 𝔇 from the left,
#     g = (Im(P + Q), Re(P - Q), Im K)        for 𝔇 from the right, and
#     g = (-Im(P + Q), Re(Q - P), -Im K)      for the harmonics.
#
# The last step, from the vector g to the cotangent R̄ of the rotor with g ⋅ v = R̄ ⋅ Ṙ, is
# R̄ = g R / ‖R‖² from the left and R̄ = R g / ‖R‖² from the right, with g taken as a
# pure-vector quaternion.  Either way R̄ is orthogonal to R, as the cotangent of a function
# of R/‖R‖ must be.
#
# The rotor, its tangents, and its cotangents are all read and written here by components,
# never by the arithmetic of `Rotor`, which assumes that a `Rotor` has unit norm: a tangent
# or cotangent stored as a `Rotor`, as Enzyme stores them, does not.


## The generators

# The vector part of Ṙ R̄ / ‖R‖² (from the left) or of R̄ Ṙ / ‖R‖² (from the right), as the
# tuple (v_x, v_y, v_z).  With 𝐯 = (X, Y, Z) and 𝐯̇ likewise, those vector parts are
# W 𝐯̇ - Ẇ 𝐯 ± 𝐯 × 𝐯̇.
@inline function rotor_generator(left::Bool, R, Ṙ)
    W, X, Y, Z = R[1], R[2], R[3], R[4]
    Ẇ, Ẋ, Ẏ, Ż = Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4]
    n² = W^2 + X^2 + Y^2 + Z^2
    cx, cy, cz = Y * Ż - Z * Ẏ, Z * Ẋ - X * Ż, X * Ẏ - Y * Ẋ
    if left
        ((W * Ẋ - Ẇ * X + cx) / n², (W * Ẏ - Ẇ * Y + cy) / n², (W * Ż - Ẇ * Z + cz) / n²)
    else
        ((W * Ẋ - Ẇ * X - cx) / n², (W * Ẏ - Ẇ * Y - cy) / n², (W * Ż - Ẇ * Z - cz) / n²)
    end
end

# The adjoint of `rotor_generator` at R: the components of g R / ‖R‖² (from the left) or of
# R g / ‖R‖² (from the right), for the vector g = (g_x, g_y, g_z) taken as a pure
# quaternion, which are -g ⋅ 𝐯 and W g ∓ 𝐯 × g.
@inline function rotor_cotangent(left::Bool, R, g)
    W, X, Y, Z = R[1], R[2], R[3], R[4]
    gx, gy, gz = g
    n² = W^2 + X^2 + Y^2 + Z^2
    cx, cy, cz = Y * gz - Z * gy, Z * gx - X * gz, X * gy - Y * gx  # 𝐯 × g
    if left
        (-(gx * X + gy * Y + gz * Z) / n², (W * gx - cx) / n², (W * gy - cy) / n², (W * gz - cz) / n²)
    else
        (-(gx * X + gy * Y + gz * Z) / n², (W * gx + cx) / n², (W * gy + cy) / n², (W * gz + cz) / n²)
    end
end

# The coefficients a_m = √((ℓ-m+1)(ℓ+m)) and b_m = √((ℓ+m+1)(ℓ-m)) of the ladder operators,
# in the real type `T`.  Both vanish exactly where the neighbor they multiply lies outside
# -ℓ:ℓ, and are then given as zero rather than as the square root of zero, whose derivative
# is infinite: when `T` is a dual number, the zero partials of the constant would be
# multiplied by that infinity, and give `NaN`.  They are computed in `float_type(T)`, the
# floating-point type underneath any dual numbers (see `src/utilities/lifting.jl`), since
# they are constants.
@inline function ladder_coefficient(n::Int, ::Type{T}) where {T}
    let F = float_type(T)
        n == 0 ? zero(F) : √F(n)
    end
end
@inline ladder_down(ℓ, m, ::Type{T}) where {T} = ladder_coefficient(Int(ℓ - m + 1) * Int(ℓ + m), T)
@inline ladder_up(ℓ, m, ::Type{T}) where {T} = ladder_coefficient(Int(ℓ + m + 1) * Int(ℓ - m), T)

# The generators of every rotor in `rotors`, in each of `N` directions, written as the
# columns of `G`: the generator in direction d of rotor iᵣ is `G[3d-2:3d, iᵣ]`.  The
# directions are given by `tangents(q)`, a tuple of `N` quaternions, each a tuple of four
# components, for the rotor `q`, and the generators are formed at the values `value(q)`.
function set_generators!(G::AbstractMatrix, left::Bool, rotors, value, tangents, ::Val{N}) where {N}
    Base.require_one_based_indexing(G, rotors)
    size(G) == (3N, length(rotors)) || throw(DimensionMismatch(
        "The generators of $(length(rotors)) rotors in $N directions need a matrix of size "
        * "$((3N, length(rotors))), not $(size(G))."
    ))
    @inbounds for iᵣ ∈ eachindex(rotors)
        q = rotors[iᵣ]
        q₀ = value(q)
        Ṙ = tangents(q)
        for d ∈ 1:N
            vx, vy, vz = rotor_generator(left, q₀, Ṙ[d])
            G[3d - 2, iᵣ], G[3d - 1, iᵣ], G[3d, iᵣ] = vx, vy, vz
        end
    end
    G
end


## The kernels

# The blocks here are 3-dimensional arrays with the rotor index first, as the calculators
# store them: [iᵣ, m′, m] for 𝔇, with rows `rows` and columns `cols`, and [iᵣ, s, m] for
# the harmonics, whose third axis is the whole of -ℓ:ℓ.  For 𝔇, the derivatives of the
# elements in the rows `outrows` ⊆ `rows` and the columns `outcols` ⊆ `cols` are written into
# `Ȧ`, laid out as [iᵣ, outrows, outcols]; for the harmonics, those of the elements in the
# spin rows `is` are written into the same positions of `Ȧ`.  Each element of `Ȧ` is
# `combine(x, ẋ)`, where `x` is the value and `ẋ` the tuple of its derivatives in the `N`
# directions whose generators are in `G`, so that a caller can assemble its own number type
# in the same pass.  Every array is indexed under `@inbounds`, so its shape is checked first.

function check_wigner_block(Ȧ, A, ℓ, rows, cols, outrows, outcols, left::Bool, Nᵣ)
    Base.require_one_based_indexing(Ȧ, A)
    within(out, r) = isempty(out) || (first(r) ≤ first(out) && last(out) ≤ last(r))
    ok = (
        size(A, 1) ≥ Nᵣ && size(A, 2) ≥ length(rows) && size(A, 3) ≥ length(cols)
        && size(Ȧ, 1) ≥ Nᵣ && size(Ȧ, 2) ≥ length(outrows) && size(Ȧ, 3) ≥ length(outcols)
        && within(outrows, rows) && within(outcols, cols)
    )
    # Every neighbor that a coefficient does not annihilate must be in the block.
    if ok && !isempty(outrows) && !isempty(outcols)
        out, r = left ? (outrows, rows) : (outcols, cols)
        ok = (first(out) == -ℓ || first(out) > first(r)) && (last(out) == ℓ || last(out) < last(r))
    end
    ok || throw(DimensionMismatch(
        "Cannot differentiate the rows $outrows and columns $outcols of a block of size "
        * "$(size(A)) with rows $rows and columns $cols at ℓ=$ℓ from the "
        * "$(left ? "left" : "right") into an array of size $(size(Ȧ)) for Nᵣ=$Nᵣ."
    ))
    nothing
end

function wigner_block_pushforward!(
    combine, Ȧ::AbstractArray{<:Any, 3}, A::AbstractArray{<:Any, 3}, ℓ, rows, cols,
    outrows, outcols, left::Bool, G::AbstractMatrix, ::Val{N}
) where {N}
    Nᵣ = size(G, 2)
    check_wigner_block(Ȧ, A, ℓ, rows, cols, outrows, outcols, left, Nᵣ)
    Base.require_one_based_indexing(G)
    size(G, 1) ≥ 3N || throw(DimensionMismatch("The generators need $(3N) rows, not $(size(G, 1))."))
    RT = real(eltype(A))
    o′ = isempty(outrows) ? 0 : Int(first(outrows) - first(rows))
    o = isempty(outcols) ? 0 : Int(first(outcols) - first(cols))
    @inbounds for (jₒ, m) ∈ enumerate(outcols), (j′ₒ, m′) ∈ enumerate(outrows)
        j′, j = j′ₒ + o′, jₒ + o
        n = left ? m′ : m
        k = Int(2n)
        a, b = ladder_down(ℓ, n, RT), ladder_up(ℓ, n, RT)
        has₋, has₊ = n > -ℓ, n < ℓ
        for iᵣ ∈ 1:Nᵣ
            x = A[iᵣ, j′, j]
            x₋ = has₋ ? (left ? A[iᵣ, j′ - 1, j] : A[iᵣ, j′, j - 1]) : zero(x)
            x₊ = has₊ ? (left ? A[iᵣ, j′ + 1, j] : A[iᵣ, j′, j + 1]) : zero(x)
            ẋ = ntuple(Val(N)) do d
                # `@inbounds` does not reach into a closure.
                w = @inbounds Complex(G[3d - 2, iᵣ], G[3d - 1, iᵣ])
                vz = @inbounds G[3d, iᵣ]
                if left
                    -im * (vz * k * x + conj(w) * a * x₋ + w * b * x₊)
                else
                    -im * (vz * k * x + w * a * x₋ + conj(w) * b * x₊)
                end
            end
            Ȧ[iᵣ, j′ₒ, jₒ] = combine(x, ẋ)
        end
    end
    Ȧ
end

# The vector g of the note above for each rotor, added into the columns of `Ḡ`, from the
# values `A` and the cotangent `Ā` of the rows `outrows` and columns `outcols`, laid out as
# `Ȧ` is above.
function wigner_block_pullback!(
    Ḡ::AbstractMatrix, A::AbstractArray{<:Any, 3}, Ā::AbstractArray{<:Any, 3}, ℓ, rows, cols,
    outrows, outcols, left::Bool
)
    Nᵣ = size(Ḡ, 2)
    check_wigner_block(Ā, A, ℓ, rows, cols, outrows, outcols, left, Nᵣ)
    Base.require_one_based_indexing(Ḡ)
    size(Ḡ, 1) ≥ 3 || throw(DimensionMismatch("The cotangents need 3 rows, not $(size(Ḡ, 1))."))
    RT = real(eltype(A))
    o′ = isempty(outrows) ? 0 : Int(first(outrows) - first(rows))
    o = isempty(outcols) ? 0 : Int(first(outcols) - first(cols))
    @inbounds for (jₒ, m) ∈ enumerate(outcols), (j′ₒ, m′) ∈ enumerate(outrows)
        j′, j = j′ₒ + o′, jₒ + o
        n = left ? m′ : m
        k = Int(2n)
        a, b = ladder_down(ℓ, n, RT), ladder_up(ℓ, n, RT)
        has₋, has₊ = n > -ℓ, n < ℓ
        for iᵣ ∈ 1:Nᵣ
            ā = conj(Ā[iᵣ, j′ₒ, jₒ])
            x = A[iᵣ, j′, j]
            p = has₋ ? ā * a * (left ? A[iᵣ, j′ - 1, j] : A[iᵣ, j′, j - 1]) : zero(ā * x)
            q = has₊ ? ā * b * (left ? A[iᵣ, j′ + 1, j] : A[iᵣ, j′, j + 1]) : zero(ā * x)
            Ḡ[1, iᵣ] += imag(p + q)
            Ḡ[2, iᵣ] += left ? real(q - p) : real(p - q)
            Ḡ[3, iᵣ] += imag(ā * k * x)
        end
    end
    Ḡ
end

function check_harmonic_block(Ȧ, A, ℓ, is, Nᵣ)
    Base.require_one_based_indexing(Ȧ, A)
    n = Int(2ℓ) + 1
    if !(
        size(A, 1) ≥ Nᵣ && size(A, 3) ≥ n && size(Ȧ, 1) ≥ Nᵣ && size(Ȧ, 3) ≥ n
        && (isempty(is) || (1 ≤ first(is) && last(is) ≤ min(size(A, 2), size(Ȧ, 2))))
    )
        throw(DimensionMismatch(
            "Cannot differentiate the spin rows $is of a block of size $(size(A)) at ℓ=$ℓ into "
            * "an array of size $(size(Ȧ)) for Nᵣ=$Nᵣ."
        ))
    end
    nothing
end

function harmonic_block_pushforward!(
    combine, Ȧ::AbstractArray{<:Any, 3}, A::AbstractArray{<:Any, 3}, ℓ, is,
    G::AbstractMatrix, ::Val{N}
) where {N}
    Nᵣ = size(G, 2)
    check_harmonic_block(Ȧ, A, ℓ, is, Nᵣ)
    Base.require_one_based_indexing(G)
    size(G, 1) ≥ 3N || throw(DimensionMismatch("The generators need $(3N) rows, not $(size(G, 1))."))
    RT = real(eltype(A))
    n = Int(2ℓ) + 1
    @inbounds for j ∈ 1:n
        m = -ℓ + (j - 1)
        k = Int(2m)
        a, b = ladder_down(ℓ, m, RT), ladder_up(ℓ, m, RT)
        for i ∈ is, iᵣ ∈ 1:Nᵣ
            y = A[iᵣ, i, j]
            y₋ = j > 1 ? A[iᵣ, i, j - 1] : zero(y)
            y₊ = j < n ? A[iᵣ, i, j + 1] : zero(y)
            ẏ = ntuple(Val(N)) do d
                w = @inbounds Complex(G[3d - 2, iᵣ], G[3d - 1, iᵣ])
                vz = @inbounds G[3d, iᵣ]
                im * (vz * k * y + w * a * y₋ + conj(w) * b * y₊)
            end
            Ȧ[iᵣ, i, j] = combine(y, ẏ)
        end
    end
    Ȧ
end

function harmonic_block_pullback!(
    Ḡ::AbstractMatrix, A::AbstractArray{<:Any, 3}, Ā::AbstractArray{<:Any, 3}, ℓ, is
)
    Nᵣ = size(Ḡ, 2)
    check_harmonic_block(Ā, A, ℓ, is, Nᵣ)
    Base.require_one_based_indexing(Ḡ)
    size(Ḡ, 1) ≥ 3 || throw(DimensionMismatch("The cotangents need 3 rows, not $(size(Ḡ, 1))."))
    RT = real(eltype(A))
    n = Int(2ℓ) + 1
    @inbounds for j ∈ 1:n
        m = -ℓ + (j - 1)
        k = Int(2m)
        a, b = ladder_down(ℓ, m, RT), ladder_up(ℓ, m, RT)
        for i ∈ is, iᵣ ∈ 1:Nᵣ
            ȳ = conj(Ā[iᵣ, i, j])
            y = A[iᵣ, i, j]
            p = j > 1 ? ȳ * a * A[iᵣ, i, j - 1] : zero(ȳ * y)
            q = j < n ? ȳ * b * A[iᵣ, i, j + 1] : zero(ȳ * y)
            Ḡ[1, iᵣ] -= imag(p + q)
            Ḡ[2, iᵣ] += real(q - p)
            Ḡ[3, iᵣ] -= imag(ȳ * k * y)
        end
    end
    Ḡ
end


## Whole arrays of harmonics
#
# The arrays that `sYlm_array` and `sYlm_matrix` return hold every ℓ, with the modes on the
# last axis in the canonical order, and with a leading rotor axis and spin axis where the
# values were computed for several of either.  Viewed as [iᵣ, s, mode], the modes of each ℓ
# are a harmonic block, so these apply the kernels above to each ℓ in turn.

# The array of harmonics `Y`, as a 3-dimensional array [iᵣ, s, mode], given whether it has a
# rotor axis and a spin axis.
function harmonic_blocks_view(Y::AbstractArray, batched::Bool)
    if batched
        ndims(Y) == 2 ? reshape(Y, size(Y, 1), 1, size(Y, 2)) : Y
    else
        reshape(Y, 1, (ndims(Y) == 1 ? 1 : size(Y, 1)), size(Y)[end])
    end
end

function harmonic_array_pushforward(
    combine, Y::AbstractArray, batched::Bool, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, G::AbstractMatrix, ::Val{N}
) where {IT<:IntegerHalf, N}
    Y₃ = harmonic_blocks_view(Y, batched)
    check_harmonic_array(Y₃, ℓₘᵢₙ, ℓₘₐₓ, size(G, 2))
    CT = promote_type(eltype(Y), Complex{eltype(G)})
    Ẏ = similar(Y, typeof(combine(zero(eltype(Y)), ntuple(_ -> zero(CT), Val(N)))))
    Ẏ₃ = harmonic_blocks_view(Ẏ, batched)
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        r = Yindex(ℓ, -ℓ, ℓₘᵢₙ):Yindex(ℓ, ℓ, ℓₘᵢₙ)
        harmonic_block_pushforward!(
            combine, view(Ẏ₃, :, :, r), view(Y₃, :, :, r), ℓ, axes(Y₃, 2), G, Val(N)
        )
    end
    Ẏ
end

function harmonic_array_pullback!(
    Ḡ::AbstractMatrix, Y::AbstractArray, Ȳ::AbstractArray, batched::Bool, ℓₘᵢₙ::IT, ℓₘₐₓ::IT
) where {IT<:IntegerHalf}
    size(Ȳ) == size(Y) || throw(DimensionMismatch(
        "The cotangent has size $(size(Ȳ)), but the harmonics $(size(Y))."
    ))
    Y₃, Ȳ₃ = harmonic_blocks_view(Y, batched), harmonic_blocks_view(Ȳ, batched)
    check_harmonic_array(Y₃, ℓₘᵢₙ, ℓₘₐₓ, size(Ḡ, 2))
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        r = Yindex(ℓ, -ℓ, ℓₘᵢₙ):Yindex(ℓ, ℓ, ℓₘᵢₙ)
        harmonic_block_pullback!(Ḡ, view(Y₃, :, :, r), view(Ȳ₃, :, :, r), ℓ, axes(Y₃, 2))
    end
    Ḡ
end

function check_harmonic_array(Y₃, ℓₘᵢₙ, ℓₘₐₓ, Nᵣ)
    if size(Y₃, 1) != Nᵣ || size(Y₃, 3) != Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        throw(DimensionMismatch(
            "An array of harmonics of size $(size(Y₃)) as [iᵣ, s, mode] does not hold "
            * "ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ for Nᵣ=$Nᵣ rotors."
        ))
    end
    nothing
end


## Lifting the blocks of a calculator of values
#
# A calculator whose rotors hold derivatives holds a calculator of their values (see
# `src/utilities/lifting.jl`), and after each of that calculator's steps it writes into its
# own block each value together with its derivatives, from the generators of its rotors'
# tangents.  Those generators depend only on the rotors, so they are computed when the
# rotors are set, and a step allocates nothing.  This is how forward-mode numbers are
# lifted; an extension for a reverse-mode tool defines methods of `set_generators!` and
# `lift!` for its own numbers, which record the step instead.

set_generators!(lift::Lift, left::Bool, rotors::AbstractVector{Quaternion{RT}}) where {RT} =
    set_generators!(lift.G, left, rotors, rotor_value, rotor_tangents, Val(ndirections(RT)))

# The stored rows and columns of this calculator's block, from those of the calculator of
# values, which stores one more on each side wherever this one needs it; see
# `stored_limits`.
function lift!(c::WignerCalculator{IT, RT, NT, ST, B, FT, <:Lift}, ℓ::IT) where {IT, RT, NT, ST, B, FT<:Real}
    inner = c.lift.inner
    rows, cols = stored_m′range(inner, ℓ), stored_mrange(inner, ℓ)
    outrows, outcols = stored_m′range(c, ℓ), stored_mrange(c, ℓ)
    wigner_block_pushforward!(
        lift_combine(RT), view(c.Wˡ, :, 1:length(outrows), 1:length(outcols)),
        view(inner.Wˡ, :, 1:length(rows), 1:length(cols)), ℓ, rows, cols, outrows, outcols,
        derivatives_from_left(c), c.lift.G, Val(ndirections(RT))
    )
    c
end
# The spin rows `is` of the block of degree ℓ, written into `Y` after its first `j₀` modes,
# from the calculator of values, which computes them into its own block.
function lift!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, <:Lift}, ℓ::IT, is, Y, j₀::Int
) where {IT, RT, NT, ST, S, B, FT<:Real}
    harmonic_block_pushforward!(
        lift_combine(RT), view(Y, :, :, (j₀ + 1):(j₀ + Int(2ℓ) + 1)), c.lift.inner.Yˡ, ℓ, is,
        c.lift.G, Val(ndirections(RT))
    )
    c
end


## Whole arrays, for the rules of the tools that cannot follow a calculator
#
# ChainRules and ReverseDiff cannot differentiate code that mutates arrays, as a calculator
# does, so their extensions define rules for the functions that return whole arrays instead:
# `D_array`, `sYlm_array`, and `sYlm_matrix_array`.  These are what those rules compute.

# The generators of the tangents `Ṙ` of the rotors `R` as the matrix `G` that the kernels
# read, for one direction, and the cotangents of the rotors from the vectors in the columns
# of `Ḡ` that the kernels return.
function rotor_generators(left::Bool, R::AbstractVector, Ṙ::AbstractVector)
    G = Matrix{float(real(eltype(eltype(R))))}(undef, 3, length(R))
    for i ∈ eachindex(R, Ṙ)
        G[1, i], G[2, i], G[3, i] = rotor_generator(left, R[i], Ṙ[i])
    end
    G
end
rotor_cotangents(left::Bool, R::AbstractVector, Ḡ::AbstractMatrix) =
    [rotor_cotangent(left, R[i], (Ḡ[1, i], Ḡ[2, i], Ḡ[3, i])) for i ∈ eachindex(R)]

# The blocks of `D_array(R, …)`, and the stored rows and columns of each, as [1, m′, m],
# with the calculator that computed them, from which the rules read the ranges of each
# block.
function D_array_with_stored(
    R, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    calc = DCalculator(R, ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    ℓs = ℓₘᵢₙ(IT):ℓₘₐₓ
    blocks = Vector{Matrix{eltype(calc.Wˡ)}}(undef, length(ℓs))
    stored = Vector{Array{eltype(calc.Wˡ), 3}}(undef, length(ℓs))
    for (i, ℓ) ∈ enumerate(ℓs)
        blocks[i] = copy(parent(recurrence!(calc, ℓ)))
        stored[i] = calc.Wˡ[1:1, 1:length(stored_m′range(calc, ℓ)), 1:length(stored_mrange(calc, ℓ))]
    end
    blocks, stored, calc
end

# The derivatives of the blocks of `D_array` along the generator `v`, from the stored rows
# and columns of each, as the matrices of the blocks themselves.
function D_array_pushforward(calc::WignerCalculator{IT}, stored, v) where {IT}
    G = reshape([v[1], v[2], v[3]], 3, 1)
    left = derivatives_from_left(calc)
    map(enumerate(ℓₘᵢₙ(IT):ℓₘₐₓ(calc))) do (i, ℓ)
        rows, cols = stored_m′range(calc, ℓ), stored_mrange(calc, ℓ)
        outrows, outcols = m′range(calc, ℓ), mrange(calc, ℓ)
        Aˢ = stored[i]
        Ȧ = similar(Aˢ, promote_type(eltype(Aˢ), Complex{eltype(G)}), 1, length(outrows), length(outcols))
        wigner_block_pushforward!(
            (x, ẋ) -> only(ẋ), Ȧ, Aˢ, ℓ, rows, cols, outrows, outcols, left, G, Val(1)
        )
        reshape(Ȧ, length(outrows), length(outcols))
    end
end

# The vector g of the note at the top of this file, from the stored rows and columns of each
# block and the cotangents `Ā` of the blocks, of which any may be `nothing` for a zero
# cotangent.
function D_array_pullback(calc::WignerCalculator{IT}, stored, Ā) where {IT}
    Ḡ = zeros(real(eltype(first(stored))), 3, 1)
    left = derivatives_from_left(calc)
    for (i, ℓ) ∈ enumerate(ℓₘᵢₙ(IT):ℓₘₐₓ(calc))
        Āᵢ = Ā[i]
        Āᵢ === nothing && continue
        rows, cols = stored_m′range(calc, ℓ), stored_mrange(calc, ℓ)
        outrows, outcols = m′range(calc, ℓ), mrange(calc, ℓ)
        wigner_block_pullback!(
            Ḡ, stored[i], reshape(Āᵢ, 1, length(outrows), length(outcols)), ℓ,
            rows, cols, outrows, outcols, left
        )
    end
    (Ḡ[1, 1], Ḡ[2, 1], Ḡ[3, 1])
end
