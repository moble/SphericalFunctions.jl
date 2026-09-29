### Derivatives with respect to the rotor
#
# Wigner's 𝔇 matrices and the spin-weighted harmonics are differentiated here by the
# angular-momentum operators, rather than through the recurrence that computes them.  The
# extensions for ChainRulesCore, Enzyme, ForwardDiff, Mooncake, and ReverseDiff attach the
# functions below, as rules, to `D_array` and `sYlm_array`, so that no tool differentiates
# the recurrence when `D` or `sYlm` is given a rotor.  The derivative along any direction is
# a combination of values of the same ℓ, which are exact wherever the values are, so the
# derivatives are as accurate as the values at every rotor, the poles included; and since
# the rule is applied again to the values it combines, nested applications give exact
# derivatives of every order.
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
#     a_{m′} = √((ℓ-m′+1)(ℓ+m′)),        b_{m′} = √((ℓ+m′+1)(ℓ-m′)),
#
# in which the terms beyond m′ = ±ℓ vanish with their coefficients.  The harmonics are
# ₛY_{ℓ,m} = (-1)^s √((2ℓ+1)/4π) conj(𝔇^ℓ_{m,-s}), a conjugated row of 𝔇 at a fixed column,
# so that
#
#     dₛY_{ℓ,m} = i [2 v_z m ₛY_{ℓ,m} + w a_m ₛY_{ℓ,m-1} + w̄ b_m ₛY_{ℓ,m+1}],
#
# at the same spin weight.  Differentiating with respect to R from the left in this way,
# rather than from the right, is what keeps the spin weight fixed; from the right the
# derivative would couple the harmonics of weight s to those of weights s ± 1.  For 𝔇 either
# would do, but the left is used for both, so that the derivatives of a block restricted in
# m′ need its values one row beyond each of its limits.
#
# The reverse-mode rules need the adjoint of that linear map.  A complex cotangent z̄ is
# taken to be ∂L/∂(Re z) + i ∂L/∂(Im z), for a real function L, which is the convention of
# ChainRules, Enzyme, Mooncake, and ReverseDiff alike, so that dL = Σ Re(conj(z̄) dz).  For
# the harmonics that sum is
#
#     Σₘ Re(conj(Ȳ_m) dY_m) = g ⋅ v,        g = (-Im(P + Q), Re(Q - P), -Im(2A)),
#
# with A = Σ conj(Ȳ_m) m Y_m, P = Σ conj(Ȳ_m) a_m Y_{m-1}, and Q = Σ conj(Ȳ_m) b_m Y_{m+1}
# summed over every ℓ, m, and spin weight; for 𝔇, with the same sums over m′, it is
#
#     g = (Im(P + Q), Re(Q - P), Im(2A)).
#
# The last step, from the vector g to the cotangent R̄ of the rotor with g ⋅ v = R̄ ⋅ Ṙ, is
# R̄ = g R / ‖R‖², with g taken as a pure-vector quaternion.  That R̄ is orthogonal to R, as
# the cotangent of a function of R/‖R‖ must be.
#
# The rotor, its tangents and its cotangents are all read and written here by components,
# never by the arithmetic of `Rotor`, which assumes that a `Rotor` has unit norm: a tangent
# or cotangent stored as a `Rotor`, as Enzyme stores them, does not.


# The vector part of Ṙ R̄ / ‖R‖², as the tuple (v_x, v_y, v_z).  With 𝐯 = (X, Y, Z) and 𝐯̇
# likewise, the vector part of Ṙ R̄ is W 𝐯̇ - Ẇ 𝐯 + 𝐯 × 𝐯̇.
@inline function rotor_generator(R, Ṙ)
    W, X, Y, Z = R[1], R[2], R[3], R[4]
    Ẇ, Ẋ, Ẏ, Ż = Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4]
    n² = W^2 + X^2 + Y^2 + Z^2
    (
        (W * Ẋ - Ẇ * X + (Y * Ż - Z * Ẏ)) / n²,
        (W * Ẏ - Ẇ * Y + (Z * Ẋ - X * Ż)) / n²,
        (W * Ż - Ẇ * Z + (X * Ẏ - Y * Ẋ)) / n²,
    )
end

# The adjoint of `rotor_generator` at R: the components of g R / ‖R‖², for the vector g =
# (g_x, g_y, g_z) taken as a pure quaternion, which is -g ⋅ 𝐯 + W g + g × 𝐯.
@inline function rotor_cotangent(R, g)
    W, X, Y, Z = R[1], R[2], R[3], R[4]
    gx, gy, gz = g
    n² = W^2 + X^2 + Y^2 + Z^2
    (
        -(gx * X + gy * Y + gz * Z) / n²,
        (W * gx + (gy * Z - gz * Y)) / n²,
        (W * gy + (gz * X - gx * Z)) / n²,
        (W * gz + (gx * Y - gy * X)) / n²,
    )
end

# The coefficients a_m = √((ℓ-m+1)(ℓ+m)) and b_m = √((ℓ+m+1)(ℓ-m)) of the ladder operators,
# in the real type `T`.  Both vanish exactly where the neighbor they multiply lies outside
# -ℓ:ℓ, and are then given as zero rather than as the square root of zero, whose derivative
# is infinite: when `T` is a dual number, as it is under nested differentiation, the zero
# partials of the constant would be multiplied by that infinity, and give `NaN`.
@inline ladder_coefficient(n::Int, ::Type{T}) where {T} = n == 0 ? zero(T) : √T(n)
@inline ladder_down(ℓ, m, ::Type{T}) where {T} = ladder_coefficient(Int(ℓ - m + 1) * Int(ℓ + m), T)
@inline ladder_up(ℓ, m, ::Type{T}) where {T} = ladder_coefficient(Int(ℓ + m + 1) * Int(ℓ - m), T)


## The harmonics

# The coefficients a_m and b_m for every m ∈ -ℓ:ℓ, in the leading entries of `a` and `b`,
# which must hold at least 2ℓ+1 of them.  They are computed once for each ℓ, rather than once
# for each value, because their square roots would otherwise cost more than the rest of the
# derivative.
function ladder_coefficients!(a, b, ℓ, ms)
    T = eltype(a)
    for (j, m) ∈ enumerate(ms)
        a[j] = ladder_down(ℓ, m, T)
        b[j] = ladder_up(ℓ, m, T)
    end
    nothing
end

# The arrays of the harmonics are indexed under `@inbounds` below, so their shapes are
# checked first.
function check_sYlm_shape(Y, Ȳ, ℓₘᵢₙ, ℓₘₐₓ)
    if size(Y)[end] != Ysize(ℓₘᵢₙ, ℓₘₐₓ) || size(Ȳ) != size(Y)
        throw(DimensionMismatch(
            "Arrays of sizes $(size(Y)) and $(size(Ȳ)) do not both hold the harmonics for "
            * "ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ, whose mode axis has length $(Ysize(ℓₘᵢₙ, ℓₘₐₓ))."
        ))
    end
    nothing
end

# The derivatives of the values `Y` returned by `sYlm_array` along each of the generators
# in the tuple `vs`, each of the form v = (v_x, v_y, v_z).  The result has the shape of `Y`,
# and each of its elements is `combine(y, ẏ)`, where `y` is the value and `ẏ` the tuple of
# its derivatives along each generator, so that a rule can assemble its own number type,
# such as a dual number, in the same pass.  The mode axis is last, and each row of a leading
# axis of spin weights, if there is one, is differentiated alike.
function sYlm_pushforward(
    combine, Y::AbstractArray, vs::NTuple{N, Any}, ℓₘᵢₙ::IT, ℓₘₐₓ::IT
) where {N, IT<:IntegerHalf}
    RT = real(eltype(Y))
    ws = map(v -> Complex(v[1], v[2]), vs)
    vzs = map(v -> v[3], vs)
    CT = promote_type(eltype(Y), map(typeof, ws)...)
    Ẏ = similar(Y, typeof(combine(zero(eltype(Y)), ntuple(_ -> zero(CT), Val(N)))))
    check_sYlm_shape(Y, Ẏ, ℓₘᵢₙ, ℓₘₐₓ)
    Yₘ = reshape(Y, :, size(Y)[end])
    Ẏₘ = reshape(Ẏ, size(Yₘ))
    a, b = Vector{RT}(undef, Int(2ℓₘₐₓ) + 1), Vector{RT}(undef, Int(2ℓₘₐₓ) + 1)
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        n = Int(2ℓ) + 1
        ladder_coefficients!(a, b, ℓ, -ℓ:ℓ)
        @inbounds for j ∈ 1:n
            twom = 2(j - 1) - Int(2ℓ)
            for k ∈ axes(Yₘ, 1)
                y = Yₘ[k, i₀ + j]
                y₋ = j > 1 ? Yₘ[k, i₀ + j - 1] : zero(y)
                y₊ = j < n ? Yₘ[k, i₀ + j + 1] : zero(y)
                ẏ = ntuple(Val(N)) do d
                    im * (vzs[d] * twom * y + ws[d] * a[j] * y₋ + conj(ws[d]) * b[j] * y₊)
                end
                Ẏₘ[k, i₀ + j] = combine(y, ẏ)
            end
        end
    end
    Ẏ
end

# The derivative along the single generator `v`, as an array of the shape of `Y`.
sYlm_pushforward(Y::AbstractArray, v, ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {IT<:IntegerHalf} =
    sYlm_pushforward((y, ẏ) -> only(ẏ), Y, (v,), ℓₘᵢₙ, ℓₘₐₓ)

# The vector g of the note above, from the values `Y` and the cotangent `Ȳ` of the same
# shape.
function sYlm_pullback(Y::AbstractArray, Ȳ::AbstractArray, ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {IT<:IntegerHalf}
    RT = real(eltype(Y))
    check_sYlm_shape(Y, Ȳ, ℓₘᵢₙ, ℓₘₐₓ)
    Yₘ = reshape(Y, :, size(Y)[end])
    Ȳₘ = reshape(Ȳ, size(Yₘ))
    a, b = Vector{RT}(undef, Int(2ℓₘₐₓ) + 1), Vector{RT}(undef, Int(2ℓₘₐₓ) + 1)
    A = P = Q = zero(promote_type(eltype(Y), eltype(Ȳ)))
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        n = Int(2ℓ) + 1
        ladder_coefficients!(a, b, ℓ, -ℓ:ℓ)
        @inbounds for j ∈ 1:n
            twom = 2(j - 1) - Int(2ℓ)
            for k ∈ axes(Yₘ, 1)
                ȳ = conj(Ȳₘ[k, i₀ + j])
                A += ȳ * twom * Yₘ[k, i₀ + j]
                if j > 1
                    P += ȳ * a[j] * Yₘ[k, i₀ + j - 1]
                end
                if j < n
                    Q += ȳ * b[j] * Yₘ[k, i₀ + j + 1]
                end
            end
        end
    end
    (-imag(P + Q), real(Q - P), -imag(A))
end


## Wigner's 𝔇

# The m′ limits of the values that the derivatives of a block need: one row beyond each
# limit, but not beyond ±ℓₘₐₓ.  Valid limits remain valid when widened; `D_array_widened`
# checks the originals, since only the widened ones reach a calculator.
@inline function D_widened_limits(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where {IT}
    (min(m′ₘₐₓ + 1, ℓₘₐₓ), max(m′ₘᵢₙ - 1, -ℓₘₐₓ))
end

# The values of `D_array` with the m′ limits widened as the derivatives need them, and
# whether that widened anything.  When it did not, the values are exactly those of
# `D_array(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)`.
function D_array_widened(
    R, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    validate_index_ranges(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    m′ₘₐₓʷ, m′ₘᵢₙʷ = D_widened_limits(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
    D_array(R, ℓₘₐₓ, m′ₘₐₓʷ, m′ₘᵢₙʷ, mₘₐₓ, mₘᵢₙ)
end
@inline function D_is_widened(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where {IT}
    D_widened_limits(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ) != (m′ₘₐₓ, m′ₘᵢₙ)
end

# The widened values and the arrays laid out as the block itself are indexed under
# `@inbounds` below, from 1, so their indexing and their lengths are checked first.
function check_D_lengths(Aʷ, A, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    Base.require_one_based_indexing(Aʷ, A)
    m′ₘₐₓʷ, m′ₘᵢₙʷ = D_widened_limits(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
    nʷ = D_offset(ℓₘₐₓ + 1, m′ₘₐₓʷ, m′ₘᵢₙʷ, mₘₐₓ, mₘᵢₙ)
    n = D_offset(ℓₘₐₓ + 1, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    if length(Aʷ) != nʷ || length(A) != n
        throw(DimensionMismatch(
            "The widened values of 𝔇 and the array laid out as its block have lengths "
            * "$(length(Aʷ)) and $(length(A)), but the limits call for $nʷ and $n."
        ))
    end
    nothing
end

# Call `f(ℓ, m′r, n′, i, iʷ, nₘ, n′ʷ)` for each ℓ, where `m′r` is the range of m′ of the
# block, `n′` its length, `nₘ` the number of its columns, and `n′ʷ` the number of rows of the
# widened block.  The offsets `i` and `iʷ` are those of the block and of the widened block,
# advanced so that the j′-th element of column jₘ is at `i + j′ + n′ (jₘ - 1)` in the one,
# and at `iʷ + j′ + n′ʷ (jₘ - 1)` in the other.
@inline function foreach_D_block(f, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where {IT}
    m′ₘₐₓʷ, m′ₘᵢₙʷ = D_widened_limits(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
    o = oʷ = 0
    for ℓ ∈ ℓₘᵢₙ(IT):ℓₘₐₓ
        m′r, mr = D_block_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        m′rʷ = first(D_block_ranges(ℓ, m′ₘₐₓʷ, m′ₘᵢₙʷ, mₘₐₓ, mₘᵢₙ))
        n′, n′ʷ, nₘ = length(m′r), length(m′rʷ), length(mr)
        f(ℓ, m′r, n′, o, oʷ + Int(first(m′r) - first(m′rʷ)), nₘ, n′ʷ)
        o += n′ * nₘ
        oʷ += n′ʷ * nₘ
    end
    nothing
end

# The derivatives of `D_array(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)` along each of the
# generators in the tuple `vs`, from the widened values `Aʷ`, laid out as `D_array` lays out
# the values themselves.  As for `sYlm_pushforward`, each element is `combine(x, ẋ)`, where
# `x` is the value of the block — which is thereby read from the widened values in the same
# pass — and `ẋ` the tuple of its derivatives.
function D_pushforward(
    combine, Aʷ::AbstractVector, vs::NTuple{N, Any},
    ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {N, IT<:IntegerHalf}
    RT = real(eltype(Aʷ))
    ws = map(v -> Complex(v[1], v[2]), vs)
    vzs = map(v -> v[3], vs)
    CT = promote_type(eltype(Aʷ), map(typeof, ws)...)
    ET = typeof(combine(zero(eltype(Aʷ)), ntuple(_ -> zero(CT), Val(N))))
    Ȧ = similar(Aʷ, ET, D_offset(ℓₘₐₓ + 1, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ))
    check_D_lengths(Aʷ, Ȧ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    a, b = Vector{RT}(undef, Int(2ℓₘₐₓ) + 1), Vector{RT}(undef, Int(2ℓₘₐₓ) + 1)
    foreach_D_block(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ) do ℓ, m′r, n′, i, iʷ, nₘ, n′ʷ
        ladder_coefficients!(a, b, ℓ, m′r)
        lower, upper = first(m′r) > -ℓ, last(m′r) < ℓ
        twom′₁ = Int(2first(m′r))
        @inbounds for jₘ ∈ 1:nₘ, j′ ∈ 1:n′
            kʷ = iʷ + j′ + n′ʷ * (jₘ - 1)
            x = Aʷ[kʷ]
            x₋ = j′ > 1 || lower ? Aʷ[kʷ - 1] : zero(x)
            x₊ = j′ < n′ || upper ? Aʷ[kʷ + 1] : zero(x)
            twom′ = twom′₁ + 2(j′ - 1)
            ẋ = ntuple(Val(N)) do d
                -im * (vzs[d] * twom′ * x + conj(ws[d]) * a[j′] * x₋ + ws[d] * b[j′] * x₊)
            end
            Ȧ[i + j′ + n′ * (jₘ - 1)] = combine(x, ẋ)
        end
    end
    Ȧ
end

# The derivative along the single generator `v`.
function D_pushforward(
    Aʷ::AbstractVector, v, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    D_pushforward((x, ẋ) -> only(ẋ), Aʷ, (v,), ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

# The vector g of the note above, from the widened values `Aʷ` and the cotangent `Ā` of
# `D_array(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)`.
function D_pullback(
    Aʷ::AbstractVector, Ā::AbstractVector, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    RT = real(eltype(Aʷ))
    check_D_lengths(Aʷ, Ā, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    a, b = Vector{RT}(undef, Int(2ℓₘₐₓ) + 1), Vector{RT}(undef, Int(2ℓₘₐₓ) + 1)
    T = promote_type(eltype(Aʷ), eltype(Ā))
    sums = Ref((zero(T), zero(T), zero(T)))  # the sums A, P, and Q
    foreach_D_block(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ) do ℓ, m′r, n′, i, iʷ, nₘ, n′ʷ
        ladder_coefficients!(a, b, ℓ, m′r)
        lower, upper = first(m′r) > -ℓ, last(m′r) < ℓ
        twom′₁ = Int(2first(m′r))
        A, P, Q = sums[]
        @inbounds for jₘ ∈ 1:nₘ, j′ ∈ 1:n′
            kʷ = iʷ + j′ + n′ʷ * (jₘ - 1)
            ā = conj(Ā[i + j′ + n′ * (jₘ - 1)])
            A += ā * (twom′₁ + 2(j′ - 1)) * Aʷ[kʷ]
            if j′ > 1 || lower
                P += ā * a[j′] * Aʷ[kʷ - 1]
            end
            if j′ < n′ || upper
                Q += ā * b[j′] * Aʷ[kʷ + 1]
            end
        end
        sums[] = (A, P, Q)
    end
    A, P, Q = sums[]
    (imag(P + Q), real(Q - P), imag(A))
end

# The values of `D_array(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)`, taken from the widened values
# `Aʷ` computed for the same rotor.
function D_narrowed(
    Aʷ::AbstractVector, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    m′ₘₐₓʷ, m′ₘᵢₙʷ = D_widened_limits(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
    A = similar(Aʷ, D_offset(ℓₘₐₓ + 1, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ))
    o = oʷ = 0
    for ℓ ∈ ℓₘᵢₙ(IT):ℓₘₐₓ
        m′r, mr = D_block_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        m′rʷ = first(D_block_ranges(ℓ, m′ₘₐₓʷ, m′ₘᵢₙʷ, mₘₐₓ, mₘᵢₙ))
        n′, n′ʷ = length(m′r), length(m′rʷ)
        for jₘ ∈ eachindex(mr), (j′, m′) ∈ enumerate(m′r)
            A[o + j′ + n′ * (jₘ - 1)] = Aʷ[oʷ + Int(m′ - first(m′rʷ)) + 1 + n′ʷ * (jₘ - 1)]
        end
        o += n′ * length(mr)
        oʷ += n′ʷ * length(mr)
    end
    A
end
