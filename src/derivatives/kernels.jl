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
# ReverseDiff use them in rules for `D_array`, `d_array`, `sYlm_array`, and `sYlm_matrix`.
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
# in which the terms beyond ±ℓ vanish with their coefficients, which are computed by
# `ladder_down` and `ladder_up` in `src/mode_weights/operators.jl`.  This is the derivative
# from the left, and it couples each element to its neighbors in the same column.  Writing
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
# The matrices `d`, and the real harmonics ₛλ_{ℓ,m}, belong to the rotations about y,
# exp(β𝐣/2), which commute with 𝐣.  So a tangent β̇ of the angle is the generator
# v = (0, β̇/2, 0) from the left and from the right alike, and the formulas above become real:
#
#     ḋ_{m′,m} = (β̇/2) (b_{m′} d_{m′+1,m} - a_{m′} d_{m′-1,m}) = (β̇/2) (a_m d_{m′,m-1} - b_m d_{m′,m+1}),
#     ₛλ̇_{ℓ,m} = (θ̇/2) (b_m ₛλ_{ℓ,m+1} - a_m ₛλ_{ℓ,m-1}).
#
# The two forms of ḋ are equal, so the side is chosen as for 𝔇, by the neighbors that a
# block holds, and otherwise for speed (see `derivatives_from_left`).  A phase e^{iβ} is
# differentiated through its argument, and a rotor through the β ∈ [0, π] of its Euler
# decomposition (see `rotation_angle`).  Each of these derivatives is a difference of two
# products that nearly cancel near the poles, and a gradient sums many of them with weights
# of either sign.  So the coefficients a and b are used to about twice the working precision
# (see `step_coefficients`), since the rounding error of a coefficient is the same in every
# element of the line that it multiplies, and would not average out in such a sum; and the
# reverse kernels add each element's term, with its exact rounding error, into a compensated
# sum.
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


# The generators of the rotors of a calculator of `d` or of ₛλₗₘ, which are the rotations
# about the y axis by the angles β or θ it is given: the generator of a tangent β̇ is (0,
# β̇/2, 0), so only its y component is held, in `vy[d, iᵣ]` for direction d of rotor iᵣ,
# where `G[3d-1, iᵣ]` would hold it.  The values of such a calculator are real, and so are
# their derivatives along these generators, which the kernels compute in real arithmetic.
# - `T` is the element type.
# - `M` is the type of the matrix.
struct AngleGenerators{T, M<:AbstractMatrix{T}}
    vy::M
end

# The cotangents of the generators of a calculator of `d` or of ₛλₗₘ, as the pullbacks below
# add them: the y component of the vector g of each rotor in `vy[1, iᵣ]`, and the rounding
# error of that sum in `c[1, iᵣ]`, so that the cotangent of the angle of rotor iᵣ is half
# their sum.  The terms of these sums have either sign and nearly cancel, and the
# compensation keeps their rounding errors from accumulating.
# - `T` is the element type.
# - `M` is the type of the matrices.
struct AngleCotangents{T, M<:AbstractMatrix{T}}
    vy::M
    c::M
    function AngleCotangents(vy::M, c::M) where {T, M<:AbstractMatrix{T}}
        size(c) == size(vy) || throw(DimensionMismatch(
            "The compensation of cotangents of size $(size(vy)) has size $(size(c))."
        ))
        new{T, M}(vy, c)
    end
end
angle_cotangent(Ḡ::AngleCotangents, iᵣ) = (Ḡ.vy[1, iᵣ] + Ḡ.c[1, iᵣ]) / 2
# The sum of `x` and `y` and its rounding error, by Knuth's algorithm, which needs no branch.
@inline function two_sum(x, y)
    s = x + y
    z = s - x
    (s, (x - (s - z)) + (y - z))
end
# The sum s + x, with the compensation c of s, to which the rounding error of the sum is
# added.
@inline function compensated_sum(s, c, x)
    t, e = two_sum(s, x)
    (t, c + e)
end
# Add x, and its own compensation c, into the cotangent of rotor iᵣ.
@inline function compensated_add!(Ḡ::AngleCotangents, iᵣ, x, c=zero(x))
    @inbounds Ḡ.vy[1, iᵣ], Ḡ.c[1, iᵣ] = compensated_sum(Ḡ.vy[1, iᵣ], Ḡ.c[1, iᵣ] + c, x)
    nothing
end

# The generators of the tangents β̇ of the angles of a calculator of `d` or of ₛλₗₘ, in one
# direction; the cotangents of zero for `Nᵣ` rotors of a calculator of elements of type
# `NT`, into which a pullback adds; and the addition of the cotangents of the angles from
# those of their generators into `β̄`.
angle_generators(β̇::AbstractVector) = AngleGenerators(reshape([b / 2 for b ∈ β̇], 1, length(β̇)))
zero_cotangents(::Type{NT}, Nᵣ::Int) where {NT<:Complex} = zeros(real(NT), 3, Nᵣ)
zero_cotangents(::Type{NT}, Nᵣ::Int) where {NT<:Real} =
    AngleCotangents(zeros(NT, 1, Nᵣ), zeros(NT, 1, Nᵣ))
function add_angle_cotangents!(β̄::AbstractVector, Ḡ::AngleCotangents)
    for i ∈ eachindex(β̄)
        β̄[i] += angle_cotangent(Ḡ, i)
    end
    β̄
end

# The type of the derivatives along the generators `G` of values of type `T`.
derivative_type(::Type{T}, G::AbstractMatrix) where {T} = promote_type(T, Complex{eltype(G)})
derivative_type(::Type{T}, G::AngleGenerators) where {T} = promote_type(T, eltype(G.vy))

# The number of rotors for which the generators or cotangents `G` are given, and the check
# that they hold `N` directions, for a kernel that reads them under `@inbounds`.
generator_columns(G::AbstractMatrix) = size(G, 2)
generator_columns(G::Union{AngleGenerators, AngleCotangents}) = size(G.vy, 2)
function check_generators(G::AbstractMatrix, N, name)
    Base.require_one_based_indexing(G)
    size(G, 1) ≥ 3N || throw(DimensionMismatch("The $name need $(3N) rows, not $(size(G, 1))."))
    nothing
end
function check_generators(G::Union{AngleGenerators, AngleCotangents}, N, name)
    Base.require_one_based_indexing(G.vy)
    size(G.vy, 1) ≥ N || throw(DimensionMismatch("The $name need $N rows, not $(size(G.vy, 1))."))
    nothing
end

# The coefficients of a derivative that steps the index n: k = 2n, a_n, and b_n.  Along the
# generators of angles, a_n and b_n are each given as a pair whose sum is the coefficient to
# about twice the working precision (see `ladder_pair`): the derivative of an element of `d`
# or of ₛλₗₘ is then the difference of two products, which nearly cancel, and the rounding
# error of a coefficient, the same in every element of a line that it multiplies, would not
# average out when the derivatives are summed, as a gradient sums them.
@inline step_coefficients(::AbstractMatrix, ℓ, n, ::Type{RT}) where {RT} =
    (Int(2n), ladder_down(ℓ, n, RT), ladder_up(ℓ, n, RT))
@inline step_coefficients(::Union{AngleGenerators, AngleCotangents}, ℓ, n, ::Type{RT}) where {RT} =
    (Int(2n), ladder_pair(Int(ℓ - n + 1) * Int(ℓ + n), RT), ladder_pair(Int(ℓ + n + 1) * Int(ℓ - n), RT))
# The coefficient √N of `ladder_coefficient`, as the pair (c, e) whose sum is √N to about
# twice the precision of `T`: e = (N - c²)/(2c), whose numerator `fma` forms exactly.
@inline function ladder_pair(N::Int, ::Type{T}) where {T}
    c = ladder_coefficient(N, T)
    (c, N == 0 ? zero(c) : fma(-c, c, oftype(c, N)) / 2c)
end
# The difference b y₊ - a y₋ for the coefficient pairs a and b, with each coefficient to that
# precision.
@inline neighbor_difference(a, b, y₋, y₊) = fma(b[1], y₊, b[2] * y₊) - fma(a[1], y₋, a[2] * y₋)

# The derivative of an element x of a block of 𝔇 or of `d` along the generator of direction
# d of rotor iᵣ, from the left or from the right, given the neighbors x₋ and x₊ that it
# couples, and the coefficients k, a, and b of the index that it steps (see
# `step_coefficients`); and the addition of that element's part of the vector g into `Ḡ`,
# given ā, the conjugate of its cotangent, and p = ā a x₋ and q = ā b x₊.  Along a generator
# about y, v = (0, v_y, 0), the formulas of the note at the top of this file are v_y (b x₊ -
# a x₋) from the left and v_y (a x₋ - b x₊) from the right.  The callers have checked the
# shapes of the arrays, which are therefore indexed under `@inbounds`.
@inline function wigner_derivative(G::AbstractMatrix, d, iᵣ, left::Bool, k, a, b, x, x₋, x₊)
    w = @inbounds Complex(G[3d - 2, iᵣ], G[3d - 1, iᵣ])
    vz = @inbounds G[3d, iᵣ]
    if left
        -im * (vz * k * x + conj(w) * a * x₋ + w * b * x₊)
    else
        -im * (vz * k * x + w * a * x₋ + conj(w) * b * x₊)
    end
end
@inline function wigner_derivative(G::AngleGenerators, d, iᵣ, left::Bool, k, a, b, x, x₋, x₊)
    vy = @inbounds G.vy[d, iᵣ]
    δ = neighbor_difference(a, b, x₋, x₊)
    vy * (left ? δ : -δ)
end
@inline function add_wigner_cotangent!(Ḡ::AbstractMatrix, iᵣ, left::Bool, ā, k, x, p, q)
    @inbounds begin
        Ḡ[1, iᵣ] += imag(p + q)
        Ḡ[2, iᵣ] += left ? real(q - p) : real(p - q)
        Ḡ[3, iᵣ] += imag(ā * k * x)
    end
    nothing
end

# The generators of angles of which some are not finite, as at a pole for `d` of a rotor,
# where the angle has no derivative (see `rotation_angle_gradient`).  An element whose
# neighbors cancel exactly is then given the derivative zero: at a pole these are the
# elements whose neighbors both vanish, which are those whose derivative as a function of
# the rotor exists and is zero.  Generators that are all finite are not wrapped in this
# type, since its test of each element keeps the kernel from vectorizing.
# - `G` is the type of the generators.
struct PoleGenerators{G<:AngleGenerators}
    generators::G
end
all_finite(::AbstractMatrix) = true
all_finite(G::AngleGenerators) = all(isfinite, G.vy)
@inline step_coefficients(G::PoleGenerators, ℓ, n, ::Type{RT}) where {RT} =
    step_coefficients(G.generators, ℓ, n, RT)
@inline function wigner_derivative(G::PoleGenerators, d, iᵣ, left::Bool, k, a, b, x, x₋, x₊)
    ẋ = wigner_derivative(G.generators, d, iᵣ, left, k, a, b, x, x₋, x₊)
    finite = isfinite(@inbounds G.generators.vy[d, iᵣ])
    finite || !iszero(neighbor_difference(a, b, x₋, x₊)) ? ẋ : zero(ẋ)
end

# The derivative of an element y of a block of the harmonics, and the addition of its part
# of g, as for 𝔇 above; it is differentiated from the left, conjugated, and along a
# generator about y its derivative is v_y (b y₊ - a y₋).
@inline function harmonic_derivative(G::AbstractMatrix, d, iᵣ, k, a, b, y, y₋, y₊)
    w = @inbounds Complex(G[3d - 2, iᵣ], G[3d - 1, iᵣ])
    vz = @inbounds G[3d, iᵣ]
    im * (vz * k * y + w * a * y₋ + conj(w) * b * y₊)
end
@inline harmonic_derivative(G::AngleGenerators, d, iᵣ, k, a, b, y, y₋, y₊) =
    (@inbounds G.vy[d, iᵣ]) * neighbor_difference(a, b, y₋, y₊)
@inline function add_harmonic_cotangent!(Ḡ::AbstractMatrix, iᵣ, ȳ, k, y, p, q)
    @inbounds begin
        Ḡ[1, iᵣ] -= imag(p + q)
        Ḡ[2, iᵣ] += real(q - p)
        Ḡ[3, iᵣ] -= imag(ȳ * k * y)
    end
    nothing
end


## The kernels

# The blocks here are 3-dimensional arrays with the rotor index first, as the calculators
# store them: [iᵣ, m′, m] for 𝔇, with rows `rows` and columns `cols`, and [iᵣ, s, m] for
# the harmonics, whose third axis is the whole of -ℓ:ℓ.  For 𝔇, the derivatives of the
# elements in the rows `outrows` ⊆ `rows` and the columns `outcols` ⊆ `cols` are written
# into `Ȧ`, laid out as [iᵣ, outrows, outcols]; for the harmonics, those of every element
# are written into the same position of `Ȧ`.  Each element of `Ȧ` is `combine(x, ẋ)`, where
# `x` is the value and `ẋ` the tuple of its derivatives in the `N` directions whose
# generators are in `G`, so that a caller can assemble its own number type in the same pass.
# Every array is indexed under `@inbounds`, so its shape is checked first.

function check_wigner_block(Ȧ, A, ℓ, rows, cols, outrows, outcols, left::Bool, Nᵣ)
    Base.require_one_based_indexing(Ȧ, A)
    within(out, r) = isempty(out) || (first(r) ≤ first(out) && last(out) ≤ last(r))
    ok = (
        size(A, 1) == Nᵣ && size(A, 2) ≥ length(rows) && size(A, 3) ≥ length(cols)
        && size(Ȧ, 1) == Nᵣ && size(Ȧ, 2) ≥ length(outrows) && size(Ȧ, 3) ≥ length(outcols)
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

# `combine` is typed by a parameter in each method, so that the methods are compiled for it
# although they only pass it on: Julia does not specialize a method on a function that it
# merely passes to another, and calls it dynamically otherwise.
function wigner_block_pushforward!(
    combine::F, Ȧ::AbstractArray{<:Any, 3}, A::AbstractArray{<:Any, 3}, ℓ, rows, cols,
    outrows, outcols, left::Bool, G::Union{AbstractMatrix, AngleGenerators}, ::Val{N}
) where {F, N}
    Nᵣ = generator_columns(G)
    check_wigner_block(Ȧ, A, ℓ, rows, cols, outrows, outcols, left, Nᵣ)
    check_generators(G, N, "generators")
    # One rotor is differentiated by a method in which `Nᵣ` is `Val(1)`, as in
    # `materialize!`, and several by one in which it is the `Int` it is.
    if !all_finite(G)
        wigner_block_pushforward!(
            combine, Ȧ, A, ℓ, rows, cols, outrows, outcols, left, PoleGenerators(G), Val(N), Nᵣ
        )
    elseif Nᵣ == 1
        wigner_block_pushforward!(combine, Ȧ, A, ℓ, rows, cols, outrows, outcols, left, G, Val(N), Val(1))
    else
        wigner_block_pushforward!(combine, Ȧ, A, ℓ, rows, cols, outrows, outcols, left, G, Val(N), Nᵣ)
    end
    Ȧ
end
# The arrays are read and written by linear index, the element (iᵣ, j′, j) of `A` at iᵣ + s′
# (j′ - 1) + s (j - 1), and likewise for `Ȧ`.  The coefficients depend only on the index
# that the derivative steps, so the loop over that index is the outer one, and they are
# computed once for each of its values.  Each element is computed by the same expression
# whatever the order of the loops.
function wigner_block_pushforward!(
    combine::F, Ȧ, A, ℓ, rows, cols, outrows, outcols, left::Bool, G, ::Val{N},
    Nᵣ::Union{Val{1}, Int}
) where {F, N}
    RT = real(eltype(A))
    s′ = rotor_count(Nᵣ)
    s, ṡ = s′ * size(A, 2), s′ * size(Ȧ, 2)
    o′ = isempty(outrows) ? 0 : Int(first(outrows) - first(rows))
    o = isempty(outcols) ? 0 : Int(first(outcols) - first(cols))
    if left
        @inbounds for (j′ₒ, m′) ∈ enumerate(outrows)
            k, a, b = step_coefficients(G, ℓ, m′, RT)
            has₋, has₊ = m′ > -ℓ, m′ < ℓ
            i, iₒ = s′ * (j′ₒ + o′ - 1) + s * o, s′ * (j′ₒ - 1)
            for jₒ ∈ 1:length(outcols)
                element_pushforward!(
                    combine, Ȧ, iₒ + ṡ * (jₒ - 1), A, i + s * (jₒ - 1), s′, has₋, has₊, G,
                    true, k, a, b, Val(N), Nᵣ
                )
            end
        end
    else
        @inbounds for (jₒ, m) ∈ enumerate(outcols)
            k, a, b = step_coefficients(G, ℓ, m, RT)
            has₋, has₊ = m > -ℓ, m < ℓ
            i, iₒ = s′ * o′ + s * (jₒ + o - 1), ṡ * (jₒ - 1)
            for j′ₒ ∈ 1:length(outrows)
                element_pushforward!(
                    combine, Ȧ, iₒ + s′ * (j′ₒ - 1), A, i + s′ * (j′ₒ - 1), s, has₋, has₊, G,
                    false, k, a, b, Val(N), Nᵣ
                )
            end
        end
    end
    nothing
end
# The derivatives of one element for every rotor, written into `Ȧ` after its first `iₒ`
# entries, from the element after the first `i` entries of `A`, whose neighbors are `step`
# before and after it, where they exist.
@inline function element_pushforward!(
    combine::F, Ȧ, iₒ, A, i, step, has₋::Bool, has₊::Bool, G, left::Bool, k, a, b, ::Val{N},
    Nᵣ::Union{Val{1}, Int}
) where {F, N}
    @inbounds for iᵣ ∈ 1:rotor_count(Nᵣ)
        x = A[i + iᵣ]
        x₋ = has₋ ? A[i + iᵣ - step] : zero(x)
        x₊ = has₊ ? A[i + iᵣ + step] : zero(x)
        ẋ = ntuple(d -> wigner_derivative(G, d, iᵣ, left, k, a, b, x, x₋, x₊), Val(N))
        Ȧ[iₒ + iᵣ] = combine(x, ẋ)
    end
    nothing
end

# The vector g of the note above for each rotor, added into the columns of `Ḡ`, from the
# values `A` and the cotangent `Ā` of the rows `outrows` and columns `outcols`, laid out as
# `Ȧ` is above.
function wigner_block_pullback!(
    Ḡ::AbstractMatrix, A::AbstractArray{<:Any, 3},
    Ā::AbstractArray{<:Any, 3}, ℓ, rows, cols, outrows, outcols, left::Bool
)
    Nᵣ = generator_columns(Ḡ)
    check_wigner_block(Ā, A, ℓ, rows, cols, outrows, outcols, left, Nᵣ)
    check_generators(Ḡ, 1, "cotangents")
    RT = real(eltype(A))
    s′, s = size(A, 1), size(A, 1) * size(A, 2)
    ṡ′, ṡ = size(Ā, 1), size(Ā, 1) * size(Ā, 2)
    o′ = isempty(outrows) ? 0 : Int(first(outrows) - first(rows))
    o = isempty(outcols) ? 0 : Int(first(outcols) - first(cols))
    # The contributions are added in the order of the storage of the elements, the rotors
    # innermost.  The coefficients a_n of the rows are tabulated once, since b_n = a_{n+1};
    # those of a column are the same throughout it.
    (isempty(outrows) || isempty(outcols)) && return Ḡ
    out = left ? outrows : outcols
    ladder = [ladder_down(ℓ, n, RT) for n ∈ first(out):(last(out) + 1)]
    @inbounds for (jₒ, m) ∈ enumerate(outcols), (j′ₒ, m′) ∈ enumerate(outrows)
        p = left ? j′ₒ : jₒ
        n = left ? m′ : m
        k, a, b = Int(2n), ladder[p], ladder[p + 1]
        has₋, has₊ = n > -ℓ, n < ℓ
        i, iₒ = s′ * (j′ₒ + o′ - 1) + s * (jₒ + o - 1), ṡ′ * (j′ₒ - 1) + ṡ * (jₒ - 1)
        step = left ? s′ : s
        for iᵣ ∈ 1:Nᵣ
            ā = conj(Ā[iₒ + iᵣ])
            x = A[i + iᵣ]
            p₋ = has₋ ? ā * a * A[i + iᵣ - step] : zero(ā * x)
            q₊ = has₊ ? ā * b * A[i + iᵣ + step] : zero(ā * x)
            add_wigner_cotangent!(Ḡ, iᵣ, left, ā, k, x, p₋, q₊)
        end
    end
    Ḡ
end

# The same for a block of `d`, whose cotangents are those of the angles of its rotors (see
# `AngleCotangents`).  The cotangent of each angle is the sum of ā δ over the elements, with
# δ the difference that the derivative of the element is v_y times, formed as the forward
# kernel forms it.  The terms have either sign, and near the poles they nearly cancel, so
# each product is added with its rounding error, which `fma` gives exactly, into a
# compensated sum.  The elements of a line, along which the coefficients are constant, are
# added into two such sums in turn, which run in parallel.
function wigner_block_pullback!(
    Ḡ::AngleCotangents, A::AbstractArray{<:Real, 3}, Ā::AbstractArray{<:Real, 3}, ℓ, rows,
    cols, outrows, outcols, left::Bool
)
    Nᵣ = generator_columns(Ḡ)
    check_wigner_block(Ā, A, ℓ, rows, cols, outrows, outcols, left, Nᵣ)
    check_generators(Ḡ, 1, "cotangents")
    (isempty(outrows) || isempty(outcols)) && return Ḡ
    RT = eltype(A)
    s′, s, ṡ = Nᵣ, Nᵣ * size(A, 2), Nᵣ * size(Ā, 2)
    o′, o = Int(first(outrows) - first(rows)), Int(first(outcols) - first(cols))
    # The lines, the number of elements in each, the steps from one element of a line to the
    # next in `A` and in `Ā`, and the step to a neighbor in `A`
    lines, len = left ? (outrows, length(outcols)) : (outcols, length(outrows))
    step, stepₒ, neighbor = left ? (s, ṡ, s′) : (s′, s′, s)
    @inbounds for iᵣ ∈ 1:Nᵣ
        t₁ = t₂ = c₁ = c₂ = zero(RT)
        for (l, n) ∈ enumerate(lines)
            # The first element of the line, in `A` and in `Ā`
            i, iₒ = if left
                (iᵣ + s′ * (l + o′ - 1) + s * o, iᵣ + s′ * (l - 1))
            else
                (iᵣ + s′ * o′ + s * (l + o - 1), iᵣ + ṡ * (l - 1))
            end
            _, a, b = step_coefficients(Ḡ, ℓ, n, RT)
            # δ from the right is the negative of that from the left
            a, b = left ? (a, b) : (.-a, .-b)
            has₋, has₊ = n > -ℓ, n < ℓ
            u = 0
            while u < len
                p, e = angle_term(Ā, A, iₒ + stepₒ * u, i + step * u, neighbor, has₋, has₊, a, b)
                t₁, c₁ = compensated_sum(t₁, c₁ + e, p)
                if u + 1 < len
                    p, e = angle_term(
                        Ā, A, iₒ + stepₒ * (u + 1), i + step * (u + 1), neighbor, has₋, has₊, a, b
                    )
                    t₂, c₂ = compensated_sum(t₂, c₂ + e, p)
                end
                u += 2
            end
        end
        t, e = two_sum(t₁, t₂)
        compensated_add!(Ḡ, iᵣ, t, (c₁ + c₂) + e)
    end
    Ḡ
end
# The term ā δ of the element of `A` at index `i`, whose cotangent is at index `iₒ` of `Ā`,
# and the rounding error of that product.
@inline function angle_term(Ā, A, iₒ, i, neighbor, has₋::Bool, has₊::Bool, a, b)
    @inbounds begin
        ā = Ā[iₒ]
        x₋ = has₋ ? A[i - neighbor] : zero(ā)
        x₊ = has₊ ? A[i + neighbor] : zero(ā)
    end
    δ = neighbor_difference(a, b, x₋, x₊)
    p = ā * δ
    (p, fma(ā, δ, -p))
end

function check_harmonic_block(Ȧ, A, ℓ, Nᵣ)
    Base.require_one_based_indexing(Ȧ, A)
    n = Int(2ℓ) + 1
    if !(
        size(A, 1) == Nᵣ && size(A, 3) ≥ n && size(Ȧ, 1) == Nᵣ && size(Ȧ, 2) ≥ size(A, 2)
        && size(Ȧ, 3) ≥ n
    )
        throw(DimensionMismatch(
            "Cannot differentiate a block of size $(size(A)) at ℓ=$ℓ into an array of size "
            * "$(size(Ȧ)) for Nᵣ=$Nᵣ."
        ))
    end
    nothing
end

function harmonic_block_pushforward!(
    combine::F, Ȧ::AbstractArray{<:Any, 3}, A::AbstractArray{<:Any, 3}, ℓ,
    G::Union{AbstractMatrix, AngleGenerators}, ::Val{N}
) where {F, N}
    Nᵣ = generator_columns(G)
    check_harmonic_block(Ȧ, A, ℓ, Nᵣ)
    check_generators(G, N, "generators")
    if Nᵣ == 1
        harmonic_block_pushforward!(combine, Ȧ, A, ℓ, G, Val(N), Val(1))
    else
        harmonic_block_pushforward!(combine, Ȧ, A, ℓ, G, Val(N), Nᵣ)
    end
    Ȧ
end
# As for 𝔇, the arrays are read and written by linear index, and the coefficients, which
# depend on m alone, are computed once for each column.
function harmonic_block_pushforward!(
    combine::F, Ȧ, A, ℓ, G, ::Val{N}, Nᵣ::Union{Val{1}, Int}
) where {F, N}
    RT = real(eltype(A))
    n = Int(2ℓ) + 1
    s′ = rotor_count(Nᵣ)
    s, ṡ = s′ * size(A, 2), s′ * size(Ȧ, 2)
    @inbounds for j ∈ 1:n
        m = -ℓ + (j - 1)
        k, a, b = step_coefficients(G, ℓ, m, RT)
        i, iₒ = s * (j - 1), ṡ * (j - 1)
        for p ∈ 1:size(A, 2)
            element_pushforward!(
                combine, Ȧ, iₒ + s′ * (p - 1), A, i + s′ * (p - 1), s, j > 1, j < n, G,
                k, a, b, Val(N), Nᵣ
            )
        end
    end
    nothing
end
@inline function element_pushforward!(
    combine::F, Ȧ, iₒ, A, i, step, has₋::Bool, has₊::Bool, G, k, a, b, ::Val{N},
    Nᵣ::Union{Val{1}, Int}
) where {F, N}
    @inbounds for iᵣ ∈ 1:rotor_count(Nᵣ)
        y = A[i + iᵣ]
        y₋ = has₋ ? A[i + iᵣ - step] : zero(y)
        y₊ = has₊ ? A[i + iᵣ + step] : zero(y)
        ẏ = ntuple(d -> harmonic_derivative(G, d, iᵣ, k, a, b, y, y₋, y₊), Val(N))
        Ȧ[iₒ + iᵣ] = combine(y, ẏ)
    end
    nothing
end

function harmonic_block_pullback!(
    Ḡ::AbstractMatrix, A::AbstractArray{<:Any, 3},
    Ā::AbstractArray{<:Any, 3}, ℓ
)
    Nᵣ = generator_columns(Ḡ)
    check_harmonic_block(Ā, A, ℓ, Nᵣ)
    check_generators(Ḡ, 1, "cotangents")
    RT = real(eltype(A))
    n = Int(2ℓ) + 1
    s′, s = size(A, 1), size(A, 1) * size(A, 2)
    ṡ′, ṡ = size(Ā, 1), size(Ā, 1) * size(Ā, 2)
    @inbounds for j ∈ 1:n
        m = -ℓ + (j - 1)
        k, a, b = Int(2m), ladder_down(ℓ, m, RT), ladder_up(ℓ, m, RT)
        for i ∈ axes(A, 2)
            o, oₒ = s′ * (i - 1) + s * (j - 1), ṡ′ * (i - 1) + ṡ * (j - 1)
            for iᵣ ∈ 1:Nᵣ
                ȳ = conj(Ā[oₒ + iᵣ])
                y = A[o + iᵣ]
                p = j > 1 ? ȳ * a * A[o + iᵣ - s] : zero(ȳ * y)
                q = j < n ? ȳ * b * A[o + iᵣ + s] : zero(ȳ * y)
                add_harmonic_cotangent!(Ḡ, iᵣ, ȳ, k, y, p, q)
            end
        end
    end
    Ḡ
end

# The same for a block of ₛλₗₘ.  Since b_m = a_{m+1}, the neighbors m and m+1 of a spin row
# contribute b_m (ȳ_m y_{m+1} - ȳ_{m+1} y_m) together.  Near the poles this difference
# nearly cancels, so it is formed exactly, as the sum of a number and its rounding error,
# multiplied by b_m to twice the working precision, and added with the rounding errors of
# these steps into a compensated sum.
function harmonic_block_pullback!(
    Ḡ::AngleCotangents, A::AbstractArray{<:Real, 3}, Ā::AbstractArray{<:Real, 3}, ℓ
)
    Nᵣ = generator_columns(Ḡ)
    check_harmonic_block(Ā, A, ℓ, Nᵣ)
    check_generators(Ḡ, 1, "cotangents")
    RT = eltype(A)
    n = Int(2ℓ) + 1
    b = [ladder_pair(Int(ℓ + m + 1) * Int(ℓ - m), RT) for m ∈ -ℓ:(ℓ - 1)]
    s, ṡ = Nᵣ * size(A, 2), Nᵣ * size(Ā, 2)
    @inbounds for iᵣ ∈ 1:Nᵣ, row ∈ axes(A, 2)
        i = iᵣ + Nᵣ * (row - 1)
        t = c = zero(RT)
        for j ∈ 1:(n - 1)
            ȳ, ȳ₊ = Ā[i + ṡ * (j - 1)], Ā[i + ṡ * j]
            y, y₊ = A[i + s * (j - 1)], A[i + s * j]
            p = ȳ * y₊
            q = ȳ₊ * y
            d, e = two_sum(p, -q)
            e += fma(ȳ, y₊, -p) - fma(ȳ₊, y, -q)
            x = b[j][1] * d
            t, c = compensated_sum(t, c + (fma(b[j][1], d, -x) + (b[j][1] * e + b[j][2] * d)), x)
        end
        compensated_add!(Ḡ, iᵣ, t, c)
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
    check_harmonic_array(Y₃, ℓₘᵢₙ, ℓₘₐₓ, generator_columns(G))
    CT = derivative_type(eltype(Y), G)
    Ẏ = similar(Y, typeof(combine(zero(eltype(Y)), ntuple(_ -> zero(CT), Val(N)))))
    Ẏ₃ = harmonic_blocks_view(Ẏ, batched)
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        r = Yindex(ℓ, -ℓ, ℓₘᵢₙ):Yindex(ℓ, ℓ, ℓₘᵢₙ)
        harmonic_block_pushforward!(
            combine, view(Ẏ₃, :, :, r), view(Y₃, :, :, r), ℓ, G, Val(N)
        )
    end
    Ẏ
end

function harmonic_array_pullback!(
    Ḡ::Union{AbstractMatrix, AngleCotangents}, Y::AbstractArray, Ȳ::AbstractArray, batched::Bool, ℓₘᵢₙ::IT, ℓₘₐₓ::IT
) where {IT<:IntegerHalf}
    size(Ȳ) == size(Y) || throw(DimensionMismatch(
        "The cotangent has size $(size(Ȳ)), but the harmonics $(size(Y))."
    ))
    Y₃, Ȳ₃ = harmonic_blocks_view(Y, batched), harmonic_blocks_view(Ȳ, batched)
    check_harmonic_array(Y₃, ℓₘᵢₙ, ℓₘₐₓ, generator_columns(Ḡ))
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        r = Yindex(ℓ, -ℓ, ℓₘᵢₙ):Yindex(ℓ, ℓ, ℓₘᵢₙ)
        harmonic_block_pullback!(Ḡ, view(Y₃, :, :, r), view(Ȳ₃, :, :, r), ℓ)
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
# `src/derivatives/lifting.jl`), and after each of that calculator's steps it writes into
# its own block each value together with its derivatives, from the generators of its rotors'
# tangents.  Those generators depend only on the rotors, so they are computed when the
# rotors are set, and a step allocates nothing.  This is how forward-mode numbers are
# lifted; an extension for a reverse-mode tool defines methods of `set_generators!` and
# `lift!` for its own numbers, which record the step instead.

set_generators!(lift::Lift, left::Bool, rotors::AbstractVector{Quaternion{RT}}) where {RT} =
    set_generators!(lift.G, left, rotors, rotor_value, rotor_tangents, Val(ndirections(RT)))
# The generators of the angles of a calculator of `d` or of ₛλₗₘ, from either side: the y
# component β̇/2 of each, in each direction (see `AngleGenerators`).
function set_generators!(lift::Lift, ::Bool, angles::AbstractVector{RT}) where {RT<:Real}
    G = lift.G
    N = ndirections(RT)
    Base.require_one_based_indexing(G, angles)
    size(G) == (N, length(angles)) || throw(DimensionMismatch(
        "The generators of $(length(angles)) angles in $N directions need a matrix of size "
        * "$((N, length(angles))), not $(size(G))."
    ))
    @inbounds for iᵣ ∈ eachindex(angles)
        β̇ = angle_tangents(angles[iᵣ])
        for d ∈ 1:N
            G[d, iᵣ] = β̇[d] / 2
        end
    end
    lift
end

# The generators of a lifting calculator's rotors, as the kernels read them.
lift_generators(
    c::Union{WignerCalculator{IT, RT, NT}, HarmonicCalculator{IT, RT, NT}}
) where {IT, RT, NT} = NT <: Complex ? c.lift.G : AngleGenerators(c.lift.G)

# The block of degree ℓ, written into `A` after its first `o` entries, from the block of the
# calculator of values, whose limits are this calculator's derivative limits, so that it
# holds one row or column more on each side wherever this one needs it; see
# `derivative_limits`.
function lift!(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, <:Lift}, ℓ::IT, A, o::Int
) where {IT, RT, NT, ST, B, FT<:Real}
    inner = c.lift.inner
    wigner_block_pushforward!(
        lift_combine(RT), block_array(c, A, ℓ, o), block_array(inner, inner.Wˡ, ℓ), ℓ,
        m′range(inner, ℓ), mrange(inner, ℓ), m′range(c, ℓ), mrange(c, ℓ),
        derivatives_from_left(c), lift_generators(c), Val(ndirections(RT))
    )
    c
end
# The values from which the derivatives of the block of degree ℓ are computed, as for a
# calculator of values (see `derivative_values` in `src/calculators/wigner.jl`).  Where they
# reach beyond the block, the wedge alone cannot give them, since it holds the values of the
# rotors and not their derivatives, so they are lifted, as the block itself is, from the
# values of the calculator of values over its own derivative ranges, which reach one row or
# column further.
function derivative_values(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, <:Lift}, ℓ::IT, A::AbstractArray=c.Wˡ,
    o::Int=0
) where {IT, RT, NT, ST, B, FT<:Real}
    rows, cols = derivative_m′range(c, ℓ), derivative_mrange(c, ℓ)
    if rows == m′range(c, ℓ) && cols == mrange(c, ℓ)
        block_array(c, A, ℓ, o)
    else
        inner = c.lift.inner
        values = Array{NT, 3}(undef, Nᵣ(c), length(rows), length(cols))
        wigner_block_pushforward!(
            lift_combine(RT), values, derivative_values(inner, ℓ), ℓ,
            derivative_m′range(inner, ℓ), derivative_mrange(inner, ℓ), rows, cols,
            derivatives_from_left(c), lift_generators(c), Val(ndirections(RT))
        )
        values
    end
end
# The spin rows `is` of the block of degree ℓ, written into `A` after its first `o` entries,
# from the calculator of values, which computes them into its own block.
function lift!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, <:Lift}, ℓ::IT, is, A, o::Int
) where {IT, RT, NT, ST, S, B, FT<:Real}
    inner = c.lift.inner
    harmonic_block_pushforward!(
        lift_combine(RT), block_array(c, A, ℓ, is, o), block_array(inner, inner.Yˡ, ℓ, is),
        ℓ, lift_generators(c), Val(ndirections(RT))
    )
    c
end


## Whole arrays, for the rules of the tools that cannot follow a calculator
#
# ChainRules and ReverseDiff cannot differentiate code that mutates arrays, as a calculator
# does, so their extensions define rules for the functions that return whole arrays instead:
# `D_array`, `d_array`, `sYlm_array`, and `sYlm_matrix_array`.  These are what those rules
# compute.

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

# The gradient of the angle β of the rotor data `x` (see `rotation_angle`) with respect to
# its components: 1 for an angle; for a phase z, the complex number g with β̇ = Re(conj(g) ż),
# which is i z / |z|²; and for a rotor, the tuple of the four partial derivatives of
# 2 atan(√(X²+Y²), √(W²+Z²)), which are `NaN` at the poles, where β has no derivative.
rotation_angle_gradient(::Real) = 1
rotation_angle_gradient(z::Complex) = im * z / abs2(z)
function rotation_angle_gradient(R::RotorLike)
    W, X, Y, Z = R[1], R[2], R[3], R[4]
    a, b = W^2 + Z^2, X^2 + Y^2
    c, s = 2sqrt(b) / (sqrt(a) * (a + b)), 2sqrt(a) / (sqrt(b) * (a + b))
    (-c * W, s * X, s * Y, -c * Z)
end

# The same gradient with respect to the components of the rotor data (see
# `rotor_data_components`), as a tuple of reals.
angle_component_gradient(β::Real) = (one(β),)
function angle_component_gradient(x::Real, y::Real)
    g = rotation_angle_gradient(Complex(x, y))
    (real(g), imag(g))
end
angle_component_gradient(w::Real, x::Real, y::Real, z::Real) =
    rotation_angle_gradient(Quaternion(w, x, y, z))

# The cotangent of the rotor data `x` from the cotangent β̄ of its angle: a number for an
# angle, a complex number for a phase, as for `rotation_angle_gradient`, and a tuple of four
# components for a rotor.  A zero β̄ gives zero even at the poles, where the gradient is not
# finite, so that the elements whose derivative there is zero (see `wigner_derivative`) give
# zero cotangents too.
function rotation_angle_cotangent(x, β̄)
    x̄ = rotation_angle_gradient(x) .* β̄
    iszero(β̄) ? map(zero, x̄) : x̄
end

# The blocks of the arrays of `D` or `d` of the rotor data `x`, with elements of type `NT`,
# as `D_array` and `d_array` compute them, and the values from which the derivatives of each
# are computed, as [1, m′, m], with the calculator that computed them, from which the rules
# read the ranges of each block.  The values are copied, since they may be the blocks
# themselves, which are handed to the caller.
function wigner_arrays_with_derivative_values(
    ::Type{NT}, x, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {NT, IT<:IntegerHalf}
    calc = series_calculator(NT, x, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    values = Array{eltype(calc.Wˡ), 3}[]
    blocks = wigner_arrays(calc) do ℓ, A, o
        push!(values, copy(derivative_values(calc, ℓ, A, o)))
    end
    blocks, values, calc
end

# The derivatives of those blocks along the generators `G` of one direction, as the matrices
# of the blocks themselves, from `values`, the values from which the derivatives of each
# block are computed.
function wigner_arrays_pushforward(calc::WignerCalculator{IT}, values, G) where {IT}
    left = derivatives_from_left(calc)
    map(enumerate(lowest_index(IT):ℓₘₐₓ(calc))) do (i, ℓ)
        rows, cols = derivative_m′range(calc, ℓ), derivative_mrange(calc, ℓ)
        outrows, outcols = m′range(calc, ℓ), mrange(calc, ℓ)
        Aᵈ = values[i]
        Ȧ = similar(Aᵈ, derivative_type(eltype(Aᵈ), G), 1, length(outrows), length(outcols))
        wigner_block_pushforward!(
            (x, ẋ) -> only(ẋ), Ȧ, Aᵈ, ℓ, rows, cols, outrows, outcols, left, G, Val(1)
        )
        reshape(Ȧ, length(outrows), length(outcols))
    end
end

# The cotangents of the generators, from the values from which the derivatives of each
# block are computed and the cotangents `Ā` of the blocks, of which any may be `nothing` for
# a zero cotangent: the vector g of the note at the top of this file for 𝔇, as a 3×1 matrix,
# and its y component for `d`, as `AngleCotangents`.
function wigner_arrays_pullback(calc::WignerCalculator{IT}, values, Ā) where {IT}
    Ḡ = zero_cotangents(eltype(first(values)), 1)
    left = derivatives_from_left(calc)
    for (i, ℓ) ∈ enumerate(lowest_index(IT):ℓₘₐₓ(calc))
        Āᵢ = Ā[i]
        Āᵢ === nothing && continue
        rows, cols = derivative_m′range(calc, ℓ), derivative_mrange(calc, ℓ)
        outrows, outcols = m′range(calc, ℓ), mrange(calc, ℓ)
        wigner_block_pullback!(
            Ḡ, values[i], reshape(Āᵢ, 1, length(outrows), length(outcols)), ℓ,
            rows, cols, outrows, outcols, left
        )
    end
    Ḡ
end
