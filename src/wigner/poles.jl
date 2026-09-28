# Rotors at and near a pole
#
# The recurrence sees a rotor through `spinor_phases`, which splits it into the half angles
# cos(β/2) and sin(β/2) and the phases z₊ = e^{i(α+γ)/2} and z₋ = e^{i(α-γ)/2}.  At the
# poles, β = 0 and β = π, that split is singular: one of X² + Y² and W² + Z² is exactly zero
# there, and its square root has an infinite derivative, so that a derivative taken through
# it by automatic differentiation is NaN, and so is every derivative of 𝔇 and of the
# harmonics computed from it — although both are smooth functions of the rotor at the poles.
# Near a pole the split is merely ill-conditioned, but that is almost as bad: at a distance
# r from it (r = sin(β/2) from β = 0, or cos(β/2) from β = π), the k-th derivatives that the
# recurrence gives are wrong by about ε r⁻ᵏ relative to their size, for machine epsilon ε,
# because they are formed from factors whose own derivatives are of order r⁻ᵏ and cancel.
# The values are unaffected.
#
# The smooth quantities are the products σ = cos(β/2) z₊ = (W + iZ)/‖R‖ and ρ = sin(β/2) z₋
# = (Y - iX)/‖R‖, which are the rotor's normalized Cayley–Klein parameters, and 𝔇 is a
# polynomial in them and their conjugates:
#
#     𝔇^ℓ_{m′,m} = Σₛ (-1)^{k+s} Cₛ σ̄^{ℓ+m-s} σ^{ℓ-m′-s} ρ̄^{k+s} ρ^s,        k = m′ - m,
#
#     Cₛ² = binomial(ℓ+m, s) binomial(ℓ-m′, s) binomial(ℓ+m′, k+s) binomial(ℓ-m, k+s),
#
# with the sum over every s for which the exponents are non-negative.  This is Wigner's
# formula for d, with the phases of the convention absorbed into σ and ρ, and
# `test/wigner/poles.jl` compares it with the recurrence.  All of the exponents are integers
# for half-odd indices too.  The full sum is no use in general — it cancels badly away from
# the poles, and its coefficients overflow for large ℓ — but near a pole it is exactly what
# is needed.
#
# At a pole one of σ and ρ — call it ζ — vanishes, and a term of degree n in ζ and ζ̄ then
# vanishes together with all of its derivatives of order less than n.  So every derivative
# of order N or less at the pole is given exactly by the terms of degree at most N, of which
# there are at most ⌊N/2⌋ + 1 in each element, and only in the elements with |m′ - m| ≤ N at
# β = 0 or |m′ + m| ≤ N at β = π; every other element vanishes to that order.  Near the pole
# those terms are accurate to the size of the first term left out, which is about ((ℓ+1)
# r)^{N+1-k} relative to the size of a k-th derivative, and their coefficients are only of
# order ℓᴺ.  The factor that does not vanish is written as its modulus κ times a phase,
# which is the power of z₊ or z₋ that `materialize!` would apply to the same element, read
# from the same table.
#
# So every rotor within the radius r_s = ε^{1/(N+1)}/(ℓₘₐₓ+1) of a pole is evaluated from
# those terms, with N = `pole_order`, in place of the recurrence.  That is the largest
# radius at which the terms left out cannot change the values, and beyond it the
# recurrence's derivatives are wrong by no more than about ε^{1-k/(N+1)} (ℓₘₐₓ+1)^k relative
# to their size.  In `Float64`, for ℓₘₐₓ = 32, the errors measured just outside the radius
# are of order 1e-13 for first derivatives and 1e-10 for second ones, depending on the
# direction, where the expansion is within a few ulps.  (A radius of √ε, for instance, would
# leave the recurrence's second derivatives just outside it wrong by 10⁻³.)  The radius is
# set by ℓₘₐₓ rather than by each ℓ so that a rotor is treated the same way at every ℓ; the
# terms left out are only smaller at the lower ℓ.
#
# The test is applied to every number type alike.  An automatic-differentiation tool that
# works on the plain floating-point code, such as Enzyme, cannot be told apart from an
# ordinary evaluation by its type, and would otherwise differentiate the singular split.
# The test must also look at the value alone, and not at any derivatives that a number
# carries with it: X² + Y² is zero at the north pole, but its second derivatives are not,
# and since version 1.0 ForwardDiff's `iszero` and `==` compare those too.  So it is written
# as a strict comparison with a positive threshold, which ForwardDiff decides by the values
# whenever they differ, and which a reverse-mode tool does not record.
#
# A rotor that is exactly at a pole — which is taken to mean that X² + Y² or W² + Z² is
# below `floatmin(Float64)` relative to ‖R‖², as nothing but zero is in `Float32`, and
# nothing but zero or a subnormal in `Float64` — also needs the engine to be given the exact
# values e^{iβ} = ±1 and the half angles, without the square root of zero that
# `spinor_phases` would take.  The engine's results for that rotor are overwritten, but a
# reverse-mode tool, such as ReverseDiff, runs its reverse pass through every operation it
# recorded, including those whose results are later overwritten, and would otherwise find 0
# × ∞ there.  A rotor near a pole but not at it is given to the engine as usual.
#
# None of this applies to `d` and `H` of a rotor, which see it only through β.  As a
# function of the rotor, β has a cone-shaped singularity at each pole, so some of their
# elements — those with |m′ - m| = 1 at β = 0, which vary as sin(β/2), for example —
# actually have no derivative there, and `spinor_phases` gives them the NaN that says so.

# The number of orders of the expansion about a pole that are kept.  Every derivative up to
# this order is exact at a rotor exactly at a pole, while a derivative of higher order there
# would be wrong, and silently so.  Eight is far beyond what is needed at the pole itself —
# a Hessian is of order two — but it is also what sets the radius `pole_radius`, and a
# larger order widens the neighbourhood in which the recurrence's inaccurate derivatives are
# replaced.  It costs almost nothing: at most seventeen elements of a column are nonzero,
# with at most five terms each.
const pole_order = 8

# The radius r_s = ε^{1/(N+1)}/(ℓₘₐₓ+1) described above, for a calculator working in type
# `RT`.
@inline function pole_radius(::Type{RT}, ℓₘₐₓ) where {RT<:Real}
    eps(RT)^(1 / (pole_order + 1)) / (Float64(ℓₘₐₓ) + 1)
end

# A rotor at or near a pole, as a calculator of 𝔇 or of the harmonics records it: its index
# among the calculator's rotors, which pole it is near, the parameter ζ that vanishes at
# that pole, and the modulus κ of the other.
struct PoleRotor{RT<:Real}
    iᵣ::Int
    north::Bool  # near β = 0, where ζ = ρ; otherwise near β = π, where ζ = σ
    ζ::Complex{RT}
    κ::RT
end

# Every recorded rotor must be one of the calculator's `n`, since `materialize!` writes its
# values at that index under `@inbounds`.  The calculators' constructors check this, as they
# check the sizes of their buffers.
function check_pole_records(poles, n::Int)
    for p ∈ poles
        if !(1 ≤ p.iᵣ ≤ n)
            throw(DimensionMismatch(
                "A rotor near a pole is recorded as rotor $(p.iᵣ), but the calculator has only "
                * "Nᵣ=$n rotors."
            ))
        end
    end
    nothing
end

# The data of the expansion about one pole, for a rotor near it: about the north pole, ζ =
# ρ, κ = cos(β/2) and the phase z₊; about the south pole, ζ = σ, κ = sin(β/2) and the phase
# z₋.  The only square roots taken are of W² + Z² about the north pole and of X² + Y² about
# the south, each of which is close to ‖R‖² there, so all three are smooth near that pole.
# They are the expressions `spinor_phases` uses for κ and the phase, and give the same
# values.
@inline function pole_data(R::AbstractQuaternion, ::Type{RT}, north::Bool) where {RT<:Real}
    W, X, Y, Z = RT(R[1]), RT(R[2]), RT(R[3]), RT(R[4])
    a = W^2 + Z^2
    b = X^2 + Y^2
    nrm = √(a + b)
    σ = Complex{RT}(W, Z)
    ρ = Complex{RT}(Y, -X)
    if north
        let sqrta = √a
            (ρ / nrm, sqrta / nrm, σ / sqrta)
        end
    else
        let sqrtb = √b
            (σ / nrm, sqrtb / nrm, ρ / sqrtb)
        end
    end
end

# Store the rotor `R` as rotor `i` of the engine `w` of a calculator of 𝔇 or of the
# harmonics, and return the phases z₊ and z₋ for its power tables, as `store_rotor!` does.
# A rotor within the radius `r` of a pole is also recorded in `poles`.  A rotor exactly at a
# pole is given to the engine as the exact values e^{iβ} = ±1 and the half angles there,
# with the phase that is undefined at that pole set to 1, as `spinor_phases` sets it, and
# the other computed by `pole_data`, which takes no square root of zero.  The recurrence
# then runs on exact constants for that rotor, and `materialize!` overwrites whatever it
# gives.
@inline function store_complex_rotor!(
    w::HCalculator{IT, RT}, poles::Vector{PoleRotor{RT}}, i::Int, R, r::Real
) where {IT, RT}
    a = RT(R[1])^2 + RT(R[4])^2
    b = RT(R[2])^2 + RT(R[3])^2
    n² = a + b
    north = b < r^2 * n²  # see the note on the test for a pole, above
    if north || a < r^2 * n²
        ζ, κ, z = pole_data(R, RT, north)
        push!(poles, PoleRotor{RT}(i, north, ζ, κ))
        if (north ? b : a) < floatmin(Float64) * n²
            let o = one(RT), n = zero(RT)
                @inbounds w.eⁱᵝ[i] = Complex{RT}(north ? o : -o, n)
                north ? set_half_angles!(w, i, o, n) : set_half_angles!(w, i, n, o)
                return north ? (z, one(Complex{RT})) : (one(Complex{RT}), z)
            end
        end
    end
    store_rotor!(w, i, R)
end

# xⁿ for a non-negative integer n, by repeated squaring.  This is `Base.power_by_squaring`
# without its conversion of `x` to the type of `x * x`, which a reverse-mode wrapper type
# such as ReverseDiff's may not survive; and it is used rather than `^` because
# `^(::Complex, ::Integer)` goes through a polar form for every element type but the floats
# and integers, which gives wrong second derivatives for `Complex{<:ForwardDiff.Dual}`.
@inline function pole_power(x, n::Int)
    y = one(x)
    while n > 0
        if isodd(n)
            y *= x
        end
        n >>= 1
        if n > 0
            x *= x
        end
    end
    y
end

# √binomial(n, j) in type `RT`.  The binomial is a product of min(j, n - j) ratios, each of
# which leaves an integer, exactly so long as it fits in `RT`.  Every binomial of the
# expansion about a pole has a small lower index, or one close to n, so this is a short
# loop.
@inline function pole_sqrtbinomial(n::Int, j::Int, ::Type{RT}) where {RT}
    let j = min(j, n - j)
        b = one(RT)
        for i ∈ 1:j
            b = b * (n - j + i) / i
        end
        √b
    end
end

# The element 𝔇^ℓ_{m′,m} from the terms of degree at most `N` in ζ of the expansion above,
# about the north pole (β = 0, ζ = ρ) or the south (β = π, ζ = σ), given the modulus κ of
# the other parameter and the `phase` conj(z₊^{m′+m}) or conj(z₋^{m′-m}) respectively.  The
# expansion is exact where ζ vanishes, and near the pole it is accurate to the size of the
# first term left out.
#
# The three need not be of one type.  ReverseDiff's number type records where a value came
# from as a type parameter, so that the conjugate of an entry of a power table, say, is not
# of the entry's own type; the arithmetic below promotes them as it goes, and `materialize!`
# converts the result into the calculator's storage.
function pole_element(
    ℓ::IT, m′::IT, m::IT, north::Bool, ζ::Complex, κ::Real, phase::Complex, N::Int=pole_order
) where {IT<:IntegerHalf}
    RT = typeof(κ)
    A, B, C, E = Int(ℓ + m), Int(ℓ - m), Int(ℓ + m′), Int(ℓ - m′)
    k = Int(m′ - m)
    # Elements outside the band are of degree N + 1 or more.  (The range of s below would be
    # empty for them anyway.)
    if north ? abs(k) > N : abs(Int(m′ + m)) > N
        return zero(Complex{RT})
    end
    # The degree of a term in ζ is k + 2s about the north pole and A + E - 2s about the south.
    slo, shi = max(0, -k), min(A, E)
    if north
        shi = min(shi, fld(N - k, 2))
    else
        slo = max(slo, cld(A + E - N, 2))
    end
    ζ̄ = conj(ζ)
    total = zero(Complex{RT})
    for s ∈ slo:shi
        # The power of ζ comes first, and the square roots of the four binomials — each
        # taken with its smaller lower index: s and k + s about the north pole, and A - s
        # and E - s about the south — are multiplied into it one at a time.  Near the pole
        # the product grows only as ((ℓ+1) r)ⁿ, and at the pole it stays exactly zero, so
        # that neither the coefficient nor any partial product of it can overflow, even in
        # `Float32` at ℓ in the hundreds of thousands.
        term = if north
            pole_power(ζ̄, k + s) * pole_power(ζ, s)
        else
            pole_power(ζ̄, A - s) * pole_power(ζ, E - s)
        end
        term *= pole_sqrtbinomial(A, s, RT)
        term *= pole_sqrtbinomial(E, s, RT)
        term *= pole_sqrtbinomial(C, k + s, RT)
        term *= pole_sqrtbinomial(B, k + s, RT)
        term *= north ? pole_power(κ, A + E - 2s) : pole_power(κ, k + 2s)
        total += iseven(k + s) ? term : -term
    end
    total * phase
end
