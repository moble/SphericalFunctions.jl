md"""
# Boyle (2016)

!!! info "Summary"
    The Wigner ``𝔇`` matrices of [Boyle (2016)](@cite Boyle_2016) — which were the
    convention of this package before version 3.0 — are the complex conjugates of the ones
    now used in the `SphericalFunctions` package, for integer *and* half-integer indices:
    ``𝔇^{(ℓ)}_{m',m}(𝐑)|_{2016} = \overline{𝔇^{(ℓ)}_{m',m}(𝐑)}``.  Because that paper
    defines the spin-weighted spherical harmonics as ``(-1)^s \sqrt{(2ℓ+1)/4π}\,
    𝔇^{(ℓ)}_{m,-s}(𝐑)|_{2016}``, its ``{}_sY_{ℓ,m}`` agree with ours.

[Boyle_2016](@citet) argued that spin-weighted spherical functions should be defined as
functions on the spin group ``\mathrm{Spin}(3)``, represented by unit quaternions, rather
than on coordinates of the 2-sphere.  That paper is the origin of most of the conventions in
this package, with one important exception: in the interim, the convention for the ``𝔇``
matrices themselves has been changed to the complex conjugate — so that they agree with
LALSuite, Wikipedia, Sakurai, and the other sources in our priority list — as explained on
the [conventions pages](@ref summary_wigner_D).

The paper defines the spin-weighted spherical harmonics as functions of a unit quaternion
``𝐑`` by [Eq. (21)]
```math
{}_sY_{ℓ,m}(𝐑) = (-1)^s \sqrt{\frac{2ℓ+1}{4π}}\; 𝔇^{(ℓ)}_{m,-s}(𝐑),
```
and gives an explicit expression for ``𝔇`` directly in terms of the components of the
quaternion [Eq. (35)].  Writing ``𝐑 = R_s + R_a`` with ``R_s = W + Z𝐤`` and ``R_a = Y𝐣 +
X𝐢`` (the parts that commute and anticommute with ``𝐤``), and treating these as complex
numbers ``R_s = r_s e^{iϕ_s}`` and ``R_a = r_a e^{iϕ_a}``, the expression is a sum over
``ρ`` of terms proportional to ``r_s^{2ℓ-(m-m')-2ρ}\, r_a^{(m-m')+2ρ}\, e^{i[(m+m')ϕ_s +
(m-m')ϕ_a]}``, with an alternative form when ``r_a > r_s``; Appendix A of the paper
describes how to evaluate the sum stably, and that algorithm is transcribed below as
`WignerDElement`.  Since ``e^{iϕ_s}`` and ``e^{iϕ_a}`` are the half-angle phases
``e^{i(α+γ)/2}`` and ``e^{i(α-γ)/2}``, the ``e^{+i(m'α + mγ)}`` dependence is the complex
conjugate of [ours](@ref summary_wigner_D).

The same paper defines left and right operators ``L`` and ``K`` [Eqs. (42)–(43)] with
``K_z\, {}_sY_{ℓ,m} = -s\, {}_sY_{ℓ,m}`` [Eq. (47)], and identifies ``\eth = -K_-`` and
``\bar{\eth} = K_+`` [Eq. (46)]; in the present conventions ``K = -R``, so this is the
statement ``R_z\, {}_sY_{ℓ,m} = s\, {}_sY_{ℓ,m}``, ``\eth = R_+``, and ``\bar{\eth} =
-R_-`` of the [summary page](@ref summary_spin_weight).

## Implementing formulas

Because the `WignerDElement` function is also needed on the [Varshalovich page](@ref
"Varshalovich et al. (1988)") — where it serves as one of two independent references for
half-integer ``ℓ`` — we define it in a test module that both pages can use.
"""

using TestItems: @testmodule, @testitem  #hide
@testmodule Boyle2016 begin  #hide

using Quaternionic
#+

# Compute a single Wigner-D matrix element for half-integer or integer ``(ℓ, m', m)``,
# following Eq. (35) and Appendix A of [Boyle (2016)](@cite Boyle_2016).  `R` is a
# `Rotor`, and `ℓ`, `m′`, and `m` are the indices of the Wigner-D matrix element.  The
# indices must all be integers or all be `Rational` with denominators of 2.
function WignerDElement(R::Rotor{T}, ℓ::I, m′::I, m::I) where {T, I}
    ## If `I` is Rational, check that the denominators are 2
    if I <: Rational
        if (denominator(ℓ) != 2) || (denominator(m′) != 2) || (denominator(m) != 2)
            error("The indices ℓ, m′, and m must all be integers or all be half-integers")
        end
    end

    ## Convert to twice the input values for half-integer support
    L  = Int(2ℓ)
    M′ = Int(2m′)
    M  = Int(2m)

    if L > 16
        error(
            "The maximum supported ℓ for this function is 8; " *
            "larger numbers become numerically unstable.\n" *
            "Consider using the `WignerD` function instead."
        )
    end

    ## Simple helper for 0+0im
    zeroCT = zero(Complex{T})

    if abs(M′) > L || abs(M) > L
        return zeroCT
    end

    let π = T(π)
        ## Split input `R` into its two complex components and extract magnitude and phase
        Rₛ = Complex(R[1], R[4])
        Rₐ = Complex(R[3], R[2])
        rₛ = abs(Rₛ)
        rₐ = abs(Rₐ)
        ϕₛ = angle(Rₛ)
        ϕₐ = angle(Rₐ)

        ## Check simple limiting cases
        if rₐ ≤ 4eps(rₛ)
            if M′ != M
                return zeroCT
            else
                return cis(M * ϕₛ)
            end

        elseif rₛ ≤ 4eps(rₐ)
            if -M′ != M
                return zeroCT
            else
                return cis(M * ϕₐ) * (((L - M) % 4 == 0) ? 1 : -1)
            end

        elseif rₐ ≤ rₛ
            λ = -(rₐ/rₛ)^2
            ρₘᵢₙ = max(0, (M′ - M)÷2)
            κ = √T(
                    (factorial((L + M)÷2) * factorial((L - M)÷2))
                    / (factorial((L + M′)÷2) * factorial((L - M′)÷2))
                ) *
                binomial((L + M′)÷2, ρₘᵢₙ) * binomial((L - M′)÷2, (L - M)÷2 - ρₘᵢₙ)
            if (ρₘᵢₙ % 2) != 0
                κ = -κ
            end
            ρₘₐₓ = min((L + M′)÷2, (L - M)÷2)
            N₁ = L + M′ + 2
            N₂ = L - M + 2
            N₃  = M - M′

            total = one(T)
            for P in reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ)
                total *= λ * ((N₁ - P)*(N₂ - P)) / (P*(N₃ + P))
                total += one(T)
            end
            return κ *
                (rₛ ^ (L - (M - M′)÷2 - 2ρₘᵢₙ)) *
                (rₐ ^ ((M - M′)÷2 + 2ρₘᵢₙ)) *
                cis((M + M′)÷2 * ϕₛ + (M - M′)÷2 * ϕₐ) *
                total

        else # rₛ < rₐ
            λ = -(rₛ/rₐ)^2
            ρₘᵢₙ = max(0, -(M′ + M)÷2)
            κ = √T(
                    (factorial((L + M)÷2) * factorial((L - M)÷2))
                    / (factorial((L + M′)÷2) * factorial((L - M′)÷2))
                ) *
                binomial((L + M′)÷2, (L - M)÷2 - ρₘᵢₙ) * binomial((L - M′)÷2, ρₘᵢₙ)
            if (((L - M)÷2 - ρₘᵢₙ) % 2) != 0
                κ = -κ
            end
            ρₘₐₓ = min((L - M′)÷2, (L - M)÷2)
            N₁ = L - M′ + 2
            N₂ = L - M + 2
            N₃  = M + M′

            total = one(T)
            for P in reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ)
                total *= λ * ((N₁ - P)*(N₂ - P)) / (P*(N₃ + P))
                total += one(T)
            end
            return κ *
                (rₐ ^ (L - (M + M′)÷2 - 2ρₘᵢₙ)) *
                (rₛ ^ ((M + M′)÷2 + 2ρₘᵢₙ)) *
                cis((M - M′)÷2 * ϕₐ + (M + M′)÷2 * ϕₛ) *
                total

        end
    end
end
#+

# Equation (21), the spin-weighted spherical harmonics as functions of a rotor:
function ₛYₗₘ(s, ℓ, m, R::Rotor{T}) where {T}
    (-1)^s * √((2ℓ+1) / (4T(π))) * WignerDElement(R, ℓ, m, -s)
end
#+

end  #hide

md"""
## Tests

We can now test the functions against the equivalent functions from the `SphericalFunctions`
package.
"""

@testitem "Boyle 2016 conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities, Boyle2016] begin  #hide
using Quaternionic
#+

# We will need to test approximate floating-point equality, so we set absolute and relative
# tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 1000eps()
#+

# The algorithm is accurate up to
ℓₘₐₓ = 8
#+
# so we test up to that point, on a set of rotors that includes the identity, the basis
# rotations, rotors near the special cases of the algorithm, and random rotors:
Rs = [
    Rotor{Float64}(1);
    [Rotor{Float64}(𝐯) for 𝐯 ∈ (imx, imy, imz)];
    [exp(3eps() * 𝐯) for 𝐯 ∈ (imx, imy, imz)];
    [Rotor{Float64}(𝐮) * exp(3eps() * 𝐯) for 𝐮 ∈ (imx, imy, imz) for 𝐯 ∈ (imx, imy, imz)];
    randn(Rotor{Float64}, 20)
]
#+

# For integer indices, the 2016 matrices are the complex conjugates of ours:
for R ∈ Rs
    for (ℓ, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Boyle2016.WignerDElement(R, ℓ, m′, m) ≈ conj(ConventionsUtilities.D(ℓ, m′, m, R)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Because Eq. (21) uses ``𝔇_{m,-s}`` without a conjugate — where [our definition](@ref
# summary_swsh) uses ``\overline{𝔇_{m,-s}}`` — the spin-weighted spherical harmonics agree:
import Quaternionic: from_spherical_coordinates
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    R = from_spherical_coordinates(θ, ϕ)
    for (s, ℓ, m) ∈ sℓmrange(4, 2)
        @test Boyle2016.ₛYₗₘ(s, ℓ, m, R) ≈ ConventionsUtilities.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# For half-integer indices, the package itself does not (yet) provide a reference, so here
# we check two properties that any representation of ``\mathrm{Spin}(3)`` must satisfy —
# the representation property ``𝔇(𝐑_1 𝐑_2) = 𝔇(𝐑_1)\, 𝔇(𝐑_2)`` and the sign change
# ``𝔇(-𝐑) = -𝔇(𝐑)`` for half-integer ``ℓ`` — and defer the comparison against an
# independent closed form to the [Varshalovich page](@ref "Varshalovich et al. (1988)").
for R₁ ∈ Rs[1:8], R₂ ∈ Rs[end-4:end]
    for J ∈ (1//2, 3//2, 5//2, 7//2)
        𝔇₁ = [Boyle2016.WignerDElement(R₁, J, M′, M) for M′ ∈ -J:J, M ∈ -J:J]
        𝔇₂ = [Boyle2016.WignerDElement(R₂, J, M′, M) for M′ ∈ -J:J, M ∈ -J:J]
        𝔇₁₂ = [Boyle2016.WignerDElement(R₁ * R₂, J, M′, M) for M′ ∈ -J:J, M ∈ -J:J]
        @test 𝔇₁₂ ≈ 𝔇₁ * 𝔇₂ atol=ϵₐ rtol=ϵᵣ
        𝔇₋ = [Boyle2016.WignerDElement(-R₁, J, M′, M) for M′ ∈ -J:J, M ∈ -J:J]
        @test 𝔇₋ ≈ -𝔇₁ atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that the ``𝔇`` matrices of Boyle (2016) are the complex
# conjugates of those defined by the `SphericalFunctions` package, and that the
# spin-weighted spherical harmonics agree.

end  #hide
