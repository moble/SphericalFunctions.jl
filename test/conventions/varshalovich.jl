raw"""
Formulas and conventions from [Varshalovich's "Quantum Theory of Angular Momentum"](@cite
Varshalovich_1988).

Note that Varshalovich labels his indices with `M` and `M′`, respectively, but if we just
plug in `m′` and `m` (note the order), we get the expected result — his formulas are the
same as this package's, except with a conjugate.

Varshalovich defines his Euler angles (scheme B, page 22) in the same way we do, except that
he specifies that this describes the rotation *of the coordinate system*.

Sec. 4.8.2 (page 92) relates the integer-index elements to the following half-integer-index
elements.  Specifically, Eqs. (14) and (15) derive the relationships from the Clebsch-Gordan
coefficients.  That is, the product of two Wigner matrices can be given as a sum over a
Wigner matrices times a pair of Clebsch-Gordan coefficients.  If one of the matrices has
spin 1/2, this gives us a series of relationships between the integer-index elements and the
half-integer-index elements, which can be combined to give the desired relationship.  Then,
given knowledge of the 1/2-spin representation (which is essentially the standard
$\mathrm{SU}(2)$ representation), we can then get any half-integer spin result from the
preceding whole-integer spin results.

Specifically, we have (from Table 4.3, page 119):
```julia
D(1//2,  1//2,  1//2, α, β, γ) =  exp(-𝒾*α/2) * cos(β/2) * exp(-𝒾*γ/2)
D(1//2,  1//2, -1//2, α, β, γ) = -exp(-𝒾*α/2) * sin(β/2) * exp( 𝒾*γ/2)
D(1//2, -1//2,  1//2, α, β, γ) =  exp( 𝒾*α/2) * sin(β/2) * exp(-𝒾*γ/2)
D(1//2, -1//2, -1//2, α, β, γ) =  exp( 𝒾*α/2) * cos(β/2) * exp( 𝒾*γ/2)
```
"""
@testmodule Varshalovich begin

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


@doc raw"""
    D(J, M, M′, α, β, γ)

Eq. 4.3(1) of [Varshalovich](@cite Varshalovich_1988), implementing
```math
    D^{J}_{M,M'}(\alpha, \beta, \gamma).
```

See also [`d`](@ref) for Varshalovich's version the Wigner d-function.
"""
function D(J, M, M′, α, β, γ)
    exp(-𝒾*M*α) * d(J, M, M′, β) * exp(-𝒾*M′*γ)
end


@doc raw"""
    d(J, M, M′, β)

Eqs. 4.3.1(2) of [Varshalovich](@cite Varshalovich_1988), implementing
```math
    d^{J}_{M,M'}(\beta).
```

See also [`D`](@ref) for Varshalovich's version the Wigner D-function.
"""
function d(J::I, M::I, M′::I, β::T) where {I, T}
    if J < 0
        throw(DomainError("J=$J must be non-negative"))
    end
    if abs(M) > J || abs(M′) > J
        if I <: Rational && abs(M) ≤ J+2 && abs(M′) ≤ J+2
            return zero(β)  # Simplify half-integer formulas by accepting this
        end
        #throw(DomainError("abs(M=$M) and abs(M=$M′) must be ≤ J=$J"))
    end
    if J ≥ 8
        throw(DomainError("J=$J≥8 will lead to overflow errors"))
    end

    # The summation index `k` ranges over all values for which the factorials are
    # non-negative.
    kₘᵢₙ = max(0, -(M+M′))
    kₘₐₓ = min(J-M, J-M′)

    # Note that Varshalovich's actual formula is reproduced here, even though it leads to
    # overflow errors for `J ≥ 8`, which could be eliminated by other means.
    return (-1)^(J-M′) * √T((J+M)❗ * (J-M)❗ * (J+M′)❗ * (J-M′)❗) *
    sum(
        k -> (
            (-1)^(k) * cos(β/2)^(M+M′+2k) * sin(β/2)^(2J-M-M′-2k)
            / T((k)❗ * (J-M-k)❗ * (J-M′-k)❗ * (M+M′+k)❗)
        ),
        kₘᵢₙ:kₘₐₓ,
        init=zero(T)
    )
end

end  # @testmodule Varshalovich


@testitem "Varshalovich conventions" setup=[Utilities, Varshalovich] begin
    using Random
    using Quaternionic: from_spherical_coordinates

    Random.seed!(1234)
    const 𝒾 = im
    const T = Float64
    const ℓₘₐₓ = 7
    ϵₐ = 8eps(T)
    ϵᵣ = 20eps(T)

    # Tests for 𝒟(j, m′, m, α, β, γ)
    let ϵₐ=√ϵᵣ, ϵᵣ=√ϵᵣ, 𝒟=Varshalovich.D
        n = 4
        for α ∈ αrange(T, n)
            for β ∈ βrange(T, n)
                if abs(sin(β)) ≤ eps(T)
                    continue
                end

                for γ ∈ γrange(T, n)
                    D = D_matrices(α, β, γ, ℓₘₐₓ)
                    i = 1
                    for j in 0:ℓₘₐₓ
                        for m′ in -j:j
                            for m in -j:j
                                @test 𝒟(j, m′, m, α, β, γ) ≈ conj(D[i]) atol=ϵₐ rtol=ϵᵣ
                                i += 1
                            end
                        end
                    end

                    # Test half-integer formula
                    for j in 1//2:ℓₘₐₓ
                        for m′ in -j:j
                            for m in -j:j
                                D1 = 𝒟(j, m, m′, α, β, γ)
                                D2 = if m′ ≠ j  # use Eq. 4.8.2(14)
                                    (
                                        √((j-m)/(j-m′)) * cos(β/2) * exp(𝒾*(α+γ)/2) *
                                        𝒟(j-1//2, m+1//2, m′+1//2, α, β, γ)
                                        -
                                        √((j+m)/(j-m′)) * sin(β/2) * exp(-𝒾*(α-γ)/2) *
                                        𝒟(j-1//2, m-1//2, m′+1//2, α, β, γ)
                                    )
                                else  # use Eq. 4.8.2(15)
                                    (
                                        √((j-m)/(j+m′)) * sin(β/2) * exp(𝒾*(α-γ)/2) *
                                        𝒟(j-1//2, m+1//2, m′-1//2, α, β, γ)
                                        +
                                        √((j+m)/(j+m′)) * cos(β/2) * exp(-𝒾*(α+γ)/2) *
                                        𝒟(j-1//2, m-1//2, m′-1//2, α, β, γ)
                                    )
                                end
                                @test D1 ≈ D2 atol=ϵₐ rtol=ϵᵣ
                            end
                        end
                    end
                end
            end
        end
    end

end
