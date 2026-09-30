@testmodule ExplicitWignerMatrices begin

    include("naive_factorial.jl")
    import .NaiveFactorials: ❗

    function d_explicit(n, m′, m, expiβ::Complex{T}) where T
        if abs(m′) < abs(m)
            return (-1)^(m-m′) * d_explicit(n, m, m′, expiβ)
        end
        if m′ < 0
            return (-1)^(m-m′) * d_explicit(n, -m′, -m, expiβ)
        end
        cosβ = expiβ.re
        sinβ = expiβ.im
        if (n,m′,m) == (0,0,0)
            T(1)
        elseif (n,m′,m) == (1,0,0)
            cosβ
        elseif (n,m′,m) == (1,1,-1)
            (1-cosβ) / 2
        elseif (n,m′,m) == (1,1,0)
            -sinβ / √T(2)
        elseif (n,m′,m) == (1,1,1)
            (1+cosβ) / 2
        elseif (n,m′,m) == (2,0,0)
            (3cosβ^2-1) / 2
        elseif (n,m′,m) == (2,1,-1)
            (1+cosβ-2cosβ^2) / 2
        elseif (n,m′,m) == (2,1,0)
            -√(T(3)/8) * 2 * sinβ * cosβ
        elseif (n,m′,m) == (2,1,1)
            (-1+cosβ+2cosβ^2) / 2
        elseif (n,m′,m) == (2,2,-2)
            (1-cosβ)^2 / 4
        elseif (n,m′,m) == (2,2,-1)
            -sinβ * (1-cosβ) / 2
        elseif (n,m′,m) == (2,2,0)
            √(T(3)/8) * sinβ^2
        elseif (n,m′,m) == (2,2,1)
            -sinβ * (1+cosβ) / 2
        elseif (n,m′,m) == (2,2,2)
            (1+cosβ)^2/4
        else
            T(NaN)
        end
    end

    function d_formula(n, m′, m, expiβ::Complex{T}) where T
        # https://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_.28small.29_d-matrix
        cosβ = expiβ.re
        sin½β = √((1-cosβ)/2)
        cos½β = √((1+cosβ)/2)
        √T((n + m′)❗ * (n - m′)❗ * (n + m)❗ * (n - m)❗) * sum(
            (-1)^(m′ - m + s)
            * cos½β ^ (2n + m - m′ - 2s)
            * sin½β ^ (m′ - m + 2s)
            / T((n + m - s)❗ * (s)❗ * (m′ - m + s)❗ * (n - m′ - s)❗)
            for s in max(0, m - m′):min(n + m, n - m′)
        )
    end

    function D_formula(n, m′, m, expiα::Complex{T}, expiβ::Complex{T}, expiγ::Complex{T}) where T
        # https://en.wikipedia.org/wiki/Wigner_D-matrix#Definition_of_the_Wigner_D-matrix
        return expiα^(-m′) * d_formula(n, m′, m, expiβ) * expiγ^(-m)
    end

    # 𝔇^ℓ_{m′,m}(R) as a polynomial in the rotor's normalized Cayley–Klein parameters σ =
    # (W + iZ)/‖R‖ and ρ = (Y - iX)/‖R‖ and their conjugates,
    #
    #     𝔇^ℓ_{m′,m} = Σₛ (-1)^{k+s} Cₛ σ̄^{ℓ+m-s} σ^{ℓ-m′-s} ρ̄^{k+s} ρ^s,        k = m′ - m,
    #
    # which is `d_formula` with σ = cos(β/2) e^{i(α+γ)/2} and ρ = sin(β/2) e^{i(α-γ)/2}: the
    # phases e^{-im′α} and e^{-imγ} of the convention are exactly what the powers of those
    # two phases give.  It shares no code with the package, which uses it only at the poles,
    # and in a truncated form.  The coefficients are exact integers until the final square
    # root, which is taken in `BigFloat`, and the powers are formed by repeated squaring
    # rather than by `^`, which for `Complex{<:ForwardDiff.Dual}` gives wrong second
    # derivatives.  `R` may have components of any real type, including dual numbers, and
    # need not be normalized; the indices may be integers or `Rational`s with denominator 2.
    #
    # The sum cancels badly away from the poles, by as much as 2^ℓ, so the sum of the
    # magnitudes of its terms is returned too, as the scale of its rounding error; that is
    # `nothing` for a type other than an `AbstractFloat`.
    function D_polynomial(ℓ, m′, m, R)
        W, X, Y, Z = R[1], R[2], R[3], R[4]
        nrm = sqrt(W^2 + X^2 + Y^2 + Z^2)
        σ = Complex(W, Z) / nrm
        ρ = Complex(Y, -X) / nrm
        T = typeof(real(σ))
        pw(z, n) = Base.power_by_squaring(z, n)
        A, B, C, E = Int(ℓ + m), Int(ℓ - m), Int(ℓ + m′), Int(ℓ - m′)
        k = Int(m′ - m)
        value = zero(σ)
        scale = zero(T)
        for s ∈ max(0, -k):min(A, E)
            C² = binomial(big(A), s) * binomial(big(E), s) * binomial(big(C), k + s) * binomial(big(B), k + s)
            Cₛ = convert(T, sqrt(big(C²)))
            value += (-1)^(k + s) * Cₛ * (pw(conj(σ), A - s) * pw(σ, E - s)) * (pw(conj(ρ), k + s) * pw(ρ, s))
            if T <: AbstractFloat
                scale += Cₛ * abs(σ)^(A + E - 2s) * abs(ρ)^(k + 2s)
            end
        end
        value, (T <: AbstractFloat ? scale : nothing)
    end

    # ₛYₗₘ(R) = (-1)^s √((2ℓ+1)/4π) conj(𝔇^ℓ_{m,-s}(R)), with (-1)^s = i^{2s}, from the
    # polynomial above; the second value is again the scale of the rounding error.
    function sYlm_polynomial(ℓ, m, s, R)
        𝔇, scale = D_polynomial(ℓ, m, -s, R)
        T = typeof(real(𝔇))
        i²ˢ = (1, im, -1, -im)[mod(Int(2s), 4) + 1]
        prefactor = sqrt((2ℓ + 1) / (4 * T(π)))
        prefactor * i²ˢ * conj(𝔇), (scale === nothing ? nothing : prefactor * scale)
    end

end  # module ExplicitWignerMatrices
