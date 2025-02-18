@testmodule NINJA begin

    include("../utilities/naive_factorial.jl")
    import .NaiveFactorials: ❗

    function Wigner_d(ι::T, ℓ, m, s) where {T<:Real}
        # Eq. II.8 of v3 of Ajith et al. (2007) 'Data formats...'
        k_min = max(0, m - s)
        k_max = min(ℓ + m, ℓ - s)
        sum(
            (-1)^k
             * cos(ι / 2) ^ (2 * ℓ + m - s - 2 * k)
             * sin(ι / 2) ^ (2 * k + s - m)
             * T(
                √((ℓ + m)❗ * (ℓ - m)❗ * (ℓ + s)❗ * (ℓ - s)❗)
                / ((ℓ + m - k)❗ * (ℓ - s - k)❗ * (k)❗ * (k + s - m)❗)
             )
            for k in k_min:k_max
        )
    end

    raw"""

    Eq. II.7 of v3 of Ajith et al. (2007) 'Data formats...' says
    ```math
    {}_sY_{\ell,m} = (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} d^\ell_{m,-s}(\iota) e^{im\phi}
    ```

    Below Eq. (2.53) of [Torres del Castillo](@cite TorresDelCastillo_2003), we see
    ```math
    {}_sY_{j,m} = (-1)^m \sqrt{\frac{2j+1}{4\pi}} d^j_{-m,s}(\iota) e^{im\phi}
    ```
    We can use identities to modify the latter as follows:
    ```math
    \begin{aligned}
    {}_sY_{j,m} &= (-1)^m \sqrt{\frac{2j+1}{4\pi}} d^j_{-m,s}(\iota) e^{im\phi} \\
               &= (-1)^m \sqrt{\frac{2j+1}{4\pi}} d^j_{-s,m}(\iota) e^{im\phi} \\
               &= (-1)^{j-s+m} \sqrt{\frac{2j+1}{4\pi}} d^j_{-s,-m}(\pi-\iota) e^{im\phi} \\
               &= (-1)^{j-s+m} \sqrt{\frac{2j+1}{4\pi}} d^j_{m,s}(\pi-\iota) e^{im\phi} \\
               &= (-1)^{2j-s+2m} \sqrt{\frac{2j+1}{4\pi}} d^j_{m,-s}(\pi-(\pi-\iota)) e^{im\phi} \\
               &= (-1)^{s} \sqrt{\frac{2j+1}{4\pi}} d^j_{m,-s}(\iota) e^{im\phi} \\
    \end{aligned}
    ```
    The last line assumes that `j`, `m`, and `s` are integers.  But in that case, the NINJA
    expression agrees with the Torres del Castillo expression.

    """

    function sYlm(s, ell, m, ι::T, ϕ::T) where {T<:Real}
        # Eq. II.7 of v3 of Ajith et al. (2007) 'Data formats...'
        # Note the weird definition w.r.t. `-s`
        if abs(s) > ell || abs(m) > ell
            return zero(complex(T))
        end
        (
            (-1)^s
            * √((2ell + 1) / (4T(π)))
            * Wigner_d(ι, ell, m, -s)
            * cis(m * ϕ)
        )
    end

    function sYlm(s, ℓ, m, ιϕ)
        sYlm(s, ℓ, m, ιϕ[1], ιϕ[2])
    end

    # Eqs. (II.9) through (II.13) of https://arxiv.org/abs/0709.0093v3 [Ajith_2007](@cite)
    m2Y22(ι::T, ϕ::T) where {T<:Real} = √(5 / (64T(π))) * (1 + cos(ι))^2 * cis(2ϕ)
    m2Y21(ι::T, ϕ::T) where {T<:Real} = √(5 / (16T(π))) * sin(ι) * (1 + cos(ι)) * cis(ϕ)
    m2Y20(ι::T, ϕ::T) where {T<:Real} = √(15 / (32T(π))) * sin(ι)^2
    m2Y2m1(ι::T, ϕ::T) where {T<:Real} = √(5 / (16T(π))) * sin(ι) * (1 - cos(ι)) * cis(-1ϕ)
    m2Y2m2(ι::T, ϕ::T) where {T<:Real} = √(5 / (64T(π))) * (1 - cos(ι))^2 * cis(-2ϕ)

    m_m2Y2m = [
        (2, m2Y22),
        (1, m2Y21),
        (0, m2Y20),
        (-1, m2Y2m1),
        (-2, m2Y2m2)
    ]

end  # module NINJA
