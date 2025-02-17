@testmodule ExplicitWignerMatrices begin

    include("../utilities/naive_factorial.jl")
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

    function D_explicit(n, m′, m, expiα::Complex{T}, expiβ::Complex{T}, expiγ::Complex{T}) where T
        return expiα^(-m′) * d_explicit(n, m′, m, expiβ) * expiγ^(-m)
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

end  # module ExplicitWignerMatrices
