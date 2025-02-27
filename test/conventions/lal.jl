@testmodule LAL begin
    # The code in this section is translated from C code in LALSuite:
    #
    #   https://git.ligo.org/lscsoft/lalsuite/-/blob/6e653c91b6e8a6728c4475729c4f967c9e09f020/lal/lib/utilities/SphericalHarmonics.c
    #
    # That code is licensed under GPLv2.  See the link for details.

    function XLALSpinWeightedSphericalHarmonic(
        theta::Float64,  # polar angle (rad)
        phi::Float64,    # azimuthal angle (rad)
        s::Int,   # spin weight
        l::Int,   # mode number l
        m::Int    # mode number m
    )
        # Sanity checks
        if l < abs(s)
            error("Invalid mode s=$s, l=$l, m=$m - require |s| <= l")
        end
        if l < abs(m)
            error("Invalid mode s=$s, l=$l, m=$m - require |m| <= l")
        end
        if s != -2
            error("Unsupported mode s=$s (only s=-2 implemented)")
        end
        if l < 2 || l > 8
            error("Unsupported mode l=$l (only l in [2,8] implemented)")
        end
        # Compute real factor
        fac = if l == 2
            if m == -2
                sqrt(5.0 / (64.0 * π)) * (1.0 - cos(theta)) * (1.0 - cos(theta))
            elseif m == -1
                sqrt(5.0 / (16.0 * π)) * sin(theta) * (1.0 - cos(theta))
            elseif m == 0
                sqrt(15.0 / (32.0 * π)) * sin(theta) * sin(theta)
            elseif m == 1
                sqrt(5.0 / (16.0 * π)) * sin(theta) * (1.0 + cos(theta))
            elseif m == 2
                sqrt(5.0 / (64.0 * π)) * (1.0 + cos(theta)) * (1.0 + cos(theta))
            end
        elseif l == 3
            if m == -3
                sqrt(21.0 / (2.0 * π)) * cos(theta/2.0) * sin(theta/2.0)^5
            elseif m == -2
                sqrt(7.0 / (4.0 * π)) * (2.0 + 3.0 * cos(theta)) * sin(theta/2.0)^4
            elseif m == -1
                sqrt(35.0 / (2.0 * π)) * (sin(theta) + 4.0 * sin(2.0 * theta) - 3.0 * sin(3.0 * theta)) / 32.0
            elseif m == 0
                sqrt(105.0 / (2.0 * π)) * cos(theta) * sin(theta)^2 / 4.0
            elseif m == 1
                -sqrt(35.0 / (2.0 * π)) * (sin(theta) - 4.0 * sin(2.0 * theta) - 3.0 * sin(3.0 * theta)) / 32.0
            elseif m == 2
                sqrt(7.0 / π) * cos(theta/2.0)^4 * (-2.0 + 3.0 * cos(theta)) / 2.0
            elseif m == 3
                -sqrt(21.0 / (2.0 * π)) * cos(theta/2.0)^5 * sin(theta/2.0)
            end
        elseif l == 4
            if m == -4
                3.0 * sqrt(7.0 / π) * cos(theta/2.0)^2 * sin(theta/2.0)^6
            elseif m == -3
                3.0 * sqrt(7.0 / (2.0 * π)) * cos(theta/2.0) * (1.0 + 2.0 * cos(theta)) * sin(theta/2.0)^5
            elseif m == -2
                3.0 * (9.0 + 14.0 * cos(theta) + 7.0 * cos(2.0 * theta)) * sin(theta/2.0)^4 / (4.0 * sqrt(π))
            elseif m == -1
                3.0 * (3.0 * sin(theta) + 2.0 * sin(2.0 * theta) + 7.0 * sin(3.0 * theta) - 7.0 * sin(4.0 * theta)) / (32.0 * sqrt(2.0 * π))
            elseif m == 0
                3.0 * sqrt(5.0 / (2.0 * π)) * (5.0 + 7.0 * cos(2.0 * theta)) * sin(theta)^2 / 16.0
            elseif m == 1
                3.0 * (3.0 * sin(theta) - 2.0 * sin(2.0 * theta) + 7.0 * sin(3.0 * theta) + 7.0 * sin(4.0 * theta)) / (32.0 * sqrt(2.0 * π))
            elseif m == 2
                3.0 * cos(theta/2.0)^4 * (9.0 - 14.0 * cos(theta) + 7.0 * cos(2.0 * theta)) / (4.0 * sqrt(π))
            elseif m == 3
                -3.0 * sqrt(7.0 / (2.0 * π)) * cos(theta/2.0)^5 * (-1.0 + 2.0 * cos(theta)) * sin(theta/2.0)
            elseif m == 4
                3.0 * sqrt(7.0 / π) * cos(theta/2.0)^6 * sin(theta/2.0)^2
            end
        elseif l == 5
            if m == -5
                sqrt(330.0 / π) * cos(theta/2.0)^3 * sin(theta/2.0)^7
            elseif m == -4
                sqrt(33.0 / π) * cos(theta/2.0)^2 * (2.0 + 5.0 * cos(theta)) * sin(theta/2.0)^6
            elseif m == -3
                sqrt(33.0 / (2.0 * π)) * cos(theta/2.0) * (17.0 + 24.0 * cos(theta) + 15.0 * cos(2.0 * theta)) * sin(theta/2.0)^5 / 4.0
            elseif m == -2
                sqrt(11.0 / π) * (32.0 + 57.0 * cos(theta) + 36.0 * cos(2.0 * theta) + 15.0 * cos(3.0 * theta)) * sin(theta/2.0)^4 / 8.0
            elseif m == -1
                sqrt(77.0 / π) * (2.0 * sin(theta) + 8.0 * sin(2.0 * theta) + 3.0 * sin(3.0 * theta) + 12.0 * sin(4.0 * theta) - 15.0 * sin(5.0 * theta)) / 256.0
            elseif m == 0
                sqrt(1155.0 / (2.0 * π)) * (5.0 * cos(theta) + 3.0 * cos(3.0 * theta)) * sin(theta)^2 / 32.0
            elseif m == 1
                sqrt(77.0 / π) * (-2.0 * sin(theta) + 8.0 * sin(2.0 * theta) - 3.0 * sin(3.0 * theta) + 12.0 * sin(4.0 * theta) + 15.0 * sin(5.0 * theta)) / 256.0
            elseif m == 2
                sqrt(11.0 / π) * cos(theta/2.0)^4 * (-32.0 + 57.0 * cos(theta) - 36.0 * cos(2.0 * theta) + 15.0 * cos(3.0 * theta)) / 8.0
            elseif m == 3
                -sqrt(33.0 / (2.0 * π)) * cos(theta/2.0)^5 * (17.0 - 24.0 * cos(theta) + 15.0 * cos(2.0 * theta)) * sin(theta/2.0) / 4.0
            elseif m == 4
                sqrt(33.0 / π) * cos(theta/2.0)^6 * (-2.0 + 5.0 * cos(theta)) * sin(theta/2.0)^2
            elseif m == 5
                -sqrt(330.0 / π) * cos(theta/2.0)^7 * sin(theta/2.0)^3
            end
        elseif l == 6
            if m == -6
                (3.0 * sqrt(715.0 / π) * cos(theta/2.0)^4 * sin(theta/2.0)^8) / 2.0
            elseif m == -5
                (sqrt(2145.0 / π) * cos(theta/2.0)^3 * (1.0 + 3.0 * cos(theta)) * sin(theta/2.0)^7) / 2.0
            elseif m == -4
                (sqrt(195.0 / (2.0 * π)) * cos(theta/2.0)^2 * (35.0 + 44.0 * cos(theta) + 33.0 * cos(2.0 * theta)) * sin(theta/2.0)^6) / 8.0
            elseif m == -3
                (3.0 * sqrt(13.0 / π) * cos(theta/2.0) * (98.0 + 185.0 * cos(theta) + 110.0 * cos(2.0 * theta) + 55.0 * cos(3.0 * theta)) * sin(theta/2.0)^5) / 32.0
            elseif m == -2
                (sqrt(13.0 / π) * (1709.0 + 3096.0 * cos(theta) + 2340.0 * cos(2.0 * theta) + 1320.0 * cos(3.0 * theta) + 495.0 * cos(4.0 * theta)) * sin(theta/2.0)^4) / 256.0
            elseif m == -1
                (sqrt(65.0 / (2.0 * π)) * cos(theta/2.0) * (161.0 + 252.0 * cos(theta) + 252.0 * cos(2.0 * theta) + 132.0 * cos(3.0 * theta) + 99.0 * cos(4.0 * theta)) * sin(theta/2.0)^3) / 64.0
            elseif m == 0
                (sqrt(1365.0 / π) * (35.0 + 60.0 * cos(2.0 * theta) + 33.0 * cos(4.0 * theta)) * sin(theta)^2) / 512.0
            elseif m == 1
                (sqrt(65.0 / (2.0 * π)) * cos(theta/2.0)^3 * (161.0 - 252.0 * cos(theta) + 252.0 * cos(2.0 * theta) - 132.0 * cos(3.0 * theta) + 99.0 * cos(4.0 * theta)) * sin(theta/2.0)) / 64.0
            elseif m == 2
                (sqrt(13.0 / π) * cos(theta/2.0)^4 * (1709.0 - 3096.0 * cos(theta) + 2340.0 * cos(2.0 * theta) - 1320.0 * cos(3.0 * theta) + 495.0 * cos(4.0 * theta))) / 256.0
            elseif m == 3
                (-3.0 * sqrt(13.0 / π) * cos(theta/2.0)^5 * (-98.0 + 185.0 * cos(theta) - 110.0 * cos(2.0 * theta) + 55.0 * cos(3.0 * theta)) * sin(theta/2.0)) / 32.0
            elseif m == 4
                (sqrt(195.0 / (2.0 * π)) * cos(theta/2.0)^6 * (35.0 - 44.0 * cos(theta) + 33.0 * cos(2.0 * theta)) * sin(theta/2.0)^2) / 8.0
            elseif m == 5
                (-sqrt(2145.0 / π) * cos(theta/2.0)^7 * (-1.0 + 3.0 * cos(theta)) * sin(theta/2.0)^3) / 2.0
            elseif m == 6
                (3.0 * sqrt(715.0 / π) * cos(theta/2.0)^8 * sin(theta/2.0)^4) / 2.0
            end
        elseif l == 7
            if m == -7
                sqrt(15015.0 / (2.0 * π)) * cos(theta/2.0)^5 * sin(theta/2.0)^9
            elseif m == -6
                (sqrt(2145.0 / π) * cos(theta/2.0)^4 * (2.0 + 7.0 * cos(theta)) * sin(theta/2.0)^8) / 2.0
            elseif m == -5
                (sqrt(165.0 / (2.0 * π)) * cos(theta/2.0)^3 * (93.0 + 104.0 * cos(theta) + 91.0 * cos(2.0 * theta)) * sin(theta/2.0)^7) / 8.0
            elseif m == -4
                (sqrt(165.0 / (2.0 * π)) * cos(theta/2.0)^2 * (140.0 + 285.0 * cos(theta) + 156.0 * cos(2.0 * theta) + 91.0 * cos(3.0 * theta)) * sin(theta/2.0)^6) / 16.0
            elseif m == -3
                (sqrt(15.0 / (2.0 * π)) * cos(theta/2.0) * (3115.0 + 5456.0 * cos(theta) + 4268.0 * cos(2.0 * theta) + 2288.0 * cos(3.0 * theta) + 1001.0 * cos(4.0 * theta)) * sin(theta/2.0)^5) / 128.0
            elseif m == -2
                (sqrt(15.0 / π) * (5220.0 + 9810.0 * cos(theta) + 7920.0 * cos(2.0 * theta) + 5445.0 * cos(3.0 * theta) + 2860.0 * cos(4.0 * theta) + 1001.0 * cos(5.0 * theta)) * sin(theta/2.0)^4) / 512.0
            elseif m == -1
                (3.0 * sqrt(5.0 / (2.0 * π)) * cos(theta/2.0) * (1890.0 + 4130.0 * cos(theta) + 3080.0 * cos(2.0 * theta) + 2805.0 * cos(3.0 * theta) + 1430.0 * cos(4.0 * theta) + 1001.0 * cos(5.0 * theta)) * sin(theta/2.0)^3) / 512.0
            elseif m == 0
                (3.0 * sqrt(35.0 / π) * cos(theta) * (109.0 + 132.0 * cos(2.0 * theta) + 143.0 * cos(4.0 * theta)) * sin(theta)^2) / 512.0
            elseif m == 1
                (3.0 * sqrt(5.0 / (2.0 * π)) * cos(theta/2.0)^3 * (-1890.0 + 4130.0 * cos(theta) - 3080.0 * cos(2.0 * theta) + 2805.0 * cos(3.0 * theta) - 1430.0 * cos(4.0 * theta) + 1001.0 * cos(5.0 * theta)) * sin(theta/2.0)) / 512.0
            elseif m == 2
                (sqrt(15.0 / π) * cos(theta/2.0)^4 * (-5220.0 + 9810.0 * cos(theta) - 7920.0 * cos(2.0 * theta) + 5445.0 * cos(3.0 * theta) - 2860.0 * cos(4.0 * theta) + 1001.0 * cos(5.0 * theta))) / 512.0
            elseif m == 3
                -(sqrt(15.0 / (2.0 * π)) * cos(theta/2.0)^5 * (3115.0 - 5456.0 * cos(theta) + 4268.0 * cos(2.0 * theta) - 2288.0 * cos(3.0 * theta) + 1001.0 * cos(4.0 * theta)) * sin(theta/2.0)) / 128.0
            elseif m == 4
                (sqrt(165.0 / (2.0 * π)) * cos(theta/2.0)^6 * (-140.0 + 285.0 * cos(theta) - 156.0 * cos(2.0 * theta) + 91.0 * cos(3.0 * theta)) * sin(theta/2.0)^2) / 16.0
            elseif m == 5
                -(sqrt(165.0 / (2.0 * π)) * cos(theta/2.0)^7 * (93.0 - 104.0 * cos(theta) + 91.0 * cos(2.0 * theta)) * sin(theta/2.0)^3) / 8.0
            elseif m == 6
                (sqrt(2145.0 / π) * cos(theta/2.0)^8 * (-2.0 + 7.0 * cos(theta)) * sin(theta/2.0)^4) / 2.0
            elseif m == 7
                -(sqrt(15015.0 / (2.0 * π)) * cos(theta/2.0)^9 * sin(theta/2.0)^5)
            end
        elseif l == 8
            if m == -8
                sqrt(34034.0 / π) * cos(theta/2.0)^6 * sin(theta/2.0)^10
            elseif m == -7
                sqrt(17017.0 / (2.0 * π)) * cos(theta/2.0)^5 * (1.0 + 4.0 * cos(theta)) * sin(theta/2.0)^9
            elseif m == -6
                sqrt(255255.0 / π) * cos(theta/2.0)^4 * (1.0 + 2.0 * cos(theta)) * sin(π/4.0 - theta/2.0) * sin(π/4.0 + theta/2.0) * sin(theta/2.0)^8
            elseif m == -5
                (sqrt(12155.0 / (2.0 * π)) * cos(theta/2.0)^3 * (19.0 + 42.0 * cos(theta) + 21.0 * cos(2.0 * theta) + 14.0 * cos(3.0 * theta)) * sin(theta/2.0)^7) / 8.0
            elseif m == -4
                (sqrt(935.0 / (2.0 * π)) * cos(theta/2.0)^2 * (265.0 + 442.0 * cos(theta) + 364.0 * cos(2.0 * theta) + 182.0 * cos(3.0 * theta) + 91.0 * cos(4.0 * theta)) * sin(theta/2.0)^6) / 32.0
            elseif m == -3
                (sqrt(561.0 / (2.0 * π)) * cos(theta/2.0) * (869.0 + 1660.0 * cos(theta) + 1300.0 * cos(2.0 * theta) + 910.0 * cos(3.0 * theta) + 455.0 * cos(4.0 * theta) + 182.0 * cos(5.0 * theta)) * sin(theta/2.0)^5) / 128.0
            elseif m == -2
                (sqrt(17.0 / π) * (7626.0 + 14454.0 * cos(theta) + 12375.0 * cos(2.0 * theta) + 9295.0 * cos(3.0 * theta) + 6006.0 * cos(4.0 * theta) + 3003.0 * cos(5.0 * theta) + 1001.0 * cos(6.0 * theta)) * sin(theta/2.0)^4) / 512.0
            elseif m == -1
                (sqrt(595.0 / (2.0 * π)) * cos(theta/2.0) * (798.0 + 1386.0 * cos(theta) + 1386.0 * cos(2.0 * theta) + 1001.0 * cos(3.0 * theta) + 858.0 * cos(4.0 * theta) + 429.0 * cos(5.0 * theta) + 286.0 * cos(6.0 * theta)) * sin(theta/2.0)^3) / 512.0
            elseif m == 0
                (3.0 * sqrt(595.0 / π) * (210.0 + 385.0 * cos(2.0 * theta) + 286.0 * cos(4.0 * theta) + 143.0 * cos(6.0 * theta)) * sin(theta)^2) / 4096.0
            elseif m == 1
                (sqrt(595.0 / (2.0 * π)) * cos(theta/2.0)^3 * (798.0 - 1386.0 * cos(theta) + 1386.0 * cos(2.0 * theta) - 1001.0 * cos(3.0 * theta) + 858.0 * cos(4.0 * theta) - 429.0 * cos(5.0 * theta) + 286.0 * cos(6.0 * theta)) * sin(theta/2.0)) / 512.0
            elseif m == 2
                (sqrt(17.0 / π) * cos(theta/2.0)^4 * (7626.0 - 14454.0 * cos(theta) + 12375.0 * cos(2.0 * theta) - 9295.0 * cos(3.0 * theta) + 6006.0 * cos(4.0 * theta) - 3003.0 * cos(5.0 * theta) + 1001.0 * cos(6.0 * theta))) / 512.0
            elseif m == 3
                -(sqrt(561.0 / (2.0 * π)) * cos(theta/2.0)^5 * (-869.0 + 1660.0 * cos(theta) - 1300.0 * cos(2.0 * theta) + 910.0 * cos(3.0 * theta) - 455.0 * cos(4.0 * theta) + 182.0 * cos(5.0 * theta)) * sin(theta/2.0)) / 128.0
            elseif m == 4
                (sqrt(935.0 / (2.0 * π)) * cos(theta/2.0)^6 * (265.0 - 442.0 * cos(theta) + 364.0 * cos(2.0 * theta) - 182.0 * cos(3.0 * theta) + 91.0 * cos(4.0 * theta)) * sin(theta/2.0)^2) / 32.0
            elseif m == 5
                -(sqrt(12155.0 / (2.0 * π)) * cos(theta/2.0)^7 * (-19.0 + 42.0 * cos(theta) - 21.0 * cos(2.0 * theta) + 14.0 * cos(3.0 * theta)) * sin(theta/2.0)^3) / 8.0
            elseif m == 6
                (sqrt(255255.0 / π) * cos(theta/2.0)^8 * (-1.0 + 2.0 * cos(theta)) * sin(theta/2.0)^4) * sin(π/4.0 - theta/2.0) * sin(π/4.0 + theta/2.0);
            elseif m == 7
                -(sqrt(17017.0 / (2.0 * π)) * cos(theta/2.0)^9 * (-1.0 + 4.0 * cos(theta)) * sin(theta/2.0)^5)
            elseif m == 8
                sqrt(34034.0 / π) * cos(theta/2.0)^10 * sin(theta/2.0)^6
            end
        end
        # Include complex phase factor
        if m ≠ 0
            ans = cis(m*phi) * fac
        else
            ans = fac
        end
    end

    function XLALJacobiPolynomial(
        n::Int,    # degree
        alpha::Int,  # alpha parameter
        beta::Int,   # beta parameter
        x::Float64  # argument
    )
        f1 = (x-1.0)/2.0
        f2 = (x+1.0)/2.0
        sum = 0.0
        if n == 0
            return 1.0
        end
        for s = 0:n
            val = 1.0
            val *= binomial(n+alpha, s)
            val *= binomial(n+beta, n-s)
            if n-s != 0
                val *= f1^(n-s)
            end
            if s != 0
                val *= f2^s
            end
            sum += val
        end
        return sum
    end

    function XLALWignerdMatrix(
        l::Int,        # mode number l
        mp::Int,       # mode number m'
        m::Int,        # mode number m
        beta::Float64  # euler angle (rad)
    )
        k = min(l+m, min(l-m, min(l+mp, l-mp)))
        a = 0
        lam = 0
        if k == l+m
            a = mp-m
            lam = mp-m
        elseif k == l-m
            a = m-mp
            lam = 0
        elseif k == l+mp
            a = m-mp
            lam = 0
        elseif k == l-mp
            a = mp-m
            lam = mp-m
        end
        b = 2*l-2*k-a
        pref = (-1)^lam * sqrt(binomial(2*l-k, k+a)) / sqrt(binomial(k+b, b))
        return pref * sin(beta/2.0)^a * cos(beta/2.0)^b * XLALJacobiPolynomial(k, a, b, cos(beta))
    end

    function XLALWignerDMatrix(
        l::Int,          # mode number l
        mp::Int,         # mode number m'
        m::Int,          # mode number m
        alpha::Float64,  # euler angle (rad)
        beta::Float64,   # euler angle (rad)
        gam::Float64     # euler angle (rad)
    )
        return cis(-1im*mp*alpha) * XLALWignerdMatrix(l, mp, m, beta) * cis(-1im*m*gam)
    end

end  # module LAL
