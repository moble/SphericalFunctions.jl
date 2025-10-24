"""
This module provides functions for encoding a (small) set of (small) integers into a single
integer, and decoding them back again.  This is useful for storing those integers in data
structures that only accept single integers or floats, such as HWedge.
"""
@testmodule EncodeDecode begin
    encode(i::Int) = i + 50
    encode(i::Rational) = numerator(i) + 50
    function encode(i...)
        output = 0
        max_integers = log10(typemax(output)) ÷ 2
        if length(i) > max_integers
            error("Can only encode up to $max_integers integers into a single integer")
        end
        for (j, v) in enumerate(i)
            if !(-50 ≤ v ≤ 49)
                error("Can only encode integers in -50:49; got $v")
            end
            output += encode(v) * 10^(2j-2)
        end
        output
    end

    function decode(i; pad=2)
        d = digits(i; pad)
        output = Int[]
        for k in 1:2:length(d)
            v = (d[k] + 10d[k+1]) - 50
            push!(output, v)
        end
        return Tuple(output)
    end
end

@testitem "EncodeDecode" setup=[EncodeDecode] begin
    using .EncodeDecode: encode, decode
    for ℓₘₐₓ ∈ (5, 7//2)
        for iᵣ ∈ 1:7
            for m′ ∈ -ℓₘₐₓ:ℓₘₐₓ
                for m ∈ -ℓₘₐₓ:ℓₘₐₓ
                    code = encode(iᵣ, m′, m)
                    (iᵣ₂, m′₂, m₂) = decode(code; pad=6)
                    @test iᵣ == iᵣ₂
                    if ℓₘₐₓ isa Int
                        @test m′ == m′₂
                        @test m == m₂
                    else
                        @test m′ == m′₂ // 2
                        @test m == m₂ // 2
                    end
                end
            end
        end
    end
end
