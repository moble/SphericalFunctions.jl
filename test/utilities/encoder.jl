"""
This module provides functions for encoding a (small) set of (small) integers into a single
integer, and decoding them back again.  This is useful for storing those integers in data
structures that only accept single integers or floats, such as HWedge.
"""
@testmodule EncodeDecode begin
    encode(i::Int) = i
    encode(i::Rational) = numerator(i)
    function encode(i...)
        output = 0
        for (j, v) in enumerate(i)
            if abs(v) >= 10
                error("Can only encode integers with absolute value less than 10")
            end
            output += (v < 0) * 10^(2j-2)
            output += encode(abs(v)) * 10^(2j-1)
        end
        output
    end

    function decode(i; pad=2)
        d = digits(i; pad)
        output = Int[]
        for k in 1:(length(d) ÷ 2)
            v = d[2k] * (d[2k-1] == 1 ? -1 : 1)
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
