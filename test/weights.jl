@testitem "weights" begin
    using DoubleFloats: Double64

    # The closed forms of the weights given by Waldvogel, evaluated in `BigFloat` once for
    # each `n`, are the reference for every type.  They share no code with the package,
    # which computes the weights by a Fourier transform.
    ϑ(k, n) = k * (big(π) / n)  # k ∈ 0:n
    b(j, n) = j==n/2 ? 1 : 2
    c(k, n) = k%n==0 ? 1 : 2
    Σ(f, r) = sum(f, r; init=zero(BigFloat))
    wᶠ¹(k, n) = (2/big(n)) * (1 - 2Σ(j -> cos(j*ϑ(2k+1, n))/big(4j^2-1), 1:(n÷2)))  # Eq. (2.3a) of Waldvogel; k ∈ 0:n-1
    wᶠ²(k, n) = (4/big(n)) * sin(ϑ(k, n)) * Σ(j -> sin((2j-1)*ϑ(k, n))/big(2j-1), 1:(n÷2))  # Eq. (2.3b) of Waldvogel; k ∈ 0:n
    wᶜᶜ(k, n) = (c(k, n)/big(n)) * (1 - Σ(j -> b(j, n) * cos(2j*ϑ(k, n))/big(4j^2-1), 1:(n÷2)))  # Eq. (4) of Waldvogel; k ∈ 0:n

    for n in [1, 2, 3, 10, 11, 12, 13, 170, 171, 1070, 1071]
        ref¹ = wᶠ¹.(0:n-1, n)
        ref² = wᶠ².(1:n, n+1)
        refᶜᶜ = n ≥ 2 ? wᶜᶜ.(0:n-1, n-1) : nothing
        @testset "$T, n=$n" for T in (Float16, Float32, Float64, Double64, BigFloat)
            # The worst errors measured over these n, in units of eps(T), are 1.2 (Float16),
            # 0.7 (Float32 and Float64), 0.8 (Double64) and 3.5 (BigFloat, where the
            # rounding of the reference itself, at the same precision, is part of the
            # difference).
            ϵ = 10eps(T)
            # That absolute bound is the natural one, since the weights sum to 2, but for
            # `Float16` it exceeds the largest weight once n is in the hundreds, where a
            # vector of zeros would satisfy it.  For the machine floats, the error is
            # therefore also bounded relative to the largest weight.  Measured in units of
            # eps(T), that error is at most 4.5, except for `fejer1` in `Float16`, where it
            # reaches 47 at n = 171.  (For `Double64` and `BigFloat`, the absolute bound is
            # already tiny compared with every weight.)
            ϵᵣ = T <: Base.IEEEFloat ? 10eps(T) : T(Inf)
            agrees(w, ref, ϵᵣ=ϵᵣ) =
                maximum(abs, w .- ref) < ϵ && maximum(abs, w .- ref) < ϵᵣ * maximum(abs, ref)
            # `fejer1` transforms a `Float16` vector with GenericFFT, whose `Float16`
            # arithmetic overflows in forming k² for n above 256, so there it is checked
            # only up to that size.
            if !(T === Float16 && n > 256)
                w = fejer1(n, T)
                @test w isa Vector{T} && length(w) == n
                @test agrees(w, ref¹, T === Float16 ? 100eps(T) : ϵᵣ)
            end
            w = fejer2(n, T)
            @test w isa Vector{T} && length(w) == n
            @test agrees(w, ref²)
            if n ≥ 2
                w = clenshaw_curtis(n, T)
                @test w isa Vector{T} && length(w) == n
                @test agrees(w, refᶜᶜ)
            end
        end
        # The default type is `Float64`
        @test fejer1(n) == fejer1(n, Float64)
        @test fejer2(n) == fejer2(n, Float64)
        n ≥ 2 && @test clenshaw_curtis(n) == clenshaw_curtis(n, Float64)
    end
end

@testitem "FastTransforms" begin
    import FastTransforms
    ϵ = eps()

    for N in [3, 10, 11, 12, 13, 170, 171, 1070, 1071]
        μ1 = FastTransforms.chebyshevmoments1(Float64, N)
        μ2 = FastTransforms.chebyshevmoments2(Float64, N)
        w1 = FastTransforms.fejerweights1(μ1)
        w2 = FastTransforms.fejerweights2(μ2)
        wc = FastTransforms.clenshawcurtisweights(μ1)

        @test w1 ≈ fejer1(N) rtol=ϵ atol=ϵ
        @test w2 ≈ fejer2(N) rtol=ϵ atol=ϵ
        @test wc ≈ clenshaw_curtis(N) rtol=ϵ atol=ϵ
    end
end

@testitem "weights: the number of nodes is validated" begin
    import DoubleFloats: Double64

    # Each rule needs at least one node, and the Clenshaw–Curtis rule, whose nodes include
    # both poles, at least two.  (The node counts for which a rule's buffer would be empty
    # are also the subject of the items in `test/bounds.jl`.)
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        for n ∈ (0, -1, -2)
            @test_throws ArgumentError fejer1(n, T)
            @test_throws "`fejer1` needs at least one node; got n=$n." fejer1(n, T)
        end
        for n ∈ (0, -1, -2)
            @test_throws ArgumentError fejer2(n, T)
            @test_throws "`fejer2` needs at least one node; got n=$n." fejer2(n, T)
        end
        for n ∈ (1, 0, -3)
            @test_throws ArgumentError clenshaw_curtis(n, T)
            @test_throws "`clenshaw_curtis` needs at least two nodes" clenshaw_curtis(n, T)
        end
    end
    @test_throws ArgumentError fejer1(0)
    @test_throws ArgumentError fejer2(0)
    @test_throws ArgumentError fejer2(-1)
    @test_throws ArgumentError clenshaw_curtis(1)
    @test_throws ArgumentError clenshaw_curtis(0)

    # The smallest rules are exact for the polynomials they can integrate: with one node at
    # the equator, and with two nodes placed symmetrically, every weight is the same, and
    # the weights sum to ∫ d(cos θ) = 2
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        ϵ = 10eps(T)
        @test fejer1(1, T) ≈ [2] atol=ϵ
        @test fejer2(1, T) ≈ [2] atol=ϵ
        @test fejer1(2, T) ≈ [1, 1] atol=ϵ
        @test fejer2(2, T) ≈ [1, 1] atol=ϵ
        @test clenshaw_curtis(2, T) ≈ [1, 1] atol=ϵ
        @test eltype(fejer1(2, T)) === eltype(fejer2(2, T)) === eltype(clenshaw_curtis(2, T)) === T
    end

    # The weights are those of ∫₀^π f(θ) sin θ dθ = ∫₋₁¹ f(x) dx, so they sum to 2, and the
    # rules integrate low-degree polynomials in x = cos θ exactly
    for (w, θ) ∈ (
        (fejer1(9), fejer1_rings(9)), (fejer2(9), fejer2_rings(9)),
        (clenshaw_curtis(9), clenshaw_curtis_rings(9)),
    )
        @test sum(w) ≈ 2 atol=10eps()
        @test sum(w .* cos.(θ).^2) ≈ 2/3 atol=10eps()
        @test sum(w .* cos.(θ).^3) ≈ 0 atol=10eps()
    end
end
