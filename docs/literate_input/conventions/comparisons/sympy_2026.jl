md"""
# SymPy (2026)

!!! info "Summary"
    SymPy's spherical harmonics `Ynm` agree with the spherical harmonics used in the
    `SphericalFunctions` package.  Its Wigner ``D`` matrix `wigner_d` is related to ours by
    ```math
    \mathtt{wigner\_d(J, α, β, γ)[J-m', J-m]}
    = 𝔇^{(J)}_{-m',-m}(α, β, γ)
    = (-1)^{m'-m}\, \overline{𝔇^{(J)}_{m',m}(α, β, γ)},
    ```
    which is [Wigner's](@ref "Wigner (1959)") convention.  SymPy documents its ``D`` as
    following Edmonds, but its transcription of Edmonds is not exact.

[SymPy](@cite SymPy_2025) is the standard symbolic-algebra library for Python.  The
statements below refer to the source code of version 1.14.0 (released 2025-04-27), accessed
2026-09-08.

## Spherical harmonics

The class `sympy.functions.special.spherical_harmonics.Ynm` is documented as
```math
Y_n^m(θ, φ)
:=
\sqrt{\frac{(2n+1)\,(n-m)!}{4π\,(n+m)!}}\; e^{imφ}\, P_n^m(\cos θ),
```
"where ``n \geq 0`` is an integer, ``m`` is an integer with ``-n \leq m \leq n``, ``θ \in [0,
π]`` is the polar angle, and ``φ \in [0, 2π]`` is the azimuthal angle"; its `expand_func`
method evaluates exactly this expression in terms of `assoc_legendre`, and the code also
records the relation `Ynm(n, -m, theta, phi) = (-1)**m*exp(-2*I*m*phi)*Ynm(n, m, theta,
phi)`, i.e., ``Y_n^{-m} = (-1)^m \overline{Y_n^m}``.  The function
`sympy.functions.special.polynomials.assoc_legendre` is documented as
```math
P_n^m(x) = (-1)^m (1-x^2)^{m/2} \frac{d^m P_n(x)}{dx^m},
```
which includes the Condon–Shortley phase, and for negative ``m`` its code evaluates
```python
S.NegativeOne**(-m) * (factorial(m + n)/factorial(n - m)) * assoc_legendre(n, -m, x)
```
that is, ``P_n^{-|m|}(x) = (-1)^{|m|} \frac{(n-|m|)!}{(n+|m|)!} P_n^{|m|}(x)``.  We therefore
expect agreement with [our spherical harmonics](@ref summary_spherical_harmonics).

## Wigner ``D`` matrix

The functions `wigner_d_small` and `wigner_d` in `sympy/physics/wigner.py` are documented as
returning the matrices of
```math
\mathcal{d}_β = \exp\left( \frac{iβ}{\hbar} J_y \right)
\qquad \text{and} \qquad
\mathcal{D}_{αβγ} =
\exp\left( \frac{iα}{\hbar} J_z\right)
\exp\left( \frac{iβ}{\hbar} J_y\right)
\exp\left( \frac{iγ}{\hbar} J_z\right),
```
"such that ``d^{(J)}_{m',m}(β) = \mathtt{wigner\_d\_small(J, beta)[J-mprime, J-m]}``" and
"``\mathcal{D}^{(J)}_{m',m}(α, β, γ) = \mathtt{wigner\_d(J, alpha, beta, gamma)[J-mprime,
J-m]}``".  The docstrings say that the components are calculated using [Edmonds](@cite
Edmonds_2016) Eq. (4.1.15) for ``d`` and Eq. (4.1.12) for ``D``, "however note that angles
alpha and gamma are swapped".  The code is
```python
def wigner_d_small(J, beta):
    M = [J-i for i in range(2*J+1)]
    d = zeros(2*J+1)
    for i, Mi in enumerate(M):
        for j, Mj in enumerate(M):
            sigmamax = min([J-Mi, J-Mj])
            sigmamin = max([0, -Mi-Mj])
            dij = sqrt(factorial(J+Mi)*factorial(J-Mi) /
                       factorial(J+Mj)/factorial(J-Mj))
            terms = [(-1)**(J-Mi-s) *
                     binomial(J+Mj, J-Mi-s) *
                     binomial(J-Mj, s) *
                     cos(beta/2)**(2*s+Mi+Mj) *
                     sin(beta/2)**(2*J-2*s-Mj-Mi)
                     for s in range(sigmamin, sigmamax+1)]
            d[i, j] = dij*Add(*terms)
    return ImmutableMatrix(d)

def wigner_d(J, alpha, beta, gamma):
    d = wigner_d_small(J, beta)
    M = [J-i for i in range(2*J+1)]
    D = [[exp(I*Mi*alpha)*d[i, j]*exp(I*Mj*gamma)
          for j, Mj in enumerate(M)] for i, Mi in enumerate(M)]
    return ImmutableMatrix(D)
```
So `wigner_d_small` is exactly [Edmonds' ``d``](@ref "Edmonds (1960)") of Eq. (4.1.15) —
which is the transpose of ours — and `wigner_d` attaches the phases ``e^{+im'α}`` and
``e^{+imγ}``.  This is *not* Edmonds' Eq. (4.1.12), which reads ``e^{im'γ} d^{(j)}_{m'm}(β)
e^{imα}``; SymPy has swapped the roles of ``α`` and ``γ`` (as its docstring admits), and the
operator it quotes, ``e^{iαJ_z} e^{iβJ_y} e^{iγJ_z}``, is likewise Edmonds' Eq. (4.1.9) with
``α`` and ``γ`` interchanged.  The result is neither Edmonds' matrix nor its conjugate.  Using
``d_{m'm}(β)|_{\text{Edmonds}} = d_{m,m'}(β) = (-1)^{m'-m} d_{m',m}(β)``, we find
```math
\mathtt{wigner\_d}_{m'm}(α, β, γ)
=
e^{im'α}\, d_{m,m'}(β)\, e^{imγ}
=
(-1)^{m'-m}\, \overline{𝔇^{(J)}_{m',m}(α, β, γ)}
=
𝔇^{(J)}_{-m',-m}(α, β, γ),
```
which happens to be [Wigner's original convention](@ref "Wigner (1959)"), and also
[Mathematica's](@ref "Mathematica (2026)").  The tests below confirm this.

## Implementing formulas

We begin by transcribing the SymPy code into Julia.  We encapsulate the formulas in a module
so that we can test them against the `SphericalFunctions` package.  A second test — which
is skipped by default because it needs a Python installation, and is run by the scheduled CI
workflow — calls the actual SymPy functions through `PythonCall`.
"""

using TestItems: @testitem  #hide
@testitem "SymPy conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module SymPy
#+

# We'll use some predefined utilities to make the code look more like the equations,
# including `∂ⁿ`, which computes derivatives symbolically so that we can transcribe the
# Legendre formulas literally.
import ..ConventionsUtilities: 𝒾, ❗, ∂ⁿ
#+

# `legendre` (Rodrigues' formula) and `assoc_legendre`, with SymPy's rule for negative
# ``m``:
function legendre(n, x)
    1 / (2^n * (n)❗) * ∂ⁿ(x -> (x^2 - 1)^n, n)(x)
end
function assoc_legendre(n, m, x)
    if m < 0
        return (-1)^(-m) * (n+m)❗ / (n-m)❗ * assoc_legendre(n, -m, x)
    end
    (-1)^m * (1 - x^2)^(m/2) * ∂ⁿ(x -> legendre(n, x), m)(x)
end
#+

# `Ynm`.  We capture the floating-point type `T` to ensure that we don't lose precision when
# converting π and the factorials to floating-point numbers.
function Ynm(n, m, θ::T, φ::T) where {T<:Real}
    √T((2n+1) * (n-m)❗ / (4big(π) * (n+m)❗)) * exp(𝒾 * m * φ) * T(assoc_legendre(n, m, cos(θ)))
end
#+

# `wigner_d_small` and `wigner_d`, transcribed element by element.  Row `i` and column `j`
# of the Python matrices correspond to ``M_i = J - i`` and ``M_j = J - j``, so the element
# `[J-mprime, J-m]` has ``M_i = m'`` and ``M_j = m``:
function wigner_d_small(J, Mᵢ, Mⱼ, β::T) where {T<:Real}
    σₘₐₓ = min(J-Mᵢ, J-Mⱼ)
    σₘᵢₙ = max(0, -Mᵢ-Mⱼ)
    dᵢⱼ = √T((J+Mᵢ)❗ * (J-Mᵢ)❗ / (J+Mⱼ)❗ / (J-Mⱼ)❗)
    dᵢⱼ * sum(
        (-1)^(J-Mᵢ-s) * T(binomial(big(J+Mⱼ), J-Mᵢ-s) * binomial(big(J-Mⱼ), s))
        * cos(β/2)^(2s+Mᵢ+Mⱼ) * sin(β/2)^(2J-2s-Mⱼ-Mᵢ)
        for s ∈ σₘᵢₙ:σₘₐₓ;
        init=zero(T)
    )
end
function wigner_d(J, Mᵢ, Mⱼ, α, β, γ)
    exp(𝒾 * Mᵢ * α) * wigner_d_small(J, Mᵢ, Mⱼ, β) * exp(𝒾 * Mⱼ * γ)
end
#+

end  # module SymPy
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 1000eps()
#+

# We only test up to
ℓₘₐₓ = 4
#+
# because the formulas are slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.  For the same reason
# we use a modest grid of Euler angles.
αβγs = αβγrange(Float64, 5)
#+

# `Ynm` agrees with ours:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test SymPy.Ynm(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# `wigner_d_small` is the transpose of our ``d``, and `wigner_d` is ``𝔇_{-m',-m}``:
for β ∈ βrange()
    for (J, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test SymPy.wigner_d_small(J, m′, m, β) ≈ ConventionsUtilities.d(J, m, m′, β) atol=ϵₐ rtol=ϵᵣ
    end
end
for (α, β, γ) ∈ αβγs
    for (J, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test SymPy.wigner_d(J, m′, m, α, β, γ) ≈ ConventionsUtilities.D(J, -m′, -m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test SymPy.wigner_d(J, m′, m, α, β, γ) ≈ (-1)^(m′-m) * conj(ConventionsUtilities.D(J, m′, m, α, β, γ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that SymPy's spherical harmonics agree with ours, and that its
# Wigner ``D`` matrix is ``(-1)^{m'-m}`` times the complex conjugate of ours.

end  #hide

md"""
## Cross-check against the actual SymPy library

The test above transcribes the SymPy source.  To make sure the transcription is faithful,
the following test calls `sympy.physics.wigner.wigner_d` and `sympy.functions.Ynm`
themselves, via `PythonCall`, evaluating the symbolic results numerically.  It is tagged
`:python` (and `:skipci`) so that it runs only when explicitly requested — e.g., with `julia
--project=. scripts/test.jl :python` — or in the scheduled CI workflow, which installs SymPy
through `CondaPkg`.
"""

@testitem "SymPy cross-check" tags=[:python, :skipci] setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide
using PythonCall
sympy = pyimport("sympy")
wigner = pyimport("sympy.physics.wigner")
sympy_version = pyconvert(String, sympy.__version__)
@info "Cross-checking against SymPy version $sympy_version"
ϵₐ = 1e-12
ϵᵣ = 1e-12
for (θ, ϕ) ∈ θϕrange(Float64, 3)
    for (ℓ, m) ∈ ℓmrange(3)
        Y_sympy = pyconvert(ComplexF64, pybuiltins.complex(sympy.N(sympy.Ynm(ℓ, m, θ, ϕ).expand(func=true), 20)))
        @test Y_sympy ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
for (α, β, γ) ∈ αβγrange(Float64, 1)
    for J ∈ 0:3
        D_sympy = wigner.wigner_d(J, α, β, γ)
        for m′ ∈ -J:J, m ∈ -J:J
            element = pyconvert(ComplexF64, pybuiltins.complex(sympy.N(D_sympy[J-m′, J-m], 20)))
            @test element ≈ ConventionsUtilities.D(J, -m′, -m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
        end
    end
end
end  #hide
