md"""
# LALSuite (2025)

!!! info "Summary"
    The LALSuite definitions of the spherical harmonics and Wigner's ``d`` and ``D``
    functions agree with the definitions used in the `SphericalFunctions` package.

"""
md"""
[LALSuite (LSC Algorithm Library Suite)](@cite LALSuite_2018) is a collection of software
routines, comprising the primary official software used by the LIGO-Virgo-KAGRA
Collaboration to detect and characterize gravitational waves.  As far as I can tell, the
ultimate source for all spin-weighted spherical harmonic values used in LALSuite is the
function
[`XLALSpinWeightedSphericalHarmonic`](https://git.ligo.org/lscsoft/lalsuite/-/blob/6e653c91b6e8a6728c4475729c4f967c9e09f020/lal/lib/utilities/SphericalHarmonics.c),
which cites the NINJA paper [AjithEtAl_2011](@cite) as its source.  Unfortunately, it cites
version *1*, which contained a serious error, using ``\tfrac{\cos\iota}{2}`` instead of
``\cos \tfrac{\iota}{2}`` and similarly for ``\sin``.  This error was corrected in version
2, but the citation was not updated.  Nonetheless, it appears that the actual code is
consistent with the *corrected* versions of the NINJA paper.

They also (quite separately) define Wigner's ``D`` matrices in terms of the ``d`` matrices,
which are — in turn — defined in terms of Jacobi polynomials.  For all of these, they cite
Wikipedia (despite the fact that the NINJA paper defined the spin-weighted spherical
harmonics in terms of the ``d`` matrices).  Nonetheless, the definitions in the code are
consistent with the definitions in the NINJA paper, which are consistent with the
definitions in the `SphericalFunctions` package.


## Implementing formulas

We will call the python module `lal` directly, but there are some minor inconveniences to
deal with first.  We have to install the `lalsuite` package, but we don't want all its
dependencies, so we run `python -m pip install --no-deps lalsuite`.  Then, we have to
translate to native Julia types, so we'll just write three quick and easy wrappers.  We
encapsulate the formulas in a module so that we can test them against the
`SphericalFunctions` package.

"""
using TestItems: @testitem  #hide
@testitem "LALSuite conventions B" setup=[ConventionsSetup, Utilities] begin  #hide

module LALSuite

include("conventions_install_lalsuite.jl")
import PythonCall

const lal = PythonCall.pyimport("lal")

function SpinWeightedSphericalHarmonic(θ, ϕ, s, ℓ, m)
    PythonCall.pyconvert(
        ComplexF64,
        lal.SpinWeightedSphericalHarmonic(θ, ϕ, s, ℓ, m),
    )
end
function WignerdMatrix(ℓ, m′, m, β)
    PythonCall.pyconvert(
        Float64,
        lal.WignerdMatrix(ℓ, m′, m, β),
    )
end
function WignerDMatrix(ℓ, m′, m, α, β, γ)
    PythonCall.pyconvert(
        ComplexF64,
        lal.WignerDMatrix(ℓ, m′, m, α, β, γ),
    )
end

end  # module LALSuite
#+


# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 100eps()
#+

# The spin-weighted spherical harmonics are defined explicitly, but only for
s = -2
#+
# and only up to
ℓₘₐₓ = 3
#+
# so we only test up to that point.
for (θ, ϕ) ∈ θϕrange()
    for (ℓ, m) ∈ ℓmrange(abs(s), ℓₘₐₓ)
        @test LALSuite.SpinWeightedSphericalHarmonic(θ, ϕ, s, ℓ, m) ≈
            SphericalFunctions.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Now, the Wigner ``d`` matrices are defined generally, but we only need to test up to
ℓₘₐₓ = 4
#+
# because the formulas are fairly inefficient and inaccurate, and this will be sufficient to
# sort out any sign or normalization differences, which are the most likely sources of
# error.
for β ∈ βrange()
    for (ℓ, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test LALSuite.WignerdMatrix(ℓ, m′, m, β) ≈
            SphericalFunctions.d(ℓ, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# We can see more-or-less by inspection that the code defines the ``D`` matrix in agreement
# with our convention, the key line being
# ```c
# cexp( -(1.0I)*mp*alpha ) * XLALWignerdMatrix( l, mp, m, beta ) * cexp( -(1.0I)*m*gam );
# ```
# And because of the higher dimensionality of the space in which to test, we want to
# restrict the range of the tests to avoid excessive computation.  We will test up to
ℓₘₐₓ = 2
#+
# because the space of options for disagreement is smaller.
for (α,β,γ) ∈ αβγrange()
    for (ℓ, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test LALSuite.WignerDMatrix(ℓ, m′, m, α, β, γ) ≈
            conj(SphericalFunctions.D(ℓ, m′, m, α, β, γ)) atol=ϵₐ rtol=ϵᵣ
    end
end
@test_broken false  # We haven't flipped the conjugation of D yet

#+

# These successful tests show that the spin-weighted spherical harmonics and the Wigner
# ``d`` and ``D`` matrices defined in LALSuite agree with the corresponding functions
# defined by the `SphericalFunctions` package.

end  #hide
