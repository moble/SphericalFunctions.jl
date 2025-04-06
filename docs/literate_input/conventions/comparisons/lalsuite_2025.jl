md"""
# LALSuite (2025)

!!! info "Summary"
    The `LALSuite` definitions of the spherical harmonics and Wigner's ``d`` and ``D``
    functions agree with the definitions used in the `SphericalFunctions` package.

[`LALSuite` (the LSC Algorithm Library Suite)](@cite LALSuite_2018) is a collection of
routines, comprising the primary official software used by the LIGO-Virgo-KAGRA
Collaboration to detect and characterize gravitational waves.  As far as I can tell, the
ultimate source for all spin-weighted spherical harmonics used in `LALSuite` is the function
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

We begin by directly translating the C code of LALSuite over to Julia code.  There are three
functions that we will want to compare with the definitions in this package:
```c
COMPLEX16 XLALSpinWeightedSphericalHarmonic( REAL8 theta, REAL8 phi, int s, int l, int m );
double XLALWignerdMatrix( int l, int mp, int m, double beta );
COMPLEX16 XLALWignerDMatrix( int l, int mp, int m, double alpha, double beta, double gam );
```
The [original source code](./lalsuite_SphericalHarmonics.md) (as of early 2025) is stored
alongside this file, so we will read it in to a `String` and then apply a series of regular
expressions to convert it to Julia code, parse it and evaluate it to turn it into runnable
Julia.  We encapsulate the formulas in a module so that we can test them against the
`SphericalFunctions` package.

We begin by setting up that module, and introducing a set of basic replacements that would
usually be defined in separate C headers.

"""
using TestItems: @testitem  #hide
@testitem "LALSuite conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module LALSuite

using Printf: @sprintf

const I = im
const LAL_PI = π
const XLAL_EINVAL = "XLAL Error: Invalid arguments"
MIN(a, b) = min(a, b)
gsl_sf_choose(a, b) = binomial(a, b)
pow(a, b) = a^b
cexp(a) = exp(a)
cpolar(a, b) = a * cis(b)
macro XLALPrError(msg, args...)
    quote
        @error @sprintf($msg, $(args...))
    end
end
#+

# Next, we simply read the source file into a string.
lalsource = read(joinpath(@__DIR__, "lalsuite_SphericalHarmonics.c"), String)
#+

# Now we define a series of replacements to apply to the C code to convert it to Julia code.
# Note that some of these will be quite specific to this particular file, and may not be
# generally applicable.
replacements = (
    ## Deal with newlines in the middle of an assignment
    r"( = .*[^;]\s*)\n" => s"\1",

    ## Remove a couple old, unused functions
    r"(?ms)XLALScalarSphericalHarmonic.*?\n}" => "# Removed",
    r"(?ms)XLALSphHarm.*?\n}" => "# Removed",

    ## Remove type annotations
    r"COMPLEX16 ?" => "",
    r"REAL8 ?" => "",
    r"INT4 ?" => "",
    r"int ?" => "",
    r"double ?" => "",

    ## Translate comments
    "/*" => "#=",
    "*/" => "=#",

    ## Brackets
    r" ?{" => "",
    r"}.*(\n *else)" => s"\1",
    r"} *else" => "else",
    r"^}" => "",
    "}" => "end",

    ## Flow control
    r"( *if.*);"=>s"\1 end",  ## one-line `if` statements
    "for( s=0; n-s >= 0; s++ )" => "for s=0:n",
    "else if" => "elseif",
    r"(?m)  break;\n *\n *case(.*?):" => s"elseif m == \1",
    r"(?m)  break;\n\s*case(.*?):" => s"elseif m == \1",
    r"(?m)  break;\n *\n *default:" => "else",
    r"(?m)  break;\n *default:" => "else",
    r"(?m)switch.*?\n *\n( *)case(.*?):" => s"\n\1if m == \2",
    r"\n *break;" => "",
    r"(?m)( *ans = fac;)\n" => s"\1\n  end\n",

    ## Deal with ugly C declarations
    "f1 = (x-1)/2.0, f2 = (x+1)/2.0" => "f1 = (x-1)/2.0; f2 = (x+1)/2.0",
    "sum=0, val=0" => "sum=0; val=0",
    "a=0, lam=0" => "a=0; lam=0",
    r"\n *fac;" => "",
    r"\n *ans;" => "",
    r"\n *gslStatus;" => "",
    r"\n *gsl_sf_result pLm;" => "",
    r"\n ?XLAL" => "\nfunction XLAL",

    ## Differences in Julia syntax
    "++" => "+=1",
    ".*" => ". *",
    "./" => ". /",
    ".+" => ". +",
    ".-" => ". -",

    ## Deal with random bad syntax
    "if (m)" => "if m != 0",
    "case 4:" => "elseif m == 4",
    "XLALPrError" => "@XLALPrError",
    "__func__" => "\"\"",
)
#+

# And we apply the replacements to the source code to convert it to Julia code.  Note that
# we apply them successively, even though `replace` can handle multiple "simultaneous"
# replacements, because the order of replacements is important.
for (pattern, replacement) in replacements
    global lalsource = replace(lalsource, pattern => replacement)
end
#+

# Finally, we just parse and evaluate the code to turn it into a runnable Julia, and we are
# done defining the module
eval(Meta.parseall(lalsource))

end  # module LALSuite
#+


# ## Tests
#
# We can now test the `LALSuite` functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 100eps()
#+

# The spin-weighted spherical harmonics are defined explicitly, but only for
s = -2
#+
# and only up to
ℓₘₐₓ = 8
#+
# so we only test up to that point.
for (θ, ϕ) ∈ θϕrange()
    for (ℓ, m) ∈ ℓmrange(abs(s), ℓₘₐₓ)
        @test LALSuite.XLALSpinWeightedSphericalHarmonic(θ, ϕ, s, ℓ, m) ≈
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
        @test LALSuite.XLALWignerdMatrix(ℓ, m′, m, β) ≈
            SphericalFunctions.d(ℓ, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# We can see more-or-less by inspection that the code defines the ``D`` matrix in agreement
# with our convention, the key line being
# ```c
# cexp( -(1.0I)*mp*alpha ) * XLALWignerdMatrix( l, mp, m, beta ) * cexp( -(1.0I)*m*gam );
# ```
for (α,β,γ) ∈ αβγrange()
    for (ℓ, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test LALSuite.XLALWignerDMatrix(ℓ, m′, m, α, β, γ) ≈
            conj(SphericalFunctions.D(ℓ, m′, m, α, β, γ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Now, just to remind ourselves, we will be changing the convention for ``D`` soon
@test_broken false  # We haven't flipped the conjugation of D yet
#+

# These successful tests show that the spin-weighted spherical harmonics and the Wigner
# ``d`` and ``D`` matrices defined in `LALSuite` agree with the corresponding functions
# defined by the `SphericalFunctions` package.

end  #hide
