md"""
# ``L_j`` and ``R_j`` with Euler angles

This package defines the angular-momentum operators ``L_j`` and ``R_j`` in terms of elements
of the Lie group and algebra:
```math
L_𝐮 f(𝐑) = \left. i \frac{d}{dϵ}\right|_{ϵ=0}
f\left(e^{-ϵ 𝐮/2}\, 𝐑\right)
\qquad \text{and} \qquad
R_𝐮 f(𝐑) = -\left. i \frac{d}{dϵ}\right|_{ϵ=0}
f\left(𝐑\, e^{-ϵ 𝐮/2}\right),
```
This is certainly the natural realm for these operators, but it is not the common one.  In
particular, virtually all textbooks and papers on the subject define these operators in
terms of the standard spherical coordinates on the 2-sphere, rather than quaternions or even
Euler angles.  In particular, the standard forms are essentially always given in terms of
the Cartesian basis, as in ``L_x``, ``L_y``, and ``L_z`` — though some times the first two
are expressed as ``L_{\pm} = L_x \pm i L_y``.

Here, we will use SymPy to just grind through the algebra of expressing the angular-momentum
operators in terms of Euler angles, including evaluating the commutators in that form, and
further reduce them to operators in terms of spherical coordinates.  We will find a couple
important results that help make contact with more standard expressions:

  1. Our results for ``L_x``, ``L_x``, and ``L_z`` in spherical coordinates agree with
     standard expressions.
  2. The commutators obey ``[L_x, L_y] = i L_z`` and cyclic permutations, in agreement with
     the standard expressions.
  3. We also find ``[R_x, R_y] = i R_z`` and cyclic permutations.
  4. We can explicitly compute ``[L_i, R_j] = 0``, as expected.
  5. Using the natural extension of Goldberg et al.'s SWSHs to include ``γ``, we can
     see that the natural spin-weight operator is ``R_z = i \partial_γ``.  Thus, we
     define ``R_z = s`` for a function with spin weight ``s``.
  6. The spin-raising operator for ``R_z`` is ``\eth = R_x + i R_y``; the spin-lowering
     operator is ``\bar{\eth} = R_x - i R_y``.

## Analytical groundwork

We start by defining a new set of Euler angles according to
```math
𝐑_{α', β', γ'}
= e^{-ϵ 𝐮 / 2} 𝐑_{α, β, γ}
\qquad \text{or} \qquad
𝐑_{α', β', γ'}
= 𝐑_{α, β, γ} e^{-ϵ 𝐮 / 2}
```
where ``𝐮`` will be each of the basis quaternions, and each of ``α'``,
``β'``, and ``γ'`` is a function of ``α``, ``β``, ``γ``, and
``ϵ``.  Then, we note that the chain rule tells us that
```math
\frac{\partial}{\partial ϵ}
=
\frac{\partial α'}{\partial ϵ} \frac{\partial}{\partial α'}
+ \frac{\partial β'}{\partial ϵ} \frac{\partial}{\partial β'}
+ \frac{\partial γ'}{\partial ϵ} \frac{\partial}{\partial γ'},
```
which we will use to convert the general expression for the angular-momentum operators in
terms of ``\partial_ϵ`` into an expression in terms of derivatives with respect to
these new Euler angles:
```math
\begin{align}
  L_j f(𝐑_{α, β, γ})
  &=
  \left. i \frac{\partial} {\partial ϵ} f \left( e^{-ϵ 𝐞_j / 2}
  𝐑_{α, β, γ} \right) \right|_{ϵ=0}
  \\
  &=
  i \left[ \left(
      \frac{\partial α'}{\partial ϵ} \frac{\partial}{\partial α'}
      + \frac{\partial β'}{\partial ϵ} \frac{\partial}{\partial β'}
      + \frac{\partial γ'}{\partial ϵ} \frac{\partial}{\partial γ'}
  \right) f \left(α', β', γ'\right) \right]_{ϵ=0}
  \\
  &=
  i \left[ \left(
      \frac{\partial α'}{\partial ϵ} \frac{\partial}{\partial α}
      + \frac{\partial β'}{\partial ϵ} \frac{\partial}{\partial β}
      + \frac{\partial γ'}{\partial ϵ} \frac{\partial}{\partial γ}
  \right) f \left(α, β, γ\right) \right]_{ϵ=0},
\end{align}
```
or for ``R_j``:
```math
\begin{align}
  R_j f(𝐑_{α, β, γ})
  &=
  -\left. i \frac{\partial} {\partial ϵ} f \left( 𝐑_{α, β, γ}
  e^{-ϵ 𝐞_j / 2} \right) \right|_{ϵ=0}
  \\
  &=
  -i \left[ \left(
      \frac{\partial α'}{\partial ϵ} \frac{\partial}{\partial α}
      + \frac{\partial β'}{\partial ϵ} \frac{\partial}{\partial β}
      + \frac{\partial γ'}{\partial ϵ} \frac{\partial}{\partial γ}
  \right) f \left(α, β, γ\right) \right]_{ϵ=0}.
\end{align}
```

So the objective is to find the new Euler angles, differentiate with respect to
``ϵ``, and then evaluate at ``ϵ = 0``.  We do this by first multiplying
``𝐑_{α, β, γ}`` and ``e^{-ϵ 𝐮 / 2}`` in the desired
order, then expanding the results in terms of its quaternion components, and then computing
the new Euler angles in terms of those components according to the usual expression.

"""

#src # Do this first just to hide stdout of the conda installation step.
#src # Note that we can't just use `#hide` because that still shows stdout.
# ````@setup euler_angular_momentum
# import SymPyPythonCall
# ````

# ## Computational infrastructure
# We'll use SymPy (via Julia) since `Symbolics.jl` isn't very good at trig yet.
import Memoization: @memoize
import LaTeXStrings: @L_str, LaTeXString
import Quaternionic: Quaternionic, Quaternion, components
import SymPyPythonCall
import SymPyPythonCall: sympy, symbols, sqrt, sin, cos, tan, acos, atan, latex
const expand_trig = sympy.expand_trig
const Derivative = sympy.Derivative
const π = sympy.pi
const I = sympy.I
nothing  #hide

# Define symbols we will use throughout
α, β, γ, θ, ϕ, ϵ = symbols("α β γ θ ϕ ϵ", real=true, positive=true)
f = symbols("f", cls=SymPyPythonCall.sympy.o.Function)
nothing  #hide

# Reinterpret the quaternion basis elements for compatibility with SymPy.  (`Quaternionic`
# defines the basis with `Bool` components, but SymPy can't handle that.)
const 𝐢 = Quaternion{Int}(Quaternionic.𝐢)
const 𝐣 = Quaternion{Int}(Quaternionic.𝐣)
const 𝐤 = Quaternion{Int}(Quaternionic.𝐤)
nothing  #hide

# Next, we define functions to compute the Euler components of the left and right operators
@memoize function 𝒪(u, side)
    ## Substitutions that sympy doesn't make but we want
    subs = Dict(
        cos(β)/sin(β) => 1/tan(β),
        sqrt(1 - cos(β))*sqrt(cos(β) + 1) => sin(β)
    )

    ## Define the essential quaternions
    e = cos(ϵ/2) + u * sin(-ϵ/2)
    R₀ = Quaternion(sympy.simplify.(sympy.expand.(components(
        (cos(α/2) + 𝐤 * sin(α/2)) * (cos(β/2) + 𝐣 * sin(β/2)) * (cos(γ/2) + 𝐤 * sin(γ/2))
    ))))

    ## Extract the (simplified) components of the product
    w, x, y, z = sympy.simplify.(sympy.expand.(components(
        side == :left ? e * R₀ : R₀ * e
    )))

    ## Convert back to Euler angles
    α′ = (atan(z/w) + atan(-x/y)).expand().simplify()
    β′ = (2*acos(sqrt(w^2 + z^2) / sqrt(w^2 + x^2 + y^2 + z^2))).expand().simplify()
    γ′ = (atan(z/w) - atan(-x/y)).expand().simplify()

    ## Differentiate with respect to ϵ, set ϵ to 0, and simplify
    ∂α′∂ϵ = expand_trig(Derivative(α′, ϵ).doit().subs(ϵ, 0).expand().simplify().subs(subs))
    ∂β′∂ϵ = expand_trig(Derivative(β′, ϵ).doit().subs(ϵ, 0).expand().simplify().subs(subs))
    ∂γ′∂ϵ = expand_trig(Derivative(γ′, ϵ).doit().subs(ϵ, 0).expand().simplify().subs(subs))

    return ∂α′∂ϵ, ∂β′∂ϵ, ∂γ′∂ϵ
end

## Note that we are not including the factor of ``i`` here; for simplicity, we will insert
## it manually when displaying the results below, and when applying these operators to
## functions (`Lx` and related definitions below).
L(u) = 𝒪(u, :left)
R(u) = 𝒪(u, :right)
nothing  #hide


# We need a couple quick helper macros to display the results.
#md # The details are boring, but you can expand the source code to see them.
#md # ```@raw html
#md # <details>
#md # <summary>
#md # Click here to expand the source code for display macros
#md # </summary>
#md # ```
macro display(expr)
    op = string(expr.args[1])
    arg = Dict(:𝐢 => "x", :𝐣 => "y", :𝐤 => "z")[expr.args[2]]
    if op == "L"
        quote
            ∂α′∂ϵ, ∂β′∂ϵ, ∂γ′∂ϵ = latex.($expr)  # Call expr; format results as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = i\left[
                %$(∂α′∂ϵ) \frac{\partial}{\partial α}
                + %$(∂β′∂ϵ) \frac{\partial}{\partial β}
                + %$(∂γ′∂ϵ) \frac{\partial}{\partial γ}
            \right]"""  # Display the result in LaTeX form
        end
    else
        quote
            ∂α′∂ϵ, ∂β′∂ϵ, ∂γ′∂ϵ = latex.($expr)  # Call expr; format results as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = -i\left[
                %$(∂α′∂ϵ) \frac{\partial}{\partial α}
                + %$(∂β′∂ϵ) \frac{\partial}{\partial β}
                + %$(∂γ′∂ϵ) \frac{\partial}{\partial γ}
            \right]"""  # Display the result in LaTeX form
        end
    end
end
nothing  #hide

# And we'll need another for the angular-momentum operators in standard ``𝕊²`` form.
conversion(∂) = ∂.subs(Dict(α => ϕ, β => θ, γ => 0)).simplify()
macro display2(expr)
    op = string(expr.args[1])
    element = expr.args[2]
    arg = Dict(:𝐢 => "x", :𝐣 => "y", :𝐤 => "z", :+ => "+", :- => "-")[element]
    @info element
    if op == "L" && arg ∈ ("+", "-")
        quote
            ∂φ′∂ϵ, ∂ϑ′∂ϵ, ∂γ′∂ϵ = (
                (
                    ($element)(I * $conversion(i), -$conversion(j)).rewrite(exp)
                    / exp(($element)(I) * ϕ)
                ).simplify()
                for (i, j) ∈ zip(L(𝐢), L(𝐣))
            )
            expr = $op * "_" * $arg  # Standard form of the operator
            expsign = ($arg=="+" ? "" : "-")
            L"""%$expr = e^{%$expsign i ϕ} \left[
                %$(∂ϑ′∂ϵ) \frac{\partial}{\partial θ}
                + %$(∂φ′∂ϵ) \frac{\partial}{\partial ϕ}
            \right]"""  # Display the result in LaTeX form
        end
    elseif op == "L"
        quote
            ∂φ′∂ϵ, ∂ϑ′∂ϵ, ∂γ′∂ϵ = latex.($conversion.($expr))  # Call expr; format as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = i\left[
                %$(∂ϑ′∂ϵ) \frac{\partial}{\partial θ}
                + %$(∂φ′∂ϵ) \frac{\partial}{\partial ϕ}
            \right]"""  # Display the result in LaTeX form
        end
    else
        quote
            ∂φ′∂ϵ, ∂ϑ′∂ϵ, ∂γ′∂ϵ = latex.($conversion.($expr))  # Call expr; format as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = -i\left[
                %$(∂ϑ′∂ϵ) \frac{\partial}{\partial θ}
                + %$(∂φ′∂ϵ) \frac{\partial}{\partial ϕ}
                + %$(∂γ′∂ϵ) \frac{\partial}{\partial γ}
            \right]"""  # Display the result in LaTeX form
        end
    end
end
nothing  #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Full expressions on ``𝕊³``
# Finally, we can actually compute the Euler components of the angular momentum operators.

#md # ### ``L`` operators in terms of Euler angles
@display L(𝐢)
#-
@display L(𝐣)
#-
@display L(𝐤)
#-
#md # ### ``R`` operators in terms of Euler angles
@display R(𝐢)
#-
@display R(𝐣)
#-
@display R(𝐤)

# In their description of the Wigner 𝔇 functions as wave functions of a rigid symmetric
# top, [Varshalovich_1988](@citet) provide equivalent expressions in Eqs. (6) and (7) of
# their Sec. 4.2 — except that ``R_x`` and ``R_z`` have the wrong signs.
# [Wikipedia](https://en.wikipedia.org/wiki/Wigner_D-matrix#Properties_of_the_Wigner_D-matrix),
# meanwhile, provides equivalent expressions, except that their ``\hat{𝒫}`` has
# (consistently) the opposite sign to ``R`` defined here.

# Note that the Wikipedia convention is actually entirely sensible — maybe more sensible
# than the one we use.  In that convention ``\hat{J}`` is in the inertial frame, whereas
# ``\hat{P}`` is exactly that operator in the body-fixed frame.  In our notation, we have
# ```math
# \begin{align}
# R_𝐮 f(𝐑)
# &=
# -\left. i \frac{d}{dϵ}\right|_{ϵ=0} f\left(𝐑\, e^{-ϵ 𝐮/2}\right) \\
# &=
# -\left. i \frac{d}{dϵ}\right|_{ϵ=0}
#   f\left(𝐑\, e^{-ϵ 𝐮/2}\, 𝐑^{-1}\, 𝐑\right) \\
# &=
# -\left. i \frac{d}{dϵ}\right|_{ϵ=0}
#   f\left(e^{-ϵ 𝐑\, 𝐮\, 𝐑^{-1}/2}\, 𝐑\right) \\
# &=
# -L_{𝐑\, 𝐮\, 𝐑^{-1}} f(𝐑),
# \end{align}
# ```
# which says that ``R`` is just the *negative of* the ``L`` operator transformed to the
# body-fixed frame.  That negative sign is slightly unnatural, but the reason we choose to
# define ``R`` in this way is for its more natural connection to the literature on
# spin-weighted spherical functions.  Also, ``R`` defined here obeys the same commutation
# relations as the standard angular-momentum operators, whereas the Wikipedia convention
# leads to "anomalous" commutation relations with an extra minus sign.

# ## Commutators

# We can also compute the commutators of the angular momentum operators, as derived above.
# First, we define the operators acting on functions of the Euler angles.  Note that this
# function differs from the one above because it explicitly takes the function ``f`` and the
# Euler angles as arguments — which will be necessary to compute the commutators.
function 𝒪(u, side, f, α, β, γ)
    let O = 𝒪(u, side)
        (side==:left ? I : -I) * (
            O[1] * f(α, β, γ).diff(α)
            + O[2] * f(α, β, γ).diff(β)
            + O[3] * f(α, β, γ).diff(γ)
        )
    end
end

Lx(f, α, β, γ) = 𝒪(𝐢, :left, f, α, β, γ)
Ly(f, α, β, γ) = 𝒪(𝐣, :left, f, α, β, γ)
Lz(f, α, β, γ) = 𝒪(𝐤, :left, f, α, β, γ)

Rx(f, α, β, γ) = 𝒪(𝐢, :right, f, α, β, γ)
Ry(f, α, β, γ) = 𝒪(𝐣, :right, f, α, β, γ)
Rz(f, α, β, γ) = 𝒪(𝐤, :right, f, α, β, γ)

nothing  #hide

# Now we define their commutator ``[O₁, O₂] = O₁O₂ - O₂O₁``:
function commutator(O₁, O₂)
    (
        O₁((α, β, γ)->O₂(f, α, β, γ), α, β, γ)
        - O₂((α, β, γ)->O₁(f, α, β, γ), α, β, γ)
    ).expand().simplify()
end

nothing  #hide

# And finally, evaluate each in turn.  We expect ``[L_x, L_y] = i L_z`` and cyclic
# permutations:

#md # ### ``L`` commutators in Euler angles
commutator(Lx, Ly)
# which equals ``i L_z``,
commutator(Ly, Lz)
# which equals ``i L_x``, and
commutator(Lz, Lx)
# which equals ``i L_y``.  Similarly, we expect ``[R_x, R_y] = i R_z`` and cyclic
# permutations:

#md # ### ``R`` commutators in Euler angles
commutator(Rx, Ry)
# which equals ``i R_z``,
commutator(Ry, Rz)
# which equals ``i R_x``, and
commutator(Rz, Rx)
# which equals ``i R_y`` — all as expected.

# Just for completeness, let's evaluate the commutators of the left and right operators,
# which should all be zero.

#md # ### ``L,R`` commutators in Euler angles
[commutator(L, R) for L ∈ (Lx, Ly, Lz), R ∈ (Rx, Ry, Rz)]
# This completes independent commutator results, which are all as we expect them to be.


# ## Standard expressions on ``𝕊²``
# We can substitute ``(α, β, γ) \to (φ, θ, 0)`` to get the standard expressions for the
# angular momentum operators on the 2-sphere.

#md # ### ``L`` operators in spherical coordinates
@display2 L(𝐢)
#-
@display2 L(𝐣)
#-
@display2 L(𝐤)

# We can also provide the usual expressions for the raising and lowering operators in terms
# of spherical coordinates with ``L_{\pm} = L_x \pm i L_y``:

#md # ### ``L_{\pm}`` operators in spherical coordinates
@display2 L(+)
#-
@display2 L(-)

# These are all indeed the standard expressions for the angular-momentum operators on the
# 2-sphere, as seen in numerous references, so we can declare compatibility between our
# unusual definition of ``L`` and more standard definitions.
#
# Now, note that including ``\partial_γ`` for an expression on the 2-sphere doesn't
# actually make any sense: ``γ`` isn't even a coordinate for the 2-sphere!  However,
# for historical reasons, we include it here when showing the results of the ``R`` operator
# in Euler angles.

#md # ### ``R`` operators in spherical coordinates
@display2 R(𝐢)
#-
@display2 R(𝐣)
#-
@display2 R(𝐤)

# We get nonzero components of ``\partial_γ``, showing that these operators really *do
# not* make sense for the 2-sphere, and therefore that it doesn't actually make sense to
# define spin-weighted spherical functions on the 2-sphere; they really only make sense on
# the 3-sphere.  Nonetheless, if we stipulate that the function ``\eta`` has a specific spin
# weight, that means that *on the 3-sphere* it is an eigenfunction with ``R_z\eta =
# i\partial_γ \eta = s\eta``.  So we could just substitute ``-i s`` for
# ``\partial_γ`` in the expressions above, and recover the standard spin-weight
# operators.  We get
# ```math
# \left[R_x + i R_y\right] \eta
# = \left[
#     -i \frac{1}{\sin θ} \frac{\partial}{\partial ϕ}
#     + \frac{s}{\tan θ}
#     - \frac{\partial}{\partial θ}
#   \right] \eta
# = -(\sin θ)^s \left\{
#     \frac{\partial}{\partial θ}
#     +i \frac{1}{\sin θ} \frac{\partial}{\partial ϕ}
#   \right\}
#   \left\{ (\sin θ)^{-s} \eta \right\}.
# ```
# And in the latter form, we can see that ``R_x + i R_y`` is exactly the spin-raising
# operator ``\eth`` as originally defined by [Newman_1966](@citet) in their Eq. (3.8).  The
# complex-conjugate of this operator is the spin-lowering operator ``\bar{\eth}`` for
# ``R_z``.  *By definition* of raising and lowering operators, this means that ``[R_z, \eth]
# = \eth`` and ``[R_z, \bar{\eth}] = -\bar{\eth}``.  We can verify these results by
# computing the commutators directly from the expressions above:
# ```math
# \begin{aligned}
# [R_z, \eth]
#   &= [R_z, R_x] + i [R_z, R_y] = i R_y - i i R_x = R_x + i R_y = \eth,
# \\
# [R_z, \bar{\eth}]
#   &= [R_z, R_x] - i [R_z, R_y] = i R_y + i i R_x = -R_x + i R_y = -\bar{\eth}.
# \end{aligned}
# ```
#
# So we see that we've reproduced precisely the standard expressions for the spin-weighted
# spherical functions (depicted as functions on the 2-sphere) from the expressions for the
# angular-momentum operators acting on general functions on the 3-sphere.  The standard
# expressions appear arbitrarily, and are not even well defined as functions on the 2-sphere
# because they also need input from tangent space of the sphere — which is not part of the
# 2-sphere proper.  On the other hand, the expressions from the 3-sphere are mathematically
# and physically well defined and intuitive.  Note that the latter is complete in itself; it
# can stand alone without reference to the 2-sphere.  Rather, what we have done here is just
# shown the connection to the inadequate standard presentation.  But it is important to
# recognize that our complete treatment on ``\mathrm{Spin}(3)`` is the more fundamental one,
# and can be used without reference to the older treatment.
