md"""
# ``L_j`` and ``R_j`` with Euler angles

This package defines the angular-momentum operators ``L_j`` and ``R_j`` in terms of elements
of the Lie group and algebra:
```math
L_𝐮 f(𝐑) = \left. i \frac{d}{d\epsilon}\right|_{\epsilon=0}
f\left(e^{-\epsilon 𝐮/2}\, 𝐑\right)
\qquad \text{and} \qquad
R_𝐮 f(𝐑) = -\left. i \frac{d}{d\epsilon}\right|_{\epsilon=0}
f\left(𝐑\, e^{-\epsilon 𝐮/2}\right),
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
  5. Using the natural extension of Goldberg et al.'s SWSHs to include ``\gamma``, we can
     see that the natural spin-weight operator is ``R_z = i \partial_\gamma``.  Thus, we
     define ``R_z = s`` for a function with spin weight ``s``.
  6. The spin-raising operator for ``R_z`` is ``\eth = R_x + i R_y``; the spin-lowering
     operator is ``\bar{\eth} = R_x - i R_y``.

## Analytical groundwork

We start by defining a new set of Euler angles according to
```math
\mathbf{R}_{\alpha', \beta', \gamma'}
= e^{-\theta \mathbf{u} / 2} \mathbf{R}_{\alpha, \beta, \gamma}
\qquad \text{or} \qquad
\mathbf{R}_{\alpha', \beta', \gamma'}
= \mathbf{R}_{\alpha, \beta, \gamma} e^{-\theta \mathbf{u} / 2}
```
where ``\mathbf{u}`` will be each of the basis quaternions, and each of ``\alpha'``,
``\beta'``, and ``\gamma'`` is a function of ``\alpha``, ``\beta``, ``\gamma``, and
``\theta``.  Then, we note that the chain rule tells us that
```math
\frac{\partial}{\partial \theta}
=
\frac{\partial \alpha'}{\partial \theta} \frac{\partial}{\partial \alpha'}
+ \frac{\partial \beta'}{\partial \theta} \frac{\partial}{\partial \beta'}
+ \frac{\partial \gamma'}{\partial \theta} \frac{\partial}{\partial \gamma'},
```
which we will use to convert the general expression for the angular-momentum operators in
terms of ``\partial_\theta`` into an expression in terms of derivatives with respect to
these new Euler angles:
```math
\begin{align}
  L_j f(\mathbf{R}_{\alpha, \beta, \gamma})
  &=
  \left. i \frac{\partial} {\partial \theta} f \left( e^{-\theta \mathbf{e}_j / 2}
  \mathbf{R}_{\alpha, \beta, \gamma} \right) \right|_{\theta=0}
  \\
  &=
  i \left[ \left(
      \frac{\partial \alpha'}{\partial \theta} \frac{\partial}{\partial \alpha'}
      + \frac{\partial \beta'}{\partial \theta} \frac{\partial}{\partial \beta'}
      + \frac{\partial \gamma'}{\partial \theta} \frac{\partial}{\partial \gamma'}
  \right) f \left(\alpha', \beta', \gamma'\right) \right]_{\theta=0}
  \\
  &=
  i \left[ \left(
      \frac{\partial \alpha'}{\partial \theta} \frac{\partial}{\partial \alpha}
      + \frac{\partial \beta'}{\partial \theta} \frac{\partial}{\partial \beta}
      + \frac{\partial \gamma'}{\partial \theta} \frac{\partial}{\partial \gamma}
  \right) f \left(\alpha, \beta, \gamma\right) \right]_{\theta=0},
\end{align}
```
and similarly for ``R_j``.

So the objective is to find the new Euler angles, differentiate with respect to ``\theta``,
and then evaluate at ``\theta = 0``.  We do this by first multiplying ``\mathbf{R}_{\alpha,
\beta, \gamma}`` and ``e^{-\theta \mathbf{u} / 2}`` in the desired order, then expanding the
results in terms of its quaternion components, and then computing the new Euler angles in
terms of those components according to the usual expression.

"""

#src # Do this first just to hide stdout of the conda installation step
# ````@setup euler_angular_momentum
# import SymPyPythonCall
# ````

# ## Computational infrastructure
# We'll use SymPy (via Julia) since `Symbolics.jl` isn't very good at trig yet.
import LaTeXStrings: @L_str, LaTeXString
import Quaternionic: Quaternionic, Quaternion, components
import SymPyPythonCall
import SymPyPythonCall: sympy, symbols, sqrt, sin, cos, tan, acos, atan, latex
const expand_trig = sympy.expand_trig
const Derivative = sympy.Derivative
const π = sympy.pi
const I = sympy.I
nothing  #hide

# Define coordinates we will use
α, β, γ, θ, ϕ = symbols("α β γ θ ϕ", real=true, positive=true)
nothing  #hide

# Reinterpret the quaternion basis elements for compatibility with SymPy.  (`Quaternionic`
# defines the basis with `Bool` components, but SymPy can't handle that.)
const 𝐢 = Quaternion{Int}(Quaternionic.𝐢)
const 𝐣 = Quaternion{Int}(Quaternionic.𝐣)
const 𝐤 = Quaternion{Int}(Quaternionic.𝐤)
nothing  #hide

# Next, we define functions to compute the Euler components of the left and right operators
function 𝒪(u, side)
    ## Substitutions that sympy doesn't make but we want
    subs = Dict(
        cos(β)/sin(β) => 1/tan(β),
        sqrt(1 - cos(β))*sqrt(cos(β) + 1) => sin(β)
    )

    ## Define the essential quaternions
    e = cos(θ/2) + u * sin(-θ/2)
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

    ## Differentiate with respect to θ, set θ to 0, and simplify
    ∂α′∂θ = expand_trig(Derivative(α′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
    ∂β′∂θ = expand_trig(Derivative(β′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
    ∂γ′∂θ = expand_trig(Derivative(γ′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))

    return ∂α′∂θ, ∂β′∂θ, ∂γ′∂θ
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
            ∂α′∂θ, ∂β′∂θ, ∂γ′∂θ = latex.($expr)  # Call expr; format results as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = i\left[
                %$(∂α′∂θ) \frac{\partial}{\partial \alpha}
                + %$(∂β′∂θ) \frac{\partial}{\partial \beta}
                + %$(∂γ′∂θ) \frac{\partial}{\partial \gamma}
            \right]"""  # Display the result in LaTeX form
        end
    else
        quote
            ∂α′∂θ, ∂β′∂θ, ∂γ′∂θ = latex.($expr)  # Call expr; format results as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = -i\left[
                %$(∂α′∂θ) \frac{\partial}{\partial \alpha}
                + %$(∂β′∂θ) \frac{\partial}{\partial \beta}
                + %$(∂γ′∂θ) \frac{\partial}{\partial \gamma}
            \right]"""  # Display the result in LaTeX form
        end
    end
end
nothing  #hide

# And we'll need another for the angular-momentum operators in standard ``S^2`` form.
conversion(∂) = latex(∂.subs(Dict(α => ϕ, β => θ, γ => 0)).simplify())
macro display2(expr)
    op = string(expr.args[1])
    arg = Dict(:𝐢 => "x", :𝐣 => "y", :𝐤 => "z")[expr.args[2]]
    if op == "L"
        quote
            ∂φ′∂θ, ∂ϑ′∂θ, ∂γ′∂θ = $conversion.($expr)  # Call expr; format results as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = i\left[
                %$(∂ϑ′∂θ) \frac{\partial}{\partial \theta}
                + %$(∂φ′∂θ) \frac{\partial}{\partial \phi}
            \right]"""  # Display the result in LaTeX form
        end
    else
        quote
            ∂φ′∂θ, ∂ϑ′∂θ, ∂γ′∂θ = $conversion.($expr)  # Call expr; format results as LaTeX
            expr = $op * "_" * $arg  # Standard form of the operator
            L"""%$expr = -i\left[
                %$(∂ϑ′∂θ) \frac{\partial}{\partial \theta}
                + %$(∂φ′∂θ) \frac{\partial}{\partial \phi}
                + %$(∂γ′∂θ) \frac{\partial}{\partial \gamma}
            \right]"""  # Display the result in LaTeX form
        end
    end
end
nothing  #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Full expressions on ``S^3``
# Finally, we can actually compute the Euler components of the angular momentum operators.
@display L(𝐢)
#-
@display L(𝐣)
#-
@display L(𝐤)
#-
@display R(𝐢)
#-
@display R(𝐣)
#-
@display R(𝐤)

# In their description of the Wigner 𝔇 functions as wave functions of a rigid symmetric
# top, [Varshalovich_1988](@citet) provide equivalent expressions in Eqs. (6) and (7) of
# their Sec. 4.2 — except that ``R_x`` and ``R_z`` have the wrong signs.
# [Wikipedia](https://en.wikipedia.org/wiki/Wigner_D-matrix#Properties_of_the_Wigner_D-matrix),
# meanwhile, provides equivalent expressions, except that their ``\hat{\mathcal{P}}`` has
# (consistently) the opposite sign to ``R`` defined here.

# Note that the Wikipedia convention is actually entirely sensible — maybe more sensible
# than the one we use.  In that convention ``\hat{J}`` is in the inertial frame, whereas
# ``\hat{P}`` is exactly that operator in the body-fixed frame.  In our notation, we have
# ```math
# \begin{align}
# R_𝐮 f(𝐑)
# &=
# -\left. i \frac{d}{d\epsilon}\right|_{\epsilon=0} f\left(𝐑\, e^{-\epsilon 𝐮/2}\right) \\
# &=
# -\left. i \frac{d}{d\epsilon}\right|_{\epsilon=0}
#   f\left(𝐑\, e^{-\epsilon 𝐮/2}\, 𝐑^{-1}\, 𝐑\right) \\
# &=
# -\left. i \frac{d}{d\epsilon}\right|_{\epsilon=0}
#   f\left(e^{-\epsilon 𝐑\, 𝐮\, 𝐑^{-1}/2}\, 𝐑\right) \\
# &=
# -L_{𝐑\, 𝐮\, 𝐑^{-1}} f(𝐑),
# \end{align}
# ```
# which says that ``R`` is just the *negative of* the ``L`` operator transformed to the
# body-fixed frame.  That negative sign is slightly unnatural, but the reason we choose to
# define ``R`` in this way is for its more natural connection to the literature on
# spin-weighted spherical functions.

# ### Commutators

# We can also compute the commutators of the angular momentum operators, as derived above.
# First, we define the operators acting on functions of the Euler angles.
f = symbols("f", cls=SymPyPythonCall.sympy.o.Function)
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
commutator(Lx, Ly)
#-
commutator(Ly, Lz)
#-
commutator(Lz, Lx)
# Similarly, we expect ``[R_x, R_y] = i R_z`` and cyclic permutations:
commutator(Rx, Ry)
#-
commutator(Ry, Rz)
#-
commutator(Rz, Rx)

# Just for completeness, let's evaluate the commutators of the left and right operators,
# which should all be zero.
commutator(Lx, Rx)
#-
commutator(Lx, Ry)
#-
commutator(Lx, Rz)
#-
commutator(Ly, Rx)
#-
commutator(Ly, Ry)
#-
commutator(Ly, Rz)
#-
commutator(Lz, Rx)
#-
commutator(Lz, Ry)
#-
commutator(Lz, Rz)


# ## Standard expressions on ``S^2``
# We can substitute ``(α, β, γ) \to (φ, θ, 0)`` to get the standard expressions for the
# angular momentum operators on the 2-sphere.
@display2 L(𝐢)
#-
@display2 L(𝐣)
#-
@display2 L(𝐤)

# Those are indeed the standard expressions for the angular-momentum operators on the
# 2-sphere, so we can declare success!
#
# Now, note that including ``\partial_\gamma`` for an expression on the 2-sphere doesn't
# actually make any sense: ``\gamma`` isn't even a coordinate for the 2-sphere!  However,
# for historical reasons, we include it here when showing the results of the ``R`` operator
# in Euler angles.

@display2 R(𝐢)
#-
@display2 R(𝐣)
#-
@display2 R(𝐤)

# We get nonzero components of ``\partial_\gamma``, showing that these operators really *do
# not* make sense for the 2-sphere, and therefore that it doesn't actually make sense to
# define spin-weighted spherical functions on the 2-sphere; they really only make sense on
# the 3-sphere.  Nonetheless, if we stipulate that the function ``\eta`` has a specific spin
# weight, that means that *on the 3-sphere* it is an eigenfunction with ``R_z\eta =
# i\partial_\gamma \eta = s\eta``.  So we could just substitute ``-i s`` for
# ``\partial_\gamma`` in the expressions above, and recover the standard spin-weight
# operators.  We get
# ```math
# \left[R_x + i R_y\right] \eta
# = \left[
#     -i \frac{1}{\sin\theta} \frac{\partial}{\partial \phi}
#     + \frac{s}{\tan \theta}
#     - \frac{\partial}{\partial \theta}
#   \right] \eta
# = -(\sin \theta)^s \left\{
#     \frac{\partial}{\partial \theta}
#     +i \frac{1}{\sin\theta} \frac{\partial}{\partial \phi}
#   \right\}
#   \left\{ (\sin \theta)^{-s} \eta \right\}.
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
