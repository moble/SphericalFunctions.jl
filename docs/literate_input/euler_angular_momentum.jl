md"""
# ``L_j`` and ``R_j`` in Euler angles
## Analytical groundwork
Here, we will use SymPy to just grind through the algebra of expressing the angular-momentum
operators in terms of Euler angles.

The plan starts by defining a new set of Euler angles according to
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
these new Euler angles.

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

# ## Computing infrastructure
# We'll use SymPy (via Julia) since `Symbolics.jl` isn't very good at trig yet.
import LaTeXStrings: @L_str, LaTeXString
import Quaternionic: Quaternionic, Quaternion, components
import SymPyPythonCall
import SymPyPythonCall: sympy, symbols, sqrt, sin, cos, tan, acos, atan, latex
const expand_trig = sympy.expand_trig
const Derivative = sympy.Derivative
const π = sympy.pi
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

L(u) = 𝒪(u, :left)
R(u) = 𝒪(u, :right)
nothing  #hide


# We need a quick helper macro to format the results.
macro display(expr)
    op = string(expr.args[1])
    arg = Dict(:𝐢 => "x", :𝐣 => "y", :𝐤 => "z")[expr.args[2]]
    quote
        ∂α′∂θ, ∂β′∂θ, ∂γ′∂θ = latex.($expr)  # Call expr; format results as LaTeX
        expr = $op * "_" * $arg  # Standard form of the operator
        L"""%$expr = i\left[
            %$(∂α′∂θ) \frac{\partial}{\partial \alpha}
            + %$(∂β′∂θ) \frac{\partial}{\partial \beta}
            + %$(∂γ′∂θ) \frac{\partial}{\partial \gamma}
        \right]"""  # Display the result in LaTeX form
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
            L"""%$expr = i\left[
                %$(∂ϑ′∂θ) \frac{\partial}{\partial \theta}
                + %$(∂φ′∂θ) \frac{\partial}{\partial \phi}
                + %$(∂γ′∂θ) \frac{\partial}{\partial \gamma}
            \right]"""  # Display the result in LaTeX form
        end
    end
end
nothing  #hide


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
# their Sec. 4.2.


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
# actually make any sense.  However, for historical reasons, we include it here when showing
# the results of the ``R`` operator in Euler angles.  These operators are really only
# relevant for spin-weighted spherical harmonics.  This nonsensicality signals the fact that
# it doesn't actually make sense to define spin-weighted spherical functions on the
# 2-sphere; they really only make sense on the 3-sphere.  Nonetheless, if we stipulate that
# the function in question has a specific spin weight, that means that it is an
# eigenfunction of ``-i\partial_\gamma`` on the 3-sphere, so we could just substitute the
# eigenvalue ``s`` for that derivative in the expression below, and recover the standard
# spin-weight operators.

@display2 R(𝐢)
#-
@display2 R(𝐣)
#-
@display2 R(𝐤)

# This last operator shows us just how little sense it makes to try to define spin-weighted
# spherical functions on the 2-sphere.  The spin eigenvalue ``s`` has to come out of
# nowhere, like some sort of deus ex machina.  Nonetheless, we can see that if we substitute
# the eigenvalue, we get
# ```math
# R_x \eta
# = i\left[ \frac{\partial}{\partial \phi} - \frac{s}{\tan \theta} \right] \eta
# = -(\sin \theta)^s \left\{\frac{\partial}{\partial \theta} \right\}
#   \left\{ (\sin \theta)^{-s} \eta \right\}.
# ```
# And in the latter form, we can see that ``R_x - i R_y`` is exactly the spin-raising
# operator ``\eth`` as originally defined by [Newman_1966](@citet) in their Eq. (3.8).
