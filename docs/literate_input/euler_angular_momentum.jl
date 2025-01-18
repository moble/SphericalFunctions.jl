md"""
# ``L_j`` in Euler angles
Here, we will use SymPy to just grind through the algebra of expressing the
angular-momentum operators in terms of Euler angles.

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

So the objective is to find the new Euler angles, differentiate with respect to
``\theta``, and then evaluate at ``\theta = 0``.  We do this by first expressing
``\mathbf{R}_{\alpha, \beta, \gamma}`` in terms of its quaternion components, then
multiplying by ``e^{-\theta \mathbf{u} / 2}`` and expanding.  We then find the new Euler
angles in terms of the components of the resulting quaternion according to the usual
expression.

"""

#src # Do this first just to hide stdout of the conda installation step
# ````@setup euler_angular_momentum
# import SymPyPythonCall
# ````


# We'll use SymPy (via Julia) since `Symbolics.jl` isn't very good at trig yet.
using LaTeXStrings
import SymPyPythonCall
import SymPyPythonCall: sympy, symbols, sqrt, sin, cos, tan, acos, atan, latex
const expand_trig = sympy.expand_trig
const Derivative = sympy.Derivative
const Quaternion =  sympy.Quaternion
const π = sympy.pi
nothing  #hide

# Define coordinates we will use
α, β, γ, θ = symbols("α β γ θ", real=true, positive=true)
uw, ux, uy, uz = symbols("u_w u_x u_y u_z", real=true)
nothing  #hide

# Define the basis quaternions
i = Quaternion(0, 1, 0, 0)
j = Quaternion(0, 0, 1, 0)
k = Quaternion(0, 0, 0, 1)
nothing  #hide

# Check that multiplication agrees with our conventions
@assert i*j == k
@assert j*k == i
@assert k*i == j
@assert i*j*k == Quaternion(-1, 0, 0, 0)
nothing  #hide

# Next, we define functions to compute the Euler components of the left and right operators
function 𝒪(u, side)
    subs = Dict(  # Substitutions that sympy doesn't make but we want
        cos(β)/sin(β) => 1/tan(β),
        sqrt(1 - cos(β))*sqrt(cos(β) + 1) => sin(β)
    )
    e = cos(θ/2) + u * sin(-θ/2)
    R₀ = ((cos(α/2) + k * sin(α/2)) * (cos(β/2) + j * sin(β/2)) * (cos(γ/2) + k * sin(γ/2))).expand().simplify()
    w, x, y, z = (
        side == :left ? e * R₀ : R₀ * e
    ).expand().simplify().to_Matrix().transpose().tolist()[1]
    α′ = (atan(z/w) + atan(-x/y)).expand().simplify()
    β′ = (2*acos(sqrt(w^2 + z^2) / sqrt(w^2 + x^2 + y^2 + z^2))).expand().simplify()
    γ′ = (atan(z/w) - atan(-x/y)).expand().simplify()
    ∂α′∂θ = expand_trig(Derivative(α′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
    ∂β′∂θ = expand_trig(Derivative(β′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
    ∂γ′∂θ = expand_trig(Derivative(γ′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
    return ∂α′∂θ, ∂β′∂θ, ∂γ′∂θ
end

L(u) = 𝒪(u, :left)
R(u) = 𝒪(u, :right)

# function L(u)
#     e = cos(θ/2) + u * sin(-θ/2)
#     R₀ = ((cos(α/2) + k * sin(α/2)) * (cos(β/2) + j * sin(β/2)) * (cos(γ/2) + k * sin(γ/2))).expand().simplify()
#     eR = (e * R₀).expand().simplify()
#     w, x, y, z = eR.to_Matrix().transpose().tolist()[1]
#     α′ = (atan(z/w) + atan(-x/y)).expand().simplify()
#     β′ = (2*acos(sqrt(w^2 + z^2) / sqrt(w^2 + x^2 + y^2 + z^2))).expand().simplify()
#     γ′ = (atan(z/w) - atan(-x/y)).expand().simplify()
#     ∂α′∂θ = expand_trig(Derivative(α′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
#     ∂β′∂θ = expand_trig(Derivative(β′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
#     ∂γ′∂θ = expand_trig(Derivative(γ′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
#     return ∂α′∂θ, ∂β′∂θ, ∂γ′∂θ
# end

# function R(u)
#     e = cos(θ/2) + u * sin(-θ/2)
#     R1 = ((cos(α/2) + k * sin(α/2)) * (cos(β/2) + j * sin(β/2)) * (cos(γ/2) + k * sin(γ/2))).expand().simplify()
#     Re = (R1 * e).expand().simplify()
#     w, x, y, z = Re.to_Matrix().transpose().tolist()[1]
#     α′ = (atan(z/w) + atan(-x/y)).expand().simplify()
#     β′ = (2*acos(sqrt(w^2 + z^2) / sqrt(w^2 + x^2 + y^2 + z^2))).expand().simplify()
#     γ′ = (atan(z/w) - atan(-x/y)).expand().simplify()
#     ∂α′∂θ = expand_trig(Derivative(α′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
#     ∂β′∂θ = expand_trig(Derivative(β′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
#     ∂γ′∂θ = expand_trig(Derivative(γ′, θ).doit().subs(θ, 0).expand().simplify().subs(subs))
#     return ∂α′∂θ, ∂β′∂θ, ∂γ′∂θ
# end

nothing  #hide


# Now we can compute the Euler components of the angular momentum operators for the three
# generators and both choices of left and right operators.
L_x = L(i)
L"""
L_x = i\left[
    %$(latex(L_x[1])) \frac{\partial}{\partial \alpha}
    + %$(latex(L_x[2])) \frac{\partial}{\partial \beta}
    + %$(latex(L_x[3])) \frac{\partial}{\partial \gamma}
\right]
"""

#
L_y = L(j)
L"""
L_y = i\left[
    %$(latex(L_y[1])) \frac{\partial}{\partial \alpha}
    + %$(latex(L_y[2])) \frac{\partial}{\partial \beta}
    + %$(latex(L_y[3])) \frac{\partial}{\partial \gamma}
\right]
"""

#
L_z = L(k)
L"""
L_z = i\left[
    %$(latex(L_z[1])) \frac{\partial}{\partial \alpha}
    + %$(latex(L_z[2])) \frac{\partial}{\partial \beta}
    + %$(latex(L_z[3])) \frac{\partial}{\partial \gamma}
\right]
"""

#
R_x = R(i)
L"""
R_x = i\left[
    %$(latex(R_x[1])) \frac{\partial}{\partial \alpha}
    + %$(latex(R_x[2])) \frac{\partial}{\partial \beta}
    + %$(latex(R_x[3])) \frac{\partial}{\partial \gamma}
\right]
"""

#
R_y = R(j)
L"""
R_y = i\left[
    %$(latex(R_y[1])) \frac{\partial}{\partial \alpha}
    + %$(latex(R_y[2])) \frac{\partial}{\partial \beta}
    + %$(latex(R_y[3])) \frac{\partial}{\partial \gamma}
\right]
"""

#
R_z = R(k)
L"""
R_z = i\left[
    %$(latex(R_z[1])) \frac{\partial}{\partial \alpha}
    + %$(latex(R_z[2])) \frac{\partial}{\partial \beta}
    + %$(latex(R_z[3])) \frac{\partial}{\partial \gamma}
\right]
"""
