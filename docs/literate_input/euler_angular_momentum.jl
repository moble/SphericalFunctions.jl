# # Expressing angular-momentum operators in Euler angles
# Here, we will use SymPy to just grind through the algebra of expressing the
# angular-momentum operators in terms of Euler angles.


# Essential imports
import SymPyPythonCall
import SymPyPythonCall: sympy
import SymPyPythonCall: symbols, sqrt, exp, sin, cos, acos, atan, Matrix, simplify
const Derivative = sympy.Derivative
const Quaternion = sympy.Quaternion
const π = sympy.pi

# Define the spherical coordinates
α, β, γ, θ = symbols("α β γ θ", real=true, positive=true)
uw, ux, uy, uz = symbols("u_w u_x u_y u_z", real=true)

# Define our basis quaternions
i = Quaternion(0, 1, 0, 0)
j = Quaternion(0, 0, 1, 0)
k = Quaternion(0, 0, 0, 1)

# And an arbitrary vector quaternion
u = Quaternion(0, ux, uy, uz)

# Check that multiplication agrees with our conventions
@assert i*j == k
@assert j*k == i
@assert k*i == j
@assert i*j*k == Quaternion(-1, 0, 0, 0)


#
function L(u)
    e = cos(θ/2) + u * sin(θ/2)
    R = ((cos(α/2) + k * sin(α/2)) * (cos(β/2) + j * sin(β/2)) * (cos(γ/2) + k * sin(γ/2))).expand().simplify()
    eR = (e * R).expand().simplify()
    w, x, y, z = eR.to_Matrix().transpose().tolist()[1]
    αp = (atan(z/w) + atan(-x/y)).expand().simplify()
    βp = (2*acos(sqrt(w^2 + z^2) / sqrt(w^2 + x^2 + y^2 + z^2))).expand().simplify()
    γp = (atan(z/w) - atan(-x/y)).expand().simplify()
    return (
        Derivative(αp, θ).doit().subs(θ, 0).expand().simplify(),
        Derivative(βp, θ).doit().subs(θ, 0).expand().simplify(),
        Derivative(γp, θ).doit().subs(θ, 0).expand().simplify()
    )
end

function R(u)
    e = cos(θ/2) + u * sin(θ/2)
    R1 = ((cos(α/2) + k * sin(α/2)) * (cos(β/2) + j * sin(β/2)) * (cos(γ/2) + k * sin(γ/2))).expand().simplify()
    Re = (R1 * e).expand().simplify()
    w, x, y, z = Re.to_Matrix().transpose().tolist()[1]
    αp = (atan(z/w) + atan(-x/y)).expand().simplify()
    βp = (2*acos(sqrt(w^2 + z^2) / sqrt(w^2 + x^2 + y^2 + z^2))).expand().simplify()
    γp = (atan(z/w) - atan(-x/y)).expand().simplify()
    return (
        Derivative(αp, θ).doit().subs(θ, 0).expand().simplify(),
        Derivative(βp, θ).doit().subs(θ, 0).expand().simplify(),
        Derivative(γp, θ).doit().subs(θ, 0).expand().simplify()
    )
end


#
L(i)

#
L(j)

#
L(k)

#
R(i)

#
R(j)

#
R(k)

#
