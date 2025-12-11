md"""
# Metrics and integration

As something of an exercise, we just work through the basic definitions of metrics and
volume forms in spherical coordinates using SymPy, to verify the volume-form factor that we
use when integrating on the sphere.  We then extend this to three-dimensional spherical
coordinates.


"""

#src # Do this first just to hide stdout of the conda installation step.
#src # Note that we can't just use `#hide` because that still shows stdout.
# ````@setup metrics_and_integration
# import PythonCall
# import SymPyPythonCall
# ````

# We'll use SymPy (via Julia) since `Symbolics.jl` isn't very good at trig yet.
import PythonCall
import SymPyPythonCall: sympy

const π = sympy.pi

## And here are some handy definitions to make manipulations look more natural:
Base.exp(x::PythonCall.Py) = x.exp()
Base.:*(x::SymPyPythonCall.Sym{PythonCall.Py}, q::PythonCall.Py) = x.o * q
(⋅)(a::PythonCall.Py, b::Int) = sympy.S(a.scalar() * b).simplify()
(⋅)(a::PythonCall.Py, b::PythonCall.Py) = sympy.S((a | b.rev()).scalar()).simplify()
nothing  #hide

# ## Three-dimensional space
# Define symbols we will use throughout, including Cartesian coordinates in terms of
# spherical
r, θ, ϕ = sympy.symbols("r θ ϕ", real=true, positive=true, zero=false)
x = r * sympy.sin(θ) * sympy.cos(ϕ)
y = r * sympy.sin(θ) * sympy.sin(ϕ)
z = r * sympy.cos(θ)
nothing  #hide

# The Cartesian metric is just the 3x3 identity matrix
metric_cartesian_3 = sympy.eye(3)

# so the Cartesian basis ``(𝐱, 𝐲, 𝐳)`` is also an *orthonormal* basis.  A point with
# coordinates ``(r, θ, ϕ)`` is separated from the origin by the vector
# ```math
# 𝐩 = r \sin θ \cos ϕ \, 𝐱 + r \sin θ \sin ϕ \, 𝐲 + r \cos θ \, 𝐳,
# ```
# If we treat a SymPy matrix as containing the components of a vector in the Cartesian
# basis, then we write this as
𝐩 = sympy.Matrix([x, y, z])
# and the coordinate basis vectors for spherical coordinates are given by differentiating
# this with respect to each of the spherical coordinates:
𝐫 = 𝐩.diff(r)
𝛉 = 𝐩.diff(θ)
𝛟 = 𝐩.diff(ϕ)
nothing  #hide

# We can compute the Jacobian of the coordinate transformation from Cartesian to spherical
# coordinates, and then transform to obtain the spherical metric — the metric *with respect
# to the spherical coordinate basis*:
jacobian_3 = sympy.Matrix([x, y, z]).jacobian([r, θ, ϕ])
metric_spherical = sympy.simplify((jacobian_3.T * metric_cartesian_3 * jacobian_3))

# Recall that there is also an associated *orthonormal* basis, which we can obtain by
# normalizing each of the coordinate basis vectors, which just involves dividing ``𝐫``,
# ``\mathbb{θ}``, and ``\mathbb{ϕ}`` by ``1``, ``r``, and ``r \sin θ`` respectively.  But
# that is not the basis we are implicitly using when — say — integrating with spherical
# coordinates, so we must use the above metric.
#
# The volume-form factor is given by the square root of the determinant of the metric:
three_volume_form_factor = sympy.simplify(sympy.sqrt(metric_spherical.det()))

# Note that SymPy correctly includes the absolute value of the ``\sin θ`` factor, but we
# can drop the absolute value since ``θ`` is in ``(0, π)``:
three_volume_form_factor = three_volume_form_factor.subs(abs(sympy.sin(θ)), sympy.sin(θ))

# Restricting to the unit sphere (``r=1``), we integrate naively to find the surface area:
S2_surface_area = sympy.integrate(
    sympy.integrate(
        three_volume_form_factor.subs(r, 1),
        (ϕ, 0, 2π),
    ),
    (θ, 0, π),
)

# Therefore, the normalized volume-form factor on the unit sphere 𝕊² is
S2_normalized_volume_form_factor = sympy.simplify(
    three_volume_form_factor.subs(r, 1) / S2_surface_area
)

# ## Four-dimensional space
#
# ### Geometric algebra
#
# Below, we will relate the extended Euler coordinates to Cartesian coordinates in four
# dimensions via the quaternions.  First, we set up the geometric algebra of
# three-dimensional space to check our conventions for the quaternion basis elements against
# the `galgebra` Python package.  We first import the package and set up a 3D Euclidean
# space with orthonormal basis vectors ``(𝐱, 𝐲, 𝐳)``:

## This is equivalent to `from galgebra.ga import Ga` in Python:
const Ga = PythonCall.pyimport("galgebra.ga" => "Ga")

o3d, 𝐱, 𝐲, 𝐳 = Ga.build("𝐱 𝐲 𝐳", g=[1, 1, 1])
nothing  #hide

# Now we define the quaternion basis elements and the pseudoscalar:
𝐢 = 𝐳 * 𝐲
𝐣 = 𝐱 * 𝐳
𝐤 = 𝐲 * 𝐱
𝐈 = 𝐱 * 𝐲 * 𝐳
nothing  #hide

# Now we check some basic identities related to these definitions:
PythonCall.pyall([
    𝐢 == 𝐳*𝐲,
    𝐣 == 𝐱*𝐳,
    𝐤 == 𝐲*𝐱,
    𝐳*𝐲 == -𝐲*𝐳,
    𝐱*𝐳 == -𝐳*𝐱,
    𝐲*𝐱 == -𝐱*𝐲,
    𝐈 == 𝐱*𝐲*𝐳,
    𝐈 == 𝐲*𝐳*𝐱,
    𝐈 == 𝐳*𝐱*𝐲,
    𝐈 == -𝐱*𝐳*𝐲,
    𝐈 == -𝐲*𝐱*𝐳,
    𝐈 == -𝐳*𝐲*𝐱,
])

# Next the basic quaternion relations:
PythonCall.pyall([
    𝐢*𝐣 == 𝐤,
    𝐣*𝐤 == 𝐢,
    𝐤*𝐢 == 𝐣,
    𝐢*𝐣*𝐤 == -1,
])

# And the duality relations:
PythonCall.pyall([
    𝐢 == 𝐈.inv() * 𝐱,
    𝐣 == 𝐈.inv() * 𝐲,
    𝐤 == 𝐈.inv() * 𝐳,
])


# ### Extended-Euler coordinates
# We first define the four-dimensional Cartesian coordinates in terms of Euler angles and a
# radius ``R``:
R, α, β, γ = sympy.symbols("R α β γ", real=true, positive=true, zero=false)
nothing  #hide

# The corresponding quaternion can be written as
𝐐 = R * exp(α * 𝐤 / 2) * exp(β * 𝐣 / 2) * exp(γ * 𝐤 / 2)
nothing  #hide

# Now we define the Cartesian coordinates in four dimensions via the quaternion
# components:
W = 𝐐 ⋅ 1
#-
X = 𝐐 ⋅ 𝐢
#-
Y = 𝐐 ⋅ 𝐣
#-
Z = 𝐐 ⋅ 𝐤

# The Cartesian metric in four dimensions is just the 4x4 identity matrix
metric_cartesian_4 = sympy.eye(4)

# and again we can compute the Jacobian of the coordinate transformation from Cartesian to
# extended-Euler coordinates, and transform to obtain the metric with respect to the
# extended-Euler coordinate basis:
jacobian_4 = sympy.Matrix([W, X, Y, Z]).jacobian([R, α, β, γ])
metric_extended_euler = sympy.simplify((jacobian_4.T * metric_cartesian_4 * jacobian_4))

# And again, the volume-form factor is given by the square root of the determinant of the
# metric:
four_volume_form_factor = sympy.simplify(sympy.sqrt(metric_extended_euler.det()))

# Again, SymPy correctly includes the absolute value of the ``\sin β`` factor, but in this
# case, it can actually be negative, since ``β`` is in ``[0, 2π]``, so it is correct for
# integration over ``\mathrm{Spin}(3)`` to include the absolute value here.
#
# Restricting to the unit sphere (``R=1``), we integrate naively to find the surface area:
S3_surface_area = sympy.integrate(
    sympy.integrate(
        sympy.integrate(
            four_volume_form_factor.subs(R, 1),
            (γ, 0, 2π),
        ),
        (β, 0, 2π),
    ),
    (α, 0, 2π)
)

# Therefore, the normalized volume-form factor on the unit sphere 𝕊³ is
S3_normalized_volume_form_factor = sympy.simplify(
    four_volume_form_factor.subs(R, 1) / S3_surface_area
)

# And finally, we can restrict back to ``\mathrm{SO}(3)`` by taking ``β ∈ [0, π]``, and
# integrate over that range:
SO3_volume = sympy.integrate(
    sympy.integrate(
        sympy.integrate(
            four_volume_form_factor.subs(R, 1),
            (γ, 0, 2π),
        ),
        (β, 0, π),
    ),
    (α, 0, 2π)
)

# And the normalized volume-form factor on ``\mathrm{SO}(3)`` is
SO3_normalized_volume_form_factor = sympy.simplify(
    four_volume_form_factor.subs(abs(sympy.sin(β)), sympy.sin(β)).subs(R, 1) / SO3_volume
)
