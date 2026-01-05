md"""
# Whittaker (1904)

!!! info "Summary"

    Whittaker consolidated conventions that are still used by physicists today.  That
    includes a right-handed orthogonal coordinate system labeled ``(x, y, z)``; "Eulerian
    angles" that he labels ``(θ, ϕ, ψ)`` but which we would label ``(ϕ, θ, ψ)``, meaning
    that we use the same symbols for the same angles, but in different orders; the resulting
    spherical coordinates as we now use them, with ``θ`` being colatitude, ``ϕ`` being the
    azimuthal angle from the ``x``-axis, and ``ψ=0``; and the standard quaternion basis
    elements ``(i, j, k)``, with rotation of a vector ``v`` by a quaternion ``q`` is given
    by ``q v q⁻¹``.

!!! warning

    Whittaker cited the wrong paper for the Euler angles, referring to Euler's 1776 paper on
    rigid-body dynamics [Euler_1776a](@cite), when in fact Euler introduced these angles
    in his 1767 paper on rigid-body dynamics [Euler_1767](@cite).

Whittaker's "Analytical Dynamics" [Whittaker_1947](@cite) was the most influential book on
classical mechanics and mathematical physics of the first half of the 20th century.  In
particular, quantum physicists found Whittaker's approach to be helpful when inventing their
new branch of physics.  It was originally published in 1904, with the 4th edition coming out
in 1947, which was reissued in 1988, with reprintings as late as 1993 — which perhaps says
something about its influence.  Sommerfeld's "Lectures on Theoretical Physics" came out
starting in 1943, and Goldstein's "Classical Mechanics" in 1950, marking the end of
Whittaker's dominance.

For better or for worse, essentially all of the basic conventions used by physicists today
(including many used by this package) were consolidated here.  Tragically, that includes a
dismissive snarkiness toward quaternions — presumably because Whittaker, like so many
physicists of his time, was poisoned by the harsh invective against quaternions coming
mostly from Gibbs and Heaviside, belying the impoverished understanding at the turn of the
20th century of the close *interdependence* between quaternions and vectors.  Whittaker
often declines to refer to them as "quaternions", referring to them simply as a set of
parameters.  He even goes so far as to make the petty and false claim that quaternion
multiplication had been discovered independently by three other people in addition to
Hamilton.  His muddled and tendentious treatment of quaternions surely sealed their fate
until computer graphics, aerospace, and robotics applications came to rescue them from
obscurity.


## Implementing expressions

We begin by writing code that implements the concepts described by Ref.
[Whittaker_1947](@cite).  We encapsulate the formulas in a module so that we can test them
against the `Quaternionic` and `SphericalFunctions` package.

"""

using TestItems: @testitem  #hide
@testitem "Whittaker conventions" setup=[ConventionsSetup, Utilities] begin  #hide

module Whittaker
#+

# ### Basis vectors and handedness
#
# We start with Sec. I.7, where Whittaker introduces the axis-angle representation:

# > Let rectangular axes ``Oxyz`` be taken, fixed in space: these will be supposed to form a
# > right-handed system, i.e. if the axes are so placed that ``Oz`` is directed vertically
# > upwards and ``Oy`` is directed to the northern horizon, then ``Ox`` will be directed to
# > the east.

# (Note that Whittaker will soon distinguish between the system ``OXYZ`` which will be fixed
# in space, and the system ``Oxyz`` which will be fixed in the rotating body.) This is a
# right-handed orthogonal system, which we will represent by the corresponding unit vectors
# as `QuatVec`s:

import Quaternionic: Quaternionic, Rotor, 𝐢, 𝐣, 𝐤, ⋅, ×̂ 

const Ox = Quaternionic.𝐢
const Oy = Quaternionic.𝐣
const Oz = Quaternionic.𝐤

const east = Ox
const north = Oy
const up = Oz

const south = -north
const west = -east
const down = -up
#+

# Whittaker continues:
#
# > Let the displacement considered be equivalent to a rotation through an angle ``ω`` about
# > a line whose direction-angles are ``(α, β, γ)``[...]
#
# The "direction-angles" are the angles that the line makes with the ``Ox``, ``Oy``, and
# ``Oz`` axes, respectively.  In modern parlance, we would typically represent this "line"
# with a unit vector, which would have components ``(\cos α, \cos β, \cos γ)``.
line(α, β, γ) = cos(α) * Ox + cos(β) * Oy + cos(γ) * Oz
#+

# He then specifies the handedness of the rotation as follows:
#
# > The angle ``ω`` must be taken with its appropriate sign, the sign being positive when
# > the line ``(α, β, γ)`` being directed vertically upwards, the rotation from the southern
# > horizon to the northern is round by the east.
#
# This supposes that ``(α, β, γ) = (π/2, π/2, 0)``, so the "line" is ``𝐳``, and we rotate
# the vector ``-𝐣`` by *increasing* ``ω`` from ``ω=0`` until it reaches ``+𝐣`` at ``ω=π``,
# then we pass through ``+𝐢`` on the way.  This is just a standard right-handed rotation, so
# we would write the rotation about this line as
axis_angle_rotation(ω, α, β, γ) = exp((ω / 2) * line(α, β, γ))
#+

# We will test below that the angles the line makes with the axes are what Whittaker
# intended, and that the behavior in this scenario describing handedness is as expected.


# ### Quaternions
#
# Section I.9 introduces quaternions in a very awkward way, evidently hobbled by Whittaker's
# spite, driving what we might presume to be his intentional obfuscation of the elegance of
# quaternions.  He starts out by arbitrarily introducing four combinations of the axis-angle
# parameters, which so happen to add in quadrature to 1.  He then works through some very
# tedious math to eventually show that these four parameters happen to obey the laws of
# quaternion multiplication.  He denotes the quaternion (after some evident indecision as to
# whether the scalar component should come first or last) as ``χ + ξi + ηj + ζk``, where
# ``i``, ``j``, ``k`` satisfy
# ```math
# i^2 = j^2 = k^2 = -1,
# \quad
# ij = -ji = k,
# \quad
# jk = -kj = i,
# \quad
# ki = -ik = j.
# ```
# These are exactly our conventions.  Whittaker goes on to say
#
# > The reader who is acquainted with quaternions will observe that the effect of the
# > rotation on any vector ``ρ`` is to convert it into the vector ``qρq⁻¹``, where ``q``
# > denotes the quaternion ``χ + ξi + ηj + ζk``; the quaternion itself is *not* the
# > rotational operator.
#
# Note the petty little dig attempting to diminish quaternions in the last clause.  This is
# exactly how we apply quaternions as well; certain other conventions use ``q⁻¹ρq``.


# ### "The Eulerian angles"
#
# Section I.10 is titled "The Eulerian angles".  Whittaker is not just saying that the three
# angles he introduces in this section are *akin* to angles that Euler used; he is
# specifically crediting Euler with the particular construction he uses:
#
# > The most practically useful of the various methods of representing parametrically the
# > displacement of a rigid body due to a rotation round a fixed point is likewise due to
# > Euler†: it has the disadvantage of being unsymmetrical, but is otherwise very simple and
# > convenient.
#
# That dagger cites [Euler_1776a](@citet), which discusses the axis-angle representation,
# the direction-cosine matrix, and a reparameterization of the direction-cosine matrix with
# just three angles — which proves that the rotation group is three-dimensional, but *does
# not* provide a parameterization at all related to the one Whittaker now gives, nor can I
# find Euler ever using a similar one:
#
# > Let ``O`` be the fixed point round which the rotation takes place, and let ``OXYZ`` be a
# > right-handed system of rectangular axes fixed in space. Let ``Oxyz`` be rectangular axes
# > fixed relatively to the body and moving with it, and such that before the displacement
# > the two sets of axes ``OXYZ`` and ``Oxyz`` are coincident in position.  Let ``OK`` be
# > perpendicular to the plane ``zOZ``, drawn so that if ``OZ`` is directed to the vertical
# > and the projection of ``Oz`` perpendicular to ``OZ`` is directed to the south, then
# > ``OK`` is directed to the east.  Denote the angles ``z\hat{O}Z``, ``Y\hat{O}K``,
# > ``y\hat{O}K`` by ``\theta``, ``\phi``, ``\psi``, respectively: these are known as the
# > three *Eulerian angles* defining the position of the axes ``Oxyz`` with reference to the
# > axes ``OXYZ``.
#
# The line ``OK`` is often called the "line of nodes", and is not a very natural object to
# define, except in the case of successive rotations — which, again, is not something that
# Euler did.  Nonetheless, Whittaker does not actually specify a sequence of rotations, so
# we have to devise one that will result in these angles between the various vectors (and we
# will test below that this works).
#
# The use of the angle between the initial and final positions of the ``OZ`` axis suggests
# that the direct rotation from one to the other should be one of our sequence, and the two
# in our sequence should be rotations about the inital ``OZ`` axis and the final ``Oz``
# axis.  Moreover, the use of the angles between ``OK`` and the initial and final positions
# of the ``OY`` axis suggest that ``OK`` is really the intermediate position of the ``OY``
# axis.
#
# So essentially, this rotation is equivalent to an initial rotation about ``OZ`` by ``ϕ``,
# which takes ``OY`` to ``OY'=OK`` (hence ``YÔK=ϕ``); then a rotation about ``OY'=OK`` by
# ``θ`` to take ``OZ`` onto ``OZ'=OZ''=Oz`` (hence ``zÔZ=θ``); finally a rotation about
# ``OZ'=Oz`` by ``ψ``, which takes ``OY'=OK`` onto ``OY''=Oy`` (hence ``yÔK=ψ``).  That is,
# in modern parlance this is a ``z``-``y'``-``z''`` rotation.  We might write this in
# quaternion notation as
# ```math
# \begin{aligned}
# \exp\left[ \frac{ψ}{2} Oz \right]\,
# \exp\left[ \frac{θ}{2} OK \right]\,
# \exp\left[ \frac{ϕ}{2} OZ \right]
# &=
# \exp\left[ \frac{ϕ}{2} OZ \right]\,
# \exp\left[ \frac{θ}{2} OY \right]\,
# \exp\left[ \frac{ψ}{2} OZ \right] \\
# &=
# \exp\left[ \frac{ϕ}{2} 𝐤 \right]\,
# \exp\left[ \frac{θ}{2} 𝐣 \right]\,
# \exp\left[ \frac{ψ}{2} 𝐤 \right].
# \end{aligned}
# ```
function eulerian_rotation(θ, ϕ, ψ)
    OX, OY, OZ = 𝐢, 𝐣, 𝐤
    R₁ = exp((ϕ/2) * OZ)
    OK = R₁(OY)
    R₂ = exp((θ/2) * OK)
    Oz = R₂(R₁(OZ))
    R₃ = exp((ψ/2) * Oz)
    R₁ * R₂ * R₃
end
#+

# We'll also need a function to find what Whittaker means by ``OK``, which is the normalized
# cross product ``OZ × Oz``:
function OK(θ, ϕ, ψ)
    OZ = 𝐤
    Oz = eulerian_rotation(θ, ϕ, ψ)(OZ)
    OZ ×̂ Oz  ## Normalized cross product
end
#+

# We also implement, for testing purposes, functions to evaluate the three angles Whittaker
# refers to:
function YÔK(θ, ϕ, ψ)
    OY = 𝐣
    let OK=OK(θ, ϕ, ψ)
        acos(OY ⋅ OK)
    end
end
function zÔZ(θ, ϕ, ψ)
    OZ = 𝐤
    let OK=OK(θ, ϕ, ψ)
        Oz = eulerian_rotation(θ, ϕ, ψ)(OZ)
        acos(Oz ⋅ OZ)
    end
end
function yÔK(θ, ϕ, ψ)
    OY = 𝐣
    let OK=OK(θ, ϕ, ψ)
        Oy = eulerian_rotation(θ, ϕ, ψ)(OY)
        acos(Oy ⋅ OK)
    end
end
#+

# Finally, Whittaker provides a table of values for the "direction-cosines" — meaning the
# projections of the rotated basis vectors onto the fixed basis vectors:
function direction_cosine(θ, ϕ, ψ)
    [
    ##                       X                                Y                      Z
    #= x =#  cos(ϕ)cos(θ)cos(ψ)-sin(ϕ)sin(ψ)  sin(ϕ)cos(θ)cos(ψ)+cos(ϕ)sin(ψ)  -sin(θ)cos(ψ);
    #= y =# -cos(ϕ)cos(θ)sin(ψ)-sin(ϕ)cos(ψ) -sin(ϕ)cos(θ)sin(ψ)+cos(ϕ)cos(ψ)   sin(θ)sin(ψ);
    #= z =#            cos(ϕ)sin(θ)                     sin(ϕ)sin(θ)                cos(θ)
    ]
end
#+

# ### Connecting Eulerian angles to quaternions
#
# In section I.11, Whittaker derives the relationship between the Eulerian angles and the
# quaternion components:
function quaternion_from_eulerian(θ, ϕ, ψ)
    ξ = sin(θ/2) * sin((ψ - ϕ)/2)
    η = sin(θ/2) * cos((ψ - ϕ)/2)
    ζ = cos(θ/2) * sin((ψ + ϕ)/2)
    χ = cos(θ/2) * cos((ψ + ϕ)/2)
    χ + ξ*𝐢 + η*𝐣 + ζ*𝐤
end
#+

# These are all the conventions we will need from Whittaker.

end  #module Whittaker
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 10eps()
ϵᵣ = 10eps()
#+

# ### Basis vectors and handedness
#
# Test that the angles the line makes with the axes are what Whittaker intended:
for (α,β,γ) ∈ αβγrange()
    l = line(α, β, γ)
    @test acos(l ⋅ Whittaker.Ox) ≈ α atol=ϵₐ rtol=ϵᵣ
    @test acos(l ⋅ Whittaker.Oy) ≈ β atol=ϵₐ rtol=ϵᵣ
    @test acos(l ⋅ Whittaker.Oz) ≈ γ atol=ϵₐ rtol=ϵᵣ
end
#+

# Now we'll test the behavior that
# > the line ``(α, β, γ)`` being directed vertically upwards, the rotation from the southern
# > horizon to the northern is round by the east.
let α=π/2, β=π/2, γ=0
    ω = range(0, π, length=21)
    ## Check that the first element takes `south` to `south`
    let R = Whittaker.axis_angle_rotation(ω[begin], α, β, γ)
        @test R(Whittaker.south) ≈ Whittaker.south atol=ϵₐ rtol=ϵᵣ
    end
    ## Check that the final element takes `south` to `north`
    let R = Whittaker.axis_angle_rotation(ω[end], α, β, γ)
        @test R(Whittaker.south) ≈ Whittaker.north atol=ϵₐ rtol=ϵᵣ
    end
    ## Check that everything in between is "round by the east" by testing that every
    ## intermediate point is closer to east than to west
    for ωᵢ ∈ ω[begin+1:end-1]
        let R = Whittaker.axis_angle_rotation(ωᵢ, α, β, γ)
            p = R(Whittaker.south)
            @test abs2(p - Whittaker.east) < abs2(p - Whittaker.west)
        end
    end
end
#+

# ### Quaternions
#
# We simply test the multiplication rules:
import Quaternionic: 𝐢, 𝐣, 𝐤

@test 𝐢^2 == 𝐣^2 == 𝐤^2 == -1
@test 𝐢 * 𝐣 == -𝐣 * 𝐢 == 𝐤
@test 𝐣 * 𝐤 == -𝐤 * 𝐣 == 𝐢
@test 𝐤 * 𝐢 == -𝐢 * 𝐤 == -𝐣
#+

# ### "The Eulerian angles"
# The discussion above shows that the angles ``ϕ,θ,ψ`` in that order correspond to what we
# would denote as ``α,β,γ`` in that order.  So we rename our utility function for generating
# a variety of angles:
const ϕθψrange = αβγrange
#+

# First, we test that the rotation as we've implemented it does correspond to the Euler
# rotation implemented by `Quaternionic`:
for (ϕ,θ,ψ) ∈ ϕθψrange()
    @test Whittaker.eulerian_rotation(θ, ϕ, ψ) ≈ Quaternionic.from_euler_angles(ϕ, θ, ψ)
end
#+

# Next, we test that it results in the angles between axes that Whittaker described:
for (ϕ,θ,ψ) ∈ ϕθψrange()
    @test YÔK(θ, ϕ, ψ) ≈ ϕ 
    @test zÔZ(θ, ϕ, ψ) ≈ θ
    @test yÔK(θ, ϕ, ψ) ≈ ψ
end
#+

# Finally, we'll test that the rotated axes project onto the fixed axes according to
# Whittaker's table of direction-cosines:
for (ϕ,θ,ψ) ∈ ϕθψrange()
    X = Whittaker.Ox
    Y = Whittaker.Oy
    Z = Whittaker.Oz
    R = Whittaker.eulerian_rotation(θ, ϕ, ψ)
    x = R(X)
    y = R(Y)
    z = R(Z)
    projections = [
        ##         X       Y       Z
        #= x =#  x ⋅ X   x ⋅ Y   x ⋅ Z ;
        #= y =#  y ⋅ X   y ⋅ Y   y ⋅ Z ;
        #= z =#  z ⋅ X   z ⋅ Y   z ⋅ Z
    ]
    dcm = Whittaker.direction_cosine(θ, ϕ, ψ)
    @test projections ≈ dcm atol=ϵₐ rtol=ϵᵣ
end
#+

# ### Connecting Eulerian angles to quaternions
#
# Finally, we can just test that the components Whittaker derived are the components we've
# been using.
for (ϕ,θ,ψ) ∈ ϕθψrange()
    R₁ = Whittaker.eulerian_rotation(θ, ϕ, ψ)
    R₂ = Whittaker.quaternion_from_eulerian(θ, ϕ, ψ)
    @test R₁ ≈ R₂ atol=ϵₐ rtol=ϵᵣ
end

end  #@testitem  #hide
