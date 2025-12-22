md"""
# Euler (1767 and 1776)

!!! warn "Change what's below"

    The following needs to talk mostly about [Euler_1767](@citet), which *did* actually use 
    Euler angles as we would recognize them now.  It can talk about incorrect citations to
    the other two, but that's secondary.


!!! info "Summary"

    Euler never actually used a succession of three simple rotations to represent a general
    rotation, so the term "Euler angles" is a misnomer — besmirching Euler's good name.
    Euler did, however, introduce the axis-angle representation and the direction-cosine
    matrix, which are important concepts in analyzing rotations.

Among Leonhard Euler's vast body of contributions to math and physics, we find his research
on rigid-body dynamics.  He was working at a time before vector (or quaternion) algebra had
been developed, and had to rely on awkward techniques reminiscent of Euclid's constructions.
Nonetheless, he found insights into the mathematics behind rotations that remain relevant
today.  In particular, in his 1776 papers [Euler_1776a](@cite) and [Euler_1776b](@cite), he
developed the axis-angle representation of rotations and the direction-cosine matrix.  But
nowhere in his writings — as far as I can find — did Euler ever use what we today call
"Euler angles".

By 1868, Tait [Tait_1868](@cite) writes down a rotation in exactly what we would now call
ZYZ Euler angles — which he calls "the usual angles ψ, θ, ϕ".  So somehow, this
representation had become normalized by that time.  It's interesting to note that Tait was
actually using the "usual" angles to demonstrate how quaternions were superior.  But at no
point does he refer to them as "Euler" angles.

The very influential book by [Whittaker](@cite Whittaker_1947) contained section I.9 on
"Euler's parametric specification of rotations round a point", which cites Euler's sections
6 and following [Euler_1776b](@cite).  That section, I would say, correctly cites Euler's
work on direction cosines, but doesn't yet get to "Euler angles".  However, his very next
section I.10 on "The Eulerian angles" cites Euler's *preceding* work [Euler_1776a](@cite)
when constructing a sequence of three rotations that we would now recognize as the
``ZY'Z''`` convention.  However, there is absolutely no support for this construction in
Euler's work.  Whittaker has besmirched Euler's good name by associating him with this
construction.

As far as I can find, Davenport is actually the first person to come along and try to
actually show *when exactly* it is possible to decompose an arbitrary rotation into three
rotations about fixed axes [Davenport_1973](@cite).

"""
