md"""
# Euler (1776)

!!! info "Summary"

    Euler never actually used a succession of three simple rotations to represent a general
    rotation, so the term "Euler angles" is a misnomer — besmirching Euler's good name.
    Euler did, however, introduce the axis-angle representation and the direction-cosine
    matrix, which are important concepts in analyzing rotations.


Euler wrote a lot about rigid-body dynamics, 

The poor man's name has been besmirched by association with "Euler angles", but his original
work never expresses rotations as a product of three rotations about fixed axes.  Rather, he
came up with the axis-angle representation, and essentially devised the direction cosine
matrix as we know it today.  His 1776 papers [Euler_1776a](@cite) and  [Euler_1776b](@cite)
derived these ideas in service of understanding rigid-body dynamics — rather than any purely
mathematical investigations.

Nonetheless, by 1868, Tait [Tait_1868](@cite) writes down a rotation in exactly what we
would now call Euler angles — which he calls "the usual angles ψ, θ, ϕ"!  So somehow, this
had become normalized by that time.

Whittaker has a section I.9 on "Euler's parametric specification of rotations round a point"
and cites sections 6 and following of [Euler_1776b](@cite).  That section, I would say,
correctly cites Euler's work on direction cosines.  However, Whittaker's section I.10 on
"The Eulerian angles" cites Euler's *preceding* work [Euler_1776a](@cite), which does not

As far as I can find, Davenport is actually the first person to come along and try to
actually show when it is possible to decompose an arbitrary rotation into three rotations
about fixed axes [Davenport_1973](@cite).

"""
