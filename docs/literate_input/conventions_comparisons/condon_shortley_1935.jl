md"""
# Condon-Shortley (1935)

[Condon and Shortley's "The Theory Of Atomic Spectra"](@cite CondonShortley_1935) is the
standard reference for the "Condon-Shortley phase convention".  Though some references are
not very clear about precisely what they mean by that phrase, it seems clear that the real
meaning includes the idea that the angular-momentum raising and lowering operators have
*real* eigenvalues when acting on the spherical harmonics.  To avoid ambiguity, we can just
look at the actual spherical harmonics they define.

The method we use here is as direct and explicit as possible.  In particular, Condon and
Shortley provide a formula for the φ=0 part in terms of iterated derivatives of a power of
sin(θ).  Rather than expressing these derivatives in terms of the Legendre polynomials —
which would subject us to another round of ambiguity — the functions in this module use
automatic differentiation to compute the derivatives explicitly.

Condon and Shortley are not very explicit about the meaning of the spherical coordinates,
but they do describe them as "spherical polar coordinates ``r, \theta, \varphi``"
immediately before equation (1) of section 4³ (page 50),
```math
L_z = -i \hbar \frac{\partial}{\partial \varphi},
```
followed by equation (8):
```math
\begin{aligned}
L_x + i L_y &= \hbar e^{i\varphi} \left(
  \frac{\partial}{\partial \theta}
  + i \cot\theta \frac{\partial}{\partial \varphi}
\right) \\
L_x - i L_y &= \hbar e^{-i\varphi} \left(
  -\frac{\partial}{\partial \theta}
  + i \cot\theta \frac{\partial}{\partial \varphi}
\right).
\end{aligned}
```

Because these expressions agree nicely with our results,


The result is that the original Condon-Shortley spherical harmonics agree perfectly with the
ones computed by this package.

(Condon and Shortley do not give an expression for the Wigner D-matrices.)

"""

using TestItems: @testmodule, @testitem  #hide

# ## Function definitions
@testmodule CondonShortley1935 begin  #hide

#+
# and so on
const 𝒾 = im

end  #hide

# ## Tests

@testitem "Condon-Shortley (1935)" setup=[CondonShortley1935] begin  #hide

#+
# Here's a test
@test 2+2 == 4

#+
# And another
@test CondonShortley1935.𝒾^2 == -1

end  #hide
