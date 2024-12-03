raw"""

The Condon-Shortley phase convention is a choice of phase factors in the definition of the
spherical harmonics that requires the coefficients in
```math
L_{\pm} |\ell,m\rangle = \alpha^{\pm}_{\ell,m} |\ell, m \pm 1\rangle
```
to be real and positive.

Eq. (44b) of [Boyle (2016)](@cite Boyle_2016) says
```math
L_{\pm} \mathfrak{D}^{(\ell)}_{m',m}(\mathbf{R})
= \sqrt{(\ell \mp m')(\ell \pm m' + 1)} \mathfrak{D}^{(\ell)}_{m' \pm 1, m}(\mathbf{R}).
```
while Eq. (21) relates the Wigner D-matrix to the spin-weighted spherical harmonics as
```math
{}_{s}Y_{\ell,m}(\mathbf{R})
= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \mathfrak{D}^{(\ell)}_{m,-s}(\mathbf{R}).
```
Plugging the latter into the former, we get
```math
L_{\pm} {}_{s}Y_{\ell,m}(\mathbf{R})
= \sqrt{(\ell \mp m)(\ell \pm m + 1)} {}_{s}Y_{\ell,m \pm 1}(\mathbf{R}).
```
That is, in our conventions we have
```math
\alpha^{\pm}_{\ell,m} = \sqrt{(\ell \mp m)(\ell \pm m + 1)},
```
which is always real and positive, and thus consistent with the Condon-Shortley phase
convention.

"""
