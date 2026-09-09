# Comparisons

Here, we compare our conventions to other sources, including
references in the literature as well as other software that implements
some of these.  Each comparison has its own page (listed in the
sidebar and in the table below), which quotes the source's actual
equations, implements them in Julia, and tests them numerically
against this package.  Those pages are written with
[Literate.jl](https://fredrikekre.github.io/Literate.jl/), so the very
same files are also part of [this package's test
suite](https://github.com/moble/SphericalFunctions.jl/tree/main/docs/literate_input/conventions/comparisons).

Among the items that would be good to compare are the following, when
actually used by any of these sources:
* Quaternions
  - Order of components
  - Basis and multiplication table
  - Operation as rotations
* Euler angles
* Spherical coordinates
* Angular momentum operators
  - Fundamental definitions
  - Expression in terms of spherical coordinates
  - Expression in terms of Euler angles
  - Right-derivative form
* Spherical harmonics
  - Condon-Shortley phase
  - Formula
* Spin-weighted spherical harmonics
  - Behavior under rotation
* Wigner D-matrices
  - Representation à la $\langle ℓ, m' | e^{-i α J_z} e^{-i β J_y} e^{-i γ J_z} | ℓ, m \rangle$
  - Rotation of spherical harmonics
  - Order of indices
  - Conjugation
  - Function of rotation or inverse rotation
  - Formula

One major result of this is that almost everyone since 1935 has used
the same exact expression for the (scalar) spherical harmonics.

When choosing conventions for this package, I intend to prioritize
consistency (to the extent that any of these references actually have
anything to say about the above items) with the following sources, in
order:

1. LALSuite
2. NINJA
3. Newman-Penrose
4. Goldberg
5. Thorne / MTW
6. Wikipedia
7. Sakurai
8. Shankar
9. Zettili

I think that should be sufficient to find a consensus on conventions
for each of the above — with the possible exception of quaternions,
for which I have my own strong opinions.



## Summary of comparisons

| Reference | Verdict |
|-----------|---------|
| [Blanchet (2024)](@ref "Blanchet (2024)") | Spin-weight ``-2`` spherical harmonics and ``d`` function agree with ours. |
| [Boyle (2016)](@ref "Boyle (2016)") | ``𝔇`` matrices (this package before v3.0) are the complex conjugates of ours, for integer and half-integer indices; spin-weighted spherical harmonics agree. |
| [Clifford (1878)](@ref "Clifford (1878)") | Introduces Clifford algebras; identifies quaternions with bivectors of the *opposite* sign to ours (``i = ι₂ι₃``), and does not use the sandwich product for rotations. |
| [Cohen-Tannoudji (1991)](@ref "Cohen-Tannoudji (1991)") | Spherical harmonics agree with ours. |
| [Condon-Shortley (1935)](@ref "Condon-Shortley (1935)") | Spherical harmonics agree with ours; this is the origin of the phase convention we use. |
| [Edmonds (1960)](@ref "Edmonds (1960)") | Euler angles, operators, and spherical harmonics agree with ours; his ``𝒟`` is our ``𝔇`` of the *inverse* rotation (the Hermitian conjugate), and his ``d`` is the transpose of ours. |
| [Euler (1767 and 1776)](@ref "Euler (1767 and 1776)") | Euler angles can be read in our convention up to an offset in ``γ``; spherical coordinates consistent with ours. |
| [Gibbs (1881)](@ref "Gibbs (1881)") | Right-handed ``(i, j, k)`` basis, as ours. |
| [Goldberg et al. (1967)](@ref "Goldberg et al. (1967)") | Spin-weighted spherical harmonics are ``(-1)^m`` times ours; ``D`` is the complex conjugate of ours with ``α`` and ``γ`` interchanged (equivalently ``(-1)^{m+m'}`` times the conjugate transpose); ``\eth``, ``\bar{\eth}`` agree. |
| [Griffiths (1995)](@ref "Griffiths (1995)") | Spherical harmonics agree with ours. |
| [Hamilton (1844 and 1853)](@ref "Hamilton (1844 and 1853)") | Quaternion multiplication table as ours; Hamilton's own axes (south, west, up) form a left-handed system. |
| [LALSuite (2025)](@ref "LALSuite (2025)") | Spin-weighted spherical harmonics, ``d``, and ``D`` all agree with ours. |
| [Le Bellac (2006)](@ref "Le Bellac (2006)") | ``D``, the rotation law, and the ``Y``–``D`` relation agree with ours. |
| [Mathematica (2026)](@ref "Mathematica (2026)") | `SphericalHarmonicY` and `EulerMatrix` agree with ours; every documented identity of `WignerD` is satisfied by ``𝔇_{-m_1,-m_2} = (-1)^{m_1-m_2}\overline{𝔇_{m_1 m_2}}`` (Wigner's convention). |
| [Newman-Penrose (1966)](@ref "Newman-Penrose (1966)") | Spin weight, ``\eth``, and ``\bar{\eth}`` agree with ours (``\eth = R_+``, ``\bar{\eth} = -R_-``); their unnormalized stereographic ``{}_sY_{ℓm}`` is ``(-1)^ℓ e^{-isϕ}\sqrt{4π/[(2ℓ+1)(ℓ+m)!(ℓ-m)!]}`` times ours. |
| [NIST DLMF (2026)](@ref "NIST DLMF (2026)") | Spherical harmonics (via Ferrers functions) agree with ours. |
| [NINJA (2011)](@ref "NINJA (2011)") | Spin-weighted spherical harmonics and ``d`` agree with ours. |
| [Sakurai (1994)](@ref "Sakurai (1994)") | Spherical harmonics, ``d``, and ``𝒟`` agree with ours. |
| [SciPy (2026)](@ref "SciPy (2026)") | `sph_harm_y` agrees with our spherical harmonics (note the argument order `(n, m, θ, ϕ)` with ``θ`` polar). |
| [Shankar (1994)](@ref "Shankar (1994)") | Spherical harmonics agree with ours; matrix elements of ``U[R(α,β,γ)]`` agree with our ``𝔇``. |
| [SymPy (2026)](@ref "SymPy (2026)") | `Ynm` agrees with ours; `wigner_d` is ``𝔇_{-m',-m} = (-1)^{m'-m}\overline{𝔇_{m'm}}`` (Wigner's convention), despite citing Edmonds. |
| [Tait (1867 and 1868)](@ref "Tait (1867 and 1868)") | Right-handed basis and multiplication table as ours; Euler angles ``(ψ, θ, ϕ)`` in the ``z``-``y'``-``z''`` form equal our ``(ϕ, θ, ψ)``, and his quaternion matches ours. |
| [Thorne (1980)](@ref "Thorne (1980)") | Scalar spherical harmonics agree with ours. |
| [Torres del Castillo (2003)](@ref "Torres del Castillo (2003)") | ``d``, ``D``, and spin-weighted spherical harmonics agree with ours. |
| [Varshalovich et al. (1988)](@ref "Varshalovich et al. (1988)") | ``D`` and ``d`` agree with ours (indices read in their order), for integer and half-integer ``J``; lab-fixed operators are our ``L``; body-fixed operators are ``(-R_x, R_y, -R_z)``. |
| [Whittaker (1904)](@ref "Whittaker (1904)") | Right-handed axes, Euler angles equal to ours with the symbols permuted, spherical coordinates and quaternion basis as ours, rotation by ``q v q⁻¹``. |
| [Wigner (1959)](@ref "Wigner (1959)") | Representation matrices are ``(-1)^{μ'-μ}`` times the complex conjugate of ours (equivalently ours with both indices negated); spherical harmonics agree. |
| [Wikipedia (2026)](@ref "Wikipedia (2026)") | ``D``, ``d``, ``Y_ℓ^m``, and ``{}_sY_{ℓm}`` agree with ours; its body-fixed operators are ``\mathcal{P} = -R``. |
| [Wilson (1929)](@ref "Wilson (1929)") | Right-handed ``(𝐢, 𝐣, 𝐤)`` basis, as ours. |
| [Zettili (2009)](@ref "Zettili (2009)") | Spherical harmonics, ``d``, ``D``, and the rotation law agree with ours. |
