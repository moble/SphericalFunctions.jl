raw"""

We can find conventions at [this
page](https://reference.wolfram.com/language/ref/WignerD.html).

> The Wolfram Language uses phase conventions where ``D^j_{m_1, m_2}(\psi, \theta, \phi) =
  \exp(i m_1 \psi + i m_2 \phi) D^j_{m_1, m_2}(0, \theta, 0)``.

> `WignerD[{1, 0, 1}, ψ, θ, ϕ]` ``-\sqrt{2} e^{i \phi} \cos\frac{\theta}{2}
> \sin\frac{\theta}{2}``

> `WignerD[{𝓁, 0, m}, θ, ϕ] == Sqrt[(4 π)/(2 𝓁 + 1)] SphericalHarmonicY[𝓁, m, θ, ϕ]`

> `WignerD[{j, m1, m2},ψ, θ, ϕ]] == (-1)^(m1 - m2) Conjugate[WignerD[{j, -m1, -m2}, ψ, θ,
> ϕ]]`

The Euler angles are defined generally such that

> `EulerMatrix[{α,β,γ},{a,b,c}]` is equivalent to ``R_{α,a} R_{β,b} R_{γ,c}``, where
> ``R_{α,a}``=`RotationMatrix[α,UnitVector[3,a]]`, etc.

and

> `EulerMatrix[{α,β,γ}]` is equivalent to `EulerMatrix[{α,β,γ},{3,2,3}]`

(representing the ``z-y-z`` convention).

Finally, we find that they say that `EulerMatrix`` corresponds to three rotations:

```mathematica
rα = RotationMatrix[α, {0, 0, 1}];
rβ = RotationMatrix[β, {0, 1, 0}];
rγ = RotationMatrix[γ, {0, 0, 1}];

Simplify[rα . rβ . rγ == EulerMatrix[{α, β, γ}]]
```

This agrees with the conventions used in this package, so we can directly compare
expressions in terms of Euler angles.

"""
