"""

We can find conventions at [this
page](https://reference.wolfram.com/language/ref/WignerD.html).

> The Wolfram Language uses phase conventions where ``D^j_{m_1, m_2}(\psi, \theta, \phi) =
  \exp(i m_1 \psi + i m_2 \phi) D^j_{m_1, m_2}(0, \theta, 0)``.

> `WignerD[{1, 0, 1}, ψ, θ, ϕ]`
> ``-\sqrt{2} e^{i \phi} \cos\frac{\theta}{2} \sin\frac{\theta}{2}``

> `WignerD[{𝓁, 0, m}, θ, ϕ] == Sqrt[(4 π)/(2 𝓁 + 1)] SphericalHarmonicY[𝓁, m, θ, ϕ]`

> `WignerD[{j, m1, m2},ψ, θ, ϕ]] == (-1)^(m1 - m2) Conjugate[WignerD[{j, -m1, -m2}, ψ, θ, ϕ]]`

"""
𝐑
