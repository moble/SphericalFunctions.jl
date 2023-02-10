@doc raw"""
    L²(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the total angular-momentum operator from `ℓₘᵢₙ` up to `ℓₘₐₓ`.  If not present,
`ℓₘᵢₙ` is assumed to be `abs(s)`.  The argument `s` is ignored; it is present only for
consistency with other operators, and is assumed to be 0 if not present.

This is the standard ``L²`` operator, familiar from basic physics, extended to work with
SWSHs.  It is also known as the Casimir operator, and is equal to
```math
L^2 = L_x^2 + L_y^2 + L_z^2 = \frac{L_+L_- + L_-L_+ + 2L_zL_z}{2}.
```
Note that these are the left Lie derivatives, but ``L^2 = R^2``, where `R` is the right Lie
derivative.  See [the documentation](/operators/) for more details.

In terms of the SWSHs, we can write the action of ``L^2`` as
```math
L^2 {}_{s}Y_{\ell,m} = \ell\,(\ell+1) {}_{s}Y_{\ell,m}
```

See also [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function L²(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    Diagonal(T[ℓ < abs(s) ? zero(m) : ℓ*(ℓ+1) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ])
end


@doc raw"""
    Lz(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the angular-momentum operator associated with the ``z`` direction.  This is the
standard ``L_z`` operator, familiar from basic physics, extended to work with SWSHs.  Note
that this is the left Lie derivative; see [`Rz`](@ref) for the equivalent right Lie
derivative.  See [the documentation](/operators/) for more details.

In terms of the SWSHs, we can write the action of ``L_z`` as
```math
L_z {}_{s}Y_{\ell,m} = m\, {}_{s}Y_{\ell,m}
```

See also [`L²`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function Lz(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    Diagonal(T[ℓ < abs(s) ? zero(m) : m for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ])
end


@doc raw"""
    L₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the angular-momentum raising operator.  This is the standard ``L_+`` operator,
familiar from basic physics, extended to work with SWSHs.  Note that this is the left Lie
derivative; see [`R₊`](@ref) for the equivalent right Lie derivative.  See [the
documentation](/operators/) for more details.

We define ``L_+`` to be the raising operator for the left Lie derivative with respect to
rotation about ``z``: ``L_z``.  By definition, this implies the commutator relation ``[L_z,
L_+] = L_+``, which allows us to derive ``L_+ = L_x + i\, L_y.``

In terms of the SWSHs, we can write the action of ``L_+`` as
```math
L_+ {}_{s}Y_{\ell,m} = \sqrt{(\ell-m)(\ell+m+1)}\, {}_{s}Y_{\ell,m+1}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{L_+(f)\right\}_{s,\ell,m} = \sqrt{(\ell+m)(\ell-m+1)}\,\left\{f\right\}_{s,\ell,m-1}.
```

See also [`L²`](@ref), [`Lz`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function L₊(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    Bidiagonal(
        zeros(T, (ℓₘₐₓ+1)^2-ℓₘᵢₙ^2),
        T[
            ℓ < abs(s) ? zero(T) : √T((ℓ+m)*(ℓ-m+1))
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ ifelse(ℓ==ℓₘᵢₙ,-ℓ+1,-ℓ):ℓ
        ],
        :L
    )
end


@doc raw"""
    L₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the angular-momentum lowering operator.  This is the standard ``L_-`` operator,
familiar from basic physics, extended to work with SWSHs.  Note that this is the left Lie
derivative; see [`R₋`](@ref) for the equivalent right Lie derivative.  See [the
documentation](/operators/) for more details.

We define ``L_-`` to be the lowering operator for the left Lie derivative with respect to
rotation about ``z``: ``L_z``.  By definition, this implies the commutator relation ``[L_z,
L_-] = -L_-``, which allows us to derive ``L_- = L_x - i\, L_y.``

In terms of the SWSHs, we can write the action of ``L_-`` as
```math
L_- {}_{s}Y_{\ell,m} = \sqrt{(\ell+m)(\ell-m+1)}\, {}_{s}Y_{\ell,m-1}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{L_-(f)\right\}_{s,\ell,m} = \sqrt{(\ell-m)(\ell+m+1)}\,\left\{f\right\}_{s,\ell,m+1}.
```

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref),
[`R₊`](@ref), [`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function L₋(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    Bidiagonal(
        zeros(T, (ℓₘₐₓ+1)^2-ℓₘᵢₙ^2),
        T[
            ℓ < abs(s) ? zero(T) : √T((ℓ-m)*(ℓ+m+1))
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ifelse(ℓ==ℓₘₐₓ,ℓ-1,ℓ)
        ],
        :U
    )
end



@doc raw"""
    R²(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the total angular-momentum operator from `ℓₘᵢₙ` up to `ℓₘₐₓ`.  If not present,
`ℓₘᵢₙ` is assumed to be `abs(s)`.  The argument `s` is ignored; it is present only for
consistency with other operators, and is assumed to be 0 if not present.

This is the ``R^2`` operator, much like the ``L^2`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  It is also known
as the Casimir operator, and is equal to
```math
R^2 = R_x^2 + R_y^2 + R_z^2 = \frac{R_+R_- + R_-R_+ + 2R_zR_z}{2}.
```
Note that these are the right Lie derivatives, but ``L^2 = R^2``, where `L` is the left Lie
derivative.  See [the documentation](/operators/) for more details.

In terms of the SWSHs, we can write the action of ``R^2`` as
```math
R^2 {}_{s}Y_{\ell,m} = \ell\,(\ell+1) {}_{s}Y_{\ell,m}
```

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`Rz`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function R²(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    Diagonal(T[ℓ < abs(s) ? zero(ℓ) : ℓ*(ℓ+1) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ])
end


@doc raw"""
    Rz(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the *right* angular-momentum operator associated with the ``z`` direction.

This is the ``R_z`` operator, much like the ``L_z`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  See [`Lz`](@ref)
for the equivalent left Lie derivative.  See [the documentation](/operators/) for more
details.

In terms of the SWSHs, we can write the action of ``R_z`` as
```math
R_z {}_{s}Y_{\ell,m} = -s\, {}_{s}Y_{\ell,m}
```
Note the unfortunate sign of ``s``, which seems to be opposite to what we expect, and arises
from the choice of definition of ``s`` in [the original paper by Newman and
Penrose](https://dx.doi.org/10.1063/1.1931221).

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`R₊`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    Diagonal(T[ℓ < abs(s) ? zero(ℓ) : -s for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ])
end


@doc raw"""
    R₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the *right* angular-momentum raising operator.

This is the ``R_+`` operator, much like the ``L_+`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  See [`L₊`](@ref)
for the equivalent left Lie derivative.  See [the documentation](/operators/) for more
details.

We define ``R_+`` to be the raising operator for the right Lie derivative with respect to
rotation about ``z``: ``R_z``.  By definition, this implies the commutator relation ``[R_z,
R_+] = R_+``, which allows us to derive ``R_+ = R_x - i\, R_y.``

In terms of the SWSHs, we can write the action of ``R_+`` as
```math
R_+ {}_{s}Y_{\ell,m} = \sqrt{(\ell+s)(\ell-s+1)}\, {}_{s-1}Y_{\ell,m}
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{R_+(f)\right\}_{s,\ell,m} = \sqrt{(\ell+s)(\ell-s+1)}\,\left\{f\right\}_{s-1,\ell,m}.
```
Because of the unfortunate sign of ``s`` arising from the choice of definition of ``s`` in
[the original paper by Newman and Penrose](https://dx.doi.org/10.1063/1.1931221), this is a
*lowering* operator for ``s``, though it really is a *raising* operator for ``R_z``, and
raises the eigenvalue of the corresponding Wigner matrix.

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref),
[`R₋`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    s′ = max(abs(s), abs(s-1))
    Diagonal(
        [
            ℓ < s′ ? zero(T) : √T((ℓ+s)*(ℓ-s+1))
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ
        ]
    )
end


@doc raw"""
    R₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute the *right* angular-momentum lowering operator.

This is the ``R_-`` operator, much like the ``L_-`` operator familiar from basic physics,
but in terms of the right Lie derivative, and extended to work with SWSHs.  See [`L₋`](@ref)
for the equivalent left Lie derivative.  See [the documentation](/operators/) for more
details.

We define ``R_-`` to be the raising operator for the right Lie derivative with respect to
rotation about ``z``: ``R_z``.  By definition, this implies the commutator relation ``[R_z,
R_-] = -R_-``, which allows us to derive ``R_- = R_x + i\, R_y.``

In terms of the SWSHs, we can write the action of ``R_-`` as
```math
R_- {}_{s}Y_{\ell,m} = \sqrt{(\ell-s)(\ell+s+1)}\, {}_{s+1}Y_{\ell,m}
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{R_-(f)\right\}_{s,\ell,m} = \sqrt{(\ell-s)(\ell+s+1)}\,\left\{f\right\}_{s+1,\ell,m}.
```
Because of the unfortunate sign of ``s`` arising from the choice of definition of ``s`` in
[the original paper by Newman and Penrose](https://dx.doi.org/10.1063/1.1931221), this is a
*raising* operator for ``s``, though it really is a *lowering* operator for ``R_z``, and
lowers the eigenvalue of the corresponding Wigner matrix - though that raises the eigenvalue
of the corresponding Wigner matrix.

See also [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref), [`Rz`](@ref),
[`R₊`](@ref), [`ð`](@ref), [`ð̄`](@ref).
"""
function R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    s′ = max(abs(s), abs(s+1))
    Diagonal(
        [
            ℓ < s′ ? zero(T) : √T((ℓ-s)*(ℓ+s+1))
            for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ
        ]
    )
end



@doc raw"""
    ð(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute coefficients for the spin-raising operator ``\eth``.

This operator was originally defined by [Newman and
Penrose](https://dx.doi.org/10.1063/1.1931221), but is more completely defined in [this
paper](https://arxiv.org/abs/1604.08140).  It is identical to [`R₋`](@ref).  Refer to that
function's documentation for more details.

By definition, the spin-raising operator satisfies the commutator relation ``[S, \eth] =
\eth`` (where ``S`` is the spin operator, which just multiplies the function by its spin).
In terms of the SWSHs, we can write the action of ``\eth`` as
```math
    \eth {}_{s}Y_{\ell,m} = \sqrt{(\ell-s)(\ell+s+1)} {}_{s+1}Y_{\ell,m}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{\eth f\right\}_{s,\ell,m} = \sqrt{(\ell-s)(\ell+s+1)}\,\left\{f\right\}_{s+1,\ell,m}.
```

See also [`ð̄`](@ref),  [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref),
[`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref).
"""
function ð(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
end


@doc raw"""
    ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, [T])

Compute coefficients for the spin-lowering operator ``\bar{\eth}``.

This operator was originally defined by [Newman and
Penrose](https://dx.doi.org/10.1063/1.1931221), but is more completely defined in [this
paper](https://arxiv.org/abs/1604.08140).  It is opposite to [`R₊`](@ref) — meaning that
``\bar{\eth} = -R₊``.  Refer to that function's documentation for more details.

By definition, the spin-lowering operator satisfies the commutator relation ``[S,
\bar{\eth}] = -\bar{\eth}`` (where ``S`` is the spin operator, which just multiplies the
function by its spin).  In terms of the SWSHs, we can write the action of ``\bar{\eth}`` as
```math
\bar{\eth} {}_{s}Y_{\ell,m} = -\sqrt{(\ell+s)(\ell-s+1)} {}_{s-1}Y_{\ell,m}.
```
Consequently, the *mode weights* of a function are affected as
```math
\left\{\bar{\eth} f\right\}_{s,\ell,m}
= -\sqrt{(\ell-s)(\ell+s+1)}\,\left\{f\right\}_{s+1,\ell,m}.
```

See also [`ð`](@ref),  [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref), [`R²`](@ref),
[`Rz`](@ref), [`R₊`](@ref).
"""
function ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T=Float64)
    -R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
end
