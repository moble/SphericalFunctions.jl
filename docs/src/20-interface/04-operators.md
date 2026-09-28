# [Differential operators](@id interface_differential_operators)

Each of the functions below returns the matrix of an angular-momentum
operator acting on mode weights in the canonical ordering described on
the [Utilities](@ref) page, for a given spin weight and range of
``ℓ``; applied to a [`ModeWeights`](@ref) instead, the same function
returns a new `ModeWeights`, with its spin weight adjusted where the
operator changes it.  The indices may be integers or half-integers,
the latter given as `Rational{Int}`s with denominator 2 or as
`HalfOddInteger`s, but all of the indices in one call must be of one
kind, as described on the [half-integer page](@ref
interface_half_integers).  See the [background page](@ref "Mode
weights") for discussion, including the explanation for why there
cannot be `Rx` or `Ry` operators.

Applied to mode weights, as `ð * w` or `ð(w)`, an operator builds no
matrix, and its result has the same range of ``ℓ`` as its input, even
where the spin weight changes; `ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` copies a
result into another range.  A product of several operators is applied
from right to left, so that `ð̄ * ð * w` is `ð̄ * (ð * w)`, which needs
no matrix either.  `mul!(w′, ð, w)` writes the result into a
destination that already has the labels of the result, and allocates
nothing.  The matrix form, `ð(s, ℓₘᵢₙ, ℓₘₐₓ)`, has no labels to check
against those of the weights, so it acts on the plain vector
`array_view(w)` rather than on the `ModeWeights` itself, whose product
with a plain matrix is refused.

Where the name of an operator is not ASCII, it has an ASCII alias,
public but not exported: `L2` for `L²`, `Lplus` and `Lminus` for `L₊`
and `L₋`, `R2`, `Rplus` and `Rminus` for `R²`, `R₊` and `R₋`, and
`eth` and `ethbar` for `ð` and `ð̄`.  [`Δspin`](@ref
SphericalFunctions.Δspin), the change in spin weight that an operator
makes, is also spelled `Deltaspin`.


## Docstrings

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["utilities/operators.jl"]
```
