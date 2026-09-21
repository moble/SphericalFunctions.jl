# [Differential operators](@id interface_differential_operators)

Each of the functions below returns the matrix of an angular-momentum
operator acting on mode weights in the canonical ordering described on
the [Utilities](@ref) page, for a given spin weight and range of
``ℓ``; applied to a [`ModeWeights`](@ref) instead, the same function
returns a new `ModeWeights`, with its spin weight adjusted where the
operator changes it.  The indices may be integers or half-integers,
the latter given as `Rational`s with denominator 2, but all of the
indices in one call must be of one type.  See the [background
page](@ref "Mode weights") for discussion, including the explanation
for why there cannot be `Rx` or `Ry` operators.


## Docstrings

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["utilities/operators.jl"]
```
