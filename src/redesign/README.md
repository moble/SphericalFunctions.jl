src/
├── SphericalFunctions.jl
├── utilities
│   ├── utilities.jl
│   ├── quadrature_weights.jl
│   ├── pixelizations.jl
│   └── complex_powers.jl
├── ssht
│   ├── ssht.jl  # pixels, rotors, mul!, ldiv!, synthesis, analysis, map2salm, salm2map
│   ├── direct.jl
│   ├── minimal.jl
│   └── reinecke_seljebotn.jl
├── wigner
│   ├── wigner.jl
│   ├── wigner_matrix.jl
│   ├── wigner_calculator.jl
│   ├── wigner_H.jl
│   ├── wigner_H_calculator.jl
│   └── recurrence.jl
└── mode_weights
    ├── mode_weights.jl
    └── operators.jl


# Single-Call API
```julia
D(R::AbstractQuaternion, ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ)::Vector{Matrix{Complex{T}}}
D(R::Vector{AbstractQuaternion}, ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ)::Matrix{Matrix{Complex{T}}}
```

# Incremental/Streaming API

```julia
# For iterative computation (saves memory, enables streaming)
calc = DCalculator(rotors, ℓₘₐₓ)

for ℓ in 0:ℓₘₐₓ
    Dˡ = next!(calc)  # Compute next ℓ, returns view
    # ... use Dˡ ...
end

# Or manual control
calc = DCalculator(rotors, ℓₘₐₓ)
D⁸ = compute!(calc, 8)  # Jump to specific ℓ
D⁹ = compute!(calc, 9)  # Incrementally compute next
```

# Reusable/Mutable API

```julia
# Allocate once, recompute for different rotors
calc = DCalculator(ℓₘₐₓ, Nᵣ)  # Pre-allocate storage

# Compute for first batch
D_matrices = compute!(calc, rotors₁)

# Recompute for different rotors (no allocation)
compute!(calc, rotors₂)  # Updates in-place

# Access results
Dˡ = calc[ℓ]  # View of Dˡ for current rotors
```


* SSHT
  - SSHTDirect
  - SSHTMinimal
  - SSHTRS
  * pixels
  * rotors
  * mul!(f, 𝒯, f̃)
  * ldiv!(f̃, 𝒯, f)
  * synthesis
  * analysis

* WignerRange
* AbstractWignerMatrix{IT, NT, ST}
  - WignerMatrix{IT, NT, ST} (rectangular or square sub-array of Dˡ)
    - WignerDMatrix{IT, RT<:Real, ST} = WignerMatrix{IT, Complex{RT}, ST}
    - WignerMatrix{IT, RT<:Real, ST} = WignerMatrix{IT, RT, ST}
  - HWedge{IT, RT<:Real, ST}
  - HAxis{IT, RT<:Real} (ST is implicitly FixedSizeVectorDefault{RT})
  - ?WignerdWedge (stores an HWedge; symmetrizes and adjusts phase on the fly)
  - ?WignerDWedge (" ; stores vectors of powers of eⁱᵅ, eⁱᵞ)
  - ?ₛYₗₘWedge (" ; stores vectors of powers of eⁱᵅ, eⁱᵞ)

* ?AbstractY:
  - ?YVector
* AbstractModeWeights:
  - ModeWeightsSymmetric: 0, 1, -1, 2, -2, 3, -3, ...
  - ModeWeightsIncreasing: ..., -3, -2, -1, 0, 1, 2, 3, ...
  - ModeWeightsDecreasing: ..., 3, 2, 1, 0, -1, -2, -3, ...
* ModeWeightOperators
  - L²
  - Lx
  - Ly
  - Lz
  - L₊
  - L₋
  - R²
  - Rx
  - Ry
  - Rz
  - R₊
  - R₋
  - ð
  - ð̄
