# SphericalFunctions.jl Package Structure

This document provides a comprehensive hierarchical overview of the SphericalFunctions.jl package, focusing on types, functions, and capabilities with an eye toward identifying gaps and potential areas for new implementations.

## Core Module: `SphericalFunctions`

### Type Definitions

#### Abstract Types
- **`SSHT{T<:Real}`** - Supertype for spin-spherical-harmonic transforms
  - Purpose: Base type for all SHT implementations
  - Type parameter: `T` (real number type for computations)

#### Concrete Types

##### Transform Types (subtypes of `SSHT`)
- **`SSHTDirect{T<:Real, Inplace, Tdecomp}`**
  - Direct method for spin-weighted spherical harmonic transforms
  - Fields:
    - `s::Integer` - Spin weight
    - `ℓₘₐₓ::Integer` - Highest ℓ value
    - `Rθϕ::Vector{Rotor{T}}` - Evaluation points
    - `ₛ𝐘::Matrix{Complex{T}}` - SH values
    - `ₛ𝐘decomposition::Tdecomp` - Decomposed matrix for inversion
  - Best for: `ℓₘₐₓ ≲ 24`
  - Supports in-place operations

- **`SSHTMinimal{T<:Real, Inplace}`**
  - Minimal sampling algorithm (Elahi et al., 2018)
  - Fields:
    - `s::Integer` - Spin weight
    - `ℓₘₐₓ::Integer` - Highest ℓ value
    - `θ::OffsetVector` - Colatitudes of sampling rings
    - `θindices::OffsetVector` - Index ranges for each θ ring
    - `plans::OffsetVector` - FFT plans
    - `bplans::OffsetVector` - Backward FFT plans
    - `ₛΛ::OffsetVector` - Matrices for mode weight computation
    - `luₛΛ::OffsetVector` - LU-decomposed ₛΛ matrices
  - Supports in-place operations
  - Uses minimal number of function samples

- **`SSHTRS{T<:Real}`**
  - Reinecke & Seljebotn (2013) algorithm
  - Fields:
    - `s` - Spin weight
    - `ℓₘₐₓ` - Highest ℓ value
    - `θ` - Ring locations (colatitude)
    - `quadrature_weight` - Quadrature weights
    - `Nϕ` - Number of azimuthal points per ring
    - `iθ` - Index range for each ring
    - `G` - Preallocated storage for FTs
    - `plan` - FFT plan
  - Best for larger `ℓₘₐₓ` values

##### Recursion Coefficient Types
- **`ALFRecursionCoefficients{T<:Real}`**
  - Stores coefficients for Associated Legendre Function recursion
  - Fields:
    - `nmax::Int` - Maximum n value
    - `a::Vector{T}` - First set of coefficients
    - `b::Vector{T}` - Second set of coefficients
    - `cde::Matrix{T}` - Combined c, d, e coefficients
  - Based on Xing et al. (2019) paper

##### Iterator Types
- **`D_iterator{VT<:Vector}`**
  - Iterates over sub-matrices of Wigner D matrices
  - Fields: `D`, `ℓₘₐₓ`, `ℓₘᵢₙ`
  - Returns views into original D data
  - ⚠️ Note: Returns transpose of documented behavior

- **`d_iterator{VT<:Vector}`**
  - Iterates over sub-matrices of reduced Wigner d matrices
  - Fields: `d`, `ℓₘₐₓ`, `ℓₘᵢₙ`
  - Returns views into original d data
  - ⚠️ Note: Returns transpose of documented behavior

- **`sYlm_iterator{VT<:Vector}`**
  - Iterates over sub-vectors of spin-weighted spherical harmonics
  - Fields: `Y`, `ℓₘₐₓ`, `ℓₘᵢₙ`, `iₘᵢₙ`
  - Returns views into original Y data

- **`λ_iterator{T}`**
  - Iterates over ₛλₗₘ values (real-valued SWSHs at ϕ=0)
  - Fields: `cosθ`, `sin½θ`, `cos½θ`, `s`, `m`
  - Infinite iterator over increasing ℓ values
  - Based on Kostelec & Rockmore (2008)

- **`AlternatingCountdown`**
  - Utility iterator: counts down with alternating signs
  - Field: `start::Int`

- **`AlternatingCountup`**
  - Utility iterator: counts up with alternating signs
  - Field: `stop::Int`

##### Utility Types
- **`OffsetVec{T,VV<:AbstractVector{T}}`**
  - Vector with index offset
  - Fields: `parent::VV`, `offset::Int`
  - Enables flexible array indexing

- **`OffsetMat{T,MM<:AbstractMatrix{T}}`**
  - Matrix with index offsets
  - Fields: `parent::MM`, `offset1::Int`, `offset2::Int`

#### Type Aliases
- **`MachineFloat = Union{Float16, Float32, Float64}`**
  - Convenience type for machine-precision floats

### Functions by Category

#### 1. Pixelization Functions
Generate sampling points on the sphere

- **`golden_ratio_spiral_pixels(s, ℓₘₐₓ, [T=Float64])`**
  - Returns: Vector of 2-SVectors with (θ, ϕ) coordinates
  - Also known as: Fibonacci sphere/lattice
  - No pixels at poles
  
- **`golden_ratio_spiral_rotors(s, ℓₘₐₓ, [T=Float64])`**
  - Returns: Vector of `Rotor`s
  - Companion to `golden_ratio_spiral_pixels`

- **`sorted_rings(s, ℓₘₐₓ, [T=Float64])`**
  - Returns: Colatitudes of rings sorted by proximity to equator
  - For condition number optimization

- **`sorted_ring_pixels(s, ℓₘₐₓ, [T=Float64])`**
  - Returns: Vector of 2-SVectors on sorted rings

- **`sorted_ring_rotors(s, ℓₘₐₓ, [T=Float64])`**
  - Returns: Vector of `Rotor`s on sorted rings

- **`fejer1_rings(N, [T=Float64])`**
  - Returns: θ values for Fejér's first quadrature rule

- **`fejer2_rings(N, [T=Float64])`**
  - Returns: θ values for Fejér's second quadrature rule

- **`clenshaw_curtis_rings(N, [T=Float64])`**
  - Returns: θ values for Clenshaw-Curtis quadrature

#### 2. Complex Powers
- **`complex_powers(z, pₘᵢₙ, pₘₐₓ)`**
  - Compute powers of complex number z from pₘᵢₙ to pₘₐₓ

- **`complex_powers!(zpowers, z, pₘᵢₙ, pₘₐₓ)`**
  - In-place version of `complex_powers`

#### 3. Indexing Functions
Map between (ℓ, m, m′) indices and array positions

##### Y (Spherical Harmonics) Indexing
- **`Ysize(ℓₘₐₓ)`** / **`Ysize(ℓₘᵢₙ, ℓₘₐₓ)`**
  - Returns: Total size of mode weight array

- **`Yrange(ℓₘₐₓ)`** / **`Yrange(ℓₘᵢₙ, ℓₘₐₓ)`**
  - Returns: Array of (ℓ, m) indices

- **`Yindex(ℓ, m, [ℓₘᵢₙ=0])`**
  - Returns: Index into mode weight array

- **`deduce_limits(ysize, [ℓmin])`**
  - Deduce (ℓₘᵢₙ, ℓₘₐₓ) from array size

##### Grid Construction
- **`theta_phi(Nθ, Nϕ, [T=Float64])`**
  - Returns: (θ, ϕ) grid in spinsfast order

- **`phi_theta(Nϕ, Nθ, [T=Float64])`**
  - Returns: (ϕ, θ) grid in standard order

##### Wigner H Indexing
- **`WignerHsize(ℓₘₐₓ, [m′ₘₐₓ=ℓₘₐₓ])`**
  - Returns: Size of H wedge array

- **`WignerHrange(ℓₘₐₓ, [m′ₘₐₓ=ℓₘₐₓ])`**
  - Returns: Array of (ℓ, m′, m) indices

- **`WignerHindex(ℓ, m′, m, [m′ₘₐₓ=ℓ])`**
  - Returns: Index into H wedge array

- **`_WignerHindex(ℓ, m′, m, m′ₘₐₓ)`**
  - Helper function assuming m ≥ |m′|

##### Wigner D Indexing
- **`WignerDsize(ℓₘₐₓ, [m′ₘₐₓ=ℓₘₐₓ])`**
  - Returns: Size of Wigner D matrix

- **`WignerDrange(ℓₘₐₓ, [m′ₘₐₓ=ℓₘₐₓ])`**
  - Returns: Array of (ℓ, m′, m) indices

- **`WignerDindex(ℓ, m′, m, [m′ₘₐₓ=ℓ])`**
  - Returns: Index into Wigner D matrix

#### 4. Associated Legendre Functions
Based on Xing et al. (2019)

- **`ALFRecursionCoefficients(nmax, [T=Float64])`**
  - Pre-compute recursion coefficients
  - 8x speedup when reused

- **`ALFrecurse!(P̄ₙ, P̄ₙ₋₁, n, t, u, aₙ, bₙ, cdeₙ)`**
  - Recursion step with pre-computed coefficients

- **`ALFrecurse!(P̄ₙ, P̄ₙ₋₁, n, t, u)`**
  - Recursion step computing coefficients on-the-fly

- **`ALFcompute(expiβ, nmax, [recursion_coefficients])`**
  - Compute fully normalized ALFs up to nmax

- **`ALFcompute!(P̄, expiβ, nmax, [recursion_coefficients])`**
  - In-place version of `ALFcompute`

#### 5. H Matrix Recursion
Based on Gumerov & Duraiswami (2015)

- **`H_recursion_coefficients(ℓₘₐₓ, T)`**
  - Pre-compute constants for H recursion
  - Returns: Tuple (aₙᵐ, bₙᵐ, dₙᵐ)

- **`H!(H, expiβ, ℓₘₐₓ, m′ₘₐₓ, H_rec_coeffs, [Hindex=WignerHindex])`**
  - Compute H matrix elements
  - Foundation for d, D, and sYlm computations

#### 6. Wigner Matrices Evaluation

##### Reduced Wigner d Matrices
- **`d_matrices(β, ℓₘₐₓ)`** / **`d_matrices(expiβ, ℓₘₐₓ)`**
  - Compute all d matrices up to ℓₘₐₓ

- **`d_matrices!(d_storage, β)`** / **`d_matrices!(d, expiβ, ℓₘₐₓ)`**
  - In-place computation

- **`d_prep(ℓₘₐₓ, [T=Float64])`**
  - Pre-allocate storage for d matrices
  - Returns: Tuple for use with `d_matrices!`

- **`d!(d, expiβ, ℓₘₐₓ, H_rec_coeffs)`** (Legacy API)
  - Direct computation with H coefficients

##### Full Wigner D Matrices
- **`D_matrices(R, ℓₘₐₓ)`**
  - Compute all D matrices from rotation R
  - R can be: Rotor, quaternion, or Euler angles

- **`D_matrices!(D_storage, R)`** / **`D_matrices!(D, R, ℓₘₐₓ)`**
  - In-place computation

- **`D_prep(ℓₘₐₓ, [T=Float64])`**
  - Pre-allocate storage for D matrices

- **`D!(D, R, ℓₘₐₓ, H_rec_coeffs)`** (Legacy API)
  - Direct computation

##### Preparation Functions
- **`d_prep(ℓₘₐₓ, [T=Float64])`**
  - Allocate storage and pre-compute coefficients for d

- **`D_prep(ℓₘₐₓ, [T=Float64])`**
  - Allocate storage and pre-compute coefficients for D

- **`dprep`, `Dprep`** (Legacy aliases)

#### 7. Spin-Weighted Spherical Harmonics

- **`sYlm_values(R, s, ℓₘₐₓ)`**
  - Evaluate SWSHs at rotation R
  - Returns: Mode weights for all (ℓ, m)

- **`sYlm_values!(sYlm_storage, R)`** / **`sYlm_values!(sYlm, R, s, ℓₘₐₓ)`**
  - In-place evaluation

- **`sYlm_prep(s, ℓₘₐₓ, [T=Float64])`**
  - Pre-allocate storage for sYlm

- **`Y!(Y, R, s, ℓₘₐₓ, H_rec_coeffs)`** (Legacy API)
  - Direct computation

- **`Yprep`** (Legacy alias)

- **`ₛ𝐘`** (Legacy alias for sYlm_values)

#### 8. Quadrature Weights

- **`fejer1(N, [T=Float64])`**
  - Fejér's first rule quadrature weights

- **`fejer2(N, [T=Float64])`**
  - Fejér's second rule quadrature weights

- **`clenshaw_curtis(N, [T=Float64])`**
  - Clenshaw-Curtis quadrature weights

#### 9. Spherical Harmonic Transforms

##### Constructor
- **`SSHT(s, ℓₘₐₓ; method="Direct", [kwargs...])`**
  - Factory function for SHT objects
  - Methods: "Direct", "Minimal", "RS"

##### Transform Operations
- **`pixels(𝒯)`**
  - Get spherical coordinates (θ, ϕ) for transform 𝒯

- **`rotors(𝒯)`**
  - Get Rotor evaluation points for transform 𝒯

##### Algebraic Operations
Transforms support algebraic semantics:
- **`f = 𝒯 * f̃`** - Synthesis (modes → values)
- **`f̃ = 𝒯 \ f`** - Analysis (values → modes)
- **`mul!(f, 𝒯, f̃)`** - In-place synthesis
- **`ldiv!(f̃, 𝒯, f)`** - In-place analysis

#### 10. Map to Mode Transform

- **`map2salm(map, spin, ℓmax)`**
  - Transform function values to mode weights
  - Follows Reinecke & Seljebotn (2013)

- **`map2salm!(salm, map, spin, ℓmax)`**
  - In-place version

- **`map2salm(map, plan)`** / **`map2salm!(salm, map, plan)`**
  - Use pre-computed plan

- **`plan_map2salm(map, spin, ℓmax)`**
  - Pre-compute transform plan

#### 11. Differential Operators

All operators return diagonal or bidiagonal matrices acting on mode weights.

##### Left Lie Derivative Operators
- **`L²(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Total angular momentum operator
  - Eigenvalue: ℓ(ℓ+1)

- **`Lz(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - z-component angular momentum
  - Eigenvalue: m

- **`L₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Raising operator for Lz
  - Returns: Bidiagonal matrix

- **`L₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Lowering operator for Lz
  - Returns: Bidiagonal matrix

##### Right Lie Derivative Operators
- **`R²(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Total angular momentum (right)
  - Eigenvalue: ℓ(ℓ+1)
  - Note: R² = L²

- **`Rz(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - z-component (right)
  - Eigenvalue: -s

- **`R₊(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Raising operator for Rz
  - Changes spin: s → s-1

- **`R₋(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Lowering operator for Rz
  - Changes spin: s → s+1

##### Spin Operators (Newman-Penrose)
- **`ð(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Spin-raising operator (eth)
  - Identical to R₋

- **`ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, [T])`**
  - Spin-lowering operator (eth-bar)
  - Equals -R₊

#### 12. Utility Functions

- **`sqrtbinomial(n, k, [T])`**
  - Square root of binomial coefficient
  - Handles large n (up to ~1026)
  - Uses logarithmic computation

- **`logbinomial(n, k, S)`**
  - Logarithm of binomial coefficient

- **`loggamma(a, T)`**
  - Logarithm of gamma function

##### Lambda Recursion Helpers
- **`λ_recursion_initialize(sin½θ, cos½θ, s, ℓ, m)`**
  - Initial values for λ recursion

- **`λ_recursion_coefficients(cosθ, s, ℓ, m)`**
  - Coefficients for λ recursion

##### Helper Constants
- **`ϵ(m)`** - Inline function: (-1)^m for odd m > 0, else 1

### File Organization

```
src/
├── SphericalFunctions.jl      # Main module file
├── utils.jl                   # Utility types and functions
├── pixelizations.jl           # Sphere sampling functions
├── complex_powers.jl          # Complex power computations
├── indexing.jl                # Index computation functions
├── iterators.jl               # Iterator types for matrices
├── associated_legendre.jl     # ALF computations
├── Hrecursion.jl              # H matrix recursion
├── evaluate.jl                # Wigner matrix evaluation
├── weights.jl                 # Quadrature weights
├── ssht.jl                    # SHT interface
├── ssht/
│   ├── direct.jl             # Direct method implementation
│   ├── minimal.jl            # Minimal sampling implementation
│   └── rs.jl                 # Reinecke-Seljebotn implementation
├── map2salm.jl               # Map to mode transforms
└── operators.jl              # Differential operators
```

## Identified Gaps and Potential Extensions

### 1. Missing Type Abstractions
- **No explicit type for mode weights**
  - Could create `ModeWeights{T,S}` type wrapping arrays with spin/ℓ metadata
  - Would enable type-safe operations and clearer APIs

- **No type for pixelization/sampling schemes**
  - Could create `SamplingScheme{T}` abstract type
  - Concrete types: `GoldenRatioSpiral`, `SortedRings`, `FejerQuadrature`, etc.
  - Would unify pixelization functions

- **No explicit grid/mesh type**
  - Could create `SphericalGrid{T}` for organized grid data
  - Would store θ, ϕ coordinates and connectivity

### 2. Missing Functionality

#### Transform Extensions
- **Inverse transforms not explicitly named**
  - Could add `synthesis` and `analysis` as named operations
  - Currently only available via `*` and `\` operators

- **No batch transform support**
  - Could add `BatchSSHT` for multiple simultaneous transforms
  - Useful for vector/tensor fields

- **No spherical convolution**
  - Natural operation given transform capabilities
  - Would enable filtering operations

#### Grid Operations
- **No interpolation functions**
  - Grid → arbitrary points interpolation
  - Useful for visualization and analysis

- **No grid refinement/coarsening**
  - Adaptive mesh refinement
  - Hierarchical representations

#### Analysis Functions
- **No power spectrum computation**
  - Natural given mode weights
  - Would add `power_spectrum(modes)`

- **No correlation functions**
  - Cross-correlation of fields
  - Angular correlation functions

### 3. Missing Conveniences

#### Constructor Helpers
- **No convenience constructors for common cases**
  - e.g., `SSHT(ℓₘₐₓ)` for scalar (s=0) harmonics
  - Auto-select method based on ℓₘₐₓ

- **No validation/checking utilities**
  - Validate mode weight array sizes
  - Check sampling adequacy

#### Conversion Functions
- **No coordinate system conversions**
  - Cartesian ↔ spherical
  - Different spherical conventions

- **No mode weight format conversions**
  - Different indexing schemes
  - Real-valued representations for s=0

### 4. Performance Extensions

#### Caching
- **No automatic caching of recursion coefficients**
  - Could implement LRU cache for common ℓₘₐₓ values
  - Would speed up repeated transforms

#### Multi-threading
- **Limited multi-threading support**
  - Could add thread-parallel batch operations
  - Parallelize over multiple fields

### 5. Specialized Types

#### Symmetry Types
- **Real-valued spherical harmonics** (s=0)
  - Optimized storage (only m≥0)
  - Faster transforms

- **Axisymmetric fields** (m=0 only)
  - 1D instead of 2D transforms
  - Significant speedup

#### Field Types
- **Vector fields on the sphere**
  - Store spin-1 components together
  - Natural basis operations

- **Tensor fields**
  - Spin-2 fields for CMB, gravitational waves
  - Decomposition into E/B modes

### 6. Validation and Testing

#### Missing Test Infrastructure
- **No property-based testing**
  - Could use QuickCheck-style tests
  - Verify mathematical properties

- **No benchmark suite exposed**
  - Benchmarks exist but not easily accessible
  - Could add performance regression testing

### 7. Documentation Gaps

- **No beginner tutorial**
  - Quick start guide needed
  - Common use cases

- **No algorithm complexity documentation**
  - Time/space complexity for each method
  - Guidance on method selection

- **No mathematical background primer**
  - Spherical harmonics introduction
  - Wigner D matrix overview

## Recommendations for New Types

### High Priority

1. **`ModeWeights{T,S<:Integer}` struct**
   ```julia
   struct ModeWeights{T<:Complex, S<:Integer}
       data::Vector{T}
       s::S
       ℓₘᵢₙ::Int
       ℓₘₐₓ::Int
   end
   ```
   - Encapsulate mode weight arrays with metadata
   - Enable safe indexing: `modes[ℓ, m]`
   - Automatic validation

2. **`SphericalGrid{T<:Real}` abstract type**
   - Subtypes for different grid types
   - Unified interface for pixelizations
   - Enable grid-specific optimizations

3. **`PowerSpectrum{T<:Real}` struct**
   ```julia
   struct PowerSpectrum{T<:Real}
       Cℓ::Vector{T}
       ℓₘᵢₙ::Int
       ℓₘₐₓ::Int
   end
   ```
   - Store angular power spectra
   - Operations: +, -, *, smoothing

### Medium Priority

4. **`SphericalVectorField` and `SphericalTensorField`**
   - Composite types for multi-component fields
   - Automatic handling of spin components
   - Natural bases (spherical, Cartesian)

5. **`SymmetricModes` for real harmonics**
   - Optimized storage for s=0 case
   - Automatic symmetry enforcement

6. **`SamplingScheme` hierarchy**
   - Abstract type with concrete implementations
   - Standardize pixelization interface

### Lower Priority

7. **`TransformPlan` unification**
   - Common interface for all transform types
   - Could replace current SSHT hierarchy

8. **`RotationSequence` type**
   - Store sequences of rotations
   - Efficient composition
   - Useful for time-series analysis

## References

The package currently references:
- Xing et al. (2019) - ALF recursion
- Gumerov & Duraiswami (2015) - H recursion
- Reinecke & Seljebotn (2013) - Fast SHT algorithm
- Elahi et al. (2018) - Minimal sampling
- Kostelec & Rockmore (2008) - Lambda recursion
- Newman & Penrose (1966) - Spin operators
- Boyle (2016) - Modern treatment of SWSHs
