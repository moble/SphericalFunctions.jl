"""Supertype of storage for spin-spherical-harmonic transforms"""
abstract type SSHT{T<:Real} end

"""
    pixels(𝒯)

Return the spherical coordinates (θ, ϕ) on which the spin-weighted spherical harmonics are
evaluated.  See also [`rotors`](@ref), which provides the actual `Rotor`s on which they are
evaluated.
"""
function pixels end


"""
    rotors(𝒯)

Return the `Rotor`s on which the spin-weighted spherical harmonics are evaluated.  See also
[`pixels`](@ref), which provides the corresponding spherical coordinates.
"""
function rotors end


# mul!, ldiv!, synthesis, analysis, map2salm, salm2map


"""
    synthesis(𝒯, f̃)
    synthesis!(f, 𝒯, f̃)
    synthesize(𝒯, f̃)
    synthesize!(f, 𝒯, f̃)

Synthesize function values `f` from spin-weighted spherical-harmonic mode weights `f̃` using
the `SSHT` object `𝒯`.
"""
function synthesis end
function synthesis! end
function synthesize end
function synthesize! end

"""
    analysis(𝒯, f)
    analysis!(f̃, 𝒯, f)
    analyze(𝒯, f)
    analyze!(f̃, 𝒯, f)

Analyze function values `f` to obtain spin-weighted spherical-harmonic mode weights `f̃`
using the `SSHT` object `𝒯`.
"""
function analysis end
function analysis! end
function analyze end
function analyze! end


function Base.show(io::IO, 𝒯::SSHT)
    print(io, typeof(𝒯), "($(𝒯.s), $(𝒯.ℓₘₐₓ))")
end


inplaceable(s, ℓₘₐₓ, Rθϕ) = (ℓₘₐₓ + 1)^2 - s^2 == length(Rθϕ)


include("ssht/direct.jl")
include("ssht/minimal.jl")
include("ssht/reinecke_seljebotn.jl")
include("ssht/huffenberger_wandelt.jl")
