# The ASCII aliases of the names that are not ASCII, in the order of the `public` list in
# `SphericalFunctions.jl`.  They are public but not exported, since names such as `L2` are
# generic enough to clash with a user's own.  Each alias is the function or value itself,
# under another name, so it shares its identity, its `nameof` and its display.  The
# docstring of each function mentions its alias.  An alias of a value, unlike one of a
# function or a type, does not lead the documentation system to the value's docstring, so
# each alias of an operator has a short docstring of its own.

# The accessors of `indices/accessors.jl`.
const ell = ℓ
const ell_min = ℓₘᵢₙ
const ell_max = ℓₘₐₓ
const mp_max = m′ₘₐₓ
const mp_min = m′ₘᵢₙ
const m_max = mₘₐₓ
const m_min = mₘᵢₙ
const s_max = sₘₐₓ
const s_min = sₘᵢₙ
const Nr = Nᵣ

# The operators of `mode_weights/operators.jl`, and the change they make to the spin weight.
const Deltaspin = Δspin

"""
    L2

An ASCII alias of the operator [`L²`](@ref), for use where the Unicode name is inconvenient.
"""
const L2 = L²

"""
    Lplus

An ASCII alias of the operator [`L₊`](@ref), for use where the Unicode name is inconvenient.
"""
const Lplus = L₊

"""
    Lminus

An ASCII alias of the operator [`L₋`](@ref), for use where the Unicode name is inconvenient.
"""
const Lminus = L₋

"""
    R2

An ASCII alias of the operator [`R²`](@ref), for use where the Unicode name is inconvenient.
"""
const R2 = R²

"""
    Rplus

An ASCII alias of the operator [`R₊`](@ref), for use where the Unicode name is inconvenient.
"""
const Rplus = R₊

"""
    Rminus

An ASCII alias of the operator [`R₋`](@ref), for use where the Unicode name is inconvenient.
"""
const Rminus = R₋

"""
    eth

An ASCII alias of the operator [`ð`](@ref), for use where the Unicode name is inconvenient.
"""
const eth = ð

"""
    ethbar

An ASCII alias of the operator [`ð̄`](@ref), for use where the Unicode name is inconvenient.
"""
const ethbar = ð̄

# The calculator of `calculators/harmonics.jl` and the setters of `calculators/setters.jl`.
const slambdalmCalculator = sλlmCalculator
const set_beta! = set_β!
const set_theta! = set_θ!
